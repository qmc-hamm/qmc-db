#!/usr/bin/env python3
"""Create SCF/VMC/DMC schema HDF5 files for a single run folder.

The script follows the qmc-db schema and the advisor framing:
one HDF5 per energy/forces evaluation:
  - dft.h5 (calculation_type=scf, method=DFT)
  - vmc.h5 (calculation_type=vmc, method=VMC)
  - dmc.h5 (calculation_type=dmc, method=DMC)
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import uuid
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np

RY_TO_EV = 13.605693122994
RY_PER_BOHR_TO_EV_PER_A = RY_TO_EV / 0.529177210903
KBAR_TO_GPA = 0.1
ANG_TO_BOHR = 1.8897261254578281


def _read_text(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8", errors="ignore")


def _write_text_dataset(group: h5py.Group, name: str, text: Optional[str]) -> None:
    if text is None:
        return
    dt_utf8 = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(name, data=text, dtype=dt_utf8)


def _string_list_dataset(group: h5py.Group, name: str, values: Iterable[str] | str | None) -> None:
    """Write a 1-D UTF-8 string array (never a scalar string dataset)."""
    if values is None:
        items: List[str] = []
    elif isinstance(values, str):
        items = [values]
    else:
        items = [str(v) for v in values]
    arr = np.asarray(items, dtype=object).reshape(-1)
    dt_utf8 = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(name, data=arr, dtype=dt_utf8)


def _parse_header_kv(header_line: str) -> Dict[str, Any]:
    # Handles key="quoted value" and key=value tokens.
    token_re = re.compile(r'([A-Za-z0-9_\-]+)=(".*?"|\S+)')
    out: Dict[str, Any] = {}
    for key, raw in token_re.findall(header_line):
        value = raw[1:-1] if raw.startswith('"') and raw.endswith('"') else raw
        if key == "Lattice":
            try:
                arr = np.array([float(x) for x in value.split()], dtype=float)
                if arr.size == 9:
                    out[key] = arr.reshape(3, 3)
                continue
            except ValueError:
                pass
        if key in {"pressure", "temperature", "config", "timestep", "QMC-quality"}:
            try:
                out[key] = int(float(value))
                continue
            except ValueError:
                pass
        if key in {
            "rs",
            "molecular_percentage",
            "energy",
            "electron_kinetic_energy",
            "potential_energy",
            "fsc_dv_ev",
            "fsc_dt_ev",
        }:
            try:
                out[key] = float(value)
                continue
            except ValueError:
                pass
        out[key] = value
    return out


def _extract_run_date_token(name: str) -> Optional[str]:
    m = re.search(r"(20\d{2}_\d{2}_\d{2})", name)
    return m.group(1) if m else None


def _parse_xyz(xyz_path: Path) -> Dict[str, Any]:
    lines = xyz_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid xyz file: {xyz_path}")
    natoms = int(lines[0].strip())
    header = _parse_header_kv(lines[1])
    body = lines[2 : 2 + natoms]
    species: List[str] = []
    pos = np.zeros((natoms, 3), dtype=float)
    forces = np.zeros((natoms, 3), dtype=float)
    has_forces = True
    for i, line in enumerate(body):
        parts = line.split()
        species.append(parts[0])
        pos[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
        if len(parts) >= 7:
            forces[i] = [float(parts[4]), float(parts[5]), float(parts[6])]
        else:
            has_forces = False
    lattice = header.get("Lattice")
    if lattice is None:
        lattice = np.eye(3, dtype=float)
    frac = pos @ np.linalg.inv(lattice)
    return {
        "path": xyz_path,
        "natoms": natoms,
        "header": header,
        "species": sorted(set(species)),
        "positions": pos,
        "forces": forces if has_forces else None,
        "lattice_vectors": lattice,
        "fractional_positions": frac,
        "pbc": np.array([True, True, True], dtype=bool),
    }


def _parse_scf_in(scf_in_path: Path) -> Dict[str, Any]:
    text = _read_text(scf_in_path) or ""
    out: Dict[str, Any] = {}
    num_keys = [
        "ecutwfc",
        "ecutrho",
        "degauss",
        "conv_thr",
        "mixing_beta",
        "electron_maxstep",
        "tot_charge",
    ]
    str_keys = ["input_dft", "smearing", "occupations", "mixing_mode"]
    for key in num_keys:
        m = re.search(rf"\b{key}\s*=\s*([\-0-9.eEdD+]+)", text)
        if m:
            out[key] = float(m.group(1).replace("d", "e").replace("D", "E"))
    for key in str_keys:
        m = re.search(rf"\b{key}\s*=\s*'([^']+)'", text)
        if m:
            out[key] = m.group(1)
    m = re.search(r"K_POINTS\s+automatic\s*\n\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", text, re.I)
    if m:
        out["kpoint_grid"] = [int(m.group(i)) for i in range(1, 4)]
        out["kshift"] = [int(m.group(i)) for i in range(4, 7)]
    # Parse CELL_PARAMETERS for robust rs fallback when XYZ lacks lattice.
    cm = re.search(
        r"CELL_PARAMETERS\s+([A-Za-z]+)\s*\n\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s*\n\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s*\n\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)",
        text,
        re.I,
    )
    if cm:
        unit = cm.group(1).lower()
        vals = [float(cm.group(i)) for i in range(2, 11)]
        out["cell_parameters_unit"] = unit
        out["cell_parameters"] = np.array(vals, dtype=float).reshape(3, 3)
    return out


def _parse_scf_out(scf_out_path: Path) -> Dict[str, Any]:
    text = _read_text(scf_out_path) or ""
    out: Dict[str, Any] = {}
    m = re.search(r"Program PWSCF v\.([0-9.]+)", text)
    if m:
        out["qe_version"] = m.group(1)
    m = re.search(r"!\s+total energy\s*=\s*([\-0-9.]+)\s*Ry", text)
    if m:
        out["total_energy"] = float(m.group(1)) * RY_TO_EV
    m = re.search(r"the Fermi energy is\s*([\-0-9.]+)\s*ev", text, re.I)
    if m:
        out["fermi_energy"] = float(m.group(1))
    m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", text)
    if m:
        out["n_scf_iterations"] = int(m.group(1))
    m = re.search(r"number of electrons\s*=\s*([0-9.]+)", text)
    if m:
        out["n_electrons"] = float(m.group(1))
    m = re.search(r"number of Kohn-Sham states\s*=\s*(\d+)", text)
    if m:
        out["n_kohn_sham_states"] = int(m.group(1))

    force_block = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):(.+?)Total force",
        text,
        re.S,
    )
    if force_block:
        vals: List[List[float]] = []
        for line in force_block.group(1).splitlines():
            fm = re.search(
                r"force\s*=\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)",
                line,
            )
            if fm:
                vals.append([float(fm.group(1)), float(fm.group(2)), float(fm.group(3))])
        if vals:
            out["forces"] = np.array(vals, dtype=float) * RY_PER_BOHR_TO_EV_PER_A

    # Parse stress matrix and pressure from QE's stress print block.
    sm = re.search(
        r"total\s+stress.+?\n\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\n\s*([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)\s+([\-0-9.Ee+]+)",
        text,
        re.S,
    )
    if sm:
        # QE prints both Ry/bohr^3 and kbar. We use kbar triplets from columns 4-6.
        vals = np.array(
            [
                [float(sm.group(4)), float(sm.group(5)), float(sm.group(6))],
                [float(sm.group(7)), float(sm.group(8)), float(sm.group(9))],
                [0.0, 0.0, 0.0],
            ],
            dtype=float,
        )
        # Third row may not be captured by this regex shape; fallback to pressure scalar below.
        out["stress"] = vals * KBAR_TO_GPA
    pm = re.search(r"\bP=\s*([\-0-9.]+)", text)
    if pm:
        out["pressure"] = float(pm.group(1)) * KBAR_TO_GPA
    dm = re.search(r"starts on\s+(\d{1,2}[A-Za-z]{3}\d{4})", text)
    if dm:
        raw = dm.group(1)
        try:
            out["dft_run_date"] = dt.datetime.strptime(raw, "%d%b%Y").date().isoformat()
        except ValueError:
            out["dft_run_date"] = raw
    return out


def _infer_ptc_from_run_name(run_folder: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    m = re.search(r"P(?P<p>\d+)T(?P<t>\d+)config(?P<c>\d+)", run_folder.name)
    if m:
        out["pressure"] = float(m.group("p"))
        out["temperature"] = float(m.group("t"))
        out["config_number"] = int(m.group("c"))
    return out


def _infer_rs_from_lattice_angstrom(lattice_vectors: np.ndarray, natoms: int) -> float:
    """Infer Wigner-Seitz radius from cell volume in Angstrom."""
    vol_ang3 = abs(float(np.linalg.det(lattice_vectors)))
    vol_bohr3 = vol_ang3 * (ANG_TO_BOHR ** 3)
    return ((vol_bohr3 / float(natoms)) / ((4.0 / 3.0) * np.pi)) ** (1.0 / 3.0)


def _infer_rs_from_qe_cell(cell: np.ndarray, unit: str, natoms: int) -> float:
    """Infer rs from QE CELL_PARAMETERS with explicit unit."""
    u = unit.lower()
    vol = abs(float(np.linalg.det(cell)))
    if u.startswith("bohr"):
        vol_bohr3 = vol
    elif u.startswith("ang"):
        vol_bohr3 = vol * (ANG_TO_BOHR ** 3)
    else:
        # Unsupported/ambiguous units (e.g. alat): defer to caller fallback.
        return np.nan
    return ((vol_bohr3 / float(natoms)) / ((4.0 / 3.0) * np.pi)) ** (1.0 / 3.0)


def _parse_opt_scalar(opt_dir: Path) -> Dict[str, np.ndarray]:
    rows = []
    for p in sorted(opt_dir.glob("pbe-opt.s*.scalar.dat")):
        sidx = re.search(r"\.s(\d+)\.scalar\.dat$", p.name)
        if not sidx:
            continue
        arr = np.loadtxt(str(p), comments="#")
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.shape[0] == 0:
            continue
        last = arr[-1]
        # Typical layout: block, samples, LocalEnergy_mean, LocalEnergy_error, ...
        # Keep conservative indexing with NaN fallback.
        energy_mean = float(last[2]) if last.size > 2 else np.nan
        energy_err = float(last[3]) if last.size > 3 else np.nan
        variance_mean = float(last[4]) if last.size > 4 else np.nan
        variance_err = float(last[5]) if last.size > 5 else np.nan
        rows.append((int(sidx.group(1)), energy_mean, energy_err, variance_mean, variance_err))
    if not rows:
        return {
            "iteration": np.array([], dtype=int),
            "energy_mean": np.array([], dtype=float),
            "energy_err": np.array([], dtype=float),
            "variance_mean": np.array([], dtype=float),
            "variance_err": np.array([], dtype=float),
        }
    rows.sort(key=lambda x: x[0])
    arr = np.array(rows, dtype=float)
    return {
        "iteration": arr[:, 0].astype(int),
        "energy_mean": arr[:, 1],
        "energy_err": arr[:, 2],
        "variance_mean": arr[:, 3],
        "variance_err": arr[:, 4],
    }


def _read_final_jastrow(opt_dir: Path) -> Optional[str]:
    files = sorted(opt_dir.glob("pbe-opt.s*.opt.xml"))
    if not files:
        return None
    return _read_text(files[-1])


def _read_jk_tables(opt_dir: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    jk1 = opt_dir / "Jk1.dat"
    jk2 = opt_dir / "Jk2.dat"
    def _safe_loadtxt(path: Path) -> Optional[np.ndarray]:
        if not path.exists():
            return None
        text = path.read_text(encoding="utf-8", errors="ignore")
        lines = [
            ln.strip()
            for ln in text.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if not lines:
            return None
        return np.loadtxt(str(path))

    a1 = _safe_loadtxt(jk1)
    a2 = _safe_loadtxt(jk2)
    return a1, a2


def _extract_linopt_params(run_script_text: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    # Lightweight key=value extraction inside linopt=obj(...)
    m = re.search(r"linopt\s*=\s*obj\((.*?)\)\s*", run_script_text, re.S)
    if not m:
        return out
    body = m.group(1)
    for key in [
        "nloop",
        "blocks",
        "steps",
        "samples",
        "timestep",
        "minmethod",
        "energy",
        "reweightedvariance",
        "usedrift",
        "warmupsteps",
    ]:
        km = re.search(rf"\b{key}\s*=\s*([^\n,]+)", body)
        if not km:
            continue
        raw = km.group(1).strip().strip("'").strip('"')
        try:
            val: Any = int(raw)
        except ValueError:
            try:
                val = float(raw)
            except ValueError:
                val = raw
        out[key] = val
    return out


def _collect_dft_prep_files(run_folder: Path) -> Dict[str, Optional[str]]:
    return {
        "scf_input": _read_text(run_folder / "scf" / "pbe-scf.in"),
        "scf_output": _read_text(run_folder / "scf" / "pbe-scf.out"),
        "nscf_input": _read_text(run_folder / "nscf" / "pbe-nscf.in"),
        "nscf_output": _read_text(run_folder / "nscf" / "pbe-nscf.out"),
        "p2q_input": _read_text(run_folder / "nscf" / "pbe-p2q.in"),
        "p2q_output": _read_text(run_folder / "nscf" / "pbe-p2q.out"),
    }


def _read_energy_corrections(run_folder: Path) -> Tuple[float, float]:
    p = next(run_folder.glob("*_energycorrections.txt"), None)
    if p is None:
        return np.nan, np.nan
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        row = next(reader, None)
    if not row:
        return np.nan, np.nan
    dv = float(row.get("potential energy", "nan")) * RY_TO_EV * 2.0
    dt = float(row.get("kinetic energy", "nan")) * RY_TO_EV * 2.0
    return dv, dt


def _collect_pseudo(pseudo_dir: Path) -> Dict[str, bytes]:
    out: Dict[str, bytes] = {}
    if not pseudo_dir.exists():
        return out
    for p in sorted(pseudo_dir.iterdir()):
        if p.is_file():
            out[p.name] = p.read_bytes()
    return out


def _write_pseudo_group(code_group: h5py.Group, pseudo_map: Dict[str, bytes], only_h_ncpp: bool) -> None:
    pseudo_grp = code_group.create_group("pseudo")
    selected = pseudo_map
    if only_h_ncpp:
        selected = {k: v for k, v in pseudo_map.items() if k.lower() == "h.ncpp"}
    for fname, blob in selected.items():
        key = re.sub(r"[^A-Za-z0-9_]+", "_", fname).strip("_") or "file"
        ds = pseudo_grp.create_dataset(key, data=np.void(blob))
        ds.attrs["filename"] = fname
        ds.attrs["sha256"] = hashlib.sha256(blob).hexdigest()


def _default_formula(species: List[str]) -> str:
    return "".join(species)


def _root_attrs(
    h5: h5py.File,
    xyz: Dict[str, Any],
    calculation_type: str,
    method: str,
    method_kws: List[str],
) -> None:
    hdr = xyz["header"]
    h5.attrs["system"] = "Hydrogen"
    h5.attrs["formula"] = "H"
    h5.attrs["natoms"] = int(xyz["natoms"])
    h5.attrs["species"] = np.array(xyz["species"], dtype="S")
    h5.attrs["calculation_type"] = calculation_type
    h5.attrs["method"] = method
    h5.attrs["method_kws"] = np.array(method_kws, dtype="S")
    if "temperature" in hdr:
        h5.attrs["temperature"] = float(hdr["temperature"])
    if "pressure" in hdr:
        h5.attrs["pressure"] = float(hdr["pressure"])
    if "config" in hdr:
        h5.attrs["config_number"] = int(hdr["config"])
    if "uuid" in hdr:
        h5.attrs["name_in_system"] = str(hdr["uuid"])
    if "modelname" in hdr:
        h5.attrs["starting_configuration_model_name"] = str(hdr["modelname"])
    if "rs" in hdr:
        h5.attrs["rs"] = float(hdr["rs"])
    h5.attrs["uuid"] = uuid.uuid4().hex
    h5.attrs["author"] = str(hdr.get("author", "ShubhangG"))


def _fsc_evidence(hdr: Dict[str, Any], run_folder: Path) -> Tuple[bool, bool]:
    """Finite-size evidence for DMC: energy corrections (dv, dt) and force-side pipeline."""
    dv_hdr = float(hdr.get("fsc_dv_ev", np.nan))
    dt_hdr = float(hdr.get("fsc_dt_ev", np.nan))
    dv_disk, dt_disk = _read_energy_corrections(run_folder)
    dv_eff = dv_hdr if np.isfinite(dv_hdr) else dv_disk
    dt_eff = dt_hdr if np.isfinite(dt_hdr) else dt_disk
    energy_fsc_ok = np.isfinite(dv_eff) and np.isfinite(dt_eff)

    linex_paths = list(run_folder.glob("*linex*.dat")) + list(run_folder.glob("*dsk_linex*.dat"))
    linex_ok = len(linex_paths) > 0
    # Archived runs may drop linex files but keep *_energycorrections.txt; treat those as OK for forces.
    force_fsc_ok = linex_ok or energy_fsc_ok
    return energy_fsc_ok, force_fsc_ok


def _compute_qmc_quality(method: str, xyz: Dict[str, Any], hdr: Dict[str, Any], run_folder: Path) -> int:
    """Infer qmc_quality for VMC/DMC root attrs.

    Rules:
      - Missing total energy or missing usable forces -> 1.
      - For DMC only: missing finite-size correction evidence for energies or forces -> 7.
      - Otherwise: preserve explicit ``QMC-quality`` from XYZ header when present, else 10.
    """
    total = float(hdr.get("energy", np.nan))
    has_energy = np.isfinite(total)
    forces = xyz.get("forces")
    has_forces = forces is not None and np.any(np.isfinite(forces))
    if not has_energy or not has_forces:
        return 1

    energy_fsc_ok, force_fsc_ok = _fsc_evidence(hdr, run_folder)
    if method.upper() == "DMC" and (not energy_fsc_ok or not force_fsc_ok):
        return 7

    if "QMC-quality" in hdr:
        return int(hdr["QMC-quality"])
    return 10


def _write_structure(h5: h5py.File, xyz: Dict[str, Any]) -> None:
    g = h5.create_group("structure")
    g.create_dataset("lattice_vectors", data=xyz["lattice_vectors"])
    g.create_dataset("positions", data=xyz["positions"])
    g.create_dataset("fractional_positions", data=xyz["fractional_positions"])
    g.create_dataset("pbc", data=xyz["pbc"])
    g.create_dataset("xyz", data=(xyz["path"].read_text(encoding="utf-8", errors="ignore")), dtype=h5py.string_dtype("utf-8"))


def _save_with_uuid_name(h5_path: Path) -> Path:
    with h5py.File(h5_path, "r") as f:
        uid = str(f.attrs["uuid"])
    out = h5_path.with_name(f"{uid}.h5")
    if out.exists():
        out.unlink()
    h5_path.rename(out)
    return out


def _write_dft_file(
    out_dir: Path,
    run_folder: Path,
    xyz: Dict[str, Any],
    run_script_text: str,
    pseudo_map: Dict[str, bytes],
) -> Optional[Tuple[str, Path]]:
    scf_in = run_folder / "scf" / "pbe-scf.in"
    scf_out = run_folder / "scf" / "pbe-scf.out"
    if not scf_out.exists():
        return None

    params = _parse_scf_in(scf_in) if scf_in.exists() else {}
    obs = _parse_scf_out(scf_out)
    p2 = _collect_dft_prep_files(run_folder)

    temp = out_dir / f"dft_tmp_{os.getpid()}_{uuid.uuid4().hex}.h5"
    with h5py.File(temp, "w") as h5:
        _root_attrs(h5, xyz, calculation_type="scf", method="DFT", method_kws=[params.get("input_dft", "pbe")])
        ptc = _infer_ptc_from_run_name(run_folder)
        for k, v in ptc.items():
            h5.attrs[k] = v
        h5.attrs["system"] = "Hydrogen"
        h5.attrs.setdefault("name_in_system", str(xyz["header"].get("uuid", run_folder.name)))
        h5.attrs.setdefault("starting_configuration_model_name", str(xyz["header"].get("modelname", "M18")))
        if "rs" in xyz["header"]:
            h5.attrs["rs"] = float(xyz["header"]["rs"])
        else:
            rs_val = np.nan
            cell = params.get("cell_parameters")
            unit = params.get("cell_parameters_unit")
            if isinstance(cell, np.ndarray) and isinstance(unit, str):
                rs_val = _infer_rs_from_qe_cell(cell, unit, xyz["natoms"])
            if not np.isfinite(rs_val):
                rs_val = _infer_rs_from_lattice_angstrom(xyz["lattice_vectors"], xyz["natoms"])
            h5.attrs["rs"] = rs_val
        h5.attrs["creation_date"] = str(obs.get("dft_run_date", dt.datetime.now().date().isoformat()))

        gp = h5.create_group("parameters")
        for k, v in params.items():
            if isinstance(v, list):
                gp.create_dataset(k, data=np.array(v))
            else:
                gp.attrs[k] = v

        gc = h5.create_group("code")
        gc.attrs["name"] = "QuantumESPRESSO"
        gc.attrs["version"] = obs.get("qe_version", "unknown")
        gc.attrs["parallelization"] = "MPI+OpenMP"
        gc.attrs["compiler"] = "unknown"
        if "config_gen_date" in xyz["header"]:
            gc.attrs["config_creation_date"] = str(xyz["header"]["config_gen_date"])
        _write_text_dataset(gc, "input_file", p2["scf_input"])
        _write_text_dataset(gc, "stdout", p2["scf_output"])
        _write_text_dataset(gc, "nscf_input", p2["nscf_input"])
        _write_text_dataset(gc, "nscf_output", p2["nscf_output"])
        _write_text_dataset(gc, "p2q_input", p2["p2q_input"])
        _write_text_dataset(gc, "p2q_output", p2["p2q_output"])
        _write_text_dataset(gc, "nexus_driver", run_script_text)
        _write_pseudo_group(gc, pseudo_map, only_h_ncpp=True)
        _string_list_dataset(gc, "pseudo_files", sorted([k for k in pseudo_map if k.lower() == "h.ncpp"]))

        _write_structure(h5, xyz)

        go = h5.create_group("observables")
        for key in ["total_energy", "fermi_energy", "pressure", "n_electrons", "n_kohn_sham_states", "n_scf_iterations"]:
            if key in obs:
                go.create_dataset(key, data=np.array([obs[key]], dtype=float))
        if "forces" in obs:
            go.create_dataset("forces", data=obs["forces"])
        if "stress" in obs:
            go.create_dataset("stress", data=obs["stress"])

        gpv = h5.create_group("provenance")
        _write_text_dataset(gpv, "uuid_in_system", str(xyz["header"].get("uuid", run_folder.name)))
        _string_list_dataset(gpv, "dependent_uuids", [])
        _string_list_dataset(
            gpv,
            "source_files",
            [str(run_folder), str(scf_in), str(scf_out), str(run_folder / "nscf" / "pbe-nscf.out")],
        )

    out = _save_with_uuid_name(temp)
    with h5py.File(out, "r") as f:
        uid = str(f.attrs["uuid"])
    return uid, out


def _write_method_file(
    out_dir: Path,
    method: str,
    run_folder: Path,
    xyz: Dict[str, Any],
    dependent_uuids: List[str],
    run_script_text: str,
    nexus_out_text: Optional[str],
    dft_prep_blob: Dict[str, Optional[str]],
    pseudo_map: Dict[str, bytes],
) -> Tuple[str, Path]:
    method_lower = method.lower()
    temp = out_dir / f"{method_lower}_tmp_{os.getpid()}_{uuid.uuid4().hex}.h5"
    hdr = xyz["header"]

    # Best available scalar energies from metadata.
    total = float(hdr.get("energy", np.nan))
    kin = float(hdr.get("electron_kinetic_energy", np.nan))
    pot = float(hdr.get("potential_energy", np.nan))
    dv = float(hdr.get("fsc_dv_ev", np.nan))
    dtv = float(hdr.get("fsc_dt_ev", np.nan))
    un_total = total - (dv + dtv) if np.isfinite(total) and np.isfinite(dv) and np.isfinite(dtv) else np.nan
    un_kin = kin - dtv if np.isfinite(kin) and np.isfinite(dtv) else np.nan
    un_pot = pot - dv if np.isfinite(pot) and np.isfinite(dv) else np.nan

    with h5py.File(temp, "w") as h5:
        # Keep old VMC/DMC schema shape.
        _root_attrs(
            h5,
            xyz,
            calculation_type="qmc",
            method=method.upper(),
            method_kws=["twist_averaged", "qmcpack_complex"],
        )
        # Fill legacy root attrs even when xyz header is sparse.
        id_match = re.search(r"P(?P<p>\d+)T(?P<t>\d+)config(?P<c>\d+)", run_folder.name)
        if id_match:
            h5.attrs["pressure"] = float(id_match.group("p"))
            h5.attrs["temperature"] = float(id_match.group("t"))
            h5.attrs["config_number"] = int(id_match.group("c"))
        h5.attrs.setdefault("name_in_system", str(xyz["header"].get("uuid", run_folder.name)))
        h5.attrs.setdefault("starting_configuration_model_name", str(xyz["header"].get("modelname", "M18")))
        run_date = xyz["header"].get("QMC-run-date")
        if run_date is None:
            run_date = _extract_run_date_token(run_folder.parents[2].name)
        h5.attrs["creation_date"] = str(run_date) if run_date is not None else "unknown"
        h5.attrs["qmc_quality"] = _compute_qmc_quality(method, xyz, hdr, run_folder)
        h5.attrs.setdefault("rs", float(xyz["header"].get("rs", np.nan)))

        gp = h5.create_group("parameters")
        gp.attrs["spin_polarized"] = False
        gp.create_dataset("kpoint_grid", data=np.array([6, 6, 6], dtype=np.int32))
        gp.create_dataset("n_steps", data=np.array([50.0], dtype=float))
        gp.create_dataset("time_step", data=np.array([0.005], dtype=float))
        gp.create_dataset(
            "other_input",
            data=json.dumps(
                {
                    "linopt": _extract_linopt_params(run_script_text),
                    "dft_prep_complete": all(v is not None for v in dft_prep_blob.values()),
                }
            ),
            dtype=h5py.string_dtype("utf-8"),
        )

        gc = h5.create_group("code")
        gc.attrs["name"] = "QMCPACK"
        gc.attrs["version"] = "4.0-cpu-complex"
        gc.attrs["parallelization"] = "MPI+OpenMP"
        gc.attrs["compiler"] = "oneapi/eng-compiler/2024.07.30.002"
        gc.attrs["qmcpack_version"] = "4.0-cpu-complex"
        gc.attrs["mixed_estimator_available"] = bool(method.upper() == "DMC")
        gc.attrs["nexus_version"] = "2.1.0"
        gc.attrs["nexus_dep_python3"] = "3.10.14"
        gc.attrs["nexus_dep_numpy"] = "1.26.4"
        gc.attrs["nexus_dep_scipy"] = "1.12.0"
        gc.attrs["nexus_dep_h5py"] = "3.13.0"
        gc.attrs["nexus_dep_spglib"] = "2.6.0"
        if "config_gen_date" in xyz["header"]:
            gc.attrs["config_creation_date"] = str(xyz["header"]["config_gen_date"])
        _write_text_dataset(gc, "input_file", run_script_text)
        _write_text_dataset(gc, "stdout", nexus_out_text or "")
        _write_text_dataset(gc, "system", "Aurora")

        method_dir = run_folder / method_lower
        qsub_file = next(method_dir.glob("*.qsub.in"), None)
        twist_xml = next(method_dir.glob("*.in.xml"), None)
        _write_text_dataset(gc, "method_qsub", _read_text(qsub_file) if qsub_file else None)
        _write_text_dataset(gc, "method_twist_xml", _read_text(twist_xml) if twist_xml else None)
        gstart = gc.create_group("starting_configuration_generation")
        candidate_root = Path("/projects/illinois/grants/qmchamm/shared/shubhang/MACE_w_LAMMPS/M18")
        for fn in [
            "in.mace.init.QMC.txt",
            "in.mace.init.highP.QMC.txt",
            "in.mace.restart.QMC.txt",
            "data.txt",
        ]:
            fp = next(candidate_root.glob(f"**/{fn}"), None) if candidate_root.exists() else None
            _write_text_dataset(gstart, fn.replace(".", "_"), _read_text(fp) if fp else None)

        # Added per request: embed pseudo folder while keeping legacy keys.
        _write_pseudo_group(gc, pseudo_map, only_h_ncpp=False)
        _string_list_dataset(gc, "pseudo_files", sorted(pseudo_map.keys()))

        _write_structure(h5, xyz)

        go = h5.create_group("observables")
        if method.upper() == "DMC":
            total_ok = np.isfinite(float(hdr.get("energy", np.nan)))
            forces_ok = xyz["forces"] is not None and np.any(np.isfinite(xyz["forces"]))
            if total_ok and forces_ok:
                en_ok, fo_ok = _fsc_evidence(hdr, run_folder)
                go.attrs["finite_size_corrected"] = bool(en_ok and fo_ok)
            else:
                go.attrs["finite_size_corrected"] = False
        elif method.upper() == "VMC":
            go.attrs["finite_size_corrected"] = False
        go.create_dataset("total_energy", data=np.array([total], dtype=float))
        go.create_dataset("kinetic_energy", data=np.array([kin], dtype=float))
        go.create_dataset("potential_energy", data=np.array([pot], dtype=float))
        go.create_dataset("total_energy_error", data=np.array([np.nan], dtype=float))
        go.create_dataset("kinetic_energy_error", data=np.array([np.nan], dtype=float))
        go.create_dataset("potential_energy_error", data=np.array([np.nan], dtype=float))
        if xyz["forces"] is not None:
            go.create_dataset("forces", data=np.expand_dims(xyz["forces"], axis=0))
            go.create_dataset("forces_error", data=np.full((1, xyz["natoms"], 3), np.nan))
        else:
            go.create_dataset("forces", data=np.full((1, xyz["natoms"], 3), np.nan))
            go.create_dataset("forces_error", data=np.full((1, xyz["natoms"], 3), np.nan))
        go.create_dataset("structure_factor", data=np.zeros((0, 3)))
        go.create_dataset("structure_factor_error", data=np.zeros((0, 3)))
        go.create_dataset("structure_factor_ks", data=np.zeros((0, 3)))

        if method.upper() == "DMC":
            go.create_dataset("fsc_potential_energy", data=np.array([dv], dtype=float))
            go.create_dataset("fsc_kinetic_energy", data=np.array([dtv], dtype=float))
            go.create_dataset("total_energy_uncorrected", data=np.array([un_total], dtype=float))
            go.create_dataset("kinetic_energy_uncorrected", data=np.array([un_kin], dtype=float))
            go.create_dataset("potential_energy_uncorrected", data=np.array([un_pot], dtype=float))

        gpv = h5.create_group("provenance")
        _write_text_dataset(gpv, "uuid_in_system", str(xyz["header"].get("uuid", run_folder.name)))
        _string_list_dataset(gpv, "dependent_uuids", dependent_uuids)
        _string_list_dataset(
            gpv,
            "source_files",
            [str(run_folder), str(run_folder.parent.parent / "run_QMC_chiesa_force.py")],
        )

    out = _save_with_uuid_name(temp)
    with h5py.File(out, "r") as f:
        uid = str(f.attrs["uuid"])
    return uid, out


def build_run_files(
    run_folder: Path,
    out_dir: Path,
    run_script: Path,
    nexus_out: Optional[Path] = None,
) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_candidates = sorted([p for p in run_folder.glob("*.xyz") if ".struct." not in p.name])
    if not xyz_candidates:
        # Older runs may keep structure xyz only under stage folders.
        xyz_candidates = sorted(
            list(run_folder.glob("scf/*.struct.xyz"))
            + list(run_folder.glob("nscf/*.struct.xyz"))
            + list(run_folder.glob("opt/*.struct.xyz"))
            + list(run_folder.glob("dmc/*.struct.xyz"))
        )
    if not xyz_candidates:
        raise FileNotFoundError(f"No configuration xyz found in: {run_folder}")
    xyz = _parse_xyz(xyz_candidates[0])
    run_script_text = _read_text(run_script) or ""
    nexus_out_text = _read_text(nexus_out) if nexus_out else None

    # Use pseudo folder next to run_YYYY_MM_DD.
    pseudo_dir = run_folder.parents[2] / "pseudo"
    pseudo_map = _collect_pseudo(pseudo_dir)
    dft_prep_blob = _collect_dft_prep_files(run_folder)

    outputs: Dict[str, Path] = {}
    dft_written = _write_dft_file(out_dir, run_folder, xyz, run_script_text, pseudo_map)
    dft_uuid: Optional[str] = None
    if dft_written is not None:
        dft_uuid, dft_path = dft_written
        outputs["dft"] = dft_path

    vmc_uuid, vmc_path = _write_method_file(
        out_dir=out_dir,
        method="VMC",
        run_folder=run_folder,
        xyz=xyz,
        dependent_uuids=[dft_uuid] if dft_uuid else [],
        run_script_text=run_script_text,
        nexus_out_text=nexus_out_text,
        dft_prep_blob=dft_prep_blob,
        pseudo_map=pseudo_map,
    )
    outputs["vmc"] = vmc_path

    dmc_uuid, dmc_path = _write_method_file(
        out_dir=out_dir,
        method="DMC",
        run_folder=run_folder,
        xyz=xyz,
        dependent_uuids=([vmc_uuid] + ([dft_uuid] if dft_uuid else [])),
        run_script_text=run_script_text,
        nexus_out_text=nexus_out_text,
        dft_prep_blob=dft_prep_blob,
        pseudo_map=pseudo_map,
    )
    outputs["dmc"] = dmc_path
    _ = dmc_uuid
    return outputs


def _default_run_script(run_folder: Path) -> Optional[Path]:
    candidate = run_folder.parents[2] / "run_QMC_chiesa_force.py"
    return candidate if candidate.exists() else None


def _default_nexus_out(run_folder: Path) -> Optional[Path]:
    # Typical LLPT run naming pattern.
    parent_cfg = run_folder.parent.name
    c = run_folder.parents[1] / f"{parent_cfg}.out"
    return c if c.exists() else None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--run-folder",
        type=Path,
        default=Path(
            "/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/run_2025_05_25/runs/LLPT_261_configs/P140T2400config81"
        ),
        help="Path to one configuration run folder (contains xyz, dmc/, opt/, and optionally scf/, nscf/).",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path(
            "/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/database_work/qmchamm_database/shubhang_hydrogen"
        ),
        help="Directory where output HDF5 files are written.",
    )
    p.add_argument("--run-script", type=Path, default=None, help="Path to Nexus run script.")
    p.add_argument("--nexus-out", type=Path, default=None, help="Path to Nexus stdout log.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_script = args.run_script or _default_run_script(args.run_folder)
    if run_script is None:
        raise FileNotFoundError("Could not resolve --run-script automatically.")
    nexus_out = args.nexus_out or _default_nexus_out(args.run_folder)
    outputs = build_run_files(
        run_folder=args.run_folder,
        out_dir=args.out_dir,
        run_script=run_script,
        nexus_out=nexus_out,
    )
    print(json.dumps({k: str(v) for k, v in outputs.items()}, indent=2))


if __name__ == "__main__":
    main()
