#!/usr/bin/env python3
"""Re-ingest Aurora/SG QMC configs whose schema HDF5 files are empty shells.

Targets configuration folders that were written from broken struct.xyz files
(no lattice / no energies in the header). For each config this script:

  1. Reads twist-averaged observables from ``qmc_scalardata_*.h5``.
  2. Reads finite-size corrections from ``*_energycorrections.txt``.
  3. Reconstructs the cell from QMCPACK ``dmc/*.in.xml`` (lattice in bohr).
  4. Patches the existing VMC and DMC ``{uuid}.h5`` files in place (UUID kept).

Typical use (the 1,640 high-T configs with missing rs / qmc_quality=1):

  python reingest_scalardata_h5.py --dry-run
  python reingest_scalardata_h5.py --workers 8
  python reingest_scalardata_h5.py --single /path/to/P125T3800config100

After a batch run, rebuild the global index and refresh coverage plots.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import math
import re
import sys
import traceback
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from create_hdf5_schema import (  # noqa: E402
    RY_TO_EV,
    _compute_qmc_quality,
    _extract_run_date_token,
    _fsc_evidence,
    _read_energy_corrections,
    _read_text,
)

BOHR_TO_ANG = 0.529177210903
HARTREE_TO_EV = 27.2114
EV_PER_BOHR_TO_EV_PER_ANG = HARTREE_TO_EV / BOHR_TO_ANG

DEFAULT_AURORA = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
LEDGER_NAME = "build_database_ledger.txt"
REINGEST_LEDGER = "reingest_scalardata_ledger.txt"
CONFIG_RE = re.compile(r"^P\d+T\d+config\d+$")


def rs_from_cell_angstrom(cell_ang: np.ndarray, natoms: int) -> float:
    v = abs(float(np.linalg.det(cell_ang)))
    r_ang = (3.0 * v / natoms / (4.0 * math.pi)) ** (1.0 / 3.0)
    return r_ang / BOHR_TO_ANG


def _parse_pt_config(name: str) -> Tuple[float, float, int]:
    m = re.match(r"P(?P<p>\d+)T(?P<t>\d+)config(?P<c>\d+)", name)
    if not m:
        raise ValueError(f"not a P*T*config folder: {name}")
    return float(m.group("p")), float(m.group("t")), int(m.group("c"))


def _force_col_key(name: str) -> Tuple[int, int]:
    m = re.match(r"force_(\d+)_(\d+)$", name)
    return (int(m.group(1)), int(m.group(2))) if m else (10**9, 10**9)


def extract_twist_avg_values(fname: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Twist-averaged block means/errors (blocks 32+), no qharv dependency."""
    dat = pd.read_hdf(fname)
    numeric_cols = dat.columns.difference(["calculation", "twist_number"])
    grp = dat.groupby(["calculation", "twist_number"])[numeric_cols]
    mdf = grp.apply(lambda g: g.iloc[32:].mean())
    edf = grp.apply(
        lambda g: g.iloc[32:].apply(
            lambda col: float(col.std(ddof=1) / math.sqrt(len(col))) if len(col) >= 2 else 0.0
        )
    )
    weights = np.ones(len(edf))
    weights_series = pd.Series(weights, index=mdf.index)
    weighted_means = (
        mdf.mul(weights_series, axis=0)
        .groupby(level="calculation")
        .sum()
        .div(weights_series.groupby("calculation").sum(), axis=0)
    )
    weighted_errors = (
        np.sqrt((edf**2).mul(weights_series**2, axis=0).groupby(level="calculation").sum())
        .div(weights_series.groupby("calculation").sum(), axis=0)
    )
    return weighted_means, weighted_errors


def _lattice_from_qmcpack_xml(run_folder: Path) -> np.ndarray:
    """Return lattice vectors in Angstrom."""
    candidates = sorted((run_folder / "dmc").glob("*.in.xml"))
    if not candidates:
        candidates = sorted((run_folder / "opt").glob("*.cont.xml"))
    if not candidates:
        raise FileNotFoundError(f"no QMCPACK xml in {run_folder}")
    root = ET.parse(candidates[0]).getroot()
    for elem in root.iter():
        if elem.tag.endswith("parameter") and elem.get("name") == "lattice":
            vals = [float(x) for x in elem.text.split()]
            cell_bohr = np.array(vals, dtype=float).reshape(3, 3)
            return cell_bohr * BOHR_TO_ANG
    raise ValueError(f"lattice parameter missing in {candidates[0]}")


def _positions_from_struct_xyz(run_folder: Path) -> np.ndarray:
    candidates = sorted(
        list(run_folder.glob("dmc/*.struct.xyz"))
        + list(run_folder.glob("opt/*.struct.xyz"))
        + [p for p in run_folder.glob("*.xyz") if ".struct." not in p.name]
    )
    if not candidates:
        raise FileNotFoundError(f"no struct xyz in {run_folder}")
    lines = candidates[0].read_text(encoding="utf-8", errors="ignore").splitlines()
    nat = int(lines[0].strip())
    pos = np.zeros((nat, 3), dtype=float)
    for i, line in enumerate(lines[2 : 2 + nat]):
        parts = line.split()
        pos[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return pos


def _build_extxyz(natoms: int, header_kv: Dict[str, Any], positions: np.ndarray) -> str:
    lattice = header_kv.get("Lattice")
    lat_str = ""
    if lattice is not None:
        lat_str = ' Lattice="' + " ".join(f"{x:.8f}" for x in np.asarray(lattice).reshape(-1)) + '"'
    props = " ".join(f'{k}="{v}"' for k, v in header_kv.items() if k != "Lattice")
    lines = [str(natoms), f"Properties=species:S:1:pos:R:3{lat_str} {props}".strip()]
    for i in range(natoms):
        lines.append(f"H {positions[i,0]:.8f} {positions[i,1]:.8f} {positions[i,2]:.8f}")
    return "\n".join(lines) + "\n"


def _scalar_h5_path(run_folder: Path) -> Path:
    matches = sorted(run_folder.glob("qmc_scalardata*.h5"))
    if not matches:
        raise FileNotFoundError(f"no qmc_scalardata_*.h5 in {run_folder}")
    return matches[0]


def _method_h5_files(run_folder: Path) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    for p in sorted(run_folder.glob("*.h5")):
        if p.stat().st_size > 5_000_000:
            continue
        try:
            with h5py.File(p, "r") as f:
                method = str(f.attrs.get("method", "")).upper()
        except Exception:
            continue
        if method in {"VMC", "DMC"}:
            out[method] = p
    if "VMC" not in out or "DMC" not in out:
        raise FileNotFoundError(f"expected VMC+DMC schema h5 in {run_folder}, found {list(out)}")
    return out


def needs_reingest(h5_path: Path) -> bool:
    with h5py.File(h5_path, "r") as f:
        if str(f.attrs.get("method", "")).upper() != "DMC":
            return False
        qq = f.attrs.get("qmc_quality")
        if qq is not None and int(qq) >= 10:
            te = f["observables/total_energy"][()]
            if np.isfinite(np.array(te)[0]):
                return False
        te = f["observables/total_energy"][()]
        return not np.isfinite(np.array(te)[0])


def gather_reingest_payload(run_folder: Path) -> Dict[str, Any]:
    pressure, temperature, config_number = _parse_pt_config(run_folder.name)
    lattice = _lattice_from_qmcpack_xml(run_folder)
    positions = _positions_from_struct_xyz(run_folder)
    natoms = positions.shape[0]
    rs = rs_from_cell_angstrom(lattice, natoms)
    frac = positions @ np.linalg.inv(lattice)

    means, errs = extract_twist_avg_values(_scalar_h5_path(run_folder))
    dv, dt = _read_energy_corrections(run_folder)
    has_fsc = np.isfinite(dv) and np.isfinite(dt)

    force_cols = sorted(
        [c for c in means.columns if str(c).startswith("force_")],
        key=_force_col_key,
    )

    def row_forces(method: str) -> np.ndarray:
        if not force_cols:
            return np.full((natoms, 3), np.nan)
        vals = means.loc[method, force_cols].to_numpy(dtype=float)
        return vals.reshape(-1, 3) * HARTREE_TO_EV / BOHR_TO_ANG

    vmc_forces = row_forces("VMC") if "VMC" in means.index else np.full((natoms, 3), np.nan)
    dmc_forces_raw = row_forces("DMC") if "DMC" in means.index else np.full((natoms, 3), np.nan)
    if np.all(np.isfinite(vmc_forces)) and np.all(np.isfinite(dmc_forces_raw)):
        dmc_forces = 2.0 * dmc_forces_raw - vmc_forces
        dmc_forces -= dmc_forces.mean(axis=0)
    else:
        dmc_forces = np.full((natoms, 3), np.nan)

    run_date = _extract_run_date_token(run_folder.parents[2].name) or "unknown"
    header_base = {
        "pressure": pressure,
        "temperature": temperature,
        "config": config_number,
        "rs": rs,
        "author": "ShubhangG",
        "modelname": "M18",
        "uuid": run_folder.name,
        "QMC-run-date": run_date,
        "fsc_dv_ev": dv,
        "fsc_dt_ev": dt,
        "Lattice": lattice,
    }

    def obs_block(method: str) -> Dict[str, Any]:
        if method not in means.index:
            raise KeyError(f"{method} missing in scalardata for {run_folder}")
        le = float(means.loc[method, "LocalEnergy"]) * HARTREE_TO_EV
        ke = float(means.loc[method, "Kinetic"]) * HARTREE_TO_EV
        pe = float(means.loc[method, "LocalPotential"]) * HARTREE_TO_EV
        le_err = float(errs.loc[method, "LocalEnergy"]) * HARTREE_TO_EV if "LocalEnergy" in errs.columns else np.nan
        if method == "DMC" and has_fsc:
            le_c = le + dv + dt
            ke_c = ke + dt
            pe_c = pe + dv
            unc_le, unc_ke, unc_pe = le, ke, pe
            forces = dmc_forces
        else:
            le_c, ke_c, pe_c = le, ke, pe
            unc_le = unc_ke = unc_pe = np.nan
            forces = vmc_forces if method == "VMC" else dmc_forces
        return {
            "total_energy": le_c,
            "kinetic_energy": ke_c,
            "potential_energy": pe_c,
            "total_energy_error": le_err,
            "kinetic_energy_error": float(errs.loc[method, "Kinetic"]) * HARTREE_TO_EV
            if "Kinetic" in errs.columns
            else np.nan,
            "potential_energy_error": float(errs.loc[method, "LocalPotential"]) * HARTREE_TO_EV
            if "LocalPotential" in errs.columns
            else np.nan,
            "forces": forces,
            "fsc_dv_ev": dv if method == "DMC" else np.nan,
            "fsc_dt_ev": dt if method == "DMC" else np.nan,
            "total_energy_uncorrected": unc_le,
            "kinetic_energy_uncorrected": unc_ke,
            "potential_energy_uncorrected": unc_pe,
            "finite_size_corrected": bool(method == "DMC" and has_fsc),
        }

    xyz_dict = {
        "natoms": natoms,
        "header": header_base,
        "species": ["H"],
        "positions": positions,
        "forces": dmc_forces,
        "lattice_vectors": lattice,
        "fractional_positions": frac,
        "pbc": np.array([True, True, True]),
    }

    return {
        "run_folder": run_folder,
        "pressure": pressure,
        "temperature": temperature,
        "config_number": config_number,
        "rs": rs,
        "lattice": lattice,
        "positions": positions,
        "fractional_positions": frac,
        "xyz_dict": xyz_dict,
        "header_base": header_base,
        "vmc": obs_block("VMC"),
        "dmc": obs_block("DMC"),
        "qmc_quality_vmc": _compute_qmc_quality("VMC", xyz_dict, {**header_base, "energy": obs_block("VMC")["total_energy"]}, run_folder),
        "qmc_quality_dmc": _compute_qmc_quality(
            "DMC",
            xyz_dict,
            {**header_base, "energy": obs_block("DMC")["total_energy"]},
            run_folder,
        ),
        "extxyz": _build_extxyz(natoms, header_base, positions),
    }


def _write_observable(group: h5py.Group, obs: Dict[str, Any], natoms: int) -> None:
    for key in list(group.keys()):
        del group[key]
    group.attrs["finite_size_corrected"] = bool(obs["finite_size_corrected"])
    for scalar in (
        "total_energy",
        "kinetic_energy",
        "potential_energy",
        "total_energy_error",
        "kinetic_energy_error",
        "potential_energy_error",
    ):
        group.create_dataset(scalar, data=np.array([float(obs.get(scalar, np.nan))], dtype=float))
    for scalar, key in (
        ("fsc_potential_energy", "fsc_dv_ev"),
        ("fsc_kinetic_energy", "fsc_dt_ev"),
        ("total_energy_uncorrected", "total_energy_uncorrected"),
        ("kinetic_energy_uncorrected", "kinetic_energy_uncorrected"),
        ("potential_energy_uncorrected", "potential_energy_uncorrected"),
    ):
        group.create_dataset(scalar, data=np.array([float(obs.get(key, np.nan))], dtype=float))
    forces = obs["forces"]
    group.create_dataset("forces", data=np.expand_dims(forces, axis=0))
    group.create_dataset("forces_error", data=np.full((1, natoms, 3), np.nan))
    group.create_dataset("structure_factor", data=np.zeros((0, 3)))
    group.create_dataset("structure_factor_error", data=np.zeros((0, 3)))
    group.create_dataset("structure_factor_ks", data=np.zeros((0, 3)))


def patch_method_h5(h5_path: Path, method: str, payload: Dict[str, Any], dry_run: bool = False) -> None:
    obs = payload["vmc" if method == "VMC" else "dmc"]
    qq = payload["qmc_quality_vmc" if method == "VMC" else "qmc_quality_dmc"]
    natoms = int(payload["positions"].shape[0])
    if dry_run:
        return
    tmp = h5_path.with_suffix(".h5.reingest_tmp")
    with h5py.File(h5_path, "r") as src, h5py.File(tmp, "w") as dst:
        # copy groups wholesale except structure/observables which we rebuild
        for key, val in src.attrs.items():
            dst.attrs[key] = val
        for name in src.keys():
            if name in {"structure", "observables"}:
                continue
            src.copy(src[name], dst)
        dst.attrs["pressure"] = float(payload["pressure"])
        dst.attrs["temperature"] = float(payload["temperature"])
        dst.attrs["config_number"] = int(payload["config_number"])
        dst.attrs["rs"] = float(payload["rs"])
        dst.attrs["qmc_quality"] = int(qq)
        g_struct = dst.create_group("structure")
        g_struct.create_dataset("lattice_vectors", data=payload["lattice"])
        g_struct.create_dataset("positions", data=payload["positions"])
        g_struct.create_dataset("fractional_positions", data=payload["fractional_positions"])
        g_struct.create_dataset("pbc", data=np.array([True, True, True], dtype=bool))
        g_struct.create_dataset("xyz", data=payload["extxyz"], dtype=h5py.string_dtype("utf-8"))
        g_obs = dst.create_group("observables")
        _write_observable(g_obs, obs, natoms)
    tmp.replace(h5_path)


def reingest_config(run_folder: Path, dry_run: bool = False) -> Dict[str, Any]:
    run_folder = run_folder.resolve()
    if not CONFIG_RE.match(run_folder.name):
        raise ValueError(f"not a config folder: {run_folder}")
    h5s = _method_h5_files(run_folder)
    if not needs_reingest(h5s["DMC"]):
        return {"status": "skipped", "reason": "dmc_h5_already_populated", "folder": str(run_folder)}
    payload = gather_reingest_payload(run_folder)
    patch_method_h5(h5s["VMC"], "VMC", payload, dry_run=dry_run)
    patch_method_h5(h5s["DMC"], "DMC", payload, dry_run=dry_run)
    return {
        "status": "ok",
        "folder": str(run_folder),
        "vmc": str(h5s["VMC"]),
        "dmc": str(h5s["DMC"]),
        "rs": payload["rs"],
        "qmc_quality_dmc": payload["qmc_quality_dmc"],
        "P_GPa": payload["pressure"],
        "T_K": payload["temperature"],
    }


def iter_target_configs(aurora_root: Path, all_bad: bool = False) -> List[Path]:
    configs: List[Path] = []
    for ledger in sorted(aurora_root.glob(f"run_*/{LEDGER_NAME}")):
        for line in ledger.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[1] != "ok":
                continue
            cfg = Path(parts[0])
            if not cfg.is_dir():
                continue
            xyz = next((p for p in sorted(cfg.glob("*.xyz")) if ".struct." not in p.name), None)
            if xyz is not None and not all_bad:
                continue
            try:
                dmc = _method_h5_files(cfg)["DMC"]
            except FileNotFoundError:
                continue
            if all_bad or needs_reingest(dmc):
                configs.append(cfg)
    return sorted(set(configs))


def append_ledger(ledger: Path, row: Dict[str, Any]) -> None:
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with ledger.open("a", encoding="utf-8") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            fh.write(json.dumps(row, sort_keys=True) + "\n")
            fh.flush()
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def _worker(args: Tuple[str, bool]) -> Dict[str, Any]:
    folder, dry_run = args
    try:
        return reingest_config(Path(folder), dry_run=dry_run)
    except Exception as exc:
        return {"status": "error", "folder": folder, "error": repr(exc), "trace": traceback.format_exc()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aurora-root", type=Path, default=DEFAULT_AURORA)
    ap.add_argument("--single", type=Path, help="Re-ingest one configuration folder.")
    ap.add_argument("--dry-run", action="store_true", help="Parse sources and report only.")
    ap.add_argument("--all-bad", action="store_true", help="Include every h5-only ledger ok config.")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--ledger", type=Path, default=None, help="Append JSONL status rows here.")
    args = ap.parse_args()

    if args.single:
        configs = [args.single.resolve()]
    else:
        configs = iter_target_configs(args.aurora_root, all_bad=args.all_bad)
    print(f"configs to process: {len(configs)}")
    if args.dry_run and configs:
        sample = gather_reingest_payload(configs[0])
        print(
            f"sample {configs[0].name}: rs={sample['rs']:.4f} "
            f"DMC E={sample['dmc']['total_energy']:.4f} eV "
            f"qmc_quality={sample['qmc_quality_dmc']}"
        )

    ledger = args.ledger or (args.aurora_root / REINGEST_LEDGER)
    if args.workers <= 1:
        results = [_worker((str(c), args.dry_run)) for c in configs]
    else:
        results = []
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = [pool.submit(_worker, (str(c), args.dry_run)) for c in configs]
            for fut in as_completed(futs):
                results.append(fut.result())

    counts: Dict[str, int] = {}
    for row in results:
        counts[row.get("status", "unknown")] = counts.get(row.get("status", "unknown"), 0) + 1
        if not args.dry_run:
            append_ledger(ledger, row)
    print("summary:", counts)
    errors = [r for r in results if r.get("status") == "error"]
    if errors:
        print("first error:", errors[0].get("error"), errors[0].get("folder"))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
