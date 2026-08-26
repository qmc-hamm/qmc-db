#!/usr/bin/env python3
"""Shared helpers for ingesting the legacy QMC-HAMM ("Ricky") database.

Source: https://qmc-hamm.hub.yt/data.html, backed by the Girder instance at
https://girder.hub.yt (collection tree melth2/02_efv/{a_dft,b_qmc}).

The legacy data is ASE Ulm trajectories. Each frame is one 96-atom hydrogen
configuration carrying a DMC (or DFT) energy/forces evaluation. We convert one
frame into one schema-compliant HDF5 file, matching the layout already used for
the Aurora/SG hydrogen files.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np

BOHR_PER_ANG = 1.0 / 0.529177210903
RY_TO_EV = 13.605693122994
EV_PER_A3_TO_GPA = 160.21766208

SOURCE_DATABASE = "qmc-hamm.hub.yt"
SOURCE_DATASET = "ricky_legacy"
AUTHOR = "Hongwei Niu (Ricky); QMC-HAMM Team"
GIRDER = "https://girder.hub.yt"

# Human-readable provenance for the seven b_qmc subfolders.
QMC_FOLDER_DESC = {
    "f54-ipq": "i-PI quantum (PIMD) NPT, n2p2-blyp configurations",
    "f54-ipq-s2": "i-PI quantum (PIMD) NPT, set 2",
    "f54-ipq-s3": "i-PI quantum (PIMD) NPT, set 3",
    "f54-ipq-s4": "i-PI quantum (PIMD) NPT, set 4",
    "f54-ipq-s5": "i-PI quantum (PIMD) NPT, set 5",
    "f54-ipc": "i-PI classical NPT",
    "f54-hs1": "ASE classical MD, NVT",
}


def utf8() -> Any:
    return h5py.string_dtype(encoding="utf-8")


def write_text(group: h5py.Group, name: str, text: Optional[str]) -> None:
    if text is None:
        return
    group.create_dataset(name, data=str(text), dtype=utf8())


def write_str_list(group: h5py.Group, name: str, values: Iterable[str] | str | None) -> None:
    """Always write a 1-D UTF-8 string array (schema requires string[])."""
    if values is None:
        items: List[str] = []
    elif isinstance(values, str):
        items = [values]
    else:
        items = [str(v) for v in values]
    arr = np.asarray(items, dtype=object).reshape(-1)
    group.create_dataset(name, data=arr, dtype=utf8())


def rs_from_volume(volume_ang3: float, n_electrons: int) -> float:
    """Wigner-Seitz radius in bohr. Hydrogen: one electron per proton."""
    if not np.isfinite(volume_ang3) or volume_ang3 <= 0 or n_electrons <= 0:
        return float("nan")
    r_ang = (3.0 * volume_ang3 / (4.0 * np.pi * n_electrons)) ** (1.0 / 3.0)
    return float(r_ang * BOHR_PER_ANG)


def voigt_to_matrix(voigt: np.ndarray) -> np.ndarray:
    """ASE stores stress as Voigt (xx, yy, zz, yz, xz, xy)."""
    v = np.asarray(voigt, dtype=float).ravel()
    if v.size == 9:
        return v.reshape(3, 3)
    if v.size != 6:
        return np.full((3, 3), np.nan)
    xx, yy, zz, yz, xz, xy = v
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=float)


def stress_to_gpa(voigt: np.ndarray) -> np.ndarray:
    """ASE stress (eV/A^3) -> stress tensor in GPa."""
    return voigt_to_matrix(voigt) * EV_PER_A3_TO_GPA


def pressure_gpa_from_stress(voigt: np.ndarray) -> float:
    """Electronic (virial) pressure in GPa. ASE stress sign is -pressure."""
    m = stress_to_gpa(voigt)
    if not np.all(np.isfinite(np.diag(m))):
        return float("nan")
    return float(-np.trace(m) / 3.0)


def qmc_item_base(name: str) -> str:
    """Strip the QMC trajectory suffix to get the parent MD trajectory base name."""
    for suffix in (".dmc_mean.traj", ".traj.dmc.traj"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return re.sub(r"\.dmc\.traj$", "", name)


def dft_name_key(name: str) -> Tuple[str, str]:
    """Split a DFT trajectory filename into (base, functional)."""
    m = re.match(r"^(.*)\.(pbe|vdw-df)\.traj$", name)
    if m:
        return m.group(1), m.group(2)
    return (name[:-5] if name.endswith(".traj") else name), "none"


def parse_pt_from_name(name: str) -> Tuple[Optional[int], Optional[int]]:
    """Pressure (GPa) and temperature (K) encoded in the trajectory filename."""
    m = re.search(r"-p0*(\d+)-t0*(\d+)", name)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.search(r"-t0*(\d+)", name)
    if m:
        return None, int(m.group(1))
    return None, None


def ensemble_from_name(name: str) -> str:
    if name.startswith("npt-"):
        return "npt"
    if name.startswith(("pbe-n096", "vdw-n096")):
        return "nvt"
    return "unknown"


def pt_dirname(pressure: Optional[float], temperature: Optional[float]) -> str:
    """Directory bucket name, e.g. P140T2400. Pressure is rounded to integer GPa."""
    p = "Pnan" if pressure is None or not np.isfinite(pressure) else f"P{int(round(pressure))}"
    t = "Tnan" if temperature is None or not np.isfinite(temperature) else f"T{int(round(temperature))}"
    return f"{p}{t}"


def extxyz_text(atoms, extra_header: Optional[Dict[str, Any]] = None) -> str:
    """Render one frame as extended XYZ, mirroring structure/xyz in the SG files."""
    lattice = np.asarray(atoms.cell, dtype=float).reshape(3, 3)
    positions = np.asarray(atoms.positions, dtype=float)
    symbols = atoms.get_chemical_symbols()
    try:
        forces = np.asarray(atoms.get_forces(), dtype=float)
    except Exception:  # noqa: BLE001
        forces = None

    fields = ["Lattice=\"" + " ".join(f"{v:.10f}" for v in lattice.ravel()) + "\""]
    props = "species:S:1:pos:R:3" + (":forces:R:3" if forces is not None else "")
    fields.append(f"Properties={props}")
    for key, val in (extra_header or {}).items():
        if val is None:
            continue
        if isinstance(val, (int, float, np.floating, np.integer)):
            fields.append(f"{key}={val}")
        else:
            fields.append(f'{key}="{val}"')
    lines = [str(len(atoms)), " ".join(fields)]
    for i, sym in enumerate(symbols):
        row = f"{sym} " + " ".join(f"{v:.10f}" for v in positions[i])
        if forces is not None:
            row += " " + " ".join(f"{v:.10f}" for v in forces[i])
        lines.append(row)
    return "\n".join(lines) + "\n"


def write_structure(h5: h5py.File, atoms, xyz_header: Optional[Dict[str, Any]] = None) -> None:
    lattice = np.asarray(atoms.cell, dtype=float).reshape(3, 3)
    positions = np.asarray(atoms.positions, dtype=float)
    g = h5.create_group("structure")
    g.create_dataset("lattice_vectors", data=lattice)
    g.create_dataset("positions", data=positions)
    g.create_dataset("fractional_positions", data=positions @ np.linalg.inv(lattice))
    g.create_dataset("pbc", data=np.asarray(atoms.pbc, dtype=bool))
    g.create_dataset("xyz", data=extxyz_text(atoms, xyz_header), dtype=utf8())


def base_root_attrs(
    h5: h5py.File,
    *,
    uuid_hex: str,
    atoms,
    calculation_type: str,
    method: str,
    method_kws: Iterable[str],
    pressure: Optional[float],
    temperature: Optional[float],
    creation_date: str,
    name_in_system: str,
    config_number: Optional[int],
    ensemble: str,
    quantum_ions: Optional[bool],
    starting_configuration_model_name: Optional[str],
    phase: Optional[str] = None,
) -> None:
    natoms = len(atoms)
    h5.attrs["uuid"] = uuid_hex
    h5.attrs["author"] = AUTHOR
    h5.attrs["system"] = "Hydrogen"
    h5.attrs["formula"] = "H"
    h5.attrs["natoms"] = int(natoms)
    h5.attrs["species"] = np.array(sorted(set(atoms.get_chemical_symbols())), dtype="S")
    h5.attrs["calculation_type"] = calculation_type
    h5.attrs["method"] = method
    h5.attrs["method_kws"] = np.array(list(method_kws), dtype="S")
    h5.attrs["rs"] = rs_from_volume(float(atoms.get_volume()), natoms)
    if temperature is not None and np.isfinite(temperature):
        h5.attrs["temperature"] = float(temperature)
    if pressure is not None and np.isfinite(pressure):
        h5.attrs["pressure"] = float(pressure)
    h5.attrs["creation_date"] = str(creation_date)
    h5.attrs["name_in_system"] = name_in_system
    if config_number is not None:
        h5.attrs["config_number"] = int(config_number)
    h5.attrs["ensemble"] = ensemble
    if phase:
        h5.attrs["phase"] = str(phase)
    if quantum_ions is not None:
        h5.attrs["quantum_ions"] = bool(quantum_ions)
    if starting_configuration_model_name:
        h5.attrs["starting_configuration_model_name"] = str(starting_configuration_model_name)
    # Legacy-database identifiers requested for provenance tracking.
    h5.attrs["source_database"] = SOURCE_DATABASE
    h5.attrs["source_dataset"] = SOURCE_DATASET


def write_provenance(
    h5: h5py.File,
    *,
    dependent_uuids: Iterable[str],
    uuid_in_system: Optional[str],
    source_files: Iterable[str],
    girder_item_id: str,
    girder_folder: str,
    girder_config_id: Optional[str],
    girder_conf_uuid: Optional[str],
    source_frame_index: int,
    source_iconf: Optional[int],
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    g = h5.create_group("provenance")
    write_str_list(g, "dependent_uuids", list(dependent_uuids))
    write_text(g, "uuid_in_system", uuid_in_system if uuid_in_system else girder_item_id)
    write_str_list(g, "source_files", list(source_files))
    write_text(g, "source_database", SOURCE_DATABASE)
    write_text(g, "source_dataset", SOURCE_DATASET)
    write_text(g, "girder_item_id", girder_item_id)
    write_text(g, "girder_folder", girder_folder)
    write_text(g, "girder_item_url", f"{GIRDER}/api/v1/item/{girder_item_id}/download")
    if girder_config_id:
        write_text(g, "girder_config_id", girder_config_id)
    if girder_conf_uuid:
        write_text(g, "girder_conf_uuid", girder_conf_uuid)
    g.create_dataset("source_frame_index", data=np.array([int(source_frame_index)], dtype=np.int64))
    if source_iconf is not None:
        g.create_dataset("source_iconf", data=np.array([int(source_iconf)], dtype=np.int64))
    if extra:
        write_text(g, "source_metadata", json.dumps(extra, sort_keys=True))


def atomic_write_h5(tmp_path: Path, final_path: Path) -> None:
    final_path.parent.mkdir(parents=True, exist_ok=True)
    os.replace(tmp_path, final_path)


def tmp_name(directory: Path, tag: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    return directory / f".{tag}_tmp_{os.getpid()}_{os.urandom(6).hex()}.h5"
