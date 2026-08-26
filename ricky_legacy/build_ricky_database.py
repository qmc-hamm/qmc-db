#!/usr/bin/env python3
"""Convert the legacy QMC-HAMM ("Ricky") database into schema-compliant HDF5 files.

Two phases, because DMC files must reference their DFT parents:

  phase "dft"  reads a_dft ASE trajectories and writes one method=DFT file per
               (frame, functional) for every frame that backs a DMC config,
               plus the standalone BOPIMC hcp trajectory.
  phase "qmc"  reads b_qmc ASE trajectories and writes one method=DMC file per
               frame, with provenance/dependent_uuids pointing at the DFT files
               written in phase "dft".

The link between the two sets is exact: a QMC frame with info['iconf'] == N is
frame index N of the parent MD trajectory named by the DFT item's meta['ftraj']
(verified on a random sample: cells and wrapped fractional coordinates agree).

Every written file is recorded in a ledger keyed by a stable source key
("ricky:<girder_item_id>:<frame>[:<functional>]") so re-runs never duplicate work.
"""

from __future__ import annotations

import argparse
import collections
import fcntl
import json
import hashlib
import re
import sys
import time
import uuid as uuidlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ricky_common import (  # noqa: E402
    AUTHOR,
    QMC_FOLDER_DESC,
    RY_TO_EV,
    base_root_attrs,
    rs_from_volume,
    dft_name_key,
    ensemble_from_name,
    parse_pt_from_name,
    pressure_gpa_from_stress,
    pt_dirname,
    qmc_item_base,
    stress_to_gpa,
    tmp_name,
    utf8,
    write_provenance,
    write_str_list,
    write_structure,
    write_text,
)

DEFAULT_META = Path("/scratch/sgoswam3/ricky_qmchamm/meta")
DEFAULT_RAW = Path("/scratch/sgoswam3/ricky_qmchamm")
DEFAULT_OUT = Path(
    "/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/ricky_legacy_database"
)


# --------------------------------------------------------------------------- #
# work plan
# --------------------------------------------------------------------------- #
def load_catalogs(meta: Path) -> Dict[str, Any]:
    qmc_items = json.loads((meta / "folder_items.json").read_text())
    dft_items = json.loads((meta / "dft_items.json").read_text())
    qmc_paths = json.loads((meta / "download_manifest.json").read_text())
    dft_paths = json.loads((meta / "download_manifest_dft.json").read_text())
    return {
        "qmc_items": qmc_items,
        "dft_items": dft_items,
        "qmc_paths": {k: v["path"] for k, v in qmc_paths.items()},
        "dft_paths": {k: v["path"] for k, v in dft_paths.items()},
    }


def _geometry_matches(a, b) -> bool:
    """Same configuration? The QMC trajectories store rounded coordinates, so
    compare cells loosely and positions with a minimum-image convention."""
    ca = np.asarray(a.cell, dtype=float)
    cb = np.asarray(b.cell, dtype=float)
    if not np.allclose(ca, cb, rtol=0, atol=1e-3):
        return False
    fa = np.asarray(a.get_scaled_positions(wrap=False), dtype=float)
    fb = np.asarray(b.get_scaled_positions(wrap=False), dtype=float)
    if fa.shape != fb.shape:
        return False
    d = fa - fb
    d -= np.round(d)
    return bool(np.abs(d).max() < 1e-3)


def build_plan(cat: Dict[str, Any], raw_root: Path) -> Dict[str, Any]:
    """Decide which DFT frames are needed and how QMC frames map onto them.

    Two trajectory names collide across the b_qmc folders (npt-p175-t1800-b0 and
    npt-p200-t1800-b0 exist in both the classical f54-ipc and the quantum
    f54-ipq* sets), so the base name alone is not a unique key. meta.conf.uuid
    identifies the parent MD run and disambiguates them. Every surviving link is
    then confirmed geometrically, so a naming surprise can never produce a wrong
    provenance edge -- it can only produce a missing one, which we report.
    """
    dft_items = cat["dft_items"]

    # base name -> list of (conf_uuid, functional, itemId)
    dft_index: Dict[str, List[Tuple[Optional[str], str, str]]] = collections.defaultdict(list)
    for iid, it in dft_items.items():
        base, func = dft_name_key(it["name"])
        conf_uuid = (it.get("meta", {}) or {}).get("conf", {}) or {}
        dft_index[base].append((conf_uuid.get("uuid"), func, iid))

    traj_cache: Dict[str, Any] = {}

    def dft_frames(iid: str):
        if iid not in traj_cache:
            traj_cache[iid] = read(str(raw_root / cat["dft_paths"][iid]), index=":")
        return traj_cache[iid]

    needed_dft: Dict[str, set] = collections.defaultdict(set)  # itemId -> {frame idx}
    qmc_links: Dict[str, Dict[int, List[Tuple[str, int, str]]]] = {}
    unlinked: List[str] = []
    stats = collections.Counter()

    for qid, it in cat["qmc_items"].items():
        base = qmc_item_base(it["name"])
        cand = dft_index.get(base) or dft_index.get(re.sub(r"\.(pbe|vdw-df)$", "", base))
        if not cand:
            unlinked.append(it["name"])
            continue
        # Prefer DFT items from the same parent MD run.
        q_conf_uuid = ((it.get("meta", {}) or {}).get("conf", {}) or {}).get("uuid")
        same_run = [c for c in cand if q_conf_uuid and c[0] == q_conf_uuid]
        if same_run:
            cand = same_run
            stats["matched_by_conf_uuid"] += 1
        else:
            stats["matched_by_name_only"] += 1

        frames = read(str(raw_root / cat["qmc_paths"][qid]), index=":")
        per_frame: Dict[int, List[Tuple[str, int, str]]] = {}
        for fi, atoms in enumerate(frames):
            iconf = atoms.info.get("iconf")
            links: List[Tuple[str, int, str]] = []
            if iconf is not None:
                ic = int(iconf)
                for _cu, func, diid in sorted(cand, key=lambda x: (x[1], x[2])):
                    dfr = dft_frames(diid)
                    if ic >= len(dfr):
                        stats["link_out_of_range"] += 1
                        continue
                    if not _geometry_matches(atoms, dfr[ic]):
                        stats["link_rejected_geometry"] += 1
                        continue
                    needed_dft[diid].add(ic)
                    links.append((diid, ic, func))
                    stats["link_accepted"] += 1
            if not links:
                stats["frames_without_any_link"] += 1
            per_frame[fi] = links
        qmc_links[qid] = per_frame

    # The BOPIMC hcp trajectory has no QMC counterpart but is a distinct
    # published DFT dataset, so ingest all of its frames.
    for iid, it in dft_items.items():
        if it["_folderName"] == "a_dft_root":
            needed_dft[iid] = set(range(len(dft_frames(iid))))

    return {
        "needed_dft": {k: sorted(v) for k, v in needed_dft.items()},
        "qmc_links": qmc_links,
        "unlinked": unlinked,
        "link_stats": dict(stats),
    }


# --------------------------------------------------------------------------- #
# ledger
# --------------------------------------------------------------------------- #
def append_ledger(ledger: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger, "a", encoding="utf-8") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            for r in rows:
                fh.write(json.dumps(r, sort_keys=True) + "\n")
            fh.flush()
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def read_ledger(ledger: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    if not ledger.exists():
        return out
    for line in ledger.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("status") == "ok" and r.get("source_key"):
            out[r["source_key"]] = r
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# DFT writer
# --------------------------------------------------------------------------- #
def write_dft_frames(job: Dict[str, Any]) -> List[Dict[str, Any]]:
    item = job["item"]
    traj_path = Path(job["traj_path"])
    out_root = Path(job["out_root"])
    frames_wanted = job["frames"]
    done = job["done"]

    meta = item.get("meta", {}) or {}
    dp = meta.get("dft_parameters", {}) or {}
    conf = meta.get("conf", {}) or {}
    name = item["name"]
    base, func = dft_name_key(name)
    if func == "none":
        func = dp.get("input_dft") or "unknown"
    p_nom, t_nom = parse_pt_from_name(name)
    ensemble = ensemble_from_name(name)
    if name.startswith("f15b_bopimc"):
        ensemble = "nvt"

    rows: List[Dict[str, Any]] = []
    frames = None
    for fi in frames_wanted:
        source_key = f"ricky:{item['_id']}:{fi}:{func}"
        if source_key in done:
            rows.append({**done[source_key], "status": "ok", "skipped": True})
            continue
        if frames is None:
            frames = read(str(traj_path), index=":")
        try:
            atoms = frames[fi]
            results = atoms.calc.results if atoms.calc is not None else {}
            energy = float(atoms.get_potential_energy())
            forces = np.asarray(atoms.get_forces(), dtype=float)
            voigt = np.asarray(results.get("stress", np.full(6, np.nan)), dtype=float)
            p_virial = pressure_gpa_from_stress(voigt)
            # The BOPIMC trajectory stores its stress with the opposite sign to
            # the ASE convention used by every other trajectory here. Dense
            # hydrogen (rs < 2) cannot sit at large negative pressure, and the
            # magnitude matches this dataset's own P(rs) curve, so flip it and
            # record that we did.
            rs_frame = rs_from_volume(float(atoms.get_volume()), len(atoms))
            sign_flipped = False
            if np.isfinite(p_virial) and p_virial < 0.0 and rs_frame < 2.0:
                voigt = -voigt
                p_virial = -p_virial
                sign_flipped = True
            pressure = float(p_nom) if p_nom is not None else p_virial

            # Legacy DMC geometries are solid-like (Bragg peaks in S(k)); the
            # sole named crystal is BOPIMC hcp. All configs labeled hcp.
            phase = "hcp"

            uuid_hex = uuidlib.uuid4().hex
            out_dir = out_root / pt_dirname(pressure, t_nom)
            tmp = tmp_name(out_dir, "dft")
            with h5py.File(tmp, "w") as h5:
                base_root_attrs(
                    h5,
                    uuid_hex=uuid_hex,
                    atoms=atoms,
                    calculation_type="scf",
                    method="DFT",
                    method_kws=[func, "plane_wave"],
                    pressure=pressure,
                    temperature=t_nom,
                    creation_date=str(item.get("created", ""))[:10],
                    name_in_system=f"ricky_legacy_{base}_{func}_i{fi}",
                    config_number=fi,
                    ensemble=ensemble,
                    quantum_ions=None,
                    starting_configuration_model_name=conf.get("config_dft"),
                    phase=phase,
                )

                gp = h5.create_group("parameters")
                gp.attrs["spin_polarized"] = False
                kgrid = dp.get("kgrid") or [0, 0, 0]
                gp.create_dataset("kpoint_grid", data=np.asarray(kgrid, dtype=np.int32))
                if dp.get("kshift") is not None:
                    gp.create_dataset("kpoint_shift", data=np.asarray(dp["kshift"], dtype=np.int32))
                for key, ry in (("ecutwfc", dp.get("ecutwfc")), ("ecutrho", dp.get("ecutrho"))):
                    if ry is not None:
                        gp.create_dataset(key, data=np.array([float(ry) * RY_TO_EV], dtype=float))
                if dp.get("degauss") is not None:
                    gp.create_dataset("sigma", data=np.array([float(dp["degauss"]) * RY_TO_EV], dtype=float))
                if dp.get("smearing"):
                    write_text(gp, "smearing", str(dp["smearing"]))
                write_text(gp, "other_input", json.dumps(dp, sort_keys=True))

                gc = h5.create_group("code")
                gc.attrs["name"] = "Quantum ESPRESSO"
                gc.attrs["version"] = str(dp.get("version", dp.get("code", "unknown")))
                gc.attrs["code_tag"] = str(dp.get("code", "unknown"))
                if dp.get("ase"):
                    gc.attrs["ase_version"] = str(dp["ase"])
                gc.attrs["parallelization"] = "unknown"
                gc.attrs["config_creation_date"] = str(item.get("created", ""))[:10]
                gc.attrs["source_girder_folder"] = item["_folderName"]
                write_text(gc, "system", "unknown (legacy QMC-HAMM production run)")
                write_text(gc, "location", f"girder.hub.yt item {item['_id']}")
                write_text(gc, "dft_parameters", json.dumps(dp, sort_keys=True))
                pseudos = (dp.get("pseudos") or {})
                if pseudos:
                    write_str_list(gc, "pseudo_files", [f"{k}: {v}" for k, v in sorted(pseudos.items())])

                info_header = {k: v for k, v in atoms.info.items()}
                info_header.update({"energy": energy, "functional": func,
                                    "source_dataset": "ricky_legacy"})
                write_structure(h5, atoms, info_header)

                go = h5.create_group("observables")
                go.create_dataset("total_energy", data=np.array([energy], dtype=float))
                go.create_dataset("forces", data=np.expand_dims(forces, axis=0))
                go.create_dataset("stress", data=stress_to_gpa(voigt))
                dsp = go.create_dataset("pressure", data=np.array([p_virial], dtype=float))
                # Root-attr `pressure` is the MD target (thermodynamic condition,
                # from the n2p2-blyp/vdw MD); this is what THIS functional's
                # virial gives for the same geometry, so the two differ.
                dsp.attrs["description"] = (
                    "electronic virial pressure (GPa) from the stress tensor of this "
                    "DFT evaluation; excludes the ionic kinetic contribution"
                )
                # f30c frames also carry MD-instantaneous thermodynamic state.
                info_num = {}
                for k, v in atoms.info.items():
                    try:
                        info_num[k] = float(v)
                    except (TypeError, ValueError):
                        pass
                for src, dst in (("Etot[eV]", "md_total_energy"), ("Epot[eV]", "md_potential_energy"),
                                 ("Ekin[eV]", "md_kinetic_energy"), ("T[K]", "md_temperature")):
                    if src in info_num:
                        go.create_dataset(dst, data=np.array([info_num[src]], dtype=float))
                pmd = [info_num.get(f"P{a}[GPa]") for a in ("xx", "yy", "zz")]
                if all(v is not None for v in pmd):
                    go.create_dataset("md_pressure", data=np.array([float(np.mean(pmd))], dtype=float))

                write_provenance(
                    h5,
                    dependent_uuids=[],
                    uuid_in_system=atoms.info.get("uuid") or conf.get("uuid"),
                    source_files=[name],
                    girder_item_id=item["_id"],
                    girder_folder=item["_folderName"],
                    girder_config_id=conf.get("configId") or meta.get("configId"),
                    girder_conf_uuid=conf.get("uuid"),
                    source_frame_index=fi,
                    source_iconf=int(atoms.info["iconf"]) if str(atoms.info.get("iconf", "")).isdigit() else None,
                    extra={"functional": func, "ftraj": meta.get("ftraj"), "path": meta.get("path"),
                           "stress_sign_convention_flipped": sign_flipped},
                )

            final = out_dir / f"{uuid_hex}.h5"
            tmp.replace(final)
            rows.append({
                "source_key": source_key,
                "status": "ok",
                "uuid": uuid_hex,
                "method": "DFT",
                "functional": func,
                "path": str(final),
                "sha256": sha256_file(final),
                "size": final.stat().st_size,
                "P": pressure,
                "T": t_nom,
                "girder_item_id": item["_id"],
                "girder_item_name": name,
                "frame": fi,
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({
                "source_key": source_key,
                "status": "error",
                "error": repr(exc),
                "girder_item_id": item["_id"],
                "girder_item_name": name,
                "frame": fi,
            })
    return rows


# --------------------------------------------------------------------------- #
# QMC writer
# --------------------------------------------------------------------------- #
def write_qmc_frames(job: Dict[str, Any]) -> List[Dict[str, Any]]:
    item = job["item"]
    traj_path = Path(job["traj_path"])
    out_root = Path(job["out_root"])
    links = {int(k): v for k, v in job["links"].items()}
    dft_uuid = job["dft_uuid"]
    done = job["done"]

    meta = item.get("meta", {}) or {}
    conf = meta.get("conf", {}) or {}
    qmc_meta = meta.get("qmc", []) or []
    name = item["name"]
    base = qmc_item_base(name)
    p_nom = conf.get("pgpa")
    t_nom = conf.get("tkelvin")
    if p_nom is None or t_nom is None:
        p_nom, t_nom = parse_pt_from_name(name)
    ensemble = conf.get("ens") or ensemble_from_name(name)

    rows: List[Dict[str, Any]] = []
    frames = None
    n_avail = None
    for fi in range(job["nframes"]):
        source_key = f"ricky:{item['_id']}:{fi}"
        if source_key in done:
            rows.append({**done[source_key], "status": "ok", "skipped": True})
            continue
        if frames is None:
            frames = read(str(traj_path), index=":")
            n_avail = len(frames)
        try:
            atoms = frames[fi]
            energy = float(atoms.get_potential_energy())
            forces = np.asarray(atoms.get_forces(), dtype=float)
            dv = atoms.info.get("fsc_dv_ev")
            dt_ = atoms.info.get("fsc_dt_ev")
            dv = float(dv) if dv is not None else np.nan
            dt_ = float(dt_) if dt_ is not None else np.nan
            # Convention (matches the Aurora/SG pipeline): stored energy is the
            # finite-size-corrected value; the corrections were added to it.
            uncorrected = energy - (dv + dt_) if np.isfinite(dv) and np.isfinite(dt_) else np.nan
            fsc_ok = bool(np.isfinite(dv) and np.isfinite(dt_))

            per_conf = qmc_meta[fi] if fi < len(qmc_meta) else {}
            det = per_conf.get("determinant", {}) or {}
            simcell = per_conf.get("simulationcell", {}) or {}
            jastrow = per_conf.get("jastrow", []) or []
            blocks = per_conf.get("qmc", []) or []
            dmc_block = next((b for b in blocks if b.get("method") == "dmc"), {})
            vmc_block = next((b for b in blocks if b.get("method") == "vmc"), {})

            deps = [dft_uuid[f"{d}:{ic}:{fn}"] for d, ic, fn in links.get(fi, [])
                    if f"{d}:{ic}:{fn}" in dft_uuid]

            uuid_hex = uuidlib.uuid4().hex
            out_dir = out_root / pt_dirname(p_nom, t_nom)
            tmp = tmp_name(out_dir, "dmc")
            with h5py.File(tmp, "w") as h5:
                base_root_attrs(
                    h5,
                    uuid_hex=uuid_hex,
                    atoms=atoms,
                    calculation_type="qmc",
                    method="DMC",
                    method_kws=[str(det.get("input_dft", "pbe")).strip("'"), "slater_jastrow", "twist_averaged"],
                    pressure=p_nom,
                    temperature=t_nom,
                    creation_date=str(item.get("created", ""))[:10],
                    name_in_system=f"ricky_legacy_{base}_i{atoms.info.get('iconf', fi)}",
                    config_number=int(atoms.info["iconf"]) if str(atoms.info.get("iconf", "")).isdigit() else fi,
                    ensemble=ensemble,
                    quantum_ions=conf.get("quantum"),
                    starting_configuration_model_name=conf.get("config_dft"),
                    phase="hcp",
                )
                h5.attrs["qmc_quality"] = 10

                gp = h5.create_group("parameters")
                gp.attrs["spin_polarized"] = False
                gp.create_dataset("kpoint_grid", data=np.asarray(det.get("kgrid", [0, 0, 0]), dtype=np.int32))
                if det.get("kshift") is not None:
                    gp.create_dataset("kpoint_shift", data=np.asarray(det["kshift"], dtype=np.int32))
                if dmc_block.get("timestep") is not None:
                    gp.create_dataset("time_step", data=np.array([float(dmc_block["timestep"])], dtype=float))
                if dmc_block.get("steps") is not None:
                    gp.create_dataset("n_steps", data=np.array([float(dmc_block["steps"])], dtype=float))
                if dmc_block.get("blocks") is not None:
                    gp.create_dataset("n_blocks", data=np.array([float(dmc_block["blocks"])], dtype=float))
                if dmc_block.get("targetwalkers") is not None:
                    gp.create_dataset("target_walkers", data=np.array([float(dmc_block["targetwalkers"])], dtype=float))
                for key in ("ecutwfc", "ecutrho"):
                    if det.get(key) is not None:
                        gp.create_dataset(key, data=np.array([float(det[key]) * RY_TO_EV], dtype=float))
                # Compact summary only; the full blobs live under /code to avoid
                # storing the same JSON twice in every one of ~17.5k files.
                write_text(gp, "other_input", json.dumps({
                    "orbital_functional": str(det.get("input_dft", "")).strip("'"),
                    "orbital_code": det.get("code"),
                    "nbnd": det.get("nbnd"),
                    "vmc": vmc_block,
                    "dmc": dmc_block,
                    "lr_handler": simcell.get("LR_handler"),
                    "lr_dim_cutoff": simcell.get("LR_dim_cutoff"),
                    "n_jastrow_terms": len(jastrow),
                }, sort_keys=True))

                gc = h5.create_group("code")
                gc.attrs["name"] = "QMCPACK"
                gc.attrs["version"] = "unknown (legacy QMC-HAMM production run)"
                gc.attrs["orbital_code"] = str(det.get("code", "unknown"))
                gc.attrs["parallelization"] = "MPI+OpenMP"
                gc.attrs["config_creation_date"] = str(item.get("created", ""))[:10]
                gc.attrs["source_girder_folder"] = item["_folderName"]
                gc.attrs["source_girder_folder_desc"] = QMC_FOLDER_DESC.get(item["_folderName"], "")
                gc.attrs["mixed_estimator_available"] = False
                write_text(gc, "system", "unknown (legacy QMC-HAMM production run)")
                write_text(gc, "location", f"girder.hub.yt item {item['_id']}")
                write_text(gc, "determinant", json.dumps(det, sort_keys=True))
                write_text(gc, "jastrow", json.dumps(jastrow, sort_keys=True))
                write_text(gc, "simulationcell", json.dumps(simcell, sort_keys=True))
                write_text(gc, "qmc_blocks", json.dumps(blocks, sort_keys=True))
                if det.get("pseudos") or det.get("pseudo_dir"):
                    write_str_list(gc, "pseudo_files", [str(det.get("pseudos") or det.get("pseudo_dir"))])

                write_structure(h5, atoms, {
                    "energy": energy,
                    "pressure": p_nom,
                    "temperature": t_nom,
                    "fsc_dv_ev": dv,
                    "fsc_dt_ev": dt_,
                    "uuid": atoms.info.get("uuid"),
                    "iconf": atoms.info.get("iconf"),
                    "source_dataset": "ricky_legacy",
                })

                go = h5.create_group("observables")
                go.attrs["finite_size_corrected"] = fsc_ok
                go.create_dataset("total_energy", data=np.array([energy], dtype=float))
                go.create_dataset("total_energy_error", data=np.array([np.nan], dtype=float))
                go.create_dataset("forces", data=np.expand_dims(forces, axis=0))
                go.create_dataset("forces_error", data=np.full((1, len(atoms), 3), np.nan))
                go.create_dataset("fsc_potential_energy", data=np.array([dv], dtype=float))
                go.create_dataset("fsc_kinetic_energy", data=np.array([dt_], dtype=float))
                go.create_dataset("total_energy_uncorrected", data=np.array([uncorrected], dtype=float))
                # The legacy dmc_mean trajectories publish only the total energy.
                # These are written as NaN so the whole database exposes one
                # uniform set of DMC observable names.
                nan1 = np.array([np.nan], dtype=float)
                for key in ("kinetic_energy", "kinetic_energy_error",
                            "potential_energy", "potential_energy_error",
                            "kinetic_energy_uncorrected", "potential_energy_uncorrected"):
                    go.create_dataset(key, data=nan1.copy())
                for key in ("structure_factor", "structure_factor_error", "structure_factor_ks"):
                    go.create_dataset(key, data=np.zeros((0, 3), dtype=float))

                write_provenance(
                    h5,
                    dependent_uuids=deps,
                    uuid_in_system=atoms.info.get("uuid"),
                    source_files=[name],
                    girder_item_id=item["_id"],
                    girder_folder=item["_folderName"],
                    girder_config_id=conf.get("configId"),
                    girder_conf_uuid=conf.get("uuid"),
                    source_frame_index=fi,
                    source_iconf=int(atoms.info["iconf"]) if str(atoms.info.get("iconf", "")).isdigit() else None,
                    extra={
                        "qmc_prefix": atoms.info.get("prefix"),
                        "input_dft": conf.get("input_dft"),
                        "config_dft": conf.get("config_dft"),
                        "vmc_block": vmc_block,
                    },
                )

            final = out_dir / f"{uuid_hex}.h5"
            tmp.replace(final)
            rows.append({
                "source_key": source_key,
                "status": "ok",
                "uuid": uuid_hex,
                "method": "DMC",
                "path": str(final),
                "sha256": sha256_file(final),
                "size": final.stat().st_size,
                "P": p_nom,
                "T": t_nom,
                "girder_item_id": item["_id"],
                "girder_item_name": name,
                "frame": fi,
                "dependent_uuids": deps,
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({
                "source_key": source_key,
                "status": "error",
                "error": repr(exc),
                "girder_item_id": item["_id"],
                "girder_item_name": name,
                "frame": fi,
            })
    return rows


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", type=Path, default=DEFAULT_META)
    ap.add_argument("--raw-root", type=Path, default=DEFAULT_RAW)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--phase", choices=["plan", "dft", "qmc", "all"], default="all")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit-items", type=int, default=0, help="debug: cap number of source items per phase")
    args = ap.parse_args()

    cat = load_catalogs(args.meta)
    out_root = args.out_root
    ledger = out_root / "ricky_ingest_ledger.jsonl"
    plan_cache = out_root / "ricky_work_plan.json"

    if plan_cache.exists() and args.phase != "plan":
        plan = json.loads(plan_cache.read_text())
        print(f"[plan] loaded cached plan from {plan_cache}")
    else:
        t0 = time.time()
        print("[plan] building work plan (reading trajectories to resolve iconf links)...")
        plan = build_plan(cat, args.raw_root)
        out_root.mkdir(parents=True, exist_ok=True)
        plan_cache.write_text(json.dumps(plan))
        print(f"[plan] built in {time.time()-t0:.0f}s -> {plan_cache}")

    n_dft = sum(len(v) for v in plan["needed_dft"].values())
    n_qmc = sum(len(v) for v in plan["qmc_links"].values())
    print(f"[plan] DFT frames to write : {n_dft} across {len(plan['needed_dft'])} items")
    print(f"[plan] DMC frames to write : {n_qmc} across {len(plan['qmc_links'])} items")
    print(f"[plan] unlinked QMC items  : {len(plan['unlinked'])}")
    print(f"[plan] link stats          : {plan.get('link_stats', {})}")
    if args.phase == "plan":
        return 0

    done = read_ledger(ledger)
    print(f"[ledger] {len(done)} entries already recorded ok")

    # ---------------- phase dft ----------------
    if args.phase in ("dft", "all"):
        jobs = []
        for iid, frames in plan["needed_dft"].items():
            jobs.append({
                "item": cat["dft_items"][iid],
                "traj_path": str(args.raw_root / cat["dft_paths"][iid]),
                "out_root": str(out_root),
                "frames": frames,
                "done": {k: v for k, v in done.items() if k.startswith(f"ricky:{iid}:")},
            })
        if args.limit_items:
            jobs = jobs[: args.limit_items]
        print(f"[dft] {len(jobs)} items -> {sum(len(j['frames']) for j in jobs)} frames")
        t0 = time.time()
        nok = nerr = 0
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(write_dft_frames, j) for j in jobs]
            for n, fut in enumerate(as_completed(futs), 1):
                rows = fut.result()
                append_ledger(ledger, [r for r in rows if not r.get("skipped")])
                nok += sum(1 for r in rows if r["status"] == "ok")
                nerr += sum(1 for r in rows if r["status"] != "ok")
                if n % 50 == 0:
                    print(f"  [dft] {n}/{len(jobs)} items  ok={nok} err={nerr}  {time.time()-t0:.0f}s", flush=True)
        print(f"[dft] done ok={nok} err={nerr} in {time.time()-t0:.0f}s")
        done = read_ledger(ledger)

    # ---------------- phase qmc ----------------
    if args.phase in ("qmc", "all"):
        dft_uuid = {}
        for key, row in done.items():
            if row.get("method") != "DFT":
                continue
            parts = key.split(":")
            if len(parts) == 4:
                dft_uuid[f"{parts[1]}:{parts[2]}:{parts[3]}"] = row["uuid"]
        print(f"[qmc] {len(dft_uuid)} DFT uuids available for provenance links")

        jobs = []
        for qid, per_frame in plan["qmc_links"].items():
            # Only ship the uuids this item actually references; sending the whole
            # 36k-entry map to every one of 1594 workers is gigabytes of pickling.
            wanted = {f"{d}:{ic}:{fn}" for links in per_frame.values() for d, ic, fn in links}
            jobs.append({
                "item": cat["qmc_items"][qid],
                "traj_path": str(args.raw_root / cat["qmc_paths"][qid]),
                "out_root": str(out_root),
                "links": per_frame,
                "nframes": len(per_frame),
                "dft_uuid": {k: dft_uuid[k] for k in wanted if k in dft_uuid},
                "done": {k: v for k, v in done.items() if k.startswith(f"ricky:{qid}:")},
            })
        if args.limit_items:
            jobs = jobs[: args.limit_items]
        print(f"[qmc] {len(jobs)} items -> {sum(j['nframes'] for j in jobs)} frames")
        t0 = time.time()
        nok = nerr = 0
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(write_qmc_frames, j) for j in jobs]
            for n, fut in enumerate(as_completed(futs), 1):
                rows = fut.result()
                append_ledger(ledger, [r for r in rows if not r.get("skipped")])
                nok += sum(1 for r in rows if r["status"] == "ok")
                nerr += sum(1 for r in rows if r["status"] != "ok")
                if n % 100 == 0:
                    print(f"  [qmc] {n}/{len(jobs)} items  ok={nok} err={nerr}  {time.time()-t0:.0f}s", flush=True)
        print(f"[qmc] done ok={nok} err={nerr} in {time.time()-t0:.0f}s")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
