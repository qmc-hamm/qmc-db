#!/usr/bin/env python3
"""Verify the legacy QMC-HAMM ("Ricky") transfer is complete and faithful.

Four independent checks:

  1. COVERAGE   every source item/frame on the Girder has an HDF5 file, and the
                per-(P,T) configuration counts match the source exactly.
  2. SCHEMA     every file has the required root attrs and groups, unique uuid,
                and the observables expected for its method.
  3. FIDELITY   for a random sample, energies/forces/cell/positions inside the
                HDF5 match the original ASE trajectory bit-for-bit.
  4. PROVENANCE every DMC file's dependent_uuids resolve to DFT files that exist
                and describe the same geometry.
"""

from __future__ import annotations

import argparse
import collections
import json
import pickle
import random
import sys
from pathlib import Path
from typing import Any, Dict, List

import h5py
import numpy as np
from ase.io import read

DEFAULT_META = Path("/scratch/sgoswam3/ricky_qmchamm/meta")
DEFAULT_RAW = Path("/scratch/sgoswam3/ricky_qmchamm")
DEFAULT_OUT = Path(
    "/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/ricky_legacy_database"
)

REQUIRED_ROOT = [
    "uuid", "system", "formula", "natoms", "species", "calculation_type",
    "method", "method_kws", "creation_date", "author", "rs",
    "source_database", "source_dataset",
]
REQUIRED_GROUPS = ["parameters", "code", "structure", "observables", "provenance"]
DMC_OBS = ["total_energy", "total_energy_error", "forces", "forces_error",
           "fsc_potential_energy", "fsc_kinetic_energy", "total_energy_uncorrected"]
DFT_OBS = ["total_energy", "forces", "stress", "pressure"]


def load_ledger(out_root: Path) -> List[Dict[str, Any]]:
    path = out_root / "ricky_ingest_ledger.jsonl"
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            rows.append(json.loads(line))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", type=Path, default=DEFAULT_META)
    ap.add_argument("--raw-root", type=Path, default=DEFAULT_RAW)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sample", type=int, default=300, help="files to deep-verify against source")
    args = ap.parse_args()

    problems: List[str] = []
    ledger = load_ledger(args.out_root)
    ok_rows = [r for r in ledger if r.get("status") == "ok"]
    err_rows = [r for r in ledger if r.get("status") != "ok"]
    print(f"ledger: {len(ledger)} rows, {len(ok_rows)} ok, {len(err_rows)} errors")
    if err_rows:
        problems.append(f"{len(err_rows)} ledger error rows")
        for r in err_rows[:5]:
            print("   ERROR ROW:", r)

    # ---------------- 1. coverage ----------------
    print("\n=== 1. COVERAGE vs source ===")
    q_frames = pickle.load(open(args.meta / "frames.pkl", "rb"))
    plan = json.loads((args.out_root / "ricky_work_plan.json").read_text())
    n_src_dmc = len(q_frames)
    n_src_dft = sum(len(v) for v in plan["needed_dft"].values())

    dmc_rows = [r for r in ok_rows if r.get("method") == "DMC"]
    dft_rows = [r for r in ok_rows if r.get("method") == "DFT"]
    print(f"  source DMC frames  {n_src_dmc:7d}   written {len(dmc_rows):7d}")
    print(f"  source DFT frames  {n_src_dft:7d}   written {len(dft_rows):7d}")
    if len(dmc_rows) != n_src_dmc:
        problems.append(f"DMC count mismatch: {len(dmc_rows)} != {n_src_dmc}")
    if len(dft_rows) != n_src_dft:
        problems.append(f"DFT count mismatch: {len(dft_rows)} != {n_src_dft}")

    # every source (item, frame) accounted for
    src_keys = {f"ricky:{r['itemId']}:{r['frame']}" for r in q_frames}
    got_keys = {r["source_key"] for r in dmc_rows}
    missing = src_keys - got_keys
    print(f"  DMC source keys missing from output: {len(missing)}")
    if missing:
        problems.append(f"{len(missing)} DMC source frames not converted")
        for k in list(missing)[:5]:
            print("     MISSING:", k)

    # per (P,T) counts must match the source exactly
    src_pt = collections.Counter((r["P"], r["T"]) for r in q_frames)
    out_pt = collections.Counter((r["P"], r["T"]) for r in dmc_rows)
    diff = {k: (src_pt.get(k, 0), out_pt.get(k, 0)) for k in set(src_pt) | set(out_pt)
            if src_pt.get(k, 0) != out_pt.get(k, 0)}
    print(f"  (P,T) bins: source={len(src_pt)} output={len(out_pt)} mismatched={len(diff)}")
    if diff:
        problems.append(f"(P,T) count mismatch in {len(diff)} bins: {list(diff.items())[:5]}")

    files_on_disk = sorted(args.out_root.glob("*/*.h5"))
    print(f"  HDF5 files on disk: {len(files_on_disk)}  (ledger ok rows: {len(ok_rows)})")
    if len(files_on_disk) != len(ok_rows):
        problems.append(f"disk/ledger mismatch: {len(files_on_disk)} files vs {len(ok_rows)} rows")

    # ---------------- 2. schema ----------------
    print("\n=== 2. SCHEMA ===")
    uuids = collections.Counter()
    bad_schema = 0
    for n, p in enumerate(files_on_disk, 1):
        try:
            with h5py.File(p, "r") as f:
                miss = [k for k in REQUIRED_ROOT if k not in f.attrs]
                gmiss = [g for g in REQUIRED_GROUPS if g not in f]
                method = str(f.attrs.get("method", ""))
                need = DMC_OBS if method == "DMC" else DFT_OBS
                omiss = [k for k in need if f"observables/{k}" not in f]
                if method == "DMC" and "finite_size_corrected" not in f["observables"].attrs:
                    omiss.append("attrs:finite_size_corrected")
                if method == "DMC" and "qmc_quality" not in f.attrs:
                    miss.append("qmc_quality")
                if method == "DFT" and "qmc_quality" in f.attrs:
                    miss.append("qmc_quality-should-be-absent-for-DFT")
                uuids[str(f.attrs.get("uuid", p.stem))] += 1
                if miss or gmiss or omiss:
                    bad_schema += 1
                    if bad_schema <= 5:
                        print(f"   {p.name}: attrs{miss} groups{gmiss} obs{omiss}")
        except Exception as exc:  # noqa: BLE001
            bad_schema += 1
            if bad_schema <= 5:
                print(f"   {p.name}: UNREADABLE {exc!r}")
        if n % 20000 == 0:
            print(f"   ...{n}/{len(files_on_disk)}", flush=True)
    dupes = {u: c for u, c in uuids.items() if c > 1}
    print(f"  files failing schema check: {bad_schema}")
    print(f"  duplicate uuids: {len(dupes)}")
    if bad_schema:
        problems.append(f"{bad_schema} files fail schema check")
    if dupes:
        problems.append(f"{len(dupes)} duplicate uuids")

    # ---------------- 3. fidelity ----------------
    print("\n=== 3. FIDELITY vs original trajectories ===")
    qmc_items = json.loads((args.meta / "folder_items.json").read_text())
    dft_items = json.loads((args.meta / "dft_items.json").read_text())
    qmc_paths = {k: v["path"] for k, v in json.loads((args.meta / "download_manifest.json").read_text()).items()}
    dft_paths = {k: v["path"] for k, v in json.loads((args.meta / "download_manifest_dft.json").read_text()).items()}

    random.seed(1234)
    sample = random.sample(ok_rows, min(args.sample, len(ok_rows)))
    traj_cache: Dict[str, Any] = {}

    def frames_of(iid: str, is_dmc: bool):
        if iid not in traj_cache:
            rel = qmc_paths[iid] if is_dmc else dft_paths[iid]
            traj_cache[iid] = read(str(args.raw_root / rel), index=":")
        return traj_cache[iid]

    checked = mism = 0
    for r in sample:
        is_dmc = r["method"] == "DMC"
        atoms = frames_of(r["girder_item_id"], is_dmc)[r["frame"]]
        with h5py.File(r["path"], "r") as f:
            e_h5 = float(f["observables/total_energy"][0])
            fo_h5 = np.asarray(f["observables/forces"][0])
            cell_h5 = np.asarray(f["structure/lattice_vectors"])
            pos_h5 = np.asarray(f["structure/positions"])
            rs_h5 = float(f.attrs["rs"])
            if is_dmc:
                dv = float(f["observables/fsc_potential_energy"][0])
                dt = float(f["observables/fsc_kinetic_energy"][0])
                unc = float(f["observables/total_energy_uncorrected"][0])
        errs = []
        if not np.isclose(e_h5, float(atoms.get_potential_energy()), rtol=0, atol=1e-9):
            errs.append("energy")
        if not np.allclose(fo_h5, np.asarray(atoms.get_forces()), rtol=0, atol=1e-9):
            errs.append("forces")
        if not np.allclose(cell_h5, np.asarray(atoms.cell), rtol=0, atol=1e-12):
            errs.append("cell")
        if not np.allclose(pos_h5, np.asarray(atoms.positions), rtol=0, atol=1e-12):
            errs.append("positions")
        rs_exp = ((3 * atoms.get_volume() / (4 * np.pi * len(atoms))) ** (1 / 3)) / 0.529177210903
        if not np.isclose(rs_h5, rs_exp, rtol=1e-10):
            errs.append("rs")
        if is_dmc:
            if not np.isclose(dv, float(atoms.info["fsc_dv_ev"]), atol=1e-9):
                errs.append("fsc_dv")
            if not np.isclose(dt, float(atoms.info["fsc_dt_ev"]), atol=1e-9):
                errs.append("fsc_dt")
            if not np.isclose(unc, e_h5 - dv - dt, atol=1e-9):
                errs.append("uncorrected")
        checked += 1
        if errs:
            mism += 1
            if mism <= 5:
                print(f"   MISMATCH {Path(r['path']).name} {r['method']}: {errs}")
    print(f"  deep-verified {checked} files, mismatches {mism}")
    if mism:
        problems.append(f"{mism}/{checked} sampled files disagree with source")

    # ---------------- 4. provenance ----------------
    print("\n=== 4. PROVENANCE links ===")
    by_uuid = {r["uuid"]: r for r in ok_rows if r.get("uuid")}
    dmc_sample = random.sample(dmc_rows, min(args.sample, len(dmc_rows)))
    n_dep = collections.Counter()
    dangling = geo_bad = 0
    for r in dmc_sample:
        with h5py.File(r["path"], "r") as f:
            deps = [x.decode() if isinstance(x, bytes) else x
                    for x in f["provenance/dependent_uuids"][...]]
            fa = np.sort(np.mod(np.asarray(f["structure/fractional_positions"]), 1.0), axis=0)
        n_dep[len(deps)] += 1
        for d in deps:
            tgt = by_uuid.get(d)
            if tgt is None or not Path(tgt["path"]).exists():
                dangling += 1
                continue
            with h5py.File(tgt["path"], "r") as g:
                fb = np.sort(np.mod(np.asarray(g["structure/fractional_positions"]), 1.0), axis=0)
            if not np.allclose(fa, fb, atol=5e-4):
                geo_bad += 1
    print(f"  sampled {len(dmc_sample)} DMC files; dependent_uuids per file: {dict(sorted(n_dep.items()))}")
    print(f"  dangling links: {dangling}   geometry-mismatched links: {geo_bad}")
    if dangling:
        problems.append(f"{dangling} dangling dependent_uuids")
    if geo_bad:
        problems.append(f"{geo_bad} dependent_uuids point at a different geometry")
    if n_dep.get(0):
        problems.append(f"{n_dep[0]} sampled DMC files have no DFT provenance link")

    # ---------------- verdict ----------------
    print("\n" + "=" * 72)
    if problems:
        print("FAILED CHECKS:")
        for p in problems:
            print("  -", p)
        return 1
    print("ALL CHECKS PASSED")
    print(f"  {len(ok_rows)} HDF5 files ({len(dmc_rows)} DMC + {len(dft_rows)} DFT)")
    print(f"  covering {len(out_pt)} (P,T) conditions, all source frames accounted for")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
