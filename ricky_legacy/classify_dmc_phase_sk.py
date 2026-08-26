#!/usr/bin/env python3
"""Classify every Ricky-legacy DMC config as solid / liquid / ambiguous via S(k).

References (known labels):
  * solid  — BOPIMC hcp DFT frames (ricky_legacy, name contains bopimc)
  * liquid — Aurora/SG HDF5 whose EXTXYZ header has phase=liquid

For each configuration we compute ionic structure-factor metrics on the
reciprocal-lattice shell |G| <= 12 / Angstrom:
  S(G) = |sum_j exp(i G · r_j)|^2 / N
then S_max, n_Sgt10, n_SgtN4, etc.

Clear solid / clear liquid thresholds are taken from the reference
distributions; everything in between is ambiguous.

Usage (typically via sbatch):
  OPENBLAS_NUM_THREADS=1 python classify_dmc_phase_sk.py
  OPENBLAS_NUM_THREADS=1 python classify_dmc_phase_sk.py --apply
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List

import h5py
import numpy as np
import pandas as pd

# Force single-threaded BLAS inside every worker process.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

ROOT = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
DEFAULT_INDEX = ROOT / "database_index.csv"
DEFAULT_OUT = Path("/scratch/sgoswam3/ricky_qmchamm/phase_classify")
KMAX = 12.0


def sk_metrics(pos: np.ndarray, cell: np.ndarray, kmax: float = KMAX) -> Dict[str, Any]:
    pos = np.asarray(pos, dtype=float)
    cell = np.asarray(cell, dtype=float)
    n_atoms = len(pos)
    b_mat = 2.0 * np.pi * np.linalg.inv(cell).T
    bmin = float(np.linalg.norm(b_mat, axis=1).min())
    nmax = int(np.ceil(kmax / bmin)) + 1
    rng = np.arange(-nmax, nmax + 1)
    ii, jj, kk = np.meshgrid(rng, rng, rng, indexing="ij")
    mask = ~((ii == 0) & (jj == 0) & (kk == 0))
    g_vecs = (
        ii[mask][:, None] * b_mat[0]
        + jj[mask][:, None] * b_mat[1]
        + kk[mask][:, None] * b_mat[2]
    )
    gnorm = np.linalg.norm(g_vecs, axis=1)
    g_vecs = g_vecs[gnorm <= kmax]
    if g_vecs.size == 0:
        return {
            "S_max": np.nan, "S_p95": np.nan, "S_mean": np.nan, "S_median": np.nan,
            "n_Sgt5": 0, "n_Sgt10": 0, "n_Sgt20": 0, "n_SgtN4": 0,
            "N": n_atoms, "n_G": 0, "S_max_over_mean": np.nan,
        }
    chunks: List[np.ndarray] = []
    for start in range(0, len(g_vecs), 2000):
        gc = g_vecs[start : start + 2000]
        phase = np.exp(1j * (pos @ gc.T))
        chunks.append((np.abs(phase.sum(axis=0)) ** 2) / n_atoms)
    s_vals = np.concatenate(chunks)
    s_max = float(s_vals.max())
    s_mean = float(s_vals.mean())
    return {
        "S_max": s_max,
        "S_p95": float(np.percentile(s_vals, 95)),
        "S_mean": s_mean,
        "S_median": float(np.median(s_vals)),
        "n_Sgt5": int((s_vals > 5).sum()),
        "n_Sgt10": int((s_vals > 10).sum()),
        "n_Sgt20": int((s_vals > 20).sum()),
        "n_SgtN4": int((s_vals > n_atoms / 4.0).sum()),
        "N": n_atoms,
        "n_G": int(len(s_vals)),
        "S_max_over_mean": float(s_max / max(s_mean, 1e-9)),
    }


def metrics_from_h5(path: str) -> Dict[str, Any]:
    with h5py.File(path, "r") as handle:
        pos = np.asarray(handle["structure/positions"])
        cell = np.asarray(handle["structure/lattice_vectors"])
        row = sk_metrics(pos, cell)
        row["path"] = path
        row["uuid"] = Path(path).stem
        for key in ("temperature", "pressure", "rs", "method", "name_in_system"):
            val = handle.attrs.get(key)
            if isinstance(val, bytes):
                val = val.decode("utf-8", "replace")
            row[key] = val
    return row


def collect_liquid_paths(index: pd.DataFrame, cache: Path) -> List[str]:
    if cache.exists():
        return json.loads(cache.read_text())
    paths: List[str] = []
    sg = index[index["source_dataset"] == "aurora_sg"]
    for i, path in enumerate(sg["path"], 1):
        try:
            with h5py.File(path, "r") as handle:
                xyz = handle["structure/xyz"][()]
                if isinstance(xyz, bytes):
                    xyz = xyz.decode("utf-8", "replace")
                header = xyz.split("\n")[1] if "\n" in xyz else ""
                match = re.search(r'\bphase=([^\s"]+)', header)
                if match and match.group(1).lower() == "liquid":
                    paths.append(path)
        except Exception:  # noqa: BLE001
            continue
        if i % 1000 == 0:
            print(f"  scanned aurora for liquid tags: {i}/{len(sg)}", flush=True)
    cache.write_text(json.dumps(paths))
    return paths


def map_metrics(paths: List[str], tag: str, workers: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, row in enumerate(pool.map(metrics_from_h5, paths, chunksize=8), 1):
            row["ref_class"] = tag
            rows.append(row)
            if i % 200 == 0 or i == len(paths):
                print(f"  {tag}: {i}/{len(paths)}  {time.time() - t0:.1f}s", flush=True)
    return rows


def pct(series: pd.Series, q: float) -> float:
    return float(np.percentile(series.to_numpy(dtype=float), q))


def build_thresholds(ref_df: pd.DataFrame) -> Dict[str, float]:
    solid = ref_df[ref_df["ref_class"] == "solid"]
    liquid = ref_df[ref_df["ref_class"] == "liquid"]
    thr = {
        "solid_S_max_p5": pct(solid["S_max"], 5),
        "solid_n10_p5": pct(solid["n_Sgt10"], 5),
        "solid_nN4_p5": pct(solid["n_SgtN4"], 5),
        "liquid_S_max_p95": pct(liquid["S_max"], 95),
        "liquid_S_max_p995": pct(liquid["S_max"], 99.5),
        "liquid_n10_p95": pct(liquid["n_Sgt10"], 95),
        "liquid_n10_p995": pct(liquid["n_Sgt10"], 99.5),
        "liquid_nN4_p995": pct(liquid["n_SgtN4"], 99.5),
    }
    # Clear solid must beat essentially all known liquids.
    thr["clear_solid_S_max"] = float(max(thr["liquid_S_max_p995"], 15.0))
    thr["clear_solid_n10"] = float(max(thr["liquid_n10_p995"], 3.0))
    # Clear liquid must sit in the liquid bulk and below the solid floor.
    thr["clear_liquid_S_max"] = float(min(thr["solid_S_max_p5"] * 0.7, thr["liquid_S_max_p95"]))
    thr["clear_liquid_n10"] = float(thr["liquid_n10_p95"])
    return thr


def classify(s_max: float, n10: int, n_n4: int, thr: Dict[str, float]) -> str:
    if (
        (s_max >= thr["clear_solid_S_max"] and n10 >= thr["clear_solid_n10"])
        or (n_n4 >= 2 and s_max >= thr["clear_solid_S_max"])
    ):
        return "solid"
    if s_max <= thr["clear_liquid_S_max"] and n10 <= thr["clear_liquid_n10"] and n_n4 == 0:
        return "liquid"
    # Secondary liquid: no strong Bragg peaks and well below solid floor.
    if (
        n_n4 == 0
        and s_max <= thr["liquid_S_max_p95"]
        and n10 <= max(2.0, thr["liquid_n10_p95"])
        and s_max < thr["solid_S_max_p5"] * 0.7
    ):
        return "liquid"
    return "ambiguous"


def stamp_one(args_tuple):
    """Top-level for ProcessPool pickling: (path, phase) -> phase."""
    path, phase = args_tuple
    with h5py.File(path, "a") as handle:
        handle.attrs["phase"] = phase
        handle.attrs["phase_source"] = "sk_vs_refs_v1"
    return phase


def apply_labels(dmc_df: pd.DataFrame, workers: int = 16) -> None:
    """Write root attr phase onto each DMC HDF5."""
    jobs = list(zip(dmc_df["path"].astype(str), dmc_df["phase_pred"].astype(str)))
    counts: Dict[str, int] = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, phase in enumerate(pool.map(stamp_one, jobs, chunksize=32), 1):
            counts[phase] = counts.get(phase, 0) + 1
            if i % 2000 == 0 or i == len(jobs):
                print(f"  stamped {i}/{len(jobs)}", flush=True)
    print("applied counts:", counts, flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--solid-sample", type=int, default=800)
    parser.add_argument("--apply", action="store_true",
                        help="also write phase attrs onto the DMC HDF5 files")
    parser.add_argument("--apply-only", action="store_true",
                        help="skip recomputation; stamp phases from an existing "
                             "dmc_Sk_metrics_classified.csv under --out")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}", flush=True)
    print(f"workers={args.workers}", flush=True)

    if args.apply_only:
        csv_path = args.out / "dmc_Sk_metrics_classified.csv"
        if not csv_path.exists():
            raise SystemExit(f"missing {csv_path}; run a full classify first")
        dmc_df = pd.read_csv(csv_path)
        print(f"apply-only: {len(dmc_df)} rows from {csv_path}", flush=True)
        print(dmc_df["phase_pred"].value_counts().to_string(), flush=True)
        apply_labels(dmc_df, workers=args.workers)
        print("DONE apply-only", flush=True)
        return 0

    index = pd.read_csv(args.index)

    bop = index[
        (index["source_dataset"] == "ricky_legacy")
        & (index["method"] == "DFT")
        & index["name_in_system"].fillna("").str.contains("bopimc")
    ]
    solid_paths = bop["path"].sample(min(args.solid_sample, len(bop)), random_state=0).tolist()
    liquid_paths = collect_liquid_paths(index, args.out / "liquid_ref_paths.json")
    print(f"solid refs: {len(solid_paths)} / {len(bop)} BOPIMC", flush=True)
    print(f"liquid refs: {len(liquid_paths)} Aurora phase=liquid", flush=True)

    solid_rows = map_metrics(solid_paths, "solid", args.workers)
    liquid_rows = map_metrics(liquid_paths, "liquid", args.workers)
    ref_df = pd.DataFrame(solid_rows + liquid_rows)
    ref_df.to_csv(args.out / "reference_Sk_metrics.csv", index=False)
    print("reference summary:", flush=True)
    print(
        ref_df.groupby("ref_class")[
            ["S_max", "S_p95", "S_mean", "n_Sgt10", "n_Sgt20", "n_SgtN4", "S_max_over_mean"]
        ].agg(["median", "mean", "min", "max"]).to_string(),
        flush=True,
    )

    thr = build_thresholds(ref_df)
    (args.out / "thresholds.json").write_text(json.dumps(thr, indent=2))
    print("thresholds:", json.dumps(thr, indent=2), flush=True)

    ref_df["pred"] = [
        classify(r.S_max, int(r.n_Sgt10), int(r.n_SgtN4), thr) for r in ref_df.itertuples()
    ]
    confusion = pd.crosstab(ref_df["ref_class"], ref_df["pred"])
    print("reference self-classification:\n", confusion.to_string(), flush=True)

    dmc = index[(index["source_dataset"] == "ricky_legacy") & (index["method"] == "DMC")]
    dmc_paths = dmc["path"].tolist()
    print(f"DMC sweep n={len(dmc_paths)}", flush=True)
    t0 = time.time()
    dmc_rows: List[Dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, row in enumerate(pool.map(metrics_from_h5, dmc_paths, chunksize=16), 1):
            row["phase_pred"] = classify(
                row["S_max"], int(row["n_Sgt10"]), int(row["n_SgtN4"]), thr
            )
            dmc_rows.append(row)
            if i % 1000 == 0 or i == len(dmc_paths):
                rate = i / max(time.time() - t0, 1e-9)
                print(
                    f"  DMC {i}/{len(dmc_paths)}  {rate:.1f}/s  "
                    f"elapsed={time.time() - t0:.0f}s",
                    flush=True,
                )

    dmc_df = pd.DataFrame(dmc_rows)
    dmc_df = dmc_df.merge(
        dmc[["uuid", "T_K", "P_GPa", "rs", "name_in_system"]],
        on="uuid",
        how="left",
        suffixes=("", "_idx"),
    )
    dmc_df.to_csv(args.out / "dmc_Sk_metrics_classified.csv", index=False)
    print("DMC phase_pred counts:\n", dmc_df["phase_pred"].value_counts().to_string(), flush=True)
    dmc_df["Tbin"] = (pd.to_numeric(dmc_df["temperature"], errors="coerce") / 200).round() * 200
    print("by T bin:\n", pd.crosstab(dmc_df["Tbin"], dmc_df["phase_pred"]).to_string(), flush=True)

    summary = {
        "n_dmc": int(len(dmc_df)),
        "counts": dmc_df["phase_pred"].value_counts().to_dict(),
        "n_solid_refs": len(solid_paths),
        "n_liquid_refs": len(liquid_paths),
        "thresholds": thr,
        "ref_confusion": confusion.to_dict(),
    }
    (args.out / "classification_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    if args.apply:
        print("writing phase attrs onto DMC HDF5 files...", flush=True)
        apply_labels(dmc_df, workers=args.workers)

    print("DONE", json.dumps(summary["counts"]), flush=True)
    print("outputs ->", args.out, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
