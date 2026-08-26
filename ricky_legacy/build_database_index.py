#!/usr/bin/env python3
"""Build one index of every schema HDF5 file we have produced, from all sources.

This is the single authoritative record of what exists in the database. It is
what the upload step consults so nothing is ever pushed twice, and what the
coverage plots read so every source is represented consistently.

Sources scanned:
  * Aurora/SG hydrogen runs  -- aurora_backup/{run_*,oldrun_*}/build_database_ledger.txt
  * Legacy QMC-HAMM (Ricky)  -- ricky_legacy_database/**/*.h5

Outputs (written next to --aurora-root):
  database_index.csv     one row per HDF5 file, with the indexable root attrs
  database_index_summary.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import h5py
import pandas as pd

DEFAULT_AURORA = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
LEDGER_NAME = "build_database_ledger.txt"
RICKY_DIR = "ricky_legacy_database"

FIELDS = [
    "uuid", "method", "calculation_type", "source_dataset", "source_database",
    "T_K", "P_GPa", "rs", "natoms", "formula", "ensemble", "quantum_ions",
    "creation_date", "qmc_quality", "finite_size_corrected", "name_in_system",
    "author", "starting_configuration_model_name", "n_dependent_uuids",
    "path", "size_bytes", "sha256",
]


def _as_str(v: Any) -> Optional[str]:
    if v is None:
        return None
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return str(v)


def _as_float(v: Any) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scan_one(args_tuple) -> Optional[Dict[str, Any]]:
    path_str, do_hash = args_tuple
    path = Path(path_str)
    try:
        st = path.stat()
        with h5py.File(path, "r") as f:
            a = f.attrs
            row: Dict[str, Any] = {
                "uuid": _as_str(a.get("uuid")) or path.stem,
                "method": _as_str(a.get("method")),
                "calculation_type": _as_str(a.get("calculation_type")),
                "source_dataset": _as_str(a.get("source_dataset")) or "aurora_sg",
                "source_database": _as_str(a.get("source_database")) or "local",
                "T_K": _as_float(a.get("temperature")),
                "P_GPa": _as_float(a.get("pressure")),
                "rs": _as_float(a.get("rs")),
                "natoms": _as_float(a.get("natoms")),
                "formula": _as_str(a.get("formula")),
                "ensemble": _as_str(a.get("ensemble")),
                "quantum_ions": _as_str(a.get("quantum_ions")),
                "creation_date": _as_str(a.get("creation_date")),
                "qmc_quality": _as_float(a.get("qmc_quality")),
                "name_in_system": _as_str(a.get("name_in_system")),
                "author": _as_str(a.get("author")),
                "starting_configuration_model_name": _as_str(a.get("starting_configuration_model_name")),
                "path": str(path),
                "size_bytes": int(st.st_size),
            }
            obs = f.get("observables")
            row["finite_size_corrected"] = (
                _as_str(obs.attrs.get("finite_size_corrected")) if obs is not None else None
            )
            dep = f.get("provenance/dependent_uuids")
            row["n_dependent_uuids"] = int(dep.shape[0]) if dep is not None and dep.shape else 0
        row["sha256"] = sha256_file(path) if do_hash else None
        return row
    except Exception:  # noqa: BLE001
        return {"uuid": path.stem, "path": str(path), "method": None, "error": True}


def collect_sg_paths(aurora_root: Path) -> List[Path]:
    """HDF5 files recorded ok in the Aurora/SG per-run ledgers."""
    out: List[Path] = []
    seen = set()
    patterns = [f"run_*/{LEDGER_NAME}", f"oldrun_*/{LEDGER_NAME}"]
    for pat in patterns:
        for ledger in sorted(aurora_root.glob(pat)):
            for line in ledger.read_text(encoding="utf-8", errors="ignore").splitlines():
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3 or parts[1] != "ok":
                    continue
                cfg = Path(parts[0])
                for fn in (x.strip() for x in parts[2].split(",")):
                    if not fn:
                        continue
                    p = cfg / fn
                    if p.exists() and str(p) not in seen:
                        seen.add(str(p))
                        out.append(p)
    return out


def collect_ricky_paths(aurora_root: Path) -> List[Path]:
    root = aurora_root / RICKY_DIR
    if not root.is_dir():
        return []
    return sorted(root.glob("*/*.h5"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aurora-root", type=Path, default=DEFAULT_AURORA)
    ap.add_argument("--out", type=Path, default=None, help="defaults to <aurora-root>/database_index.csv")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--hash", action="store_true", help="also compute sha256 of every file (slower)")
    args = ap.parse_args()

    out_csv = args.out or (args.aurora_root / "database_index.csv")

    sg = collect_sg_paths(args.aurora_root)
    ricky = collect_ricky_paths(args.aurora_root)
    print(f"Aurora/SG HDF5 files from ledgers : {len(sg)}")
    print(f"Ricky legacy HDF5 files on disk   : {len(ricky)}")
    all_paths = sg + ricky
    if not all_paths:
        raise SystemExit("no HDF5 files found")

    jobs = [(str(p), args.hash) for p in all_paths]
    rows: List[Dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for n, row in enumerate(ex.map(scan_one, jobs, chunksize=64), 1):
            if row is not None:
                rows.append(row)
            if n % 10000 == 0:
                print(f"  scanned {n}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows)
    bad = int(df.get("error", pd.Series(dtype=bool)).fillna(False).sum()) if "error" in df else 0
    df = df[[c for c in FIELDS if c in df.columns]]
    df.to_csv(out_csv, index=False)

    summary = {
        "n_files": int(len(df)),
        "n_unreadable": bad,
        "by_source_dataset": {str(k): int(v) for k, v in df["source_dataset"].value_counts(dropna=False).items()},
        "by_method": {str(k): int(v) for k, v in df["method"].value_counts(dropna=False).items()},
        "total_size_GB": round(float(df["size_bytes"].sum()) / 1e9, 3),
        "n_unique_uuid": int(df["uuid"].nunique()),
        "T_K_range": [_as_float(df["T_K"].min()), _as_float(df["T_K"].max())],
        "P_GPa_range": [_as_float(df["P_GPa"].min()), _as_float(df["P_GPa"].max())],
        "rs_range": [_as_float(df["rs"].min()), _as_float(df["rs"].max())],
        "index_csv": str(out_csv),
    }
    (out_csv.parent / "database_index_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    dupes = int(len(df) - df["uuid"].nunique())
    if dupes:
        print(f"WARNING: {dupes} duplicate uuids in index")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
