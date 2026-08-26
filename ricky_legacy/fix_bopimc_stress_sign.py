#!/usr/bin/env python3
"""Correct the stress/pressure sign for the BOPIMC hcp DFT files.

The a_dft root trajectory f15b_bopimc-hcp-n96-rs1.49-t800.traj stores its stress
tensor with the opposite sign to the ASE convention used by every other
trajectory in the legacy QMC-HAMM set. Taken literally it puts dense hydrogen at
rs = 1.49 near -150 GPa, which is unphysical; the magnitude instead matches this
dataset's own P(rs) curve (rs = 1.49 -> about 160 GPa from the npt DMC set).

This negates observables/stress and observables/pressure, sets the root
`pressure` attribute, moves each file into the correct P<P>T<T> bucket, and
records provenance/stress_sign_convention_flipped so the change is auditable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

DEFAULT_OUT = Path(
    "/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/ricky_legacy_database"
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--apply", action="store_true", help="write changes (default is a dry run)")
    args = ap.parse_args()

    ledger_path = args.out_root / "ricky_ingest_ledger.jsonl"
    rows = [json.loads(l) for l in ledger_path.read_text(encoding="utf-8").splitlines() if l.strip()]

    targets = []
    for r in rows:
        if r.get("status") != "ok" or r.get("method") != "DFT":
            continue
        if r.get("P") is not None and float(r["P"]) < 0:
            targets.append(r)
    print(f"files with negative pressure: {len(targets)}")
    if not targets:
        return 0

    names = {r.get("girder_item_name") for r in targets}
    print(f"source trajectories involved: {names}")
    if not args.apply:
        print("DRY RUN -- rerun with --apply to modify files")
        return 0

    moved = []
    for n, r in enumerate(targets, 1):
        path = Path(r["path"])
        with h5py.File(path, "r+") as f:
            f["observables/stress"][...] = -np.asarray(f["observables/stress"][...])
            f["observables/pressure"][...] = -np.asarray(f["observables/pressure"][...])
            new_p = float(f["observables/pressure"][0])
            f.attrs["pressure"] = new_p
            g = f["provenance"]
            if "stress_sign_convention_flipped" in g:
                del g["stress_sign_convention_flipped"]
            g.create_dataset("stress_sign_convention_flipped", data=np.array([True]))
            temp = float(f.attrs.get("temperature", float("nan")))
        # rebucket into the right P<P>T<T> directory
        p_tag = f"P{int(round(new_p))}" if np.isfinite(new_p) else "Pnan"
        t_tag = f"T{int(round(temp))}" if np.isfinite(temp) else "Tnan"
        dest_dir = args.out_root / f"{p_tag}{t_tag}"
        dest = dest_dir / path.name
        if dest != path:
            dest_dir.mkdir(parents=True, exist_ok=True)
            path.replace(dest)
            r["path"] = str(dest)
            moved.append(str(dest))
        r["P"] = new_p
        if n % 500 == 0:
            print(f"  patched {n}/{len(targets)}", flush=True)

    # rewrite the ledger with corrected paths/pressures
    tmp = ledger_path.with_suffix(".jsonl.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    tmp.replace(ledger_path)

    print(f"patched {len(targets)} files; {len(moved)} relocated")
    if moved:
        print(f"  new bucket: {Path(moved[0]).parent.name}")
    # clean up any now-empty source buckets
    for d in sorted(args.out_root.glob("P-*T*")):
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
            print(f"  removed empty dir {d.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
