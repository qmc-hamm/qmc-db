#!/usr/bin/env python3
"""Summarize build_database ledgers so you can monitor progress at a glance.

Layout:
  Each run_YYYY_MM_DD folder owns its own ledger and log:
      <aurora_root>/run_YYYY_MM_DD/build_database_ledger.txt
      <aurora_root>/run_YYYY_MM_DD/build_database.log
  HDF5 outputs live inside the configuration folder.

Usage:
    ./verify_ledger.py                              # aggregate all per-run ledgers
    ./verify_ledger.py --show-failed                # also print failing rows
    ./verify_ledger.py --include-runs run_2025_05_25,run_2025_06_01
    ./verify_ledger.py --ledger /path/to/build_database_ledger.txt  # inspect one ledger
    ./verify_ledger.py --no-aurora                  # skip discovery (faster)

Reports:
  - Counts per status (ok / skipped / error / invalid / upload_failed)
  - Per-run-root breakdown
  - Per-config HDF5 audit: expected files (from `ok` rows) vs. files on disk
    in each configuration folder
  - Discovery coverage: pending vs. total configs under <aurora_root>
  - With --show-failed, every non-ok row with its message
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DEFAULT_AURORA_ROOT = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
LEDGER_NAME = "build_database_ledger.txt"

ALL_STATUSES = ["ok", "skipped", "error", "invalid", "upload_failed"]
CONFIG_FOLDER_RE = re.compile(r"^P\d+T\d+config\d+$|^rs\d+(?:\.\d+)?T\d+config\d+$")


def parse_ledger(ledger_path: Path) -> Dict[str, Tuple[str, str]]:
    """Return latest (status, message) per config path."""
    if not ledger_path.exists():
        return {}
    latest: Dict[str, Tuple[str, str]] = {}
    for raw in ledger_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 2:
            continue
        latest[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else "")
    return latest


def discover_configs(aurora_root: Path) -> List[Path]:
    out: List[Path] = []
    for run_dir in sorted(aurora_root.glob("run_*/runs")):
        for cfg_root in sorted(p for p in run_dir.iterdir() if p.is_dir()):
            for config in sorted(p for p in cfg_root.iterdir() if p.is_dir()):
                if CONFIG_FOLDER_RE.match(config.name):
                    out.append(config)
    return out


def run_root_of(path_str: str) -> str:
    p = Path(path_str)
    return p.parents[2].name if len(p.parents) >= 3 else "<unknown>"


def gather_rows(
    aurora_root: Path,
    explicit_ledger: Optional[Path],
    include_runs: set,
) -> List[Tuple[str, str, str]]:
    """Return list of (path, status, message) from one ledger or all per-run ledgers."""
    if explicit_ledger is not None:
        return [(p, s, m) for p, (s, m) in parse_ledger(explicit_ledger).items()]
    rows: List[Tuple[str, str, str]] = []
    for ledger in sorted(aurora_root.glob(f"run_*/{LEDGER_NAME}")):
        run_name = ledger.parent.name
        if include_runs and run_name not in include_runs:
            continue
        for p, (s, m) in parse_ledger(ledger).items():
            rows.append((p, s, m))
    return rows


def fmt_table(rows: List[List[str]], headers: List[str]) -> str:
    widths = [len(h) for h in headers]
    for r in rows:
        for i, cell in enumerate(r):
            widths[i] = max(widths[i], len(str(cell)))
    sep = "  "
    lines = [sep.join(h.ljust(widths[i]) for i, h in enumerate(headers))]
    lines.append(sep.join("-" * widths[i] for i in range(len(headers))))
    for r in rows:
        lines.append(sep.join(str(c).ljust(widths[i]) for i, c in enumerate(r)))
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aurora-root", type=Path, default=DEFAULT_AURORA_ROOT)
    ap.add_argument("--ledger", type=Path, default=None, help="Inspect a single specific ledger file.")
    ap.add_argument("--include-runs", type=str, default="", help="Comma list of run_YYYY_MM_DD names to include.")
    ap.add_argument("--no-aurora", action="store_true", help="Skip aurora_root coverage scan (faster).")
    ap.add_argument("--show-failed", action="store_true", help="Print every non-ok row with its message.")
    ap.add_argument("--show-files", action="store_true", help="Print sample missing / unexpected files.")
    args = ap.parse_args()

    include_runs = {s.strip() for s in args.include_runs.split(",") if s.strip()}

    if args.ledger is not None:
        print(f"Ledger:  {args.ledger}")
    else:
        print(f"Aurora:  {args.aurora_root}")
        print(f"Ledgers: {args.aurora_root}/run_*/{LEDGER_NAME}")
        if include_runs:
            print(f"Filter:  {sorted(include_runs)}")
    print()

    rows = gather_rows(args.aurora_root, args.ledger, include_runs)
    if not rows:
        print("No ledger entries found.")
        return 0

    by_status: Counter[str] = Counter()
    by_run_status: Dict[str, Counter[str]] = defaultdict(Counter)
    expected_files_per_cfg: Dict[str, List[str]] = {}
    failed: List[Tuple[str, str, str]] = []
    for path_str, status, message in rows:
        by_status[status] += 1
        by_run_status[run_root_of(path_str)][status] += 1
        if status == "ok":
            expected_files_per_cfg[path_str] = [t.strip() for t in message.split(",") if t.strip()]
        if status != "ok":
            failed.append((path_str, status, message))

    # Status summary
    total = sum(by_status.values())
    status_rows = [[s, by_status.get(s, 0)] for s in ALL_STATUSES if by_status.get(s, 0) or s == "ok"]
    print("== Status summary ==")
    print(fmt_table(status_rows, ["status", "count"]))
    print(f"total: {total}")
    print()

    # Per-run breakdown
    print("== Per-run breakdown ==")
    headers = ["run_root"] + ALL_STATUSES + ["total"]
    run_rows: List[List[str]] = []
    for run, counts in sorted(by_run_status.items()):
        row = [run] + [counts.get(s, 0) for s in ALL_STATUSES] + [sum(counts.values())]
        run_rows.append([str(x) for x in row])
    print(fmt_table(run_rows, headers))
    print()

    # Aurora discovery coverage
    if not args.no_aurora and args.aurora_root.exists():
        discovered = discover_configs(args.aurora_root)
        if include_runs:
            discovered = [c for c in discovered if c.parents[2].name in include_runs]
        ledger_paths = {r[0] for r in rows}
        seen_ok = sum(1 for p, s, _ in rows if s == "ok")
        coverage = (seen_ok / len(discovered) * 100.0) if discovered else 0.0
        not_in_ledger = [p for p in discovered if str(p) not in ledger_paths]
        print("== Discovery coverage ==")
        print(f"discovered configs:  {len(discovered)}")
        print(f"ledger rows (any):   {len(rows)}")
        print(f"ledger rows ok:      {seen_ok}  ({coverage:.1f}%)")
        print(f"missing from ledger: {len(not_in_ledger)}")
        if args.show_files and not_in_ledger:
            print("\nfirst 10 not-yet-processed configs:")
            for p in not_in_ledger[:10]:
                print(f"  {p}")
        print()

    # On-disk HDF5 audit: each config folder should contain the files
    # named in its `ok` ledger row.
    n_expected = sum(len(v) for v in expected_files_per_cfg.values())
    missing_from_disk: List[Tuple[str, str]] = []  # (config_path, filename)
    for cfg_path, fnames in expected_files_per_cfg.items():
        cfg = Path(cfg_path)
        on_disk = {p.name for p in cfg.glob("*.h5")} if cfg.exists() else set()
        for f in fnames:
            if f not in on_disk:
                missing_from_disk.append((cfg_path, f))
    print("== HDF5 file audit (per-config) ==")
    print(f"expected (from ledger): {n_expected}")
    print(f"missing from disk:      {len(missing_from_disk)}")
    if args.show_files and missing_from_disk:
        print("\nmissing (first 10):")
        for cfg, f in missing_from_disk[:10]:
            print(f"  {f}  <- expected in {cfg}")
    print()

    # Failed rows
    if failed:
        print(f"== Failed/non-ok rows ({len(failed)}) ==")
        if args.show_failed:
            for path_str, status, message in failed:
                print(f"[{status}] {path_str}")
                if message:
                    print(f"    {message}")
        else:
            print("(re-run with --show-failed to list each)")
        print()

    print(f"SUMMARY status_counts={dict(by_status)} missing_from_disk={len(missing_from_disk)}")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
