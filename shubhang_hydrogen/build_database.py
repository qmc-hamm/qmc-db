#!/usr/bin/env python3
"""Build database HDF5 files for one or many configuration folders.

Artifact placement (no files are written to the script directory):
  * HDF5 outputs go into the configuration folder itself
    e.g. ``<aurora_root>/run_YYYY_MM_DD/runs/<group>/<Pxxx...>/`` .
  * The ledger and log for each ``run_YYYY_MM_DD`` live alongside that root
    as ``<run_root>/build_database_ledger.txt`` and ``<run_root>/build_database.log``.

These defaults can be overridden via ``--out-dir`` / ``--ledger`` / ``--log``
if you really want a flat layout instead.

Modes:
  * ``--single <path>``  Process one configuration folder, append to its
    per-run-root ledger, and exit. Safe for use under ``xargs -P`` / job
    arrays since the ledger and log are file-locked.
  * ``--list-only``      Print pending configuration folders (one per line),
    honoring all per-run-root ledgers and any ``--include-runs`` filter.
  * default              Discover and process all pending configurations
    serially.

For each configuration we:
  1. Call ``build_run_files`` from ``create_hdf5_schema`` to write DFT/VMC/DMC HDF5s
     directly inside the configuration folder.
  2. Validate the resulting files against the project schema.
  3. Optionally upload each validated file to the OSN/S3 bucket.
  4. Atomically append a row to the per-run-root ledger.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import json
import os
import re
import shlex
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from create_hdf5_schema import build_run_files  # type: ignore  # noqa: E402

DEFAULT_AURORA_ROOT = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
LEDGER_NAME = "build_database_ledger.txt"
LOG_NAME = "build_database.log"

CONFIG_FOLDER_RE = re.compile(r"^P\d+T\d+config\d+$|^rs\d+(?:\.\d+)?T\d+config\d+$")
RETRYABLE = {"error", "invalid", "upload_failed"}

REQUIRED_ROOT_ATTRS = {
    "system",
    "formula",
    "natoms",
    "species",
    "calculation_type",
    "method",
    "uuid",
    "creation_date",
}
REQUIRED_GROUPS = {"code", "structure", "observables", "provenance"}


def _now() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


def log(line: str, log_path: Path) -> None:
    """Append a single line to a log file under fcntl lock.

    No stdout/stderr emission so that bash dispatchers that redirect
    child stdout/stderr into the same file don't get duplicate lines.
    """
    msg = f"[{_now()}] pid={os.getpid()} {line}"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            fh.write(msg + "\n")
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def load_ledger(ledger_path: Path) -> Dict[str, str]:
    done: Dict[str, str] = {}
    if not ledger_path.exists():
        return done
    for raw in ledger_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        done[parts[0]] = parts[1] if len(parts) > 1 else "ok"
    return done


def append_ledger(ledger_path: Path, key: str, status: str, extra: str = "") -> None:
    """Atomically append a row to a ledger; safe across concurrent processes."""
    line = "\t".join([key, status, extra]).rstrip() + "\n"
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with ledger_path.open("a", encoding="utf-8") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            fh.write(line)
            fh.flush()
            os.fsync(fh.fileno())
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def run_root_of(config_folder: Path) -> Path:
    """Return the run_YYYY_MM_DD ancestor for a configuration folder."""
    # config_folder = <aurora_root>/run_YYYY_MM_DD/runs/<group>/<config>
    return config_folder.parents[2]


def ledger_for(config_folder: Path) -> Path:
    return run_root_of(config_folder) / LEDGER_NAME


def log_for(config_folder: Path) -> Path:
    return run_root_of(config_folder) / LOG_NAME


def find_run_script(run_folder: Path) -> Optional[Path]:
    candidate = run_root_of(run_folder) / "run_QMC_chiesa_force.py"
    return candidate if candidate.exists() else None


def find_nexus_out(run_folder: Path) -> Optional[Path]:
    parent_cfg = run_folder.parent.name
    candidate = run_folder.parents[1] / f"{parent_cfg}.out"
    return candidate if candidate.exists() else None


def iter_config_folders(aurora_root: Path) -> List[Path]:
    out: List[Path] = []
    for run_dir in sorted(aurora_root.glob("run_*/runs")):
        for cfg_root in sorted(p for p in run_dir.iterdir() if p.is_dir()):
            for config in sorted(p for p in cfg_root.iterdir() if p.is_dir()):
                if CONFIG_FOLDER_RE.match(config.name):
                    out.append(config)
    return out


def gather_status_per_config(
    configs: List[Path],
    explicit_ledger: Optional[Path] = None,
) -> Dict[str, str]:
    """Read either one explicit ledger or one ledger per run_root."""
    if explicit_ledger is not None:
        return load_ledger(explicit_ledger)
    cache: Dict[Path, Dict[str, str]] = {}
    result: Dict[str, str] = {}
    for cfg in configs:
        rl = ledger_for(cfg)
        if rl not in cache:
            cache[rl] = load_ledger(rl)
        if str(cfg) in cache[rl]:
            result[str(cfg)] = cache[rl][str(cfg)]
    return result


def validate_h5(path: Path) -> Tuple[bool, str]:
    try:
        with h5py.File(path, "r") as f:
            missing_attrs = [a for a in REQUIRED_ROOT_ATTRS if a not in f.attrs]
            missing_groups = [g for g in REQUIRED_GROUPS if g not in f]
            if missing_attrs or missing_groups:
                return False, f"missing_attrs={missing_attrs} missing_groups={missing_groups}"

            method = str(f.attrs.get("method", "")).upper()
            obs = f["observables"]
            if method == "DFT":
                if "total_energy" not in obs:
                    return False, "dft observables missing total_energy"
            elif method in {"VMC", "DMC"}:
                for key in ("total_energy", "kinetic_energy", "potential_energy", "forces"):
                    if key not in obs:
                        return False, f"qmc observables missing {key}"
                if "finite_size_corrected" not in obs.attrs:
                    return False, "qmc observables.attrs missing finite_size_corrected"
                if "qmc_quality" not in f.attrs:
                    return False, "qmc root attrs missing qmc_quality"

            prov = f["provenance"]
            for key in ("uuid_in_system", "dependent_uuids"):
                if key not in prov:
                    return False, f"provenance missing {key}"
    except Exception as exc:  # noqa: BLE001
        return False, f"open_error={exc!r}"
    return True, "ok"


def upload_one(
    local_path: Path,
    s3_bucket: str,
    endpoint: str,
    dry_run: bool,
    log_path: Path,
) -> bool:
    target = f"{s3_bucket.rstrip('/')}/{local_path.name}"
    cmd = [
        sys.executable, "-m", "awscli",
        "--endpoint-url", endpoint,
        "s3", "cp", str(local_path), target, "--no-progress",
    ]
    if dry_run:
        cmd.append("--dryrun")
    log(f"upload cmd: {' '.join(shlex.quote(c) for c in cmd)}", log_path)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        log(f"upload TIMEOUT: {local_path.name}", log_path)
        return False
    if proc.returncode != 0:
        log(f"upload FAIL rc={proc.returncode}: {proc.stderr.strip()[:500]}", log_path)
        return False
    if proc.stdout.strip():
        log(f"upload stdout: {proc.stdout.strip()[:500]}", log_path)
    return True


def process_one(
    run_folder: Path,
    out_dir: Path,
    log_path: Path,
    args: argparse.Namespace,
) -> Tuple[str, str, List[Path]]:
    run_script = args.run_script or find_run_script(run_folder)
    if run_script is None:
        return "skipped", "no run_QMC_chiesa_force.py", []
    nexus_out = args.nexus_out or find_nexus_out(run_folder)

    try:
        outputs = build_run_files(
            run_folder=run_folder,
            out_dir=out_dir,
            run_script=run_script,
            nexus_out=nexus_out,
        )
    except FileNotFoundError as exc:
        return "skipped", f"build skipped: {exc}", []
    except Exception:  # noqa: BLE001
        return "error", "build failed: " + traceback.format_exc(limit=2).strip(), []

    produced = list(outputs.values())
    invalid: List[str] = []
    for path in produced:
        ok, reason = validate_h5(path)
        if not ok:
            invalid.append(f"{path.name}: {reason}")
    if invalid:
        return "invalid", "; ".join(invalid), produced

    if args.upload:
        for path in produced:
            if not upload_one(path, args.s3_bucket, args.endpoint, args.dry_run_upload, log_path):
                return "upload_failed", path.name, produced

    return "ok", ",".join(p.name for p in produced), produced


def filter_pending(configs: List[Path], status_map: Dict[str, str], retry_failed: bool) -> List[Path]:
    pending: List[Path] = []
    for cfg in configs:
        status = status_map.get(str(cfg))
        if status == "ok":
            continue
        if status in RETRYABLE and not retry_failed:
            continue
        pending.append(cfg)
    return pending


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--aurora-root", type=Path, default=DEFAULT_AURORA_ROOT)
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Override HDF5 output directory. Default: write into each configuration's own folder.")
    p.add_argument("--ledger", type=Path, default=None,
                   help="Override ledger file. Default: <run_root>/" + LEDGER_NAME)
    p.add_argument("--log", type=Path, default=None,
                   help="Override log file. Default: <run_root>/" + LOG_NAME)
    p.add_argument("--limit", type=int, default=0, help="If > 0, process at most N configurations (orchestration mode).")
    p.add_argument("--run-script", type=Path, default=None)
    p.add_argument("--nexus-out", type=Path, default=None)
    p.add_argument("--include-runs", type=str, default="", help="Comma list of run_YYYY_MM_DD names to include (empty = all).")
    p.add_argument("--retry-failed", action="store_true", help="Re-process entries marked error/invalid/upload_failed in the ledger.")
    p.add_argument("--upload", action="store_true", help="Upload validated HDF5 files to S3/OSN.")
    p.add_argument("--dry-run-upload", action="store_true", help="Dry-run uploads only.")
    p.add_argument("--s3-bucket", type=str, default=os.environ.get("S3_BUCKET", "s3://phy240060/QMCHAMM"))
    p.add_argument("--endpoint", type=str, default=os.environ.get("OSN_ENDPOINT", "https://uri.osn.mghpcc.org"))
    p.add_argument("--single", type=Path, default=None, help="Process one configuration folder and exit.")
    p.add_argument("--list-only", action="store_true", help="Print pending configuration folders, one per line, and exit.")
    return p.parse_args()


def _validate_upload_env(log_path: Path, upload: bool) -> bool:
    if not upload:
        return True
    if not (os.environ.get("AWS_ACCESS_KEY_ID") and os.environ.get("AWS_SECRET_ACCESS_KEY")):
        log("--upload requested but AWS_ACCESS_KEY_ID/AWS_SECRET_ACCESS_KEY are not set", log_path)
        return False
    return True


def cmd_single(args: argparse.Namespace) -> int:
    cfg = args.single.resolve()
    if not cfg.exists():
        # Fall back to a stderr message since we may not yet have a log path.
        print(f"--single path missing: {cfg}", file=sys.stderr)
        return 2
    out_dir = args.out_dir if args.out_dir is not None else cfg
    log_path = args.log if args.log is not None else log_for(cfg)
    ledger_path = args.ledger if args.ledger is not None else ledger_for(cfg)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not _validate_upload_env(log_path, args.upload):
        return 2

    log(f"--> {cfg}", log_path)
    started = time.time()
    try:
        status, message, _ = process_one(cfg, out_dir, log_path, args)
    except KeyboardInterrupt:
        log(f"interrupted while processing {cfg}", log_path)
        return 130
    elapsed = time.time() - started
    log(f"<-- {cfg} :: {status} :: {message} :: {elapsed:.2f}s", log_path)
    append_ledger(ledger_path, str(cfg), status, message)
    return 0 if status in {"ok", "skipped"} else 1


def cmd_list(args: argparse.Namespace) -> int:
    include_runs = {s.strip() for s in args.include_runs.split(",") if s.strip()}
    configs = iter_config_folders(args.aurora_root)
    if include_runs:
        configs = [c for c in configs if run_root_of(c).name in include_runs]
    status_map = gather_status_per_config(configs, args.ledger)
    pending = filter_pending(configs, status_map, args.retry_failed)
    if args.limit and len(pending) > args.limit:
        pending = pending[: args.limit]
    for cfg in pending:
        print(cfg)
    return 0


def cmd_run_serial(args: argparse.Namespace) -> int:
    include_runs = {s.strip() for s in args.include_runs.split(",") if s.strip()}
    configs = iter_config_folders(args.aurora_root)
    if include_runs:
        configs = [c for c in configs if run_root_of(c).name in include_runs]

    status_map = gather_status_per_config(configs, args.ledger)
    pending = filter_pending(configs, status_map, args.retry_failed)
    if args.limit and len(pending) > args.limit:
        pending = pending[: args.limit]

    # Group counters for the final SUMMARY, written to each touched log.
    n_ok = n_skip = n_err = 0
    started = time.time()
    touched_logs: set[Path] = set()

    for cfg in pending:
        out_dir = args.out_dir if args.out_dir is not None else cfg
        log_path = args.log if args.log is not None else log_for(cfg)
        ledger_path = args.ledger if args.ledger is not None else ledger_for(cfg)
        out_dir.mkdir(parents=True, exist_ok=True)
        touched_logs.add(log_path)

        if not _validate_upload_env(log_path, args.upload):
            return 2

        log(f"--> {cfg}", log_path)
        try:
            status, message, _ = process_one(cfg, out_dir, log_path, args)
        except KeyboardInterrupt:
            log("interrupted by user; exiting", log_path)
            return 130
        log(f"<-- {cfg} :: {status} :: {message}", log_path)
        append_ledger(ledger_path, str(cfg), status, message)
        if status == "ok":
            n_ok += 1
        elif status == "skipped":
            n_skip += 1
        else:
            n_err += 1

    elapsed = time.time() - started
    summary = {
        "processed": len(pending), "ok": n_ok, "skipped": n_skip,
        "errors": n_err, "elapsed_seconds": round(elapsed, 2),
    }
    for lp in touched_logs:
        log(f"SUMMARY {json.dumps(summary)}", lp)
    return 0 if n_err == 0 else 1


def main() -> int:
    args = parse_args()
    if args.single is not None and args.list_only:
        print("--single and --list-only are mutually exclusive", file=sys.stderr)
        return 2
    if args.single is not None:
        return cmd_single(args)
    if args.list_only:
        return cmd_list(args)
    return cmd_run_serial(args)


if __name__ == "__main__":
    sys.exit(main())
