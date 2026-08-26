#!/usr/bin/env python3
"""Push the schema HDF5 database to OSN, uploading each file at most once.

The bucket is the ground truth: it is listed once at the start of a run and the
upload set is the difference between the index and that listing. On top of that
a persistent registry keyed by content hash catches the case where the same
physical file was written twice under different uuids, so identical bytes are
never sent twice even from different source trees.

Typical use:
  python3 sync_database_to_osn.py --dry-run       # show what would go up
  python3 sync_database_to_osn.py                 # do it
  python3 sync_database_to_osn.py --verify        # re-check sizes in bucket

Credentials come from the environment (AWS_ACCESS_KEY_ID /
AWS_SECRET_ACCESS_KEY); use submit_database_sync.sh to be prompted for them.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DEFAULT_ROOT = Path("/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup")
DEFAULT_BUCKET = "s3://phy240060/QMCHAMM"
DEFAULT_ENDPOINT = "https://uri.osn.mghpcc.org"

# awscli need not live in the same interpreter that runs this script, so it is
# invoked out-of-process through whichever python actually provides it.
AWS_PYTHON = os.environ.get("AWS_CLI_PYTHON") or shutil.which("python3") or sys.executable

REGISTRY_FIELDS = [
    "sha256", "uuid", "s3_key", "size_bytes", "uploaded_at",
    "source_dataset", "method", "local_path",
]


def aws(endpoint: str, *args: str, timeout: int = 900,
        unsigned: bool = False) -> Tuple[int, str, str]:
    cmd = [AWS_PYTHON, "-m", "awscli", "--endpoint-url", endpoint]
    if unsigned:
        cmd.append("--no-sign-request")
    cmd += ["s3", *args]
    env = dict(os.environ)
    env.setdefault("AWS_DEFAULT_REGION", "us-east-1")
    env.setdefault("AWS_REQUEST_CHECKSUM_CALCULATION", "when_required")
    env.setdefault("AWS_RESPONSE_CHECKSUM_VALIDATION", "when_required")
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    return p.returncode, p.stdout, p.stderr


def list_bucket(bucket: str, endpoint: str, cache: Path,
                refresh: bool, unsigned: bool = False) -> Dict[str, int]:
    """Map of object key basename -> size. Cached so repeat runs are cheap."""
    if cache.exists() and not refresh:
        out: Dict[str, int] = {}
        with open(cache, newline="") as fh:
            for key, size in csv.reader(fh, delimiter="\t"):
                out[key] = int(size)
        print(f"remote listing (cached {cache.name}): {len(out)} objects")
        return out

    print(f"listing {bucket} ...", flush=True)
    rc, out, err = aws(endpoint, "ls", bucket.rstrip("/") + "/", "--recursive",
                       timeout=3600, unsigned=unsigned)
    if rc != 0:
        raise SystemExit(f"could not list bucket:\n{err.strip()}")

    prefix = bucket.split("//", 1)[1].split("/", 1)
    sub = (prefix[1].rstrip("/") + "/") if len(prefix) > 1 else ""
    remote: Dict[str, int] = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        size, key = parts[2], parts[3]
        if sub and key.startswith(sub):
            key = key[len(sub):]
        if key:
            remote[key] = int(size)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for k, v in sorted(remote.items()):
            w.writerow([k, v])
    print(f"remote listing: {len(remote)} objects  -> {cache}")
    return remote


def load_registry(path: Path) -> Dict[str, dict]:
    if not path.exists():
        return {}
    out: Dict[str, dict] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("sha256"):
                out[row["sha256"]] = row
    return out


def append_registry(path: Path, rows: List[dict]) -> None:
    if not rows:
        return
    new = not path.exists()
    with open(path, "a", newline="") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        w = csv.DictWriter(fh, fieldnames=REGISTRY_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)
        fh.flush()
        os.fsync(fh.fileno())
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def load_index(path: Path) -> List[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    missing = [r for r in rows if not r.get("sha256")]
    if missing:
        raise SystemExit(
            f"{len(missing)} index rows have no sha256; rebuild the index with "
            "build_database_index.py --hash before syncing"
        )
    return rows


def plan(rows: List[dict], remote: Dict[str, int],
         registry: Dict[str, dict]) -> Tuple[List[dict], dict]:
    """Split the index into what still needs uploading and why the rest doesn't."""
    todo: List[dict] = []
    stats = {"in_bucket": 0, "size_mismatch": 0, "dup_content": 0,
             "missing_local": 0, "to_upload": 0}
    seen_hash: Dict[str, str] = {}

    for r in rows:
        key = f"{r['uuid']}.h5"
        size = int(r["size_bytes"])
        sha = r["sha256"]

        if not Path(r["path"]).exists():
            stats["missing_local"] += 1
            continue

        if key in remote:
            if remote[key] == size:
                stats["in_bucket"] += 1
                continue
            stats["size_mismatch"] += 1
            todo.append(r)
            continue

        # Same bytes already up under another key, or duplicated inside the
        # index itself: leave the existing object as the single copy.
        prior = registry.get(sha)
        if prior and prior.get("s3_key") in remote:
            stats["dup_content"] += 1
            continue
        if sha in seen_hash:
            stats["dup_content"] += 1
            continue

        seen_hash[sha] = r["uuid"]
        todo.append(r)

    stats["to_upload"] = len(todo)
    return todo, stats


def upload_one(row: dict, bucket: str, endpoint: str,
               retries: int = 3) -> Tuple[dict, Optional[str]]:
    key = f"{row['uuid']}.h5"
    dest = f"{bucket.rstrip('/')}/{key}"
    for attempt in range(1, retries + 1):
        rc, _, err = aws(endpoint, "cp", row["path"], dest, "--no-progress",
                         timeout=600)
        if rc == 0:
            return {
                "sha256": row["sha256"],
                "uuid": row["uuid"],
                "s3_key": key,
                "size_bytes": row["size_bytes"],
                "uploaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "source_dataset": row.get("source_dataset", ""),
                "method": row.get("method", ""),
                "local_path": row["path"],
            }, None
        if attempt < retries:
            time.sleep(2 * attempt)
    return {}, f"{key}: {err.strip().splitlines()[-1] if err.strip() else 'failed'}"


def upload_batch(todo: List[dict], bucket: str, endpoint: str, staging: Path,
                 log: Path) -> Tuple[List[dict], List[str]]:
    """Upload via one recursive copy over a flat tree of hardlinks.

    The bucket is flat and named by uuid while the database is laid out in
    per-condition directories, so the upload set is first mirrored as hardlinks
    (no extra bytes, same filesystem). One awscli invocation then transfers
    everything, which avoids paying process startup 54k times.
    """
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    by_key: Dict[str, dict] = {}
    for r in todo:
        key = f"{r['uuid']}.h5"
        dst = staging / key
        try:
            os.link(r["path"], dst)
        except OSError:
            shutil.copy2(r["path"], dst)
        by_key[key] = r
    print(f"staged {len(by_key)} files under {staging}", flush=True)

    rc, out, err = aws(endpoint, "cp", str(staging) + "/",
                       bucket.rstrip("/") + "/", "--recursive", "--no-progress",
                       timeout=48 * 3600)
    log.write_text(out + ("\n" + err if err else ""))

    done: List[dict] = []
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for line in out.splitlines():
        if not line.startswith("upload:"):
            continue
        key = line.rsplit("/", 1)[-1].strip()
        r = by_key.pop(key, None)
        if r is None:
            continue
        done.append({
            "sha256": r["sha256"], "uuid": r["uuid"], "s3_key": key,
            "size_bytes": r["size_bytes"], "uploaded_at": now,
            "source_dataset": r.get("source_dataset", ""),
            "method": r.get("method", ""), "local_path": r["path"],
        })

    shutil.rmtree(staging, ignore_errors=True)
    errors = [f"{k}: not confirmed uploaded" for k in sorted(by_key)]
    if rc != 0 and not errors:
        errors.append(f"awscli exited {rc}; see {log}")
    return done, errors


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--index", type=Path, default=None,
                    help="defaults to <root>/database_index.csv")
    ap.add_argument("--registry", type=Path, default=None,
                    help="defaults to <root>/upload_registry.csv")
    ap.add_argument("--bucket", default=DEFAULT_BUCKET)
    ap.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--refresh-remote", action="store_true",
                    help="re-list the bucket instead of using the cached listing")
    ap.add_argument("--unsigned-list", action="store_true",
                    help="list the bucket anonymously (it is public-read); lets "
                         "you plan a sync without write credentials")
    ap.add_argument("--only-source", default=None,
                    help="restrict to one source_dataset, e.g. ricky_legacy")
    ap.add_argument("--limit", type=int, default=0, help="upload at most N files")
    ap.add_argument("--per-file", action="store_true",
                    help="one awscli call per file instead of the batched "
                         "hardlink-staging copy; slower but retries per file")
    args = ap.parse_args()

    index = args.index or (args.root / "database_index.csv")
    registry_path = args.registry or (args.root / "upload_registry.csv")
    cache = args.root / "osn_remote_listing.tsv"

    probe = subprocess.run([AWS_PYTHON, "-m", "awscli", "--version"],
                           capture_output=True, text=True)
    if probe.returncode != 0:
        raise SystemExit(
            f"awscli is not importable from {AWS_PYTHON}; set AWS_CLI_PYTHON to "
            "an interpreter that has it (e.g. /usr/bin/python3)"
        )
    print(f"awscli via {AWS_PYTHON}: {(probe.stdout or probe.stderr).strip()}")

    rows = load_index(index)
    if args.only_source:
        rows = [r for r in rows if r.get("source_dataset") == args.only_source]
    print(f"index: {len(rows)} files, "
          f"{sum(int(r['size_bytes']) for r in rows) / 1e9:.2f} GB")

    remote = list_bucket(args.bucket, args.endpoint, cache, args.refresh_remote,
                         unsigned=args.unsigned_list)
    registry = load_registry(registry_path)
    print(f"registry: {len(registry)} content hashes already recorded")

    todo, stats = plan(rows, remote, registry)
    print(json.dumps(stats, indent=2))

    if args.limit:
        todo = todo[: args.limit]
    if not todo:
        print("nothing to upload; bucket is in sync with the index")
        return 0
    if args.dry_run:
        print(f"DRY RUN: would upload {len(todo)} files "
              f"({sum(int(r['size_bytes']) for r in todo) / 1e9:.2f} GB)")
        for r in todo[:10]:
            print(f"  {r['path']} -> {args.bucket.rstrip('/')}/{r['uuid']}.h5")
        if len(todo) > 10:
            print(f"  ... and {len(todo) - 10} more")
        return 0

    t0 = time.time()
    if not args.per_file:
        done, errors = upload_batch(
            todo, args.bucket, args.endpoint,
            args.root / ".upload_staging",
            args.root / f"osn_upload_awscli_{datetime.now():%Y%m%d_%H%M%S}.log",
        )
        append_registry(registry_path, done)
        print(f"uploaded {len(done)}/{len(todo)} in {time.time() - t0:.0f}s; "
              f"registry -> {registry_path}")
        if errors:
            errlog = args.root / "osn_upload_errors.txt"
            errlog.write_text("\n".join(errors))
            print(f"{len(errors)} not confirmed; listed in {errlog}")
            return 1
        return 0

    print(f"uploading {len(todo)} files with {args.workers} workers", flush=True)
    done: List[dict] = []
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(upload_one, r, args.bucket, args.endpoint): r for r in todo}
        for n, fut in enumerate(as_completed(futs), 1):
            rec, err = fut.result()
            if err:
                errors.append(err)
            else:
                done.append(rec)
            if n % 500 == 0 or n == len(todo):
                append_registry(registry_path, done)
                done = []
                rate = n / max(time.time() - t0, 1e-9)
                print(f"  {n}/{len(todo)}  {rate:.1f} files/s  "
                      f"errors={len(errors)}", flush=True)
    append_registry(registry_path, done)

    print(f"uploaded {len(todo) - len(errors)}/{len(todo)} in "
          f"{time.time() - t0:.0f}s; registry -> {registry_path}")
    if errors:
        errlog = args.root / "osn_upload_errors.txt"
        errlog.write_text("\n".join(errors))
        print(f"{len(errors)} failures listed in {errlog}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
