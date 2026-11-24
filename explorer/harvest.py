"""Mirror the OSN bucket into a local cache and build the explorer index."""
import argparse
import collections
import concurrent.futures as cf
import datetime
import gzip
import json
import math
import multiprocessing
import os
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from explorer import bucket, eos as eosmod, links, normalize  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_LOCAL = HERE.parent / "Li7P3S11_QMC_DFT_dataset"
EXTRA_FIELDS = ["config_id", "linked_dft_uuids", "linked_qmc_uuids", "missing_links",
                "geom_mismatch", "P_used", "P_estimated", "P_extrapolated"]
FIELDS = normalize.ROW_FIELDS + EXTRA_FIELDS


def _clean(v):
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    return v


def _dumps(obj):
    return json.dumps(_clean(obj), allow_nan=False)


def _atomic(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def build_index(rows, eos="fit", harvested_at=None):
    return {"fields": FIELDS,
            "columns": {f: [r.get(f) for r in rows] for f in FIELDS},
            "harvested_at": harvested_at or datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "eos": eos, "normalize_version": normalize.NORMALIZE_VERSION}


def _load_rows(path):
    try:
        with gzip.open(path, "rt") as fh:
            idx = json.load(fh)
    except (OSError, ValueError):
        return {}, None
    cols = idx.get("columns", {})
    n = len(cols.get("key", []))
    return ({cols["key"][i]: {f: cols[f][i] if f in cols else None for f in normalize.ROW_FIELDS}
             for i in range(n)}, idx.get("normalize_version"))


def _local_objs(dirs):
    out = []
    for d in dirs:
        for p in sorted(Path(d).glob("*.h5")):
            st = p.stat()
            out.append({"key": "local:" + str(p.resolve()), "etag": f"{st.st_mtime_ns}-{st.st_size}",
                        "size": st.st_size, "last_modified": None})
    return out


def _path(key, cache):
    return Path(key[6:]) if key.startswith("local:") else Path(cache) / os.path.basename(key)


def _parse(args):
    path, key, etag, size = args
    return normalize.parse_file(path, key, etag, size)


def _counts(rows, field):
    return dict(collections.Counter(str(r.get(field)) for r in rows))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default=os.environ.get("QMCDB_CACHE", "qmc-db-cache/h5"))
    ap.add_argument("--data", default=str(HERE / "data"))
    ap.add_argument("--limit", type=int)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--local-dir", action="append", default=None)
    ap.add_argument("--no-download", action="store_true")
    a = ap.parse_args(argv)
    if a.limit is not None and Path(a.data).resolve() == (HERE / "data").resolve():
        ap.error("--limit would overwrite the real index in explorer/data; pass --data <scratch dir>")
    data, cache = Path(a.data), Path(a.cache)
    data.mkdir(parents=True, exist_ok=True)
    local_dirs = a.local_dir if a.local_dir is not None else [str(DEFAULT_LOCAL)]

    remote = sorted((o for o in bucket.list_bucket() if o["key"].endswith(".h5")),
                    key=lambda o: o["key"])
    if a.limit is not None:
        remote = remote[:a.limit]
    remote = remote + _local_objs(local_dirs)
    objs = {o["key"]: o for o in remote}
    try:
        manifest = json.loads((data / "manifest.json").read_text())
    except (OSError, ValueError):
        manifest = {}
    plan = bucket.plan_sync(remote, manifest)
    prev, prev_ver = _load_rows(data / "index.json.gz")
    reparse = prev_ver != normalize.NORMALIZE_VERSION
    if reparse:
        prev = {}  # parser changed: re-parse every cached file instead of reusing old rows

    errors_dl, not_cached = [], 0
    new_s, changed_s, unch_s = set(plan["new"]), set(plan["changed"]), set(plan["unchanged"])
    stale = {k for k in changed_s if a.no_download and not k.startswith("local:")}
    todo = [k for k in plan["new"] + plan["changed"] if k not in stale]
    todo += [k for k in plan["unchanged"] if k not in prev]
    if reparse:
        todo += [k for k in stale if _path(k, cache).exists()]
    need_dl = []
    for k in todo:
        if k.startswith("local:"):
            continue
        p = _path(k, cache)
        fresh = k in new_s or k not in changed_s
        if p.exists() and (fresh and p.stat().st_size == objs[k]["size"] or a.no_download):
            continue
        if a.no_download:
            not_cached += 1
        else:
            need_dl.append(k)
    failed = set()
    if need_dl:
        def dl(k):
            try:
                bucket.download(objs[k], cache)
            except Exception as e:
                return k, f"{type(e).__name__}: {e}"
            return k, None
        with cf.ThreadPoolExecutor(a.workers) as ex:
            for k, err in ex.map(dl, need_dl):
                if err:
                    failed.add(k)
                    errors_dl.append({"key": k, "error": err})
    skip = failed | {k for k in todo if not k.startswith("local:") and not _path(k, cache).exists()}

    jobs = [(str(_path(k, cache)), k, manifest[k]["etag"] if k in stale else objs[k]["etag"],
             manifest[k]["size"] if k in stale else objs[k]["size"]) for k in todo if k not in skip]
    if len(jobs) < 50:
        parsed = [_parse(j) for j in jobs]
    else:
        with cf.ProcessPoolExecutor(a.workers, mp_context=multiprocessing.get_context("spawn")) as ex:
            parsed = list(ex.map(_parse, jobs, chunksize=16))
    new_rows = {r["key"]: r for r in parsed}

    rows = []
    for o in remote:
        k = o["key"]
        r = new_rows.get(k) or (prev.get(k) if k in unch_s or k in stale else None)
        if r is not None:
            rows.append(r)
    rows.sort(key=lambda r: r["key"])
    manifest = {r["key"]: dict(manifest[r["key"]]) if r["key"] in stale
                else {"etag": objs[r["key"]]["etag"], "size": objs[r["key"]]["size"]} for r in rows}

    path_of = {r["uuid"]: _path(r["key"], cache) for r in rows}
    res, lsum = links.resolve(rows, lambda u: links.geom_from_file(path_of[u]))
    for r in rows:
        r.update(res[r["uuid"]])

    report = {"harvested_at": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        fit = eosmod.FittedEos.fit(rows)
    except ValueError as e:
        fit = None
        for r in rows:
            r["P_used"], r["P_estimated"], r["P_extrapolated"] = r["P_recorded"], False, False
        report["eos"] = {"name": "none", "reason": str(e)}
    if fit is not None:
        eosmod.apply_eos(rows, fit)
        report["eos"] = {"name": fit.name, "a": fit.a, "b": fit.b, "cv_rms": fit.cv_rms}
        _atomic(data / "eos_fit.json", _dumps(fit.to_json()).encode())
        _atomic(data / "eos_residuals.json",
                _dumps(fit.residuals(rows).to_dict("records")).encode())

    idx = build_index(rows, "fit" if fit else "none", report["harvested_at"])
    tmp = data / "index.json.gz.tmp"
    with gzip.open(tmp, "wt") as fh:
        fh.write(_dumps(idx))
    os.replace(tmp, data / "index.json.gz")
    _atomic(data / "manifest.json", json.dumps(manifest).encode())

    report.update(
        total=len(rows), by_tab=_counts(rows, "tab"), by_source=_counts(rows, "source"),
        by_method=_counts(rows, "method"),
        new=len(plan["new"]), changed=len(plan["changed"]), deleted=len(plan["deleted"]),
        unchanged=len(plan["unchanged"]), parsed=len(parsed),
        parse_errors=[{"key": r["key"], "error": r["parse_error"]} for r in rows if r.get("parse_error")],
        download_errors=errors_dl, not_cached=not_cached, stale=len(stale), links=lsum,
        rs_suspect=sum(1 for r in rows if r.get("tab") == "hydrogen" and r.get("rs_suspect") is True))
    _atomic(data / "harvest_report.json", _dumps(report).encode())
    print(f"harvest: {len(rows)} rows, parsed {len(parsed)}, new {len(plan['new'])}, "
          f"changed {len(plan['changed'])}, deleted {len(plan['deleted'])}, "
          f"parse_errors {len(report['parse_errors'])}, download_errors {len(errors_dl)}, "
          f"not_cached {not_cached}, eos {report['eos']['name']}, rs_suspect {report['rs_suspect']}")
    return report


if __name__ == "__main__":
    main()
