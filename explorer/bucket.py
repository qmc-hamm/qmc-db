"""List and mirror the public OSN bucket (unsigned GETs, read-only)."""
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = "https://uri.osn.mghpcc.org/phy240060"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def _get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def parse_listing(xml):
    root = ET.fromstring(xml)
    objs = []
    for c in root.iter(NS + "Contents"):
        objs.append({
            "key": c.findtext(NS + "Key"),
            "etag": (c.findtext(NS + "ETag") or "").strip('"'),
            "size": int(c.findtext(NS + "Size")),
            "last_modified": c.findtext(NS + "LastModified"),
        })
    tok = None
    if (root.findtext(NS + "IsTruncated") or "").lower() == "true":
        tok = root.findtext(NS + "NextContinuationToken")
        if not tok:
            raise RuntimeError("listing truncated but no continuation token")
    return objs, tok


def list_bucket(base=BASE, prefix="QMCHAMM/"):
    out, tok, seen = [], None, set()
    while True:
        q = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if tok:
            q["continuation-token"] = tok
        objs, tok = parse_listing(_get(base + "?" + urllib.parse.urlencode(q)).decode())
        out += objs
        if not tok:
            return out
        if tok in seen:
            raise RuntimeError(f"repeated continuation token {tok}")
        seen.add(tok)


def plan_sync(remote, manifest):
    rem = {o["key"]: o["etag"] for o in remote}
    plan = {"new": [], "changed": [], "deleted": [], "unchanged": []}
    for k, e in rem.items():
        if k not in manifest:
            plan["new"].append(k)
        elif manifest[k].get("etag") != e:
            plan["changed"].append(k)
        else:
            plan["unchanged"].append(k)
    plan["deleted"] = [k for k in manifest if k not in rem]
    return plan


def download(obj, cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / os.path.basename(obj["key"])
    part = dest.with_name(dest.name + ".part")
    data = _get(BASE + "/" + urllib.parse.quote(obj["key"]))
    part.write_bytes(data)
    if len(data) != obj["size"]:
        part.unlink()
        raise IOError(f"size mismatch for {obj['key']}: {len(data)} != {obj['size']}")
    os.replace(part, dest)
    return dest
