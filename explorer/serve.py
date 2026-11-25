"""Local read-only web server for the QMC-HAMM explorer (127.0.0.1 only)."""
import argparse
import gzip
import json
import math
import os
import sys
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import h5py  # noqa: E402
import numpy as np  # noqa: E402

from explorer import export  # noqa: E402

HERE = Path(__file__).resolve().parent
STATIC = HERE / "static"
LINK_FIELDS = ["config_id", "linked_dft_uuids", "linked_qmc_uuids", "missing_links", "geom_mismatch"]
DATA_FILES = {"index.json.gz": None, "eos_residuals.json": "application/json",
              "harvest_report.json": "application/json", "eos_fit.json": "application/json"}
TYPES = {".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json",
         ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon"}


def _py(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        return _py(v.tolist())
    if isinstance(v, np.generic):
        return _py(v.item())
    if isinstance(v, (list, tuple)):
        return [_py(x) for x in v]
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def _attrs(g):
    return {k: _py(v) for k, v in g.attrs.items()}


def _scalars(group):
    out = {}
    for k, d in group.items():
        if isinstance(d, h5py.Dataset) and d.size == 1:
            out[k] = _py(np.ravel(d[()])[0] if d.ndim else d[()])
    return out


def _observables(g):
    out = {}
    for k, d in g.items():
        if not isinstance(d, h5py.Dataset):
            continue
        v = d[()]
        if k.startswith("structure_factor"):
            out[k] = _py(np.asarray(v))
        elif d.size == 1:
            out[k] = _py(np.ravel(v)[0] if d.ndim else v)
        else:
            s = {"shape": list(d.shape), "min": None, "max": None}
            fin = np.asarray(v, float)[np.isfinite(v)] if d.dtype.kind in "fiu" else []
            if len(fin):
                s["min"], s["max"] = float(fin.min()), float(fin.max())
            out[k] = s
    return out


def file_detail(path, row):
    with h5py.File(path, "r") as f:
        d = {"uuid": row["uuid"], "attrs": _attrs(f)}
        for name in ("observables", "parameters", "provenance"):
            g = f.get(name)
            if g is None:
                d[name] = {}
            elif name == "observables":
                d[name] = _observables(g)
            else:
                d[name] = {"attrs": _attrs(g), "values": _scalars(g)}
    d["links"] = {k: row.get(k) for k in LINK_FIELDS}
    return d


def load_index(path):
    with gzip.open(path, "rt") as fh:
        idx = json.load(fh)
    cols, n = idx["columns"], len(idx["columns"]["uuid"])
    return [{f: cols[f][i] for f in idx["fields"]} for i in range(n)]


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def _send(self, code, body=b"", ctype="text/plain", extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, allow_nan=False).encode(), "application/json")

    def _file(self, path, ctype, extra=None):
        size = path.stat().st_size
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(size))
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        with open(path, "rb") as fh:
            while chunk := fh.read(1 << 20):
                self.wfile.write(chunk)

    def row_path(self, row):
        k = row["key"]
        return Path(k[6:]) if k.startswith("local:") else Path(self.server.cache) / Path(k).name

    def do_GET(self):
        p = urllib.parse.unquote(urllib.parse.urlsplit(self.path).path)
        srv = self.server
        try:
            if p == "/":
                return self._static("index.html")
            if p.startswith("/static/"):
                return self._static(p[8:])
            if p.startswith("/data/") and p[6:] in DATA_FILES:
                f = Path(srv.data) / p[6:]
                if not f.is_file():
                    return self._send(404, b"not found")
                if p.endswith(".gz"):
                    return self._file(f, "application/json", {"Content-Encoding": "gzip"})
                return self._file(f, "application/json")
            if p.startswith("/api/file/"):
                row = srv.by_uuid.get(p[10:])
                if row is None or not self.row_path(row).exists():
                    return self._send(404, b"unknown uuid")
                return self._json(file_detail(self.row_path(row), row))
            if p.startswith("/h5/") and p.endswith(".h5"):
                row = srv.by_uuid.get(p[4:-3])
                if row is None or not self.row_path(row).exists():
                    return self._send(404, b"unknown uuid")
                return self._file(self.row_path(row), "application/x-hdf5",
                                  {"Content-Disposition": f'attachment; filename="{row["uuid"]}.h5"'})
        except (BrokenPipeError, ConnectionResetError):
            return
        except Exception:
            self.close_connection = True
            return self._send(500, b"internal error")
        self._send(404, b"not found")

    def _static(self, rel):
        f = (STATIC / rel).resolve()
        if not f.is_relative_to(STATIC.resolve()) or not f.is_file():
            return self._send(404, b"not found")
        self._file(f, TYPES.get(f.suffix, "application/octet-stream"))

    def do_POST(self):
        if urllib.parse.urlsplit(self.path).path != "/api/export.xyz":
            return self._send(404, b"not found")
        try:
            n = int(self.headers.get("Content-Length") or "")
            if n < 0:
                raise ValueError
            body = self.rfile.read(n).decode()
            payload = json.loads(urllib.parse.parse_qs(body)["payload"][0])
            uuids = [str(u) for u in payload["uuids"]]
        except (KeyError, ValueError, TypeError, IndexError, AttributeError):
            self.close_connection = True
            return self._send(400, b"bad payload")
        by = self.server.by_uuid
        rows = [by[u] for u in dict.fromkeys(uuids) if u in by]
        if not rows:
            return self._send(404, b"no known uuids")
        summary = {"unknown_uuids": len(set(uuids)) - len(rows), "missing_links": 0}
        if payload.get("include_linked_dft"):
            summary["missing_links"] = (
                sum(u not in by for r in rows for u in r.get("linked_dft_uuids") or [])
                + sum(len(r.get("missing_links") or []) for r in rows))
            rows = export.expand_with_linked_dft(rows, by)
        summary["suspect_skipped"] = 0
        if not payload.get("include_suspect"):
            keep = [r for r in rows if r.get("rs_suspect") is not True]
            summary["suspect_skipped"] = len(rows) - len(keep)
            rows = keep
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Disposition", f'attachment; filename="qmchamm_H_{len(rows)}frames.xyz"')
        self.send_header("X-Export-Summary", json.dumps(
            {"planned_frames": len(rows), "unknown_uuids": summary["unknown_uuids"],
             "missing_links": summary["missing_links"],
             "suspect_skipped": summary["suspect_skipped"]}, separators=(",", ":")))
        self.send_header("Transfer-Encoding", "chunked")
        self.end_headers()

        def chunk(s):
            b = s.encode()
            if b:
                self.wfile.write(b"%X\r\n%s\r\n" % (len(b), b))

        try:
            for frame in export.iter_xyz(rows, self.server.cache, summary):
                chunk(frame)
            self.wfile.write(b"0\r\n\r\n")
            self.log_message("export summary %s", json.dumps(summary))
        except (BrokenPipeError, ConnectionResetError):
            self.close_connection = True


def make_server(port=8765, data=None, cache=None):
    data = Path(data or HERE / "data")
    cache = Path(cache or os.environ.get("QMCDB_CACHE", str(HERE / "qmc-db-cache" / "h5")))
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    srv.daemon_threads = True
    srv.data, srv.cache = data, cache
    idx = data / "index.json.gz"
    srv.rows = load_index(idx) if idx.exists() else []
    srv.by_uuid = {r["uuid"]: r for r in srv.rows}
    return srv


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8771)
    ap.add_argument("--data", default=str(HERE / "data"))
    ap.add_argument("--cache", default=os.environ.get("QMCDB_CACHE", str(HERE / "qmc-db-cache" / "h5")))
    a = ap.parse_args(argv)
    srv = make_server(a.port, a.data, a.cache)
    print(f"http://127.0.0.1:{srv.server_address[1]}/", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
