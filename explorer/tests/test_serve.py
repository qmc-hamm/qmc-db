import gzip
import inspect
import io
import json
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import ase.io
import pytest

from explorer import bucket, export, harvest, serve

FIXD = Path(__file__).parent / "fixtures"
LI = Path(__file__).resolve().parents[2] / "Li7P3S11_QMC_DFT_dataset"
RD = "f3beda31a0284b3ab3a64e54432960ce"


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("data")
    remote = [{"key": "QMCHAMM/" + p.name, "etag": "e", "size": p.stat().st_size,
               "last_modified": "x"} for p in sorted(FIXD.glob("*.h5"))]
    mp = pytest.MonkeyPatch()
    mp.setattr(bucket, "list_bucket", lambda *a, **k: remote)
    harvest.main(["--no-download", "--cache", str(FIXD), "--data", str(tmp),
                  "--local-dir", str(LI)])
    mp.undo()
    srv = serve.make_server(port=0, data=tmp, cache=FIXD)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv
    srv.shutdown()
    srv.server_close()


def url(srv, path):
    return f"http://127.0.0.1:{srv.server_address[1]}{path}"


def post(srv, payload):
    body = urllib.parse.urlencode({"payload": json.dumps(payload)}).encode()
    return urllib.request.urlopen(urllib.request.Request(url(srv, "/api/export.xyz"), body))


def test_serve_routes(server):
    assert server.server_address[0] == "127.0.0.1"
    d = json.load(urllib.request.urlopen(url(server, "/api/file/" + RD)))
    assert "method" in d["attrs"] and "observables" in d and "config_id" in d["links"]
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(url(server, "/api/file/nope"))
    assert e.value.code == 404
    r = post(server, {"uuids": [RD, "bogus"], "include_linked_dft": True})
    assert "attachment" in r.headers["Content-Disposition"]
    assert r.headers["Transfer-Encoding"] == "chunked"
    text = r.read().decode()
    assert len(ase.io.read(io.StringIO(text), ":", format="extxyz")) == 3
    hdr = json.loads(r.headers["X-Export-Summary"])
    assert hdr["planned_frames"] == 3 and hdr["unknown_uuids"] == 1
    with pytest.raises(urllib.error.HTTPError) as e:
        post(server, {"uuids": ["bogus"]})
    assert e.value.code == 404


def test_other_routes(server):
    r = urllib.request.urlopen(url(server, "/h5/" + RD + ".h5"))
    assert r.headers["Content-Type"] == "application/x-hdf5" and r.read(4) == b"\x89HDF"
    r = urllib.request.urlopen(url(server, "/data/index.json.gz"))
    assert r.headers["Content-Encoding"] == "gzip"
    assert "fields" in json.loads(gzip.decompress(r.read()))
    for p in ("/static/../serve.py", "/static/nope.js", "/h5/nope.h5"):
        with pytest.raises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(url(server, p))
        assert e.value.code == 404


def test_export_streams():
    assert inspect.isgeneratorfunction(export.iter_xyz)
    assert inspect.isgenerator(export.iter_xyz([], FIXD, {}))
    src = inspect.getsource(serve.Handler.do_POST)
    assert "join" not in src


def test_bad_post_is_400(server):
    for body in (b"payload=notjson", b"x=1", b"payload=%FF"):
        req = urllib.request.Request(url(server, "/api/export.xyz"), body)
        with pytest.raises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(req)
        assert e.value.code == 400


def test_suspect_skipped_by_default(server):
    SUS = "aeafdf12efae4e76ba34a5f5b94957c9"
    r = post(server, {"uuids": [RD, SUS]})
    hdr = json.loads(r.headers["X-Export-Summary"])
    assert hdr["suspect_skipped"] == 1 and hdr["planned_frames"] == 1
    assert len(ase.io.read(io.StringIO(r.read().decode()), ":", format="extxyz")) == 1
    r = post(server, {"uuids": [RD, SUS], "include_suspect": True})
    hdr = json.loads(r.headers["X-Export-Summary"])
    assert hdr["suspect_skipped"] == 0 and hdr["planned_frames"] == 2
    fr = ase.io.read(io.StringIO(r.read().decode()), ":", format="extxyz")
    assert [bool(a.info.get("rs_suspect")) for a in fr] == [False, True]
