import gzip
import json
from pathlib import Path

import pytest

from explorer import bucket, harvest, normalize

FIXD = Path(__file__).parent / "fixtures"
LI = Path(__file__).resolve().parents[2] / "Li7P3S11_QMC_DFT_dataset"


def fake_remote():
    return [{"key": "QMCHAMM/", "etag": "d", "size": 0, "last_modified": "x"}] + [{"key": "QMCHAMM/" + p.name, "etag": "e-" + p.name, "size": p.stat().st_size,
             "last_modified": "x"} for p in sorted(FIXD.glob("*.h5"))]


def run(tmp, monkeypatch, remote, *extra):
    monkeypatch.setattr(bucket, "list_bucket", lambda *a, **k: remote)
    harvest.main(["--no-download", "--cache", str(FIXD), "--data", str(tmp),
                  "--local-dir", str(LI), *extra])
    idx = json.loads(gzip.open(tmp / "index.json.gz").read())
    return idx, json.loads((tmp / "harvest_report.json").read_text())


def test_harvest_offline(tmp_path, monkeypatch):
    idx, rep = run(tmp_path, monkeypatch, fake_remote())
    cols = idx["columns"]
    assert "P_estimated" in idx["fields"] and len(cols["P_estimated"]) == len(cols["uuid"])
    tabs = [t for t, k in zip(cols["tab"], cols["key"]) if k.startswith("local:")]
    assert len(tabs) == 200 and set(tabs) == {"li7p3s11"}
    assert rep["rs_suspect"] == 1
    assert rep["parsed"] == len(cols["uuid"])


def test_harvest_sync_deleted(tmp_path, monkeypatch):
    remote = fake_remote()
    idx1, _ = run(tmp_path, monkeypatch, remote)
    idx2, rep = run(tmp_path, monkeypatch, remote[2:])
    assert rep["deleted"] == 1 and rep["parsed"] == 0
    assert len(idx2["columns"]["uuid"]) == len(idx1["columns"]["uuid"]) - 1
    assert remote[1]["key"] not in idx2["columns"]["key"]


def _run_dl(tmp, monkeypatch, remote, *extra):
    monkeypatch.setattr(bucket, "list_bucket", lambda *a, **k: remote)
    monkeypatch.setattr(bucket, "download", lambda obj, cache: Path(cache) / Path(obj["key"]).name)
    harvest.main(["--cache", str(FIXD), "--data", str(tmp), "--local-dir", str(LI), *extra])
    return json.loads((tmp / "harvest_report.json").read_text())


def test_no_download_changed_keeps_old(tmp_path, monkeypatch):
    remote = fake_remote()
    run(tmp_path, monkeypatch, remote)
    k = remote[1]["key"]
    remote2 = [dict(o) for o in remote]
    remote2[1]["etag"] = "NEW"
    idx, rep = run(tmp_path, monkeypatch, remote2)
    assert rep["stale"] == 1 and rep["parsed"] == 0
    assert k in idx["columns"]["key"]
    man = json.loads((tmp_path / "manifest.json").read_text())
    assert man[k]["etag"] == remote[1]["etag"]


def test_changed_etag_reparsed(tmp_path, monkeypatch):
    remote = fake_remote()
    run(tmp_path, monkeypatch, remote)
    remote2 = [dict(o) for o in remote]
    remote2[1]["etag"] = "NEW"
    rep = _run_dl(tmp_path, monkeypatch, remote2)
    assert rep["changed"] == 1 and rep["parsed"] == 1 and rep["stale"] == 0
    man = json.loads((tmp_path / "manifest.json").read_text())
    assert man[remote[1]["key"]]["etag"] == "NEW"


def test_version_mismatch_reparses(tmp_path, monkeypatch):
    remote = fake_remote()
    idx, _ = run(tmp_path, monkeypatch, remote)
    assert idx["normalize_version"] == normalize.NORMALIZE_VERSION
    _, rep = run(tmp_path, monkeypatch, remote)
    assert rep["parsed"] == 0
    old = json.loads(gzip.open(tmp_path / "index.json.gz").read())
    old["normalize_version"] = normalize.NORMALIZE_VERSION - 1
    with gzip.open(tmp_path / "index.json.gz", "wt") as fh:
        json.dump(old, fh)
    idx, rep = run(tmp_path, monkeypatch, remote)
    assert rep["parsed"] == len(idx["columns"]["uuid"]) and rep["download_errors"] == []
    assert idx["normalize_version"] == normalize.NORMALIZE_VERSION


def test_limit_refused_on_default_data():
    with pytest.raises(SystemExit):
        harvest.main(["--limit", "1"])


def test_failed_changed_download_keeps_old(tmp_path, monkeypatch):
    remote = fake_remote()
    run(tmp_path, monkeypatch, remote)
    remote2 = [dict(o) for o in remote]
    remote2[1]["etag"] = "NEW"
    monkeypatch.setattr(bucket, "list_bucket", lambda *a, **k: remote2)
    monkeypatch.setattr(bucket, "download", lambda obj, cache: (_ for _ in ()).throw(IOError("boom")))
    harvest.main(["--cache", str(FIXD), "--data", str(tmp_path), "--local-dir", str(LI)])
    idx = json.loads(gzip.open(tmp_path / "index.json.gz").read())
    man = json.loads((tmp_path / "manifest.json").read_text())
    assert remote[1]["key"] in idx["columns"]["key"]
    assert man[remote[1]["key"]]["etag"] == remote[1]["etag"]  # old ETag kept, so the next run retries


def test_unchanged_missing_cache_redownloaded(tmp_path, monkeypatch):
    import shutil
    cache = tmp_path / "cache"
    shutil.copytree(FIXD, cache)
    remote = fake_remote()
    monkeypatch.setattr(bucket, "list_bucket", lambda *a, **k: remote)
    got = []
    monkeypatch.setattr(bucket, "download", lambda obj, c: got.append(obj["key"]) or shutil.copy(FIXD / Path(obj["key"]).name, c))
    args = ["--cache", str(cache), "--data", str(tmp_path / "d"), "--local-dir", str(LI)]
    harvest.main(args)
    got.clear()
    (cache / Path(remote[1]["key"]).name).unlink()
    harvest.main(args)
    assert got == [remote[1]["key"]]


def test_no_fit_removes_eos_files(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch, fake_remote())
    assert (tmp_path / "eos_residuals.json").exists()
    monkeypatch.setattr(harvest.eosmod.FittedEos, "fit", classmethod(lambda cls, rows: (_ for _ in ()).throw(ValueError("no training rows"))))
    _, rep = run(tmp_path, monkeypatch, fake_remote())
    assert rep["eos"]["name"] == "none"
    assert not (tmp_path / "eos_residuals.json").exists() and not (tmp_path / "eos_fit.json").exists()
