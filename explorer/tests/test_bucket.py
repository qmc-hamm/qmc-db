from explorer import bucket

XML = """<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
<IsTruncated>true</IsTruncated><NextContinuationToken>abc</NextContinuationToken>
<Contents><Key>QMCHAMM/a.h5</Key><LastModified>2026-01-01T00:00:00.000Z</LastModified><ETag>&quot;e1&quot;</ETag><Size>10</Size></Contents>
<Contents><Key>QMCHAMM/b.h5</Key><LastModified>2026-01-02T00:00:00.000Z</LastModified><ETag>&quot;e2&quot;</ETag><Size>20</Size></Contents>
</ListBucketResult>"""


def test_parse_listing():
    objs, tok = bucket.parse_listing(XML)
    assert len(objs) == 2 and tok == "abc"
    assert objs[0]["key"] == "QMCHAMM/a.h5" and objs[0]["etag"] == "e1"
    assert objs[1]["size"] == 20 and isinstance(objs[1]["size"], int)


def test_plan_sync():
    man = {"a": {"etag": "e1"}, "b": {"etag": "e2"}, "c": {"etag": "e3"}}
    rem = [{"key": "a", "etag": "e1"}, {"key": "b", "etag": "e9"}, {"key": "d", "etag": "e4"}]
    p = bucket.plan_sync(rem, man)
    assert (p["new"], p["changed"], p["deleted"], p["unchanged"]) == (["d"], ["b"], ["c"], ["a"])


def _page(trunc, tok=None):
    t = f"<NextContinuationToken>{tok}</NextContinuationToken>" if tok else ""
    return ('<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">'
            f"<IsTruncated>{trunc}</IsTruncated>{t}</ListBucketResult>").encode()


def test_list_bucket_truncated_without_token(monkeypatch):
    import pytest
    monkeypatch.setattr(bucket, "_get", lambda url: _page("true"))
    with pytest.raises(RuntimeError, match="token"):
        bucket.list_bucket("http://x")


def test_list_bucket_repeated_token(monkeypatch):
    import pytest
    monkeypatch.setattr(bucket, "_get", lambda url: _page("true", "same"))
    with pytest.raises(RuntimeError, match="repeat"):
        bucket.list_bucket("http://x")
