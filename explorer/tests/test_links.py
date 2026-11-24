import numpy as np

from explorer import links
from explorer.normalize import parse_file

RICKY = ["f3beda31a0284b3ab3a64e54432960ce", "0f75d82394714b668a5cd639f48866df",
         "3b7edad5edcc4925a6b456ca1c28dd3f"]
ILNUR = ["1da8a80343085ab99b7a8e334ad43e0e", "0326c6f68ad05777a371f5886ae5f136",
         "1e08265a6ee85e3d8e0fb67646563aa9"]
AURORA = ["b525945ff3354d32b7c4e143547cccbf", "4c9372432d36491ea1b1d7b41a21bb9b",
          "146bb6afc508452ebd5ebc0c8b50ae0b"]


def load(fx, uuids):
    paths = {u: fx(u + ".h5") for u in uuids}
    rows = [parse_file(p, u + ".h5", "e", 1) for u, p in paths.items()]
    return rows, lambda u: links.geom_from_file(paths[u]) if u in paths else None


def test_links_ricky(fx):
    rows, geom = load(fx, RICKY)
    res, s = links.resolve(rows, geom)
    assert res[RICKY[0]]["linked_dft_uuids"] == sorted(RICKY[1:])
    for d in RICKY[1:]:
        assert res[d]["linked_qmc_uuids"] == [RICKY[0]]
    assert {res[u]["config_id"] for u in RICKY} == {min(RICKY)}


def test_links_chain_aurora(fx):
    rows, geom = load(fx, AURORA)
    res, _ = links.resolve(rows, geom)
    assert res[AURORA[0]]["linked_dft_uuids"] == [AURORA[2]]
    assert AURORA[1] not in res[AURORA[0]]["linked_dft_uuids"]
    assert len({res[u]["config_id"] for u in AURORA}) == 1


def test_links_ilnur(fx):
    rows, geom = load(fx, ILNUR)
    res, _ = links.resolve(rows, geom)
    assert res[ILNUR[0]]["linked_dft_uuids"] == [ILNUR[1]]


def test_links_missing_and_mismatch(fx):
    rows, geom = load(fx, ILNUR)
    rows[0]["dependent_uuids"] = list(rows[0]["dependent_uuids"]) + ["deadbeef"]
    res, s = links.resolve(rows, geom)
    assert res[ILNUR[0]]["missing_links"] == ["deadbeef"] and s["missing"] == 1

    def geom2(u):
        g = geom(u)
        if u == ILNUR[1]:
            return g[0] * 1.1, g[1]
        return g
    rows, _ = load(fx, ILNUR)
    res, _ = links.resolve(rows, geom2)
    assert ILNUR[1] in res[ILNUR[0]]["geom_mismatch"]
    assert res[ILNUR[0]]["config_id"] != res[ILNUR[1]]["config_id"]


def test_rqmc_counts_as_qmc(fx):
    rows, geom = load(fx, RICKY)
    rows[0]["method"] = "RQMC"
    res, _ = links.resolve(rows, geom)
    assert res[RICKY[1]]["linked_qmc_uuids"] == [RICKY[0]]
