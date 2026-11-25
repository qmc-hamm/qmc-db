import pytest

from explorer.normalize import ROW_FIELDS, parse_file


def row(fx, name):
    p = fx(name)
    return parse_file(p, "QMCHAMM/" + name, "etag", p.stat().st_size)


def test_ricky_dmc(fx):
    r = row(fx, "f3beda31a0284b3ab3a64e54432960ce.h5")
    assert set(ROW_FIELDS) <= set(r)
    assert r["source"] == "ricky_legacy" and r["source_inferred"] is False
    assert r["method"] == "DMC" and r["T"] == 600.0 and r["P_recorded"] == 150.0
    assert abs(r["rs_computed"] - 1.5053945) < 1e-5
    assert r["rs_suspect"] is False
    assert r["E"] == pytest.approx(-1430.71269847)
    assert r["E_unc"] == pytest.approx(-1440.15100647)
    assert r["E_err"] is None
    assert r["fsc_applied"] is True and r["qmc_quality"] == 10
    assert len(r["dependent_uuids"]) == 2
    assert r["tab"] == "hydrogen"
    assert r["uuid"] == "f3beda31a0284b3ab3a64e54432960ce"
    assert r["parse_error"] is None


def test_ilnur_dmc(fx):
    r = row(fx, "1da8a80343085ab99b7a8e334ad43e0e.h5")
    assert r["source"] == "ilnur_leonardo"
    assert r["T"] == 200.0 and r["P_recorded"] == 600.0
    assert r["rs_recorded"] is None and abs(r["rs_computed"] - 1.2) < 1e-4
    assert r["phase"] == "CSIV" and r["creation_path"] == "DFT"
    assert r["has_sk"] is True and r["natoms"] == 128


def test_aurora_early_shell(fx):
    r = row(fx, "aeafdf12efae4e76ba34a5f5b94957c9.h5")
    assert r["source"] == "aurora_sg" and r["source_inferred"] is True
    assert r["rs_suspect"] is True and abs(r["rs_computed"] - 0.256) < 0.01
    assert r["P_recorded"] == 100.0 and r["T"] == 3400.0


def test_aurora_deepmd_dft_forces(fx):
    r = row(fx, "146bb6afc508452ebd5ebc0c8b50ae0b.h5")
    assert r["P_recorded"] is None and abs(r["rs_computed"] - 1.29) < 1e-4
    assert r["forces_shape_ok"] is False


def test_kevin_ly(fx):
    r = row(fx, "9683695738e55068a34a19f3a1dac17c.h5")
    assert r["T"] is None and r["P_recorded"] is None
    assert r["source"] == "kevin_ly_atomic_h"


def test_other_tabs(fx):
    assert row(fx, "6ba5fa7a97fa4caab4999e9cdb28bdad.h5")["tab"] == "graphene"
    r = row(fx, "FeAlN_32_atoms_CASCI.h5")
    assert r["tab"] == "defects" and r["E"] is None


def test_unreadable_file(tmp_path):
    p = tmp_path / "junk.h5"
    p.write_bytes(b"0123456789")
    r = parse_file(p, "QMCHAMM/junk.h5", "e", 10)
    assert r["parse_error"] and r["tab"] == "other"


def test_aurora_vmc_not_fsc(fx):
    r = row(fx, "4c9372432d36491ea1b1d7b41a21bb9b.h5")  # quality 10, observables.attrs finite_size_corrected=False
    assert r["method"] == "VMC" and r["qmc_quality"] == 10
    assert r["fsc_applied"] is False
