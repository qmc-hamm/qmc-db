import io
import shutil

import ase.io
import h5py
import numpy as np
from pytest import approx

from explorer import export, links
from explorer.normalize import parse_file

RD, RP1, RP2 = ("f3beda31a0284b3ab3a64e54432960ce", "0f75d82394714b668a5cd639f48866df",
                "3b7edad5edcc4925a6b456ca1c28dd3f")
IL = "1da8a80343085ab99b7a8e334ad43e0e"
AU = "146bb6afc508452ebd5ebc0c8b50ae0b"


def row(fx, u, key=None):
    r = parse_file(fx(u + ".h5"), key or u + ".h5", "e", 1)
    r.update(P_used=r.get("P_recorded"), P_estimated=False, P_extrapolated=False)
    return r


def frames(rows, cache, summary=None):
    s = {} if summary is None else summary
    text = "".join(export.iter_xyz(rows, cache, s))
    return ase.io.read(io.StringIO(text), ":", format="extxyz"), s


def test_energy_is_corrected(fx):
    (a,), _ = frames([row(fx, RD)], fx("").parent / "fixtures")
    assert a.get_potential_energy() == approx(-1430.71269847)
    assert a.info["energy_uncorrected"] == approx(-1440.15100647)
    assert a.info["fsc_applied"] is True
    assert "total_energy" not in a.info or a.info["total_energy"] == approx(-1430.71269847)


def test_no_substitution(fx, tmp_path):
    shutil.copy(fx(RD + ".h5"), tmp_path / (RD + ".h5"))
    with h5py.File(tmp_path / (RD + ".h5"), "r+") as f:
        f["observables/total_energy"][...] = np.nan
        assert np.isfinite(f["observables/total_energy_uncorrected"][0])
    (a,), _ = frames([parse_file(tmp_path / (RD + ".h5"), RD + ".h5", "e", 1)], tmp_path)
    assert np.isnan(a.get_potential_energy())


def test_roundtrip_geometry(fx):
    r = row(fx, IL)
    r.update(config_id="c1", linked_dft_uuids=["0326c6f68ad05777a371f5886ae5f136"])
    (a,), _ = frames([r], fx("").parent / "fixtures")
    assert a.info["linked_dft_uuids"] == "0326c6f68ad05777a371f5886ae5f136"
    with h5py.File(fx(IL + ".h5")) as f:
        cell = np.reshape(f["structure/lattice_vectors"][()], (3, 3))
        pos = np.reshape(f["structure/positions"][()], (-1, 3))
        frc = np.reshape(f["observables/forces"][()], (-1, 3))
    assert np.abs(a.cell[:] - cell).max() < 1e-10
    assert np.abs(a.positions - pos).max() < 1e-10
    assert np.abs(a.get_forces() - frc).max() < 1e-10
    assert a.pbc.all()
    for k in ("uuid", "source", "method", "T", "P", "P_estimated", "rs", "config_id",
              "phase", "creation_path"):
        assert k in a.info, k
    assert ("qmc_quality" in a.info) == (r["qmc_quality"] is not None)
    assert a.info["uuid"] == IL


def test_linked_dft_pairs(fx):
    cache = fx("").parent / "fixtures"
    rows = [row(fx, u) for u in (RD, RP1, RP2)]
    res, _ = links.resolve(rows, lambda u: links.geom_from_file(fx(u + ".h5")))
    for r in rows:
        r.update(res[r["uuid"]])
    by = {r["uuid"]: r for r in rows}
    out = export.expand_with_linked_dft([by[RD]], by)
    assert [r["uuid"] for r in out] == [RD, RP1, RP2]
    fr, _ = frames(out, cache)
    assert len(fr) == 3
    assert fr[0].info["linked_dft_uuids"] == ",".join(sorted([RP1, RP2]))
    for d in fr[1:]:
        assert np.allclose(d.cell[:], fr[0].cell[:], atol=1e-6)
        assert d.info["linked_qmc_uuids"] == RD
        assert d.info["config_id"] == fr[0].info["config_id"]
        assert np.allclose(d.positions, fr[0].positions, atol=1e-5)


def test_forces_shape_mismatch(fx):
    r = row(fx, AU)
    (a,), s = frames([r], fx("").parent / "fixtures")
    assert "forces" not in a.arrays
    assert a.info["forces_shape_mismatch"] is True
    assert s["forces_omitted"] == 1


def test_missing_cache_file(fx, tmp_path):
    r = row(fx, IL)
    fr, s = frames([r], tmp_path)
    assert fr == [] and s["missing_files"] == 1 and s["frames"] == 0


def test_multispecies_without_symbols_skipped(fx):
    n = "FeAlN_32_atoms_CASCI.h5"
    r = parse_file(fx(n), n, "e", 1)
    fr, s = frames([r], fx("").parent / "fixtures")
    assert fr == [] and s["symbols_unknown"] == 1 and s["frames"] == 0


def test_malformed_xyz_falls_back_to_species(fx, tmp_path):
    shutil.copy(fx(RD + ".h5"), tmp_path / (RD + ".h5"))
    with h5py.File(tmp_path / (RD + ".h5"), "r+") as f:
        del f["structure/xyz"]
        f.create_group("structure/xyz")
    r = parse_file(tmp_path / (RD + ".h5"), RD + ".h5", "e", 1)
    (a,), s = frames([r], tmp_path)
    assert set(a.get_chemical_symbols()) == {"H"} and s["symbols_unknown"] == 0


def test_unreadable_file_skipped(fx, tmp_path):
    (tmp_path / (IL + ".h5")).write_bytes(b"not hdf5")
    fr, s = frames([row(fx, IL)], tmp_path)
    assert fr == [] and s["unreadable_files"] == 1 and s["frames"] == 0


def test_stress_and_suspect_info(fx):
    cache = fx("").parent / "fixtures"
    (a,), _ = frames([row(fx, RP1)], cache)
    with h5py.File(fx(RP1 + ".h5")) as f:
        p_rec = float(f["observables/pressure"][0])
    # schema stress is GPa with ASE sign; exported in eV/A^3, so -tr/3 * 160.2 is the recorded virial P
    assert -np.trace(a.get_stress(voigt=False)) / 3 * export.EV_PER_A3_TO_GPA == approx(p_rec)
    r = row(fx, IL)
    r.update(rs_suspect=True, geom_mismatch=["aa", "bb"])
    (b,), s = frames([r], cache)
    assert b.info["rs_suspect"] and b.info["geom_mismatch"] == "aa,bb" and s["suspect_frames"] == 1
    (c,), s = frames([row(fx, IL)], cache)
    assert "rs_suspect" not in c.info and "geom_mismatch" not in c.info and s["suspect_frames"] == 0


def test_bad_frame_does_not_truncate(fx, tmp_path):
    shutil.copy(fx(RD + ".h5"), tmp_path / (RD + ".h5"))
    shutil.copy(fx(IL + ".h5"), tmp_path / (IL + ".h5"))
    with h5py.File(tmp_path / (RD + ".h5"), "r+") as f:
        del f["structure/lattice_vectors"]
    fr, s = frames([row(fx, RD), row(fx, IL)], tmp_path)
    assert len(fr) == 1 and s["unreadable_files"] == 1 and s["frames"] == 1


def test_stress_only_for_ricky(fx):
    cache = fx("").parent / "fixtures"
    (a,), s = frames([row(fx, RP1)], cache)
    assert a.calc is not None and "stress" in a.calc.results and "stress_omitted" not in a.info
    assert s["stress_omitted"] == 0
    (b,), s = frames([row(fx, "0326c6f68ad05777a371f5886ae5f136")], cache)
    assert "stress" not in b.calc.results and b.info["stress_omitted"] == "unverified units"
    assert s["stress_omitted"] == 1
