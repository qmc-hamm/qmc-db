"""Parse one QMC-HAMM HDF5 file into a flat, JSON-serializable index row."""
import ast
import math
from pathlib import Path

import h5py
import numpy as np

NORMALIZE_VERSION = 2
BOHR = 0.529177210903
DEFECT_SYSTEMS = {"FeAlN", "CrAlN", "NVdiamond", "SiVdiamond"}

ROW_FIELDS = [
    "key", "uuid", "etag", "size", "tab",
    "system", "formula", "method", "natoms", "author",
    "source", "source_inferred",
    "T", "P_recorded", "rs_computed", "rs_recorded", "rs_suspect",
    "phase", "creation_path", "qmc_quality",
    "E", "E_err", "E_unc", "E_per_atom", "fsc_applied",
    "dependent_uuids", "forces_shape_ok", "has_sk", "parse_error",
]


def _py(v):
    """Decode bytes / numpy scalars / arrays into plain Python."""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        if v.ndim == 0:
            return _py(v[()])
        return [_py(x) for x in v.tolist()]
    if isinstance(v, np.generic):
        return _py(v.item())
    if isinstance(v, list):
        return [_py(x) for x in v]
    return v


def _num(v):
    """Finite float or None."""
    try:
        v = _py(v)
        if isinstance(v, list):
            v = v[0] if len(v) == 1 else None
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def _str(v):
    v = _py(v)
    if isinstance(v, list):
        v = v[0] if len(v) == 1 else None
    return None if v is None else str(v)


def _list(v):
    v = _py(v)
    if isinstance(v, str):
        try:
            v = ast.literal_eval(v)
        except Exception:
            return [v]
    if not isinstance(v, list):
        v = [v]
    return [str(x) for x in v]


def _ds(f, path):
    try:
        return f[path][()] if path in f else None
    except Exception:
        return None


def _safe(fn, default=None):
    try:
        return fn()
    except Exception:
        return default


def parse_file(path, key, etag, size):
    row = {k: None for k in ROW_FIELDS}
    row.update(key=key, etag=etag, size=size, tab="other", dependent_uuids=[])
    base = str(key).rsplit("/", 1)[-1]
    row["uuid"] = base[:-3] if base.endswith(".h5") else base
    try:
        f = h5py.File(Path(path), "r")
    except Exception as e:
        row["parse_error"] = f"{type(e).__name__}: {e}"
        return row
    try:
        _fill(f, row)
    except Exception as e:
        row["parse_error"] = f"{type(e).__name__}: {e}"
    finally:
        f.close()
    return row


def _fill(f, row):
    a = f.attrs
    g = lambda k: a[k] if k in a else None
    u = _safe(lambda: _str(g("uuid")))
    if u:
        row["uuid"] = u
    row["system"] = _safe(lambda: _str(g("system")))
    row["formula"] = _safe(lambda: _str(g("formula")))
    row["method"] = _safe(lambda: _str(g("method")))
    row["natoms"] = _safe(lambda: int(_py(g("natoms"))))
    row["author"] = _safe(lambda: ", ".join(_list(g("author"))) if g("author") is not None else None)
    method, natoms = row["method"], row["natoms"]

    # source
    src = _safe(lambda: _str(g("source_dataset"))) or _safe(lambda: _str(g("dataset")))
    row["source_inferred"] = False
    if not src:
        if row["author"] == "ShubhangG":
            src, row["source_inferred"] = "aurora_sg", True
        else:
            src = "unknown"
    row["source"] = src

    # conditions
    row["T"] = _num(g("temperature"))
    if row["T"] is None:
        row["T"] = _num(g("structure_temperature"))
    row["P_recorded"] = _num(g("pressure"))
    if row["P_recorded"] is None:
        row["P_recorded"] = _num(g("structure_pressure"))
    row["rs_recorded"] = _num(g("rs"))
    species = _safe(lambda: _list(g("species")), [])
    row["rs_computed"] = _safe(lambda: _rs(f, natoms))

    # tab
    sysn, form = row["system"], row["formula"]
    if sysn == "Hydrogen" or form == "H":
        row["tab"] = "hydrogen"
    elif method == "CASCI" or sysn in DEFECT_SYSTEMS:
        row["tab"] = "defects"
    elif "Li" in species:
        row["tab"] = "li7p3s11"
    elif species == ["C"]:
        row["tab"] = "graphene"
    if row["tab"] == "hydrogen":
        rc, rr = row["rs_computed"], row["rs_recorded"]
        row["rs_suspect"] = bool(
            rc is None or not (0.8 <= rc <= 3.0)
            or (rr is not None and abs(rr - rc) > 0.01 * rc)
        )

    # labels
    phase = _safe(lambda: _str(g("phase"))) or _safe(lambda: _str(_ds(f, "provenance/phase_label_raw")))
    row["phase"] = phase or "unknown"
    cp = _safe(lambda: _str(_ds(f, "provenance/config_creation_method")))
    cp = cp or _safe(lambda: _str(_ds(f, "provenance/starting_configuration_model_name")))
    cp = cp or _safe(lambda: _str(g("starting_configuration_model_name")))
    cp = cp or _safe(lambda: _str(g("ensemble")))
    row["creation_path"] = cp
    row["qmc_quality"] = _safe(lambda: int(_py(g("qmc_quality"))) if g("qmc_quality") is not None else None)

    # energies
    te = _ds(f, "observables/total_energy")
    if te is not None and np.size(te) == 1:
        row["E"] = _num(np.ravel(te)[0])
    row["E_err"] = _num(_ds(f, "observables/total_energy_error")) if row["E"] is not None else None
    row["E_unc"] = _num(_ds(f, "observables/total_energy_uncorrected"))
    if row["E"] is not None and natoms:
        row["E_per_atom"] = row["E"] / natoms
    if method in ("VMC", "DMC"):
        fsc = _num(_ds(f, "observables/fsc_kinetic_energy"))
        row["fsc_applied"] = bool(fsc is not None or row["qmc_quality"] == 10)

    # links and flags
    deps = _ds(f, "provenance/dependent_uuids")
    row["dependent_uuids"] = _safe(lambda: [str(x) for x in _py(deps)] if deps is not None else [], [])
    row["forces_shape_ok"] = _safe(lambda: _forces_ok(f, natoms))
    sk = _ds(f, "observables/structure_factor")
    row["has_sk"] = bool(sk is not None and np.size(sk) > 0)


def _rs(f, natoms):
    L = np.asarray(f["structure/lattice_vectors"][()], dtype=float)
    vol = abs(np.linalg.det(L))
    return float((3 * vol / (4 * math.pi * natoms)) ** (1 / 3) / BOHR)


def _forces_ok(f, natoms):
    if "observables/forces" not in f:
        return None
    s = f["observables/forces"].shape
    return bool(s in ((1, natoms, 3), (natoms, 3)))
