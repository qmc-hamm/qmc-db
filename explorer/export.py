"""Export index rows as extxyz frames (one frame per file)."""
import io
from pathlib import Path

import ase.io
import h5py
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

INFO = {"uuid": "uuid", "source": "source", "method": "method", "T": "T", "P": "P_used",
        "P_estimated": "P_estimated", "P_extrapolated": "P_extrapolated",
        "rs": "rs_computed", "config_id": "config_id", "qmc_quality": "qmc_quality",
        "phase": "phase", "creation_path": "creation_path", "fsc_applied": "fsc_applied",
        "energy_uncorrected": "E_unc", "energy_error": "E_err"}


EV_PER_A3_TO_GPA = 160.21766208


def _path(row, cache):
    key = row["key"]
    return Path(key[6:]) if key.startswith("local:") else Path(cache) / Path(key).name


def _symbols(f, n):
    try:
        lines = f["structure/xyz"][()].decode().splitlines()[2:2 + n]
        sym = [ln.split()[0] for ln in lines]
        if len(sym) == n:
            return sym
    except (KeyError, AttributeError, IndexError, TypeError, ValueError):
        pass
    sp = [s.decode() if isinstance(s, bytes) else str(s) for s in np.ravel(f.attrs.get("species", []))]
    return sp * n if len(sp) == 1 else None


def _frame(row, path, summary):
    with h5py.File(path, "r") as f:
        cell = np.asarray(f["structure/lattice_vectors"][()], float).reshape(-1, 3, 3)[0]
        pos = np.asarray(f["structure/positions"][()], float)
        pos = pos.reshape(-1, pos.shape[-2], 3)[0]
        pbc = np.asarray(f["structure/pbc"][()], bool).reshape(3)
        n = len(pos)
        sym = _symbols(f, n)
        if sym is None:
            summary["symbols_unknown"] += 1
            return None
        atoms = Atoms(sym, positions=pos, cell=cell, pbc=pbc)
        forces = ferr = stress = None
        has_stress = "observables/stress" in f
        if has_stress and row.get("source") != "ricky_legacy":
            atoms.info["stress_omitted"] = "unverified units"
            summary["stress_omitted"] += 1
        elif has_stress:
            st = np.asarray(f["observables/stress"][()], float)
            if st.size == 9:
                stress = st.reshape(3, 3)
            elif st.size == 6:
                stress = np.array([[st[0], st[5], st[4]], [st[5], st[1], st[3]], [st[4], st[3], st[2]]])
            if stress is not None and np.isfinite(stress).all():
                stress = stress / EV_PER_A3_TO_GPA  # schema GPa -> ASE eV/A^3, ASE sign (P = -tr/3)
            else:
                stress = None
        if "observables/forces" in f:
            fr = np.asarray(f["observables/forces"][()], float)
            if fr.size == 3 * n:
                forces = fr.reshape(n, 3)
                if "observables/forces_error" in f:
                    fe = np.asarray(f["observables/forces_error"][()], float)
                    if fe.size == 3 * n:
                        ferr = fe.reshape(n, 3)
            else:
                atoms.info["forces_shape_mismatch"] = True
                summary["forces_omitted"] += 1
    e = row.get("E")
    atoms.calc = SinglePointCalculator(
        atoms, energy=float("nan") if e is None else e,
        **({} if forces is None else {"forces": forces}),
        **({} if stress is None else {"stress": stress}))
    if ferr is not None:
        atoms.arrays["forces_error"] = ferr
    for k, src in INFO.items():
        if row.get(src) is not None:
            atoms.info[k] = row[src]
    if row.get("linked_dft_uuids"):
        atoms.info["linked_dft_uuids"] = ",".join(row["linked_dft_uuids"])
    if row.get("linked_qmc_uuids"):
        atoms.info["linked_qmc_uuids"] = ",".join(row["linked_qmc_uuids"])
    if row.get("geom_mismatch"):
        atoms.info["geom_mismatch"] = ",".join(row["geom_mismatch"])
    if row.get("rs_suspect") is True:
        atoms.info["rs_suspect"] = True
    return atoms


def iter_xyz(rows, cache, summary):
    for k, v in (("frames", 0), ("missing_files", 0), ("forces_omitted", 0), ("symbols_unknown", 0), ("unreadable_files", 0), ("suspect_frames", 0), ("stress_omitted", 0), ("methods", {})):
        summary.setdefault(k, v)
    for row in rows:
        path = _path(row, cache)
        if not path.exists():
            summary["missing_files"] += 1
            continue
        try:
            atoms = _frame(row, path, summary)
            if atoms is None:
                continue
            buf = io.StringIO()
            ase.io.write(buf, atoms, format="extxyz")
        except Exception:
            summary["unreadable_files"] += 1
            continue
        summary["frames"] += 1
        summary["suspect_frames"] += bool(atoms.info.get("rs_suspect"))
        m = row.get("method")
        summary["methods"][m] = summary["methods"].get(m, 0) + 1
        yield buf.getvalue()


def expand_with_linked_dft(rows, by_uuid):
    out = list(rows)
    seen = {r["uuid"] for r in out}
    extra = sorted({u for r in rows for u in r.get("linked_dft_uuids") or []} - seen)
    return out + [by_uuid[u] for u in extra if u in by_uuid]
