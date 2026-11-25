"""Group DMC/VMC/DFT files of one geometry into config_id link groups."""
import h5py
import numpy as np


def geom_from_file(path):
    try:
        with h5py.File(path, "r") as f:
            cell = np.asarray(f["structure/lattice_vectors"][()], float).reshape(-1, 3, 3)[0]
            pos = np.asarray(f["structure/positions"][()], float)
            pos = pos.reshape(-1, pos.shape[-2], 3)[0]
        return cell, pos
    except Exception:
        return None


def _same(ga, gb):
    (ca, pa), (cb, pb) = ga, gb
    if len(pa) != len(pb) or not np.allclose(ca, cb, atol=1e-6):
        return False
    try:
        d = (pa - pb) @ np.linalg.inv(ca)
    except np.linalg.LinAlgError:
        return np.allclose(pa, pb, atol=1e-5)
    d -= np.round(d)
    return bool(np.abs(d).max() < 1e-5)


def resolve(rows, geom):
    by = {r["uuid"]: r for r in rows}
    parent = {u: u for u in by}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    missing = {u: [] for u in by}
    mism = {u: [] for u in by}
    total = linked = 0
    for u, r in by.items():
        for v in r.get("dependent_uuids") or []:
            total += 1
            if v not in by:
                missing[u].append(v)
                continue
            ga, gb = geom(u), geom(v)
            if ga is None or gb is None or not _same(ga, gb):
                mism[u].append(v)
                continue
            linked += 1
            parent[find(u)] = find(v)
    comps = {}
    for u in by:
        comps.setdefault(find(u), []).append(u)
    res = {}
    for members in comps.values():
        cid = min(members)
        for u in members:
            others = [m for m in members if m != u]
            res[u] = {
                "config_id": cid,
                "linked_dft_uuids": sorted(m for m in others if by[m]["method"] == "DFT"),
                "linked_qmc_uuids": sorted(m for m in others if by[m]["method"] != "DFT")
                if by[u]["method"] == "DFT" else [],
                "missing_links": sorted(missing[u]),
                "geom_mismatch": sorted(mism[u]),
            }
    summary = {"links_total": total, "linked": linked,
               "missing": sum(len(v) for v in missing.values()),
               "geom_mismatch": sum(len(v) for v in mism.values()),
               "components": len(comps)}
    return res, summary
