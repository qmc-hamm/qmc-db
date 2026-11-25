"""Pressure equation of state P(rs, T) fitted to rows with a recorded pressure."""
import math
from typing import Protocol

import numpy as np
import pandas as pd

NB = 20
MIN_PTS = 5


class EosBackend(Protocol):
    name: str

    def pressure(self, rs: np.ndarray, T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


def _finite(v):
    return v is not None and math.isfinite(v)


def _training(rows):
    seen, out = set(), []
    for r in rows:
        if r.get("tab", "hydrogen") != "hydrogen" or r.get("rs_suspect") is True:
            continue
        P, T, rs = r.get("P_recorded"), r.get("T"), r.get("rs_computed")
        if not (_finite(P) and _finite(T) and _finite(rs)) or P <= 0 or rs <= 0:
            continue
        key = r.get("config_id") or r["uuid"]
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def _design(x, t, a, b):
    return np.stack([x**i * t**j for i in range(a + 1) for j in range(b + 1)], axis=1)


def _bin(v, edges):
    return np.clip(np.searchsorted(edges, v, side="right") - 1, 0, len(edges) - 2)


class FittedEos:
    name = "fit"

    def __init__(self, coef, a, b, xedges, tedges, mask, cv_rms):
        self.coef = np.asarray(coef, float)
        self.a, self.b = int(a), int(b)
        self.xedges = np.asarray(xedges, float)
        self.tedges = np.asarray(tedges, float)
        self.mask = np.asarray(mask, bool)
        self.cv_rms = float(cv_rms)

    @classmethod
    def fit(cls, rows):
        tr = _training(rows)
        if not tr:
            raise ValueError("no training rows")
        x = np.log([r["rs_computed"] for r in tr])
        t = np.array([r["T"] for r in tr], float) / 1000
        P = np.array([r["P_recorded"] for r in tr], float)
        y = np.log(P)
        folds = np.array_split(np.random.default_rng(0).permutation(len(tr)), 5)
        best = None
        for a in (2, 3, 4):
            for b in (1, 2):
                A = _design(x, t, a, b)
                err = []
                for f in folds:
                    m = np.ones(len(tr), bool)
                    m[f] = False
                    if m.sum() == 0 or len(f) == 0:
                        continue
                    c = np.linalg.lstsq(A[m], y[m], rcond=None)[0]
                    err.append((np.exp(A[f] @ c) - P[f]) / P[f])
                rms = float(np.sqrt(np.mean(np.concatenate(err) ** 2))) if err else np.inf
                if best is None or rms < best[0] - 1e-12:
                    best = (rms, a, b)
        rms, a, b = best
        coef = np.linalg.lstsq(_design(x, t, a, b), y, rcond=None)[0]
        xe = np.linspace(x.min(), x.max(), NB + 1)
        te = np.linspace(t.min(), t.max(), NB + 1)
        cnt = np.zeros((NB, NB), int)
        np.add.at(cnt, (_bin(x, xe), _bin(t, te)), 1)
        return cls(coef, a, b, xe, te, cnt >= MIN_PTS, rms)

    def pressure(self, rs, T):
        x = np.log(np.asarray(rs, float))
        t = np.asarray(T, float) / 1000
        P = np.exp(_design(x, t, self.a, self.b) @ self.coef)
        inr = ((x >= self.xedges[0]) & (x <= self.xedges[-1])
               & (t >= self.tedges[0]) & (t <= self.tedges[-1]))
        ind = inr & self.mask[_bin(x, self.xedges), _bin(t, self.tedges)]
        return P, ind

    def to_json(self):
        return dict(name=self.name, coef=self.coef.tolist(), a=self.a, b=self.b,
                    xedges=self.xedges.tolist(), tedges=self.tedges.tolist(),
                    mask=self.mask.tolist(), cv_rms=self.cv_rms)

    @classmethod
    def from_json(cls, d):
        return cls(d["coef"], d["a"], d["b"], d["xedges"], d["tedges"], d["mask"], d["cv_rms"])

    def residuals(self, rows):
        tr = _training(rows)
        rs = np.array([r["rs_computed"] for r in tr], float)
        T = np.array([r["T"] for r in tr], float)
        P = np.array([r["P_recorded"] for r in tr], float)
        Pf, _ = self.pressure(rs, T)
        return pd.DataFrame(dict(
            source=[r.get("source") for r in tr],
            creation_path=[r.get("creation_path") for r in tr],
            rs=rs, T=T, P_recorded=P, P_fit=Pf, rel_resid=(Pf - P) / P))


def apply_eos(rows, eos):
    for r in rows:
        if r.get("P_recorded") is not None:
            r["P_used"], r["P_estimated"], r["P_extrapolated"] = r["P_recorded"], False, False
        elif (r.get("tab", "hydrogen") == "hydrogen" and r.get("T") is not None
              and r.get("rs_computed") is not None and r.get("rs_suspect") is not True):
            P, ind = eos.pressure(np.array([r["rs_computed"]]), np.array([r["T"]]))
            r["P_used"], r["P_estimated"], r["P_extrapolated"] = float(P[0]), True, not bool(ind[0])
        else:
            r["P_used"], r["P_estimated"], r["P_extrapolated"] = None, False, False
