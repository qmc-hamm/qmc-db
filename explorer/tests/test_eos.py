import numpy as np
import pytest

from explorer.eos import FittedEos, apply_eos


def truth(rs, T):
    return 100 * (1.5 / rs) ** 4.5 * (1 + 2e-5 * T)


def synth_rows(n=2000):
    rng = np.random.default_rng(0)
    rs = rng.uniform(1.2, 1.8, n)
    T = rng.uniform(500, 4000, n)
    P = truth(rs, T) * (1 + 0.01 * rng.standard_normal(n))
    return [dict(uuid=f"u{i}", tab="hydrogen", source="s", creation_path="p",
                 T=float(T[i]), P_recorded=float(P[i]), rs_computed=float(rs[i]),
                 rs_suspect=False) for i in range(n)]


def test_fit_recovers_synthetic():
    eos = FittedEos.fit(synth_rows())
    P, ind = eos.pressure(np.array([1.4]), np.array([2000.0]))
    assert abs(P[0] / truth(1.4, 2000) - 1) < 0.02
    assert ind[0]


def test_out_of_domain_flag():
    eos = FittedEos.fit(synth_rows())
    _, ind = eos.pressure(np.array([2.5]), np.array([2000.0]))
    assert not ind[0]


def test_training_filter():
    rows = synth_rows(200)
    bad = [dict(r, rs_suspect=True) for r in rows]
    with pytest.raises(ValueError, match="no training rows"):
        FittedEos.fit(bad)
    mixed = rows[:100] + [dict(r, P_recorded=None) for r in rows[100:150]] \
        + [dict(r, T=None) for r in rows[150:]]
    assert len(FittedEos.fit(mixed).residuals(mixed)) == 100


def test_apply_eos():
    eos = FittedEos.fit(synth_rows())
    base = dict(uuid="x", tab="hydrogen", rs_computed=1.4, T=2000.0, rs_suspect=False)
    a = dict(base, P_recorded=123.0)
    b = dict(base, P_recorded=None)
    c = dict(base, P_recorded=None, rs_suspect=True)
    d = dict(base, P_recorded=None, T=None)
    apply_eos([a, b, c, d], eos)
    assert a["P_used"] == 123.0 and a["P_estimated"] is False
    assert isinstance(b["P_used"], float) and b["P_estimated"] is True
    assert c["P_used"] is None and d["P_used"] is None


def test_json_roundtrip():
    eos = FittedEos.fit(synth_rows())
    e2 = FittedEos.from_json(eos.to_json())
    rs = np.linspace(1.1, 2.0, 20)
    T = np.linspace(300, 5000, 20)
    P1, d1 = eos.pressure(rs, T)
    P2, d2 = e2.pressure(rs, T)
    assert np.allclose(P1, P2, rtol=1e-12, atol=0) and (d1 == d2).all()
