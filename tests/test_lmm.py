"""Synthetic-data checks for scripts/lmm.py. No real data is touched. Run: python -I tests/test_lmm.py"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
from lmm import lmm_fit, lmm_fit_fallback, satterthwaite, ols, tost_partial_r, bh


def sim(rng, n=47, k=6, visits=(3, 10), beta3=0.0, sd_u=(4.0, 1.0), corr=0.3, sd_e=3.0, balanced=False):
    P = rng.normal(size=n); C = rng.normal(size=(n, k))
    rows = []
    cov = np.array([[sd_u[0] ** 2, corr * sd_u[0] * sd_u[1]], [corr * sd_u[0] * sd_u[1], sd_u[1] ** 2]])
    U = rng.multivariate_normal([0, 0], cov, size=n)
    for i in range(n):
        m = visits[1] if balanced else rng.integers(visits[0], visits[1] + 1)
        tt = np.arange(m, dtype=float) if balanced else np.sort(rng.uniform(0, 10, m))
        for tj in tt:
            mu = 20 + 2.0 * tj + 1.0 * P[i] + beta3 * P[i] * tj + C[i] @ np.full(k, 0.5) + C[i] @ np.full(k, 0.2) * tj
            rows.append((i, tj, mu + U[i, 0] + U[i, 1] * tj + rng.normal(0, sd_e), P[i], *C[i]))
    a = np.array(rows)
    g, t, y, Pv, Cv = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4:]
    X = np.c_[np.ones(len(y)), t, Pv, Pv * t, Cv, Cv * t[:, None]]
    return y, X, g, t


def test_vs_statsmodels():
    from statsmodels.regression.mixed_linear_model import MixedLM
    rng = np.random.default_rng(1)
    y, X, g, t = sim(rng)
    f = lmm_fit(y, X, g, t, "slope_corr")
    m = MixedLM(y, X, groups=g, exog_re=np.c_[np.ones(len(t)), t]).fit(reml=True, method=["lbfgs"])
    db = np.max(np.abs(f["beta"] - m.fe_params))
    dc = np.max(np.abs(f["re_cov"] - np.asarray(m.cov_re)))
    ds = abs(f["sigma"] ** 2 - m.scale)
    # SEs: lme4/lmerTest use sigma^2 (X'V^-1X)^-1 (implemented here); statsmodels inverts the joint
    # observed Hessian, so small differences (<1% relative) are expected.
    dse = np.max(np.abs(np.sqrt(np.diag(f["cov_beta"])) / m.bse_fe - 1))
    print(f"vs statsmodels: max|dbeta|={db:.2e} max|dRE cov|={dc:.2e} |dsigma2|={ds:.2e} max rel dSE={dse:.2e}")
    assert db < 1e-3 and dc < 1e-2 and ds < 1e-2 and dse < 0.01


def test_balanced_df():
    # Balanced design, between-donor slope terms: intercept-slope (t), P*t, 6 C*t = 8.
    # Exact df for P*t is n_groups - 8 = 39.
    rng = np.random.default_rng(2)
    y, X, g, t = sim(rng, balanced=True, visits=(6, 6))
    f = lmm_fit(y, X, g, t, "slope_corr")
    L = np.zeros(X.shape[1]); L[3] = 1
    s = satterthwaite(f, L)
    print(f"balanced df for P*t: {s['df']:.3f} (expected 39); converged={f['converged']}")
    assert abs(s["df"] - 39) < 0.05


def test_type1_and_power(nsim=300):
    rng = np.random.default_rng(3)
    L = np.zeros(4 + 12); L[3] = 1
    rej0 = rej1 = 0; fb = 0
    for i in range(nsim):
        y, X, g, t = sim(rng)
        f = lmm_fit_fallback(y, X, g, t); fb += f["re"] != "slope_corr"
        rej0 += satterthwaite(f, L)["p"] < 0.05
    for i in range(100):
        y, X, g, t = sim(rng, beta3=0.6)
        rej1 += satterthwaite(lmm_fit_fallback(y, X, g, t), L)["p"] < 0.05
    print(f"type I error {rej0/nsim:.3f} (n={nsim}; fallback used {fb}); power at beta3=0.6: {rej1/100:.2f}")
    assert 0.02 < rej0 / nsim < 0.085 and rej1 / 100 > 0.5


def test_fallback_singular():
    rng = np.random.default_rng(4)
    y, X, g, t = sim(rng, sd_u=(4.0, 0.0), corr=0.0)
    f = lmm_fit_fallback(y, X, g, t)
    print("no slope variance -> fallback log:", [(r["re"], r["converged"]) for r in f["fallback_log"]])
    assert f["re"] != "slope_corr"


def test_misc():
    rng = np.random.default_rng(5)
    X = np.c_[np.ones(100), rng.normal(size=(100, 2))]
    y = X @ [1, 2, 0] + rng.normal(size=100)
    o = ols(y, X)
    ref = np.linalg.lstsq(X, y, rcond=None)[0]
    assert np.allclose(o["beta"], ref)
    e = tost_partial_r(0.0, 39); assert e["equivalent"]
    e = tost_partial_r(3.0, 39); assert not e["equivalent"]
    q = bh([0.01, 0.04, 0.03, 0.5]); assert np.allclose(q, [0.04, 0.04 * 4 / 3 * 0 + 0.0533333, 0.0533333, 0.5], atol=1e-4)
    print("ols / tost / bh ok")


if __name__ == "__main__":
    test_misc(); test_vs_statsmodels(); test_balanced_df(); test_fallback_singular(); test_type1_and_power()
    print("ALL LMM TESTS PASSED")
