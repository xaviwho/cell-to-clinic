"""Statistics for the pre-registered models (OSF osf.io/s6w9c, Part 5), numpy/scipy only.

statsmodels.api cannot be imported on the analysis machine (Windows Application Control blocks
one of its DLLs), and statsmodels has no Satterthwaite degrees of freedom anyway. So this module
implements, following lme4/lmerTest:

- lmm_fit: linear mixed model fitted by REML. Per-donor random effects use the relative-covariance
  Cholesky parametrisation of lme4 (theta), with sigma^2 profiled out.
    re="slope_corr": random intercept + slope, correlated (registered model)
    re="slope_diag": random intercept + slope, uncorrelated (fallback 1)
    re="intercept" : random intercept only (fallback 2)
- lmm_fit_fallback: applies the registered fallback order. "Fails to converge" is defined here:
  the optimiser reports failure, the fit is singular (a random-effect SD < 1e-4 * residual SD, or
  |correlation| > 0.999, as in lme4::isSingular), or the deviance Hessian over the variance
  parameters is not positive definite.
- satterthwaite: t test for one contrast L'beta, with Satterthwaite df as in lmerTest. The
  asymptotic covariance of (theta, sigma) is 2 * inverse Hessian of the REML deviance; gradients
  and Hessians are by central finite differences.
- ols: ordinary least squares t test (H3).
- partial_r, tost_partial_r: effect size and equivalence test. The 90% CI of beta (t distribution,
  same df) maps monotonically to partial r via r = t / sqrt(t^2 + df). Equivalence means both
  ends lie inside +/- bound.
- bh: Benjamini-Hochberg q-values.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, stats

SINGULAR_REL_SD = 1e-4
SINGULAR_CORR = 0.999
RE_TYPES = ("slope_corr", "slope_diag", "intercept")


def _lambda(theta, re):
    if re == "slope_corr":
        return np.array([[theta[0], 0.0], [theta[1], theta[2]]])
    if re == "slope_diag":
        return np.array([[theta[0], 0.0], [0.0, theta[1]]])
    return np.array([[theta[0]]])


def _n_theta(re):
    return {"slope_corr": 3, "slope_diag": 2, "intercept": 1}[re]


def _bounds(re):
    return {"slope_corr": [(0, None), (None, None), (0, None)],
            "slope_diag": [(0, None), (0, None)],
            "intercept": [(0, None)]}[re]


class _Data:
    def __init__(self, y, X, groups, t, re):
        y = np.asarray(y, float); X = np.asarray(X, float); t = np.asarray(t, float)
        g = np.asarray(groups)
        ok = np.isfinite(y) & np.all(np.isfinite(X), 1) & np.isfinite(t)
        y, X, t, g = y[ok], X[ok], t[ok], g[ok]
        self.levels, idx = np.unique(g, return_inverse=True)
        self.blocks = []
        for k in range(len(self.levels)):
            m = idx == k
            Z = np.c_[np.ones(m.sum()), t[m]] if re != "intercept" else np.ones((m.sum(), 1))
            self.blocks.append((y[m], X[m], Z))
        self.N, self.p = X.shape
        self.re = re


def _pieces(theta, d: _Data):
    """Sum_i logdet(M_i), X'M^-1X, X'M^-1y, y'M^-1y with M_i = Z Lam Lam' Z' + I."""
    Lam = _lambda(theta, d.re)
    D = Lam @ Lam.T
    ld = 0.0
    XtX = np.zeros((d.p, d.p)); Xty = np.zeros(d.p); yty = 0.0
    for y, X, Z in d.blocks:
        M = Z @ D @ Z.T + np.eye(len(y))
        c = np.linalg.cholesky(M)
        ld += 2 * np.log(np.diag(c)).sum()
        Xs = np.linalg.solve(c, X); ys = np.linalg.solve(c, y)
        XtX += Xs.T @ Xs; Xty += Xs.T @ ys; yty += ys @ ys
    return ld, XtX, Xty, yty


def _profiled_dev(theta, d):
    if not np.all(np.isfinite(theta)) or np.max(np.abs(theta)) > 1e6:
        return 1e30
    try:
        ld, XtX, Xty, yty = _pieces(theta, d)
        cx = np.linalg.cholesky(XtX)
    except np.linalg.LinAlgError:
        return 1e30
    beta = np.linalg.solve(XtX, Xty)
    q = yty - beta @ Xty
    nmp = d.N - d.p
    if q <= 0:
        return 1e30
    return ld + 2 * np.log(np.diag(cx)).sum() + nmp * (1 + np.log(2 * np.pi * q / nmp))


def _full_dev(phi, d):
    """REML deviance in (theta, sigma), beta profiled by GLS (used for the Satterthwaite Hessian)."""
    theta, sigma = phi[:-1], phi[-1]
    ld, XtX, Xty, yty = _pieces(theta, d)
    beta = np.linalg.solve(XtX, Xty)
    q = yty - beta @ Xty
    s2 = sigma ** 2
    _, ldx = np.linalg.slogdet(XtX)
    return (d.N - d.p) * np.log(s2) + ld + ldx + q / s2 + (d.N - d.p) * np.log(2 * np.pi)


def _num_grad(f, x, h=None):
    x = np.asarray(x, float)
    g = np.zeros_like(x)
    for i in range(len(x)):
        hi = (h if h is not None else 1e-5 * max(1.0, abs(x[i])))
        e = np.zeros_like(x); e[i] = hi
        g[i] = (f(x + e) - f(x - e)) / (2 * hi)
    return g


def _num_hess(f, x):
    x = np.asarray(x, float)
    n = len(x)
    H = np.zeros((n, n))
    hs = [1e-4 * max(1.0, abs(v)) for v in x]
    f0 = f(x)
    for i in range(n):
        ei = np.zeros(n); ei[i] = hs[i]
        H[i, i] = (f(x + ei) - 2 * f0 + f(x - ei)) / hs[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n); ej[j] = hs[j]
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * hs[i] * hs[j])
    return H


def lmm_fit(y, X, groups, t, re="slope_corr", starts=None):
    """REML fit. X must contain the intercept column. t is the random-slope variable."""
    d = _Data(y, X, groups, t, re)
    if starts is None:
        starts = {"slope_corr": [[1.0, 0.0, 0.5], [0.5, 0.0, 0.1], [2.0, 0.0, 1.0]],
                  "slope_diag": [[1.0, 0.5], [0.5, 0.1], [2.0, 1.0]],
                  "intercept": [[1.0], [0.3], [3.0]]}[re]
    best = None
    for s0 in starts:
        r = optimize.minimize(_profiled_dev, np.asarray(s0, float), args=(d,), method="L-BFGS-B",
                              bounds=_bounds(re), options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    theta = best.x
    ld, XtX, Xty, yty = _pieces(theta, d)
    beta = np.linalg.solve(XtX, Xty)
    q = yty - beta @ Xty
    sigma = np.sqrt(q / (d.N - d.p))
    Lam = _lambda(theta, re)
    D = sigma ** 2 * Lam @ Lam.T
    sds = np.sqrt(np.diag(D))
    corr = D[0, 1] / (sds[0] * sds[1]) if (re == "slope_corr" and sds.min() > 0) else 0.0
    singular = bool(np.any(sds < SINGULAR_REL_SD * sigma) or abs(corr) > SINGULAR_CORR)
    phi = np.r_[theta, sigma]
    H = _num_hess(lambda p: _full_dev(p, d), phi)
    try:
        hess_pd = bool(np.all(np.linalg.eigvalsh((H + H.T) / 2) > 0))
    except np.linalg.LinAlgError:
        hess_pd = False
    converged = bool(best.success) and not singular and hess_pd
    return {"re": re, "beta": beta, "cov_beta": sigma ** 2 * np.linalg.inv(XtX), "theta": theta,
            "sigma": float(sigma), "re_cov": D, "re_corr": float(corr), "singular": singular,
            "hess_pd": hess_pd, "opt_success": bool(best.success), "converged": converged,
            "reml_dev": float(best.fun), "n_obs": d.N, "n_groups": len(d.levels), "_d": d, "_phi": phi,
            "_H": H}


def lmm_fit_fallback(y, X, groups, t):
    """Registered fallback order: correlated slope -> uncorrelated slope -> intercept only."""
    tried = []
    for re in RE_TYPES:
        f = lmm_fit(y, X, groups, t, re)
        tried.append({"re": re, "converged": f["converged"], "singular": f["singular"],
                      "hess_pd": f["hess_pd"], "opt_success": f["opt_success"]})
        if f["converged"] or re == "intercept":
            f["fallback_log"] = tried
            return f


def satterthwaite(fit, L):
    """t test of L'beta with Satterthwaite df (lmerTest method)."""
    L = np.asarray(L, float)
    d, phi = fit["_d"], fit["_phi"]
    est = float(L @ fit["beta"])

    def v_of(p):
        _, XtX, _, _ = _pieces(p[:-1], d)
        return float(p[-1] ** 2 * L @ np.linalg.solve(XtX, L))

    v = v_of(phi)
    g = _num_grad(v_of, phi)
    H = fit["_H"]
    try:
        A = 2 * np.linalg.inv(H)
        denom = float(g @ A @ g)
        df = 2 * v ** 2 / denom if denom > 0 else np.nan
    except np.linalg.LinAlgError:
        df = np.nan
    if not np.isfinite(df):  # boundary fit: fall back to residual-type df, flagged
        df, df_method = float(fit["n_groups"] - np.linalg.matrix_rank(np.atleast_2d(L))), "fallback_groups"
    else:
        df_method = "satterthwaite"
    se = np.sqrt(v)
    tval = est / se
    p = 2 * stats.t.sf(abs(tval), df)
    return {"estimate": est, "se": float(se), "t": float(tval), "df": float(df), "p": float(p),
            "df_method": df_method}


def ols(y, X):
    y = np.asarray(y, float); X = np.asarray(X, float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), 1)
    y, X = y[ok], X[ok]
    n, p = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    r = y - X @ beta
    s2 = r @ r / (n - p)
    se = np.sqrt(np.diag(XtX_inv) * s2)
    tv = beta / se
    return {"beta": beta, "se": se, "t": tv, "df": float(n - p), "p": 2 * stats.t.sf(np.abs(tv), n - p), "n": n}


def partial_r(t, df):
    return float(t / np.sqrt(t ** 2 + df))


def tost_partial_r(t, df, bound=0.40, level=0.90):
    """90% CI of beta mapped to partial r; equivalent if the whole CI lies inside +/- bound."""
    q = stats.t.ppf(1 - (1 - level) / 2, df)
    lo, hi = partial_r(t - q, df), partial_r(t + q, df)
    return {"r": partial_r(t, df), "ci90": [lo, hi], "equivalent": bool(lo > -bound and hi < bound)}


def classify(q, tost):
    if q < 0.05:
        return "supported"
    if tost["equivalent"]:
        return "equivalent"
    return "inconclusive"


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)
