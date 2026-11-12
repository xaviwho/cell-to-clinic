#!/usr/bin/env python
"""
Positive-control gate: does day-0 (iPSC) bulk RNA predict day-65 DA yield (%TH+)?

Benchmark: Bressan et al. 2023 report day-0 expression predicts DA yield
(AUC 0.93; CV mean 0.85). This gate re-tests it with OUR per-line %TH+.

PRE-SPECIFIED (frozen by git tag before fitting; any change must be logged in ACTIVITY.md)
  Lines     donor lines with both day-0 RNAB and day-65 scRNA. The control line (the
            only line with several scRNA runs) is excluded: it is one donor and would
            leak across folds. One day-0 sample per line; if several versions exist,
            the highest version (README: v2 = repeat after a failed first attempt).
  Outcome   primary: high vs low DA yield = per-line 100*frac_TH_pos_all (>=1 UMI)
            >= median of the analysis lines (ties go high).
            sensitivity: same with frac_TH_pos_all_sens (>=2 UMI).
  Features  day-0 VST (vstTable.tsv). Within each TRAINING fold: keep the N_TOP_GENES
            most variable genes, then z-score (fold-fitted).
  Model     L2 logistic regression, C tuned by inner stratified 5-fold CV (AUC) over
            C_GRID; class_weight balanced.
  CV        outer stratified 5-fold, N_REPEATS repeats, seed SEED. Per repeat, AUC of
            the pooled out-of-fold probabilities; report the mean (and SD) over repeats.
  Null      N_PERM label permutations, each scored with one CV repeat;
            p = (1 + #{perm AUC >= observed mean AUC}) / (1 + N_PERM).
  Baseline  batch-only logistic model (one-hot BATCH, public Table S1), same CV.
            Reported, not gating (batch/media is a known driver of DA yield).
  Continuous (secondary) ridge regression on logit(%TH+) with the same features;
            Spearman rho of pooled out-of-fold predictions vs observed, mean over repeats.
  Verdict   pass:    mean AUC >= AUC_PASS and p < 0.05
            partial: mean AUC >= AUC_PARTIAL and p < 0.05
            fail:    otherwise

Prints aggregates only (n, AUCs, rho, p). No IDs, no per-line values.

Usage:
  python -I scripts/day0_gate.py --vst <vstTable.tsv> --scrn-summary data/derived/scrn/summary.csv \
         --table-s1 <mmc2.xlsx> --out data/derived/day0_gate.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import rankdata

# scikit-learn is NOT used: Windows Application Control on the analysis machine blocks
# some of its compiled modules. Minimal numpy/scipy equivalents follow (same conventions).


def roc_auc_score(y: np.ndarray, s: np.ndarray) -> float:
    """AUC via Mann-Whitney U (ties count 1/2)."""
    y = np.asarray(y).astype(bool)
    r = rankdata(s)
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


class StratifiedKFold:
    def __init__(self, n_splits, shuffle=True, random_state=0):
        self.k, self.seed = n_splits, random_state

    def split(self, X, y):
        rng = np.random.default_rng(self.seed)
        fold = np.empty(len(y), dtype=int)
        for cls in np.unique(y):
            idx = rng.permutation(np.flatnonzero(y == cls))
            fold[idx] = (np.arange(len(idx)) + rng.integers(self.k)) % self.k
        for f in range(self.k):
            yield np.flatnonzero(fold != f), np.flatnonzero(fold == f)


class KFold:
    def __init__(self, n_splits, shuffle=True, random_state=0):
        self.k, self.seed = n_splits, random_state

    def split(self, X):
        n = len(X)
        fold = np.random.default_rng(self.seed).permutation(np.arange(n) % self.k)
        for f in range(self.k):
            yield np.flatnonzero(fold != f), np.flatnonzero(fold == f)


class LogisticRegression:
    """L2 logistic regression, sklearn convention: min 0.5||w||^2 + C*sum_i s_i*logloss_i,
    unpenalised intercept, class_weight='balanced' -> s_i = n/(2*n_class)."""

    def __init__(self, C=1.0, class_weight="balanced", max_iter=5000):
        self.C, self.max_iter = C, max_iter

    def fit(self, X, y):
        y = np.asarray(y, dtype=float)
        n, p = X.shape
        n1 = y.sum()
        sw = np.where(y == 1, n / (2 * n1), n / (2 * (n - n1)))
        ys = 2 * y - 1

        def f(theta):
            w, b = theta[:p], theta[p]
            z = ys * (X @ w + b)
            loss = np.logaddexp(0, -z)
            g = -ys * sw * self.C * (1 / (1 + np.exp(z)))
            return 0.5 * w @ w + self.C * (sw * loss).sum(), np.r_[w + X.T @ g, g.sum()]

        r = minimize(f, np.zeros(p + 1), jac=True, method="L-BFGS-B",
                     options={"maxiter": self.max_iter})
        self.w, self.b = r.x[:p], r.x[p]
        return self

    def predict_proba(self, X):
        p1 = 1 / (1 + np.exp(-(X @ self.w + self.b)))
        return np.c_[1 - p1, p1]


class Ridge:
    """Ridge with unpenalised intercept (centred), dual closed form (fast for n << p)."""

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, X, t):
        self.mx, self.mt = X.mean(0), t.mean()
        Xc = X - self.mx
        K = Xc @ Xc.T
        self.w = Xc.T @ np.linalg.solve(K + self.alpha * np.eye(len(t)), t - self.mt)
        return self

    def predict(self, X):
        return (X - self.mx) @ self.w + self.mt


sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrn_gate import read_table  # noqa: E402

GATE_VERSION = "0.1.0"
N_TOP_GENES = 1000
C_GRID = [1e-3, 1e-2, 1e-1, 1.0]
RIDGE_ALPHAS = [1e1, 1e2, 1e3, 1e4]
N_SPLITS = 5
N_REPEATS = 20
N_PERM = 200
SEED = 0
AUC_PASS = 0.75
AUC_PARTIAL = 0.65
P_MAX = 0.05

warnings.filterwarnings("ignore", category=UserWarning)


def top_var(Xtr: np.ndarray, k: int) -> np.ndarray:
    return np.argsort(Xtr.var(axis=0))[::-1][:k]


def zfit(Xtr: np.ndarray):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    return mu, sd


def fit_logit(Xtr, ytr, seed):
    """Inner-CV tuned L2 logistic regression; returns predict_proba function."""
    inner = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
    best_c, best = C_GRID[0], -np.inf
    for c in C_GRID:
        oof = np.zeros(len(ytr))
        for itr, ite in inner.split(Xtr, ytr):
            m = LogisticRegression(C=c, class_weight="balanced", max_iter=5000)
            m.fit(Xtr[itr], ytr[itr])
            oof[ite] = m.predict_proba(Xtr[ite])[:, 1]
        auc = roc_auc_score(ytr, oof)
        if auc > best:
            best, best_c = auc, c
    m = LogisticRegression(C=best_c, class_weight="balanced", max_iter=5000).fit(Xtr, ytr)
    return m


def cv_auc(X: np.ndarray, y: np.ndarray, seed: int, select: bool = True) -> float:
    outer = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in outer.split(X, y):
        Xtr, Xte = X[tr], X[te]
        if select:
            idx = top_var(Xtr, N_TOP_GENES)
            Xtr, Xte = Xtr[:, idx], Xte[:, idx]
        mu, sd = zfit(Xtr)
        m = fit_logit((Xtr - mu) / sd, y[tr], seed)
        oof[te] = m.predict_proba((Xte - mu) / sd)[:, 1]
    return roc_auc_score(y, oof)


def cv_rho(X: np.ndarray, t: np.ndarray, seed: int) -> float:
    outer = KFold(N_SPLITS, shuffle=True, random_state=seed)
    oof = np.zeros(len(t))
    for tr, te in outer.split(X):
        idx = top_var(X[tr], N_TOP_GENES)
        Xtr, Xte = X[tr][:, idx], X[te][:, idx]
        mu, sd = zfit(Xtr)
        Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
        inner = KFold(N_SPLITS, shuffle=True, random_state=seed)
        best_a, best = RIDGE_ALPHAS[0], np.inf
        for a in RIDGE_ALPHAS:
            err = 0.0
            for itr, ite in inner.split(Xtr):
                r = Ridge(alpha=a).fit(Xtr[itr], t[tr][itr])
                err += ((r.predict(Xtr[ite]) - t[tr][ite]) ** 2).sum()
            if err < best:
                best, best_a = err, a
        oof[te] = Ridge(alpha=best_a).fit(Xtr, t[tr]).predict(Xte)
    return float(pd.Series(oof).corr(pd.Series(t), method="spearman"))


def evaluate(X, y, label, seed0=SEED, perms=True):
    aucs = [cv_auc(X, y, seed0 + r) for r in range(N_REPEATS)]
    res = {"outcome": label, "n": int(len(y)), "n_high": int(y.sum()),
           "auc_mean": round(float(np.mean(aucs)), 3), "auc_sd": round(float(np.std(aucs)), 3)}
    if perms:
        rng = np.random.default_rng(seed0)
        null = [cv_auc(X, rng.permutation(y), seed0 + 1000 + k) for k in range(N_PERM)]
        res["perm_p"] = round((1 + sum(a >= res["auc_mean"] for a in null)) / (1 + N_PERM), 4)
        res["perm_null_auc_mean"] = round(float(np.mean(null)), 3)
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--vst", type=Path, required=True)
    p.add_argument("--scrn-summary", type=Path, required=True)
    p.add_argument("--table-s1", type=Path, required=True)
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args()

    # scRNA outcome per donor line (control line = the one with several runs -> excluded)
    s = pd.read_csv(a.scrn_summary)
    s["patno"] = s["run"].str.extract(r"^SCRN_(PPMI\d+)")[0]
    n_runs = s.groupby("patno").size()
    s = s[s["patno"].map(n_runs) == 1].set_index("patno")

    # day-0 VST, one sample per line (highest version)
    vst = pd.read_csv(a.vst, sep="\t", index_col=0)
    meta = pd.Series(vst.columns).str.extract(r"^RNAB_(PPMI\d+)(B\w+?)?_\d+_da0_v(\d+)$")
    meta.columns = ["patno", "batchsfx", "ver"]
    meta["col"] = vst.columns
    meta = meta.dropna(subset=["patno"])
    meta = meta[meta["batchsfx"].isna()]
    meta["ver"] = meta["ver"].astype(int)
    meta = meta.sort_values("ver").groupby("patno").tail(1).set_index("patno")

    lines = sorted(set(s.index) & set(meta.index))
    X = vst[meta.loc[lines, "col"]].T.to_numpy(dtype=float)
    X = X[:, np.isfinite(X).all(0)]

    th1 = 100 * s.loc[lines, "frac_TH_pos_all"].to_numpy()
    th2 = 100 * s.loc[lines, "frac_TH_pos_all_sens"].to_numpy()
    y1 = (th1 >= np.median(th1)).astype(int)
    y2 = (th2 >= np.median(th2)).astype(int)

    res = {"gate_version": GATE_VERSION, "lines": len(lines), "genes_total": int(X.shape[1]),
           "params": {"n_top_genes": N_TOP_GENES, "C_grid": C_GRID, "ridge_alphas": RIDGE_ALPHAS,
                      "n_splits": N_SPLITS, "n_repeats": N_REPEATS, "n_perm": N_PERM, "seed": SEED,
                      "auc_pass": AUC_PASS, "auc_partial": AUC_PARTIAL, "p_max": P_MAX}}
    res["primary"] = evaluate(X, y1, "TH>=1UMI median split")
    res["sensitivity_2umi"] = evaluate(X, y2, "TH>=2UMI median split")

    # batch-only baseline
    s1 = read_table(a.table_s1)
    s1["patno"] = s1["PATNO"].astype(str).str.extract(r"(PPMI\d+)")[0]
    batch = s1.drop_duplicates("patno").set_index("patno")["BATCH"].astype(str)
    B = pd.get_dummies(batch.reindex(lines)).to_numpy(dtype=float)
    res["batch_only_baseline"] = {
        "auc_mean": round(float(np.mean([
            roc_auc_score(y1, _oof_simple(B, y1, SEED + r)) for r in range(N_REPEATS)])), 3),
        "n_batches": int(B.shape[1])}

    # continuous secondary
    t = np.log((th1 / 100 + 1e-3) / (1 - th1 / 100 + 1e-3))
    res["continuous_spearman_rho_mean"] = round(float(np.mean([cv_rho(X, t, SEED + r) for r in range(N_REPEATS)])), 3)

    pr = res["primary"]
    ok_p = pr["perm_p"] < P_MAX
    res["verdict"] = ("pass" if pr["auc_mean"] >= AUC_PASS and ok_p
                      else "partial" if pr["auc_mean"] >= AUC_PARTIAL and ok_p else "fail")
    txt = json.dumps(res, indent=2)
    print(txt)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(txt)


def _oof_simple(B: np.ndarray, y: np.ndarray, seed: int) -> np.ndarray:
    outer = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in outer.split(B, y):
        m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000).fit(B[tr], y[tr])
        oof[te] = m.predict_proba(B[te])[:, 1]
    return oof


if __name__ == "__main__":
    main()
