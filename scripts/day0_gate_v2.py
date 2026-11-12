#!/usr/bin/env python
"""
A0-4 positive-control gate, v0.2.0: EXACT REPLICATION of the published day-0 classifier.

v0.1.0 (scripts/day0_gate.py) FAILED (CV AUC 0.573). v0.2.0 replicates the published method
(Bressan et al. 2023, Cell Genomics, STAR Methods "Prediction of neuronal differentiation
efficiency using bulk RNA-seq data at day 0"; verified verbatim from the Europe PMC full text):
  - "supervised machine learning (logistic regression) on the DESeq2 v1.26.0 normalised count
     expression values for genes at day 0"
  - outcome: fraction of cells in the 'Dopaminergic Neurons' cluster at day 65; "high (n = 62)
    and low (n = 21) ... first quartile value of cell percentages (Q25 = 15.7%)"
  - "full set of reliably expressed genes (threshold of inclusion mean normalised count >=50)"
  - "elastic net ... (alpha = 0.5) with a penalty lambda equal to 0.22"
  - "repeated (100 times) 5-fold cross-validation ... cv.glmnet"; published CV mean AUC 0.85 (SD 0.03)

PRE-SPECIFIED (frozen by git tag day0-gate-v0.2.0 before fitting)
  Lines     donor lines with day-0 RNAB and a Table S3 row, matched by CDI suffix (79 expected).
            DEVIATION from the paper (83 rows): the control line is excluded, because Table S3
            cannot link its 4 rows to its batch-specific day-0 samples, and it would leak
            across folds. One day-0 sample per line (highest version).
  Outcome   Table S3 'Dopaminergic Neurons' (%) >= 15.7 -> high (1), else low (0)  [published cutoff]
  Features  normalizedCountTable.tsv, day-0 columns of the analysis lines, all released rows;
            keep rows with mean normalised count >= 50 over those samples (unsupervised, outside CV).
  Model     binomial elastic net, glmnet objective
              -(1/n) sum[y*eta - log(1+exp(eta))] + lambda*((1-alpha)/2*||b||^2 + alpha*||b||_1),
            alpha 0.5, lambda 0.22, x standardised (glmnet default, 1/n SD) on each training set,
            unpenalised intercept, no class weights. Solved by FISTA, KKT-checked.
  CV        100 repeats x 5-fold, random (unstratified) folds as cv.glmnet; per repeat AUC =
            fold-size-weighted mean of per-fold AUCs (cv.glmnet type.measure="auc"); folds lacking
            a class are skipped. Report mean and SD over repeats. Secondary: pooled out-of-fold AUC.
  Null      200 label permutations, one CV repeat each; p = (1 + #{null >= observed}) / 201.
  Verdict   pass: mean CV AUC >= 0.75 and p < 0.05; partial: >= 0.65 and p < 0.05; else fail.
  Also reported (not gating):
    apparent (in-sample) AUC + number of non-zero coefficients (paper: 0.93, 10 genes) and
      whether the 5 published genes are among them;
    S2 the same model with outcome = OUR %TH+ (>=1 UMI) split at its own Q25 (bridges v0.1.0);
    S3 batch-only baseline (Table S1 BATCH) for the published outcome, same CV;
    S4 Spearman rho of the 5 published genes' day-0 normalised counts vs Table S3 DA %
       (published direction: HNRNPH3 +, SRSF5 +, HSD17B6 +, ZSWIM8 -, ARSA -), BH-FDR.

Prints aggregates only (n, AUCs, rho, p). No IDs, no per-line values.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrn_gate import num, read_table  # noqa: E402
from day0_gate import roc_auc_score  # noqa: E402  (numpy Mann-Whitney AUC)

GATE_VERSION = "0.2.0"
DA_CUTOFF = 15.7
MIN_MEAN_NORM = 50.0
ALPHA = 0.5
LAMBDA = 0.22
N_SPLITS = 5
N_REPEATS = 100
N_PERM = 200
SEED = 0
AUC_PASS, AUC_PARTIAL, P_MAX = 0.75, 0.65, 0.05
PUBLISHED_GENES = {  # symbol: (Ensembl ID, published direction)
    "HNRNPH3": ("ENSG00000096746", +1), "SRSF5": ("ENSG00000100650", +1),
    "HSD17B6": ("ENSG00000025423", +1), "ZSWIM8": ("ENSG00000214655", -1),
    "ARSA": ("ENSG00000100299", -1),
}


def enet_logistic(X: np.ndarray, y: np.ndarray, lam=LAMBDA, alpha=ALPHA, ws0=200, tol=1e-8):
    """glmnet-equivalent binomial elastic net at a single lambda. X standardised internally
    (1/n SD). Exact solution via a working set: solve on a subset of features, check KKT on
    all features, add violators, repeat. Returns (b0, b, mu, sd)."""
    n, p = X.shape
    mu = X.mean(0)
    sd = X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    l1, l2 = lam * alpha, lam * (1 - alpha)
    p0 = np.full(n, y.mean())
    g0 = np.abs(Z.T @ (p0 - y) / n)
    ws = np.argsort(g0)[::-1][: min(ws0, p)]
    b = np.zeros(p)
    b0 = float(np.log(y.mean() / (1 - y.mean())))
    for _ in range(50):
        b0, bw = _fista(Z[:, ws], y, l1, l2, tol, init=np.r_[b0, b[ws]])
        b = np.zeros(p)
        b[ws] = bw
        pr = 1 / (1 + np.exp(-(b0 + Z @ b)))
        g = Z.T @ (pr - y) / n
        viol = np.flatnonzero((b == 0) & (np.abs(g) > l1 * (1 + 1e-6)))
        if len(viol) == 0:
            break
        add = viol[np.argsort(np.abs(g[viol]))[::-1][:max(50, len(ws) // 2)]]
        ws = np.union1d(ws, add)
    return b0, b, mu, sd


def _fista(Z, y, l1, l2, tol=1e-8, max_iter=200000, init=None):
    n = len(y)
    # Lipschitz constant of the smooth part: 0.25 * sigma_max([1 Z])^2 / n + l2
    A = np.c_[np.ones(n), Z]
    L = 0.25 * np.linalg.norm(A @ A.T, 2) / n + l2
    t = 1.0 / L
    th = np.zeros(Z.shape[1] + 1) if init is None else init.astype(float).copy()
    zk, tk = th.copy(), 1.0

    def grad(v):
        eta = A @ v
        p = 1 / (1 + np.exp(-eta))
        g = A.T @ (p - y) / n
        g[1:] += l2 * v[1:]
        return g

    for _ in range(max_iter):
        g = grad(zk)
        new = zk - t * g
        new[1:] = np.sign(new[1:]) * np.maximum(np.abs(new[1:]) - t * l1, 0)
        if np.max(np.abs(new - th)) < tol:
            th = new
            break
        if g @ (new - th) > 0:  # adaptive restart (O'Donoghue & Candes 2015)
            tk = 1.0
        tk1 = (1 + np.sqrt(1 + 4 * tk * tk)) / 2
        zk = new + ((tk - 1) / tk1) * (new - th)
        th, tk = new, tk1
    return th[0], th[1:]


def kkt_gap(X, y, b0, b, mu, sd, lam=LAMBDA, alpha=ALPHA) -> float:
    Z = (X - mu) / sd
    p = 1 / (1 + np.exp(-(b0 + Z @ b)))
    g = Z.T @ (p - y) / len(y) + lam * (1 - alpha) * b
    l1 = lam * alpha
    nz = b != 0
    viol = np.r_[np.abs(g[nz] + l1 * np.sign(b[nz])), np.maximum(np.abs(g[~nz]) - l1, 0),
                 abs((p - y).mean())]
    return float(viol.max())


def predict(X, b0, b, mu, sd):
    return 1 / (1 + np.exp(-(b0 + ((X - mu) / sd) @ b)))


def folds(n, seed):
    return np.random.default_rng(seed).permutation(np.arange(n) % N_SPLITS)


def cv_once(X, y, seed):
    f = folds(len(y), seed)
    aucs, w, oof = [], [], np.full(len(y), np.nan)
    for k in range(N_SPLITS):
        te, tr = f == k, f != k
        if y[tr].min() == y[tr].max():
            continue
        b0, b, mu, sd = enet_logistic(X[tr], y[tr])
        oof[te] = predict(X[te], b0, b, mu, sd)
        if 0 < y[te].sum() < te.sum():
            aucs.append(roc_auc_score(y[te], oof[te]))
            w.append(te.sum())
    fold_auc = float(np.average(aucs, weights=w)) if aucs else np.nan
    m = ~np.isnan(oof)
    return fold_auc, roc_auc_score(y[m], oof[m])


def evaluate(X, y, perms=True):
    r = [cv_once(X, y, SEED + i) for i in range(N_REPEATS)]
    fa, pa = np.array([a for a, _ in r]), np.array([b for _, b in r])
    res = {"n": int(len(y)), "n_high": int(y.sum()), "n_low": int(len(y) - y.sum()),
           "cv_auc_mean": round(float(np.nanmean(fa)), 3), "cv_auc_sd": round(float(np.nanstd(fa)), 3),
           "pooled_oof_auc_mean": round(float(np.mean(pa)), 3)}
    if perms:
        rng = np.random.default_rng(SEED)
        null = np.array([cv_once(X, rng.permutation(y), SEED + 10_000 + k)[0] for k in range(N_PERM)])
        res["perm_p"] = round(float((1 + np.sum(null >= res["cv_auc_mean"])) / (1 + N_PERM)), 4)
        res["perm_null_mean"] = round(float(np.nanmean(null)), 3)
    return res


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--norm", type=Path, required=True, help="normalizedCountTable.tsv")
    ap.add_argument("--scrn-summary", type=Path, required=True)
    ap.add_argument("--table-s1", type=Path, required=True)
    ap.add_argument("--table-s3", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    # day-0 samples: donor lines only (no batch suffix), highest version per line
    cols = pd.read_csv(a.norm, sep="\t", nrows=0).columns
    meta = pd.Series(cols[1:]).str.extract(r"^RNAB_(PPMI\d+)(B\w+?)?_(\d+)_da0_v(\d+)$")
    meta.columns = ["patno", "bsfx", "cdi", "ver"]
    meta["col"] = cols[1:]
    meta = meta.dropna(subset=["patno"])
    meta = meta[meta["bsfx"].isna()]
    meta["ver"] = meta["ver"].astype(int)
    meta = meta.sort_values("ver").groupby("patno").tail(1)

    # control line = the PATNO with several scRNA runs; exclude
    s = pd.read_csv(a.scrn_summary)
    s["patno"] = s["run"].str.extract(r"^SCRN_(PPMI\d+)")[0]
    ctrl = set(s.groupby("patno").size()[lambda x: x > 1].index)
    meta = meta[~meta["patno"].isin(ctrl)]

    s3 = read_table(a.table_s3)
    s3["cdi"] = s3["FOUNDIN_ID"].astype(str).str.extract(r"(\d+)")[0]
    s3["da"] = num(s3["Dopaminergic Neurons"])
    s3 = s3.groupby("cdi").filter(lambda g: len(g) == 1).set_index("cdi")
    meta = meta[meta["cdi"].isin(s3.index)].set_index("patno")

    norm = pd.read_csv(a.norm, sep="\t", index_col=0, usecols=[cols[0]] + meta["col"].tolist())
    norm = norm[meta["col"]]
    keep = norm.mean(axis=1) >= MIN_MEAN_NORM
    Xdf = norm[keep]
    X = Xdf.T.to_numpy(float)
    da = s3.loc[meta["cdi"], "da"].to_numpy(float)
    y = (da >= DA_CUTOFF).astype(float)

    res = {"gate_version": GATE_VERSION, "lines": int(len(y)), "rows_total": int(len(norm)),
           "rows_kept_mean_ge_50": int(keep.sum()),
           "params": {"da_cutoff": DA_CUTOFF, "min_mean_norm": MIN_MEAN_NORM, "alpha": ALPHA,
                      "lambda": LAMBDA, "n_splits": N_SPLITS, "n_repeats": N_REPEATS, "n_perm": N_PERM,
                      "seed": SEED, "auc_pass": AUC_PASS, "auc_partial": AUC_PARTIAL, "p_max": P_MAX}}

    # apparent fit on all lines
    b0, b, mu, sd = enet_logistic(X, y)
    nz = np.flatnonzero(b)
    ids = Xdf.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    res["apparent"] = {"auc": round(roc_auc_score(y, predict(X, b0, b, mu, sd)), 3),
                       "n_nonzero": int(len(nz)), "kkt_max_violation": kkt_gap(X, y, b0, b, mu, sd),
                       "published_genes_selected": [g for g, (e, _) in PUBLISHED_GENES.items()
                                                    if e in set(ids[nz])]}

    res["primary_published_outcome"] = evaluate(X, y)

    # S2: our %TH+ (>=1 UMI), Q25 split (R type-7 quantile = numpy 'linear')
    s1run = s.set_index("patno")
    th = 100 * s1run.loc[meta.index, "frac_TH_pos_all"].to_numpy(float)
    y2 = (th >= np.percentile(th, 25)).astype(float)
    res["S2_our_TH_Q25"] = evaluate(X, y2, perms=False)

    # S3: batch-only baseline (one-hot), same enet + CV
    t1 = read_table(a.table_s1)
    t1["patno"] = t1["PATNO"].astype(str).str.extract(r"(PPMI\d+)")[0]
    batch = t1.drop_duplicates("patno").set_index("patno")["BATCH"].astype(str).reindex(meta.index)
    B = pd.get_dummies(batch).to_numpy(float)
    res["S3_batch_only"] = evaluate(B, y, perms=False)
    res["S3_batch_only"]["n_batches"] = int(B.shape[1])

    # S4: published genes, Spearman vs DA %
    rows = []
    for g, (e, sign) in PUBLISHED_GENES.items():
        hit = np.flatnonzero(ids == e)
        if len(hit) == 0:
            rows.append((g, None, None, sign))
            continue
        v = norm.iloc[np.flatnonzero(norm.index.astype(str).str.replace(r"\.\d+$", "", regex=True) == e)[0]]
        rho, p = spearmanr(v.to_numpy(float), da)
        rows.append((g, float(rho), float(p), sign))
    ps = [r[2] for r in rows if r[2] is not None]
    qs = iter(bh(ps)) if ps else iter([])
    res["S4_published_genes"] = [
        {"gene": g, "rho": None if r is None else round(r, 3), "p": None if p is None else float(f"{p:.3g}"),
         "fdr": None if p is None else float(f"{next(qs):.3g}"), "published_sign": sgn,
         "same_sign": None if r is None else bool(np.sign(r) == sgn)}
        for g, r, p, sgn in rows]

    pr = res["primary_published_outcome"]
    ok = pr["perm_p"] < P_MAX
    res["verdict"] = ("pass" if pr["cv_auc_mean"] >= AUC_PASS and ok
                      else "partial" if pr["cv_auc_mean"] >= AUC_PARTIAL and ok else "fail")
    txt = json.dumps(res, indent=2)
    print(txt)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(txt)


if __name__ == "__main__":
    main()
