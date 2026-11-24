"""A3: split Jerber across-batch within-line variance into cell-sampling noise and true batch variance.
EXPLORATORY, public data only. Tag jerber-a3-v0.1.0 before running.

For each DA-only phenotype of A2 (P2da SNCA, P3da OXPHOS, P4da lysosome; samples with >= 50 DA cells,
identical scaling to A2) and for a random set of 1,000 expressed genes:
  v_cell(s)  = bootstrap variance (B = 200, cells resampled with replacement within sample s)
  c          = median over samples of v_cell(s) * n_s       (noise per 1 / DA cell)
  S_w^2, S_d^2 from A2's split (one sample per line; replicated lines in >= 2 pools)
  sampling share of S_w^2 = mean over replicated samples of c / n_s, divided by S_w^2
  sigma_batch^2 = S_w^2 - mean_rep(c / n_s)   (floored at 0)
  sigma_line^2  = S_d^2 - sigma_batch^2 - mean_donor(c / n_s)
  ICC(n)        = sigma_line^2 / (sigma_line^2 + sigma_batch^2 + c / n) for n DA cells per run.
P1 (logit DA fraction of all cells) uses the analytic binomial variance 1 / (n p (1 - p)).
Output: data/derived/jerber/a3_noise.json (aggregates only).
"""
import json, subprocess, sys
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from jerber_inspect import read_col           # noqa: E402
import jerber_reliability as jr               # noqa: E402

H5 = "D:/ST-DT4/data/external/jerber/x/day52.h5"
OUT = ROOT / "data/derived/jerber"
B, MIN_DA, NGENE, SEED = 200, 50, 1000, 0
NS = (50, 100, 200, 500, 1000, 1e9)


def collect(sel_genes, keys):
    """per-DA-cell counts for sel_genes, total counts, sample index; untreated DA cells of `keys`."""
    f = h5py.File(H5, "r"); g = f["obs"]
    donor, pool, treat, ct = (read_col(g, k) for k in ("donor_id", "pool_id", "treatment", "celltype"))
    kid = {k: i for i, k in enumerate(keys)}
    sid = np.array([kid.get((d, p), -1) if (t == "NONE" and c == "DA") else -1
                    for d, p, t, c in zip(donor, pool, treat, ct)])
    X = f["raw/X"]; ip = X["indptr"][()]; n = len(f["raw/var/index"])
    rows, tots, S = [], [], []
    for r0 in range(0, len(sid), 20000):
        r1 = min(r0 + 20000, len(sid)); idx = np.flatnonzero(sid[r0:r1] >= 0)
        if not len(idx):
            continue
        a, b = ip[r0], ip[r1]
        M = sparse.csr_matrix((X["data"][a:b], X["indices"][a:b], ip[r0:r1 + 1] - a), shape=(r1 - r0, n))[idx]
        tots.append(np.asarray(M.sum(1)).ravel()); rows.append(M[:, sel_genes].tocsr()); S.append(sid[r0 + idx])
        print(f"\r{r1}/{len(sid)}", end="", flush=True)
    print()
    return sparse.vstack(rows).tocsr(), np.concatenate(tots), np.concatenate(S)


def decompose(y_boot_var, n_cells, values, lines, pools):
    dv_idx, reps = jr.split(lines, pools, np.arange(len(lines)))
    c = np.median(y_boot_var * n_cells)
    sd2 = np.var(values[dv_idx], ddof=1)
    ss = sum(((values[r] - values[r].mean()) ** 2).sum() for r in reps); df = sum(len(r) - 1 for r in reps)
    sw2 = ss / df
    rep_idx = np.concatenate(reps)
    samp_rep = np.mean(c / n_cells[rep_idx]); samp_don = np.mean(c / n_cells[dv_idx])
    sb2 = max(sw2 - samp_rep, 0.0)
    sl2 = sd2 - sb2 - samp_don
    return {"S_d2": sd2, "S_w2": sw2, "sampling_share_of_S_w2": min(samp_rep / sw2, 1.0),
            "sigma_batch2": sb2, "sigma_line2": sl2, "c_per_cell": c,
            "median_cells_replicated": float(np.median(n_cells[rep_idx])),
            "icc_by_cells": {str(int(m)) if m < 1e8 else "inf": float(sl2 / (sl2 + sb2 + c / m)) for m in NS}}


def main():
    pb = dict(np.load(OUT / "pseudobulk.npz"))
    D, genes, lines, pools, comp, cts = pb["D"], pb["genes"], pb["lines"], pb["pools"], pb["comp"], list(pb["celltypes"])
    nda = comp[:, cts.index("DA")]; ncell = comp.sum(1)
    s = nda >= MIN_DA
    Ld = jr.log2cpm(D[s]); kd = (D[s] / D[s].sum(1, keepdims=True) * 1e6).mean(0) >= 10
    cfg = json.load(open(ROOT / "config/phenotypes.json"))
    gsets = {k: [x.strip() for x in open(ROOT / v) if x.strip() and not x.startswith("#")] for k, v in cfg["gene_sets"].items()}
    m3, m4 = kd & np.isin(genes, gsets["P3"]), kd & np.isin(genes, gsets["P4"])
    rng = np.random.default_rng(SEED)
    rand = rng.choice(np.flatnonzero(kd & ~m3 & ~m4 & (genes != "SNCA")), NGENE, replace=False)
    sel = np.unique(np.r_[np.flatnonzero(genes == "SNCA"), np.flatnonzero(m3), np.flatnonzero(m4), rand])
    pos = {g: i for i, g in enumerate(sel)}
    mu, sd = Ld.mean(0), Ld.std(0, ddof=1)           # A2 scaling over the samples in s
    keys = list(zip(lines[s], pools[s]))
    Xc, tot, sid = collect(sel, keys)

    def pheno(cnt, lib):                               # cnt: (..., len(sel)) counts; lib: library size
        L = np.log2(cnt / lib[..., None] * 1e6 + 1)
        full = lambda m: (L[..., [pos[i] for i in np.flatnonzero(m)]] - mu[m]) / sd[m]
        return {"P2da": L[..., pos[np.flatnonzero(genes == "SNCA")[0]]],
                "P3da": full(m3).mean(-1), "P4da": full(m4).mean(-1),
                "genes": L[..., [pos[i] for i in rand]]}

    nS = len(keys); boot = {k: [] for k in ("P2da", "P3da", "P4da", "genes")}
    point = {k: [] for k in boot}; ncells = np.zeros(nS)
    for i in range(nS):
        rows = np.flatnonzero(sid == i); ncells[i] = len(rows)
        Xi = Xc[rows]; ti = tot[rows]
        W = rng.multinomial(len(rows), np.full(len(rows), 1 / len(rows)), size=B)
        pb_b = W @ Xi.toarray()
        ph = pheno(pb_b, W @ ti); p0 = pheno(np.asarray(Xi.sum(0)).ravel(), np.array(ti.sum()))
        for k in boot:
            boot[k].append(np.var(ph[k], axis=0, ddof=1)); point[k].append(p0[k])
    ls, ps = lines[s], pools[s]
    res = {"n_samples": nS, "B": B, "min_da": MIN_DA}
    for k in ("P2da", "P3da", "P4da"):
        res[k] = decompose(np.array(boot[k]), ncells, np.array(point[k]), ls, ps)
    G = [decompose(np.array(boot["genes"])[:, j], ncells, np.array(point["genes"])[:, j], ls, ps) for j in range(NGENE)]
    q = lambda key: np.nanpercentile([x[key] for x in G], [25, 50, 75]).round(3).tolist()
    res["genes_random1000"] = {"sampling_share_of_S_w2_q25_50_75": q("sampling_share_of_S_w2"),
                               "icc_by_cells_median": {m: float(np.nanmedian([x["icc_by_cells"][m] for x in G])) for m in G[0]["icc_by_cells"]}}
    # P1 analytic binomial
    sa = ncell >= 200; p = nda[sa] / ncell[sa]; pc = np.clip(p, 0.005, 0.995)
    res["P1"] = decompose(1 / (ncell[sa] * pc * (1 - pc)), ncell[sa], jr.logit(p), lines[sa], pools[sa])
    res["P1"]["note"] = "cells = all cells per run; binomial noise analytic"
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump(res, open(OUT / "a3_noise.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != "git_commit"}, indent=1, default=float))


if __name__ == "__main__":
    main()
