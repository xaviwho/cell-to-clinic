"""A2: across-batch reliability atlas on Jerber et al. 2021 day-52 (public, CC-BY). EXPLORATORY.

Sample = (line, pool), untreated cells only. Each pool is an independent differentiation, so
replicate samples of a line in different pools give the ACROSS-BATCH within-line variance.

Pass 1 (cached): pseudobulk raw counts per sample, all cells and DA cells (Jerber's "DA" label).
Pass 2: phenotypes mirroring FOUNDIN (config/phenotypes.json):
  P1    logit DA fraction (clip 0.005-0.995)
  P2    SNCA log2CPM, all-cell pseudobulk          } residualised on [1, logit DA, logit progenitor]
  P3/P4 mean gene-z of OXPHOS / lysosome sets       } (progenitor = FPP + P_FPP + P_Sert)
  P2da/P3da/P4da  the same on the DA-only pseudobulk (no composition adjustment)
ICC = 1 - S_w^2 / S_d^2 (cell2clinic.gates.stability_icc): S_d^2 from one sample per line (for
replicated lines, the alphabetically first pool), S_w^2 pooled within replicated lines.
Pre-set minimum cells per sample: 200 all-cell, 50 DA (sensitivity: 100/500 and 20/100).
Also: genome-wide gene ICC distribution (genes with mean CPM >= 10).

Output: data/derived/jerber/{pseudobulk.npz, a2_reliability.json}
"""
import argparse, json, subprocess, sys
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, "D:/cell2clinic/src")
from jerber_inspect import read_col          # noqa: E402
from cell2clinic.gates import stability_icc  # noqa: E402

PROG = ("FPP", "P_FPP", "P_Sert")


def pseudobulk(h5, out, chunk=20000):
    f = h5py.File(h5, "r"); g = f["obs"]
    donor, pool, treat, ct = (read_col(g, k) for k in ("donor_id", "pool_id", "treatment", "celltype"))
    keep = treat == "NONE"
    keys = sorted(set(zip(donor[keep], pool[keep])))
    kid = {k: i for i, k in enumerate(keys)}
    gid = np.array([kid.get((d, p), -1) if k else -1 for d, p, k in zip(donor, pool, keep)])
    isda = ct == "DA"
    genes = np.array([x.decode() for x in f["raw/var/index"][()]])
    X = f["raw/X"]; ip = X["indptr"][()]
    G, n = len(keys), len(genes)
    A = np.zeros((G, n)); D = np.zeros((G, n))
    for r0 in range(0, len(gid), chunk):
        r1 = min(r0 + chunk, len(gid)); a, b = ip[r0], ip[r1]
        M = sparse.csr_matrix((X["data"][a:b], X["indices"][a:b], ip[r0:r1 + 1] - a), shape=(r1 - r0, n))
        gg = gid[r0:r1]
        for mat, sel in ((A, gg >= 0), (D, (gg >= 0) & isda[r0:r1])):
            idx = np.flatnonzero(sel)
            if len(idx):
                I = sparse.csr_matrix((np.ones(len(idx)), (gg[idx], idx)), shape=(G, r1 - r0))
                mat += (I @ M).toarray()
        print(f"\r{r1}/{len(gid)}", end="", flush=True)
    print()
    cts = sorted(np.unique(ct))
    comp = np.array([[np.sum((gid == i) & (ct == c)) for c in cts] for i in range(G)])
    np.savez_compressed(out, A=A, D=D, genes=genes, lines=np.array([k[0] for k in keys]),
                        pools=np.array([k[1] for k in keys]), comp=comp, celltypes=np.array(cts))


def logit(p, lo=0.005, hi=0.995):
    p = np.clip(p, lo, hi); return np.log(p / (1 - p))


def log2cpm(C):
    return np.log2(C / C.sum(1, keepdims=True) * 1e6 + 1)


def set_score(L, genes, gset):
    m = np.isin(genes, gset)
    Z = (L[:, m] - L[:, m].mean(0)) / L[:, m].std(0, ddof=1)
    return Z.mean(1), int(m.sum())


def resid(y, Z):
    Z1 = np.c_[np.ones(len(y)), Z]
    return y - Z1 @ np.linalg.lstsq(Z1, y, rcond=None)[0]


def split(lines, pools, y):
    """one value per line (first pool) and replicate groups for lines in >= 2 pools."""
    donor_vals, reps = [], []
    for l in np.unique(lines):
        idx = np.flatnonzero(lines == l)
        idx = idx[np.argsort(pools[idx])]
        donor_vals.append(y[idx[0]])
        if len(idx) >= 2:
            reps.append(y[idx])
    return np.array(donor_vals), reps


def icc_genes(L, lines, pools):
    """vectorised ICC for every column of L."""
    dv, reps = split(lines, pools, np.arange(len(lines)))
    sd2 = L[dv].var(0, ddof=1)
    ss = sum(((L[r] - L[r].mean(0)) ** 2).sum(0) for r in reps); df = sum(len(r) - 1 for r in reps)
    return 1 - (ss / df) / sd2, len(dv), df


def analyse(pb, min_all, min_da, cfg):
    A, D, genes = pb["A"], pb["D"], pb["genes"]; lines, pools = pb["lines"], pb["pools"]
    comp, cts = pb["comp"], list(pb["celltypes"])
    ncell = comp.sum(1); nda = comp[:, cts.index("DA")]
    prog = comp[:, [cts.index(c) for c in PROG]].sum(1)
    gsets = {k: [x.strip() for x in open(ROOT / v) if x.strip() and not x.startswith("#")]
             for k, v in cfg["gene_sets"].items()}
    res = {}
    # all-cell phenotypes
    s = ncell >= min_all
    Ls = log2cpm(A[s]); cpm = A[s] / A[s].sum(1, keepdims=True) * 1e6
    kg = cpm.mean(0) >= cfg["min_mean_cpm"]
    p1 = logit(nda[s] / ncell[s]); Zc = np.c_[p1, logit(prog[s] / ncell[s])]
    ph = {"P1": p1, "P2": resid(Ls[:, genes == "SNCA"][:, 0], Zc)}
    for k, lab in (("P3", "P3"), ("P4", "P4")):
        sc, ng = set_score(Ls[:, kg], genes[kg], gsets[k]); ph[lab] = resid(sc, Zc); res[f"{lab}_n_genes"] = ng
    ph["P2_unadj"] = Ls[:, genes == "SNCA"][:, 0]
    for k, y in ph.items():
        dv, reps = split(lines[s], pools[s], y); res[k] = stability_icc(dv, reps)
    gi, nd, dfw = icc_genes(Ls[:, kg], lines[s], pools[s])
    gr = icc_genes(resid(Ls[:, kg], Zc), lines[s], pools[s])[0]
    res["genes_all"] = {"n_genes": int(kg.sum()), "n_donor": nd, "df_within": dfw,
                        "icc_quantiles_10_25_50_75_90": np.nanpercentile(gi, [10, 25, 50, 75, 90]).round(3).tolist(),
                        "share_icc_ge_0.5": float(np.mean(gi >= 0.5)),
                        "resid_icc_quantiles_10_25_50_75_90": np.nanpercentile(gr, [10, 25, 50, 75, 90]).round(3).tolist(),
                        "resid_share_icc_ge_0.5": float(np.mean(gr >= 0.5))}
    res["n_samples_all"] = int(s.sum())
    # DA-only phenotypes
    s = nda >= min_da
    Ld = log2cpm(D[s]); kd = (D[s] / D[s].sum(1, keepdims=True) * 1e6).mean(0) >= cfg["min_mean_cpm"]
    phd = {"P2da": Ld[:, genes == "SNCA"][:, 0]}
    for k in ("P3", "P4"):
        phd[k + "da"] = set_score(Ld[:, kd], genes[kd], gsets[k])[0]
    for k, y in phd.items():
        dv, reps = split(lines[s], pools[s], y); res[k] = stability_icc(dv, reps)
    gi, nd, dfw = icc_genes(Ld[:, kd], lines[s], pools[s])
    res["genes_da"] = {"n_genes": int(kd.sum()), "n_donor": nd, "df_within": dfw,
                       "icc_quantiles_10_25_50_75_90": np.nanpercentile(gi, [10, 25, 50, 75, 90]).round(3).tolist(),
                       "share_icc_ge_0.5": float(np.mean(gi >= 0.5))}
    res["n_samples_da"] = int(s.sum())
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", default="D:/ST-DT4/data/external/jerber/x/day52.h5")
    ap.add_argument("--out", default=str(ROOT / "data/derived/jerber"))
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pbf = out / "pseudobulk.npz"
    if not pbf.exists():
        pseudobulk(a.h5, pbf)
    pb = dict(np.load(pbf, allow_pickle=False))
    cfg = json.load(open(ROOT / "config/phenotypes.json"))
    res = {"primary_min_cells": [200, 50]}
    for mn_all, mn_da in ((200, 50), (100, 20), (500, 100)):
        res[f"min{mn_all}_da{mn_da}"] = analyse(pb, mn_all, mn_da, cfg)
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump(res, open(out / "a2_reliability.json", "w"), indent=1, default=float)
    print(json.dumps(res["min200_da50"], indent=1, default=float))


if __name__ == "__main__":
    main()
