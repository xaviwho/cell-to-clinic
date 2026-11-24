"""Cell-to-Clinic II, Workstream B1: DA-neuron-only phenotypes and their pre-link gates. NO clinical data.
Tag da-phen-v1.0.0 before running. Rules below are fixed before any value is seen.

Inputs: da_relabel pseudobulk (da-labels-v1.0.0: DA_neuron_v2 per run), summary.csv; scrn-v0.1.0
pseudobulk (DA_neuron column, old TH>=1 rule) for the sensitivity analysis.

Phenotypes (per run):
  D1  logit(new-rule DA fraction of QC cells), clip 0.005-0.995
  D2  SNCA log2(CPM+1) in the DA-neuron pseudobulk
  D3  mean gene-z (population SD over donor lines) of HALLMARK OXPHOS genes, DA-neuron pseudobulk
  D4  as D3 for KEGG LYSOSOME
  Genes kept: mean CPM >= 10 over donor lines in the DA pseudobulk.
  No composition adjustment for D2-D4 (they are measured inside one cell type).
Variants (report only):
  *_old  D2-D4 from the old-rule DA_neuron pseudobulk (scrn-v0.1.0)
  *_adj  D2-D4 residualised on [1, D1, logit(marker progenitor + proliferating fraction)]
Gates (no clinical data):
  G-cells   a run needs >= 100 new-rule DA cells; runs below are excluded from D2-D4.
  G-icc     ICC = 1 - S_w^2/S_d^2 (S_d^2: donor lines, one run each; S_w^2: the control line's
            4 runs, 3 df), F-based 95% CI - the same estimator as registered gate A0-2.
            Reported next to the Jerber across-batch ICC (jerber-a2-v0.1.0) for the same phenotype,
            because FOUNDIN's control runs may share a batch (E5c, A2).
Outputs (data/, never committed): data/derived/phenotypes_da/{lines_da.csv, genes_da.csv, gates_da.json}
Only aggregates go to gates_da.json.
"""
from __future__ import annotations

import argparse, json, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenotypes import load_config, read_geneset, logit_clip, parse_runs  # noqa: E402
from gates_prelink import icc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DAR = ROOT / "data" / "derived" / "da_relabel"
SCRN = ROOT / "data" / "derived" / "scrn"
MIN_DA = 100
VERSION = "1.0.0"


def load_counts(runs, kind):
    cols = {}
    for r in runs:
        if kind == "new":
            s = pd.read_csv(DAR / "pseudobulk" / f"{r}.csv.gz", usecols=["gene", "DA_neuron_v2"]).set_index("gene")["DA_neuron_v2"]
        else:
            s = pd.read_csv(SCRN / "runs" / r / "pseudobulk.csv.gz", usecols=["gene", "DA_neuron"]).set_index("gene")["DA_neuron"]
        cols[r] = s.groupby(level=0).sum()
    return pd.DataFrame(cols).fillna(0.0).T.loc[list(runs)]


def scores(C, isd, ok, gs, cfg):
    """D2-D4 from a runs x genes count matrix; rows not ok -> NaN."""
    cpm = C.div(C.sum(1).replace(0, np.nan), axis=0) * 1e6
    keep = cpm.loc[isd & ok].mean(0) >= cfg["min_mean_cpm"]
    L = np.log2(cpm.loc[:, keep] + 1)
    mu, sd = L.loc[isd & ok].mean(0), L.loc[isd & ok].std(0, ddof=0).replace(0, np.nan)
    Z = (L - mu) / sd
    out = pd.DataFrame(index=C.index)
    out["D2"] = L[cfg["P2_gene"]] if cfg["P2_gene"] in L else np.nan
    n = {}
    for k, lab in (("P3", "D3"), ("P4", "D4")):
        g = [x for x in gs[k] if x in Z.columns]; n[lab] = len(g)
        out[lab] = Z[g].mean(1)
    out.loc[~ok] = np.nan
    return out, Z, n, int(keep.sum())


def compute(cfg):
    summ = pd.read_csv(DAR / "summary.csv")
    old = pd.read_csv(SCRN / "summary.csv")[["run", "frac_Progenitor", "frac_Proliferating", "frac_DA_neuron", "n_cells_qc"]]
    runs = parse_runs(summ["run"]).merge(summ, on="run").merge(old.drop(columns="n_cells_qc"), on="run")
    lo, hi = cfg["logit_clip"]
    isd = (~runs.is_control).to_numpy()
    ok = (runs.n_da_v2 >= MIN_DA).to_numpy()
    gs = {k: read_geneset(v) for k, v in cfg["gene_sets"].items()}
    idx = runs.run.tolist()
    new, Znew, ngs, nkept = scores(load_counts(idx, "new"), pd.Series(isd, idx), pd.Series(ok, idx), gs, cfg)
    n_old = (runs.frac_DA_neuron * runs.n_qc).to_numpy()
    ok_old = n_old >= MIN_DA
    oldp, _, _, _ = scores(load_counts(idx, "old"), pd.Series(isd, idx), pd.Series(ok_old, idx), gs, cfg)
    out = runs[["run", "line", "is_control", "n_qc", "n_da_v2"]].copy()
    out["D1"] = logit_clip(runs.frac_DA_v2, lo, hi)
    for k in ("D2", "D3", "D4"):
        out[k] = new[k].to_numpy(); out[k + "_old"] = oldp[k].to_numpy()
    comp = np.c_[out["D1"], logit_clip(runs.frac_Progenitor + runs.frac_Proliferating, lo, hi)]
    fit = isd & ok
    for k in ("D2", "D3", "D4"):
        y = out[k].to_numpy(); A = np.c_[np.ones(len(y)), comp]
        b = np.linalg.lstsq(A[fit], y[fit], rcond=None)[0]
        out[k + "_adj"] = np.where(ok, y - A @ b, np.nan)
    info = {"n_runs": int(len(out)), "n_donor_lines": int(isd.sum()), "n_control_runs": int((~isd).sum()),
            "n_genes_kept_da": nkept, "geneset_sizes_kept": ngs,
            "G_cells": {"rule": f">= {MIN_DA} new-rule DA cells", "donor_lines_excluded": int((isd & ~ok).sum()),
                        "control_runs_excluded": int((~isd & ~ok).sum()),
                        "old_rule_donor_lines_excluded": int((isd & ~ok_old).sum())}}
    gene_cols = [g for k in ("P3", "P4") for g in gs[k] if g in Znew.columns]
    return out, Znew[sorted(set(gene_cols))], info


def gates(out, jerber):
    don, ctl = out[~out.is_control], out[out.is_control]
    res = {}
    jmap = {"D1": "P1", "D2": "P2da", "D3": "P3da", "D4": "P4da"}
    for k in ["D1", "D2", "D3", "D4"] + [f"D{i}_{s}" for s in ("old", "adj") for i in (2, 3, 4)]:
        r = icc(don[k], ctl[k]); r["gating"] = k in ("D1", "D2", "D3", "D4")
        if r["icc"] is not None:
            r["pass_foundin_0.5"] = bool(r["icc"] >= 0.5)
        if k in jmap and jerber:
            j = jerber["min200_da50"][jmap[k]]
            r["jerber_across_batch_icc"] = j["icc"]; r["jerber_ci95"] = j["ci95"]
        res[k] = r
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "phenotypes_da"))
    a = ap.parse_args(argv)
    cfg = load_config()
    out, genes, info = compute(cfg)
    d = Path(a.out); d.mkdir(parents=True, exist_ok=True)
    out.to_csv(d / "lines_da.csv", index=False); genes.to_csv(d / "genes_da.csv")
    jf = ROOT / "data" / "derived" / "jerber" / "a2_reliability.json"
    res = {"version": VERSION, "info": info, "G_icc": gates(out, json.load(open(jf)) if jf.exists() else None),
           "git_commit": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()}
    json.dump(res, open(d / "gates_da.json", "w"), indent=1)
    print(json.dumps(info))
    for k, r in res["G_icc"].items():
        if r["icc"] is not None:
            j = f"  Jerber {r['jerber_across_batch_icc']:+.2f}" if "jerber_across_batch_icc" in r else ""
            print(f"{k:8s} ICC {r['icc']:+.2f} ({r['ci95'][0]:+.2f}, {r['ci95'][1]:+.2f}) n_d {r['n_donor']} n_w {r['n_within']}{j}")


if __name__ == "__main__":
    main()
