"""Phase 2, Workstream B0: refined DA-neuron labels (no clinical data). Tag da-labels-v1.0.0 before running.

Why: under scrn-v0.1.0, "DA_neuron" = Neuron-compartment cell with >= 1 TH UMI. B0 found that this
label is not DA-specific: apart from TH, canonical DA markers are not enriched against other
neurons (median log2 ratio: SLC6A3 0.00, KCNJ6 -0.26, SLC18A2 -0.11, FOXA2 -0.21), while GAD1 and
SLC17A6 are (0.38, 0.49).

Rule (fixed here):
  Among QC-passing cells labelled Neuron by scrn-v0.1.0 (DA_neuron or Other_neuron), compute a
  DA module score with the scrn-v0.1.0 score_genes algorithm (log1p of 1e4-normalised counts;
  25 bins; 50 control genes per bin; seed 0 per run) over
      DA_MARKERS = TH, DDC, SLC6A3, KCNJ6, SLC18A2, NR4A2, LMX1B, FOXA2, EN1
  (computed on all QC cells, so the controls are drawn genome-wide). A neuron is DA if score > t.
Calibration of the single threshold t (public data only):
  t is chosen on the grid -0.10, -0.075, ..., 0.60 to maximise the Pearson r between
  logit(per-line DA fraction of all QC cells) and logit(Table S3 "Dopaminergic Neurons"), over the
  79 donor lines (clip 0.005-0.995).
Checks reported: r at the chosen t (vs 0.57 for the old rule); marker enrichment of new DA vs other
neurons (median log2 CPM ratio and the share of runs > 0) for the 9 markers and for GAD1, SLC17A6;
the r curve across the grid (robustness).
Pass 1: per-run DA scores of neurons -> data/derived/da_relabel/scores/<run>.npz
Pass 2 (after t is fixed): per-run pseudobulk of new DA neurons and new other neurons ->
        data/derived/da_relabel/pseudobulk/<run>.csv.gz; summary.csv (DA fraction, n DA cells).
Output summary: data/derived/da_relabel/da_relabel.json (aggregates only).
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrn_per_run import find_runs, read_10x_h5, score_genes, QC  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "derived" / "da_relabel"
SCRN = ROOT / "data" / "derived" / "scrn"
DA_MARKERS = ["TH", "DDC", "SLC6A3", "KCNJ6", "SLC18A2", "NR4A2", "LMX1B", "FOXA2", "EN1"]
CHECK = DA_MARKERS + ["GAD1", "SLC17A6"]
GRID = np.round(np.arange(-0.10, 0.6001, 0.025), 3)
SEED = 0


def qc_and_norm(h5):
    X, barcodes, genes = read_10x_h5(h5)
    total = np.asarray(X.sum(axis=1)).ravel(); n_genes = np.diff(X.indptr)
    mt = np.char.startswith(genes.astype(str), "MT-")
    pct_mt = np.divide(np.asarray(X[:, np.flatnonzero(mt)].sum(axis=1)).ravel() * 100, total,
                       out=np.zeros_like(total, dtype=float), where=total > 0)
    q = ((n_genes >= QC["min_genes"]) & (n_genes <= QC["max_genes"]) & (total >= QC["min_counts"])
         & (pct_mt <= QC["max_pct_mt"]))
    return X[q].tocsr(), barcodes[q], genes, total[q]


def neuron_mask(run, barcodes):
    lab = pd.read_csv(SCRN / "runs" / run / "cell_labels.csv.gz", usecols=["barcode", "qc_pass", "cell_type"])
    lab = lab[lab.qc_pass.astype(bool)].set_index("barcode")["cell_type"]
    ct = lab.reindex(barcodes.astype(str))
    if ct.isna().any():
        raise ValueError("barcode mismatch with scrn-v0.1.0 labels")
    return ct.isin(["DA_neuron", "Other_neuron"]).to_numpy()


def pass1(runs):
    (OUT / "scores").mkdir(parents=True, exist_ok=True)
    for i, (run, h5) in enumerate(sorted(runs.items())):
        counts, bc, genes, tot = qc_and_norm(h5)
        gidx = {g: k for k, g in enumerate(genes)}
        Xn = (sp.diags((1e4 / tot).astype(np.float32)) @ counts).tocsr(); Xn.data = np.log1p(Xn.data)
        present = np.array([gidx[g] for g in DA_MARKERS if g in gidx], dtype=int)
        s = score_genes(Xn, present, np.random.default_rng(SEED))
        neu = neuron_mask(run, bc)
        np.savez_compressed(OUT / "scores" / f"{run}.npz", score=s.astype(np.float32), neuron=neu, n_qc=len(bc))
        print(f"pass1 {i+1}/{len(runs)}", flush=True)


def calibrate(s3_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    from phenotypes import read_table_s3, load_config, parse_runs, logit_clip
    cfg = load_config()
    s3 = read_table_s3(s3_path, cfg)
    files = sorted((OUT / "scores").glob("*.npz"))
    runs = parse_runs([f.stem for f in files])
    frac = np.zeros((len(files), len(GRID)))
    for i, f in enumerate(files):
        z = np.load(f)
        frac[i] = [(z["neuron"] & (z["score"] > t)).sum() / int(z["n_qc"]) for t in GRID]
    don = (~runs.is_control).to_numpy()
    m = runs[don].merge(s3, on="cdi", how="left")
    y = logit_clip(m["s3_da"].to_numpy(), .005, .995)
    r = [float(np.corrcoef(logit_clip(frac[don, j], .005, .995), y)[0, 1]) for j in range(len(GRID))]
    j = int(np.nanargmax(r))
    old = pd.read_csv(SCRN / "summary.csv").set_index("run").loc[runs.run[don], "frac_DA_neuron"].to_numpy()
    r_old = float(np.corrcoef(logit_clip(old, .005, .995), y)[0, 1])
    return {"grid": GRID.tolist(), "r_curve": r, "t": float(GRID[j]), "r_at_t": r[j], "r_old_rule": r_old,
            "n_donor_lines": int(don.sum()),
            "median_new_da_frac": float(np.median(frac[:, j])), "median_old_da_frac": float(np.median(old))}


def pass2(runs, t):
    (OUT / "pseudobulk").mkdir(parents=True, exist_ok=True)
    rows, enr = [], []
    for i, (run, h5) in enumerate(sorted(runs.items())):
        counts, bc, genes, tot = qc_and_norm(h5)
        z = np.load(OUT / "scores" / f"{run}.npz")
        da = z["neuron"] & (z["score"] > t); oth = z["neuron"] & ~da
        pb = pd.DataFrame({"gene": genes,
                           "DA_neuron_v2": np.asarray(counts[da].sum(0)).ravel(),
                           "Other_neuron_v2": np.asarray(counts[oth].sum(0)).ravel()})
        pb.to_csv(OUT / "pseudobulk" / f"{run}.csv.gz", index=False)
        rows.append({"run": run, "n_qc": int(len(bc)), "n_da_v2": int(da.sum()), "frac_DA_v2": float(da.mean())})
        if da.sum() and oth.sum():
            p = pb.set_index("gene"); cpm = p / p.sum() * 1e6
            enr.append({g: float(np.log2((cpm.loc[g, "DA_neuron_v2"] + 1) / (cpm.loc[g, "Other_neuron_v2"] + 1)))
                        for g in CHECK if g in cpm.index})
        print(f"pass2 {i+1}/{len(runs)}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    e = pd.DataFrame(enr)
    return {"marker_log2_ratio_median": e.median().round(3).to_dict(), "marker_share_runs_gt0": (e > 0).mean().round(3).to_dict(),
            "n_da_cells_median": float(np.median([r["n_da_v2"] for r in rows])),
            "runs_lt_100_da_cells": int(sum(r["n_da_v2"] < 100 for r in rows))}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default=str(ROOT / "data" / "unpacked"))
    ap.add_argument("--s3", required=True)
    ap.add_argument("--step", choices=["pass1", "calibrate", "pass2", "all"], default="all")
    a = ap.parse_args(argv)
    runs = find_runs(Path(a.in_dir))
    res_path = OUT / "da_relabel.json"
    res = json.load(open(res_path)) if res_path.exists() else {"version": "1.0.0", "markers": DA_MARKERS}
    if a.step in ("pass1", "all"):
        pass1(runs)
    if a.step in ("calibrate", "all"):
        res["calibration"] = calibrate(a.s3)
    if a.step in ("pass2", "all"):
        res["pass2"] = pass2(runs, res["calibration"]["t"])
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(res_path, "w"), indent=1)
    print("written da_relabel.json")


if __name__ == "__main__":
    main()
