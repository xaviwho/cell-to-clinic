#!/usr/bin/env python
"""
FOUNDIN-PD day-65 scRNA: per-run QC, marker-based annotation and pseudobulk.

Processes Cell Ranger `filtered_feature_bc_matrix.h5` files ONE RUN AT A TIME, so it
fits in modest RAM (no need to load the 42.7 GB Seurat object).

For each run it writes (under --out-dir, which must be inside the git-ignored data/ area):
  runs/<run>/metrics.json        QC + phenotype metrics (one row of the summary)
  runs/<run>/cell_labels.parquet barcode -> QC pass, cell type, TH UMI count
  runs/<run>/pseudobulk.parquet  raw-count sums, genes x cell type
and finally:
  summary.csv                    all runs' metrics combined
  run_info.json                  parameters, package versions, git commit

Governance: the console only reports counts of runs processed and timing; it never
prints per-run values or sample IDs. All outputs stay local (data/ is git-ignored).
Re-running skips runs whose outputs already exist (use --force to redo).

IMPORTANT: the QC thresholds and marker sets below are analysis decisions. Freeze them
(commit + tag) before running on the data, and describe them in the OSF registration.
Any change after seeing outputs must be logged in ACTIVITY.md.

Usage:
  python scripts/scrn_per_run.py --in-dir data/unpacked/SCRN_CellRanger \
         --out-dir data/derived/scrn [--dry-run] [--force] [--limit N]
"""
from __future__ import annotations

import argparse
import gc
import json
import subprocess
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

PIPELINE_VERSION = "0.1.0"

# ---------------------------------------------------------------- frozen parameters
QC = {
    "min_genes": 500,       # genes detected per cell
    "max_genes": 8000,      # crude doublet guard
    "min_counts": 1000,     # UMIs per cell
    "max_pct_mt": 20.0,     # % mitochondrial UMIs
}

# Broad compartments, scored with sc.tl.score_genes (fixed seed). Highest score wins.
MARKERS = {
    "Neuron": ["STMN2", "SNAP25", "SYT1", "GAP43", "DCX", "ELAVL4", "INA", "RBFOX3"],
    "Progenitor": ["SOX2", "VIM", "NES", "FOXA2", "LMX1A", "CORIN", "HES1", "OTX2"],
    "Proliferating": ["MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C", "BIRC5"],
    "Ependymal_like": ["FOXJ1", "PIFO", "RSPH1", "CCDC153", "TPPP3", "DYNLRB2"],
}
# Neurons are then split into DA vs other neurons by TH expression.
DA_GENE = "TH"
DA_MIN_UMI = 1
# Extra genes whose detection rate is reported per run (phenotype candidates).
REPORT_GENES = ["TH", "DDC", "SLC6A3", "KCNJ6", "NR4A2", "SNCA", "LRRK2", "GBA1"]
SCORE_SEED = 0
MIN_ASSIGN_SCORE = 0.0      # below this for all compartments -> "Unassigned"
# ----------------------------------------------------------------------------------


def find_runs(in_dir: Path) -> dict[str, Path]:
    """Map run name -> filtered h5. Handles <run>/outs/filtered_feature_bc_matrix.h5
    and flat files such as <run>_filtered_feature_bc_matrix.h5."""
    runs: dict[str, Path] = {}
    for h5 in sorted(in_dir.rglob("*.h5")):
        name = h5.name
        if "raw" in name.lower() or "molecule_info" in name:
            continue
        if "filtered" not in name.lower():
            continue
        if name == "filtered_feature_bc_matrix.h5":
            parent = h5.parent
            run = parent.parent.name if parent.name == "outs" else parent.name
        else:
            run = name.replace("_filtered_feature_bc_matrix.h5", "").replace(".h5", "")
        if run in runs:
            raise SystemExit(f"Duplicate run name found ({len(runs)} runs so far); check layout.")
        runs[run] = h5
    return runs


def process_run(run: str, h5: Path, out: Path) -> None:
    adata = sc.read_10x_h5(h5)
    adata.var_names_make_unique()
    n_raw = adata.n_obs

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    obs = adata.obs
    qc_pass = (
        (obs["n_genes_by_counts"] >= QC["min_genes"])
        & (obs["n_genes_by_counts"] <= QC["max_genes"])
        & (obs["total_counts"] >= QC["min_counts"])
        & (obs["pct_counts_mt"] <= QC["max_pct_mt"])
    ).to_numpy()

    ad = adata[qc_pass].copy()
    counts = ad.X.tocsr() if sp.issparse(ad.X) else sp.csr_matrix(ad.X)

    # TH UMI per cell (raw counts)
    if DA_GENE in ad.var_names:
        th_umi = np.asarray(counts[:, ad.var_names.get_loc(DA_GENE)].todense()).ravel()
    else:
        th_umi = np.zeros(ad.n_obs)

    # Normalised copy for scoring only
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    score_cols = []
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g in ad.var_names]
        col = f"score_{ct}"
        if len(present) >= 2:
            sc.tl.score_genes(ad, present, score_name=col, random_state=SCORE_SEED, use_raw=False)
        else:
            ad.obs[col] = -np.inf
        score_cols.append(col)
    scores = ad.obs[score_cols].to_numpy()
    best = scores.argmax(axis=1)
    labels = np.array(list(MARKERS))[best].astype(object)
    labels[scores.max(axis=1) <= MIN_ASSIGN_SCORE] = "Unassigned"
    is_neuron = labels == "Neuron"
    labels[is_neuron & (th_umi >= DA_MIN_UMI)] = "DA_neuron"
    labels[is_neuron & (th_umi < DA_MIN_UMI)] = "Other_neuron"

    # Per-cell table (all cells, QC-failed ones labelled as such)
    full_label = np.full(n_raw, "QC_fail", dtype=object)
    full_label[qc_pass] = labels
    full_th = np.zeros(n_raw)
    full_th[qc_pass] = th_umi
    cells = pd.DataFrame(
        {
            "barcode": adata.obs_names.astype(str),
            "qc_pass": qc_pass,
            "cell_type": full_label,
            "th_umi": full_th.astype(np.int32),
            "n_genes": obs["n_genes_by_counts"].to_numpy(),
            "total_counts": obs["total_counts"].to_numpy(),
            "pct_mt": obs["pct_counts_mt"].to_numpy(),
        }
    )

    # Pseudobulk: raw counts summed per cell type (genes x cell type)
    cts = sorted(set(labels))
    pb = {}
    for ct in cts:
        m = labels == ct
        pb[ct] = np.asarray(counts[m].sum(axis=0)).ravel()
    pb["All_QC_pass"] = np.asarray(counts.sum(axis=0)).ravel()
    pseudobulk = pd.DataFrame(pb, index=ad.var_names.astype(str)).astype(np.int64)

    # Metrics
    n = int(qc_pass.sum())
    vc = pd.Series(labels).value_counts()
    neurons = int(vc.get("DA_neuron", 0) + vc.get("Other_neuron", 0))
    metrics = {
        "run": run,
        "pipeline_version": PIPELINE_VERSION,
        "n_cells_raw": int(n_raw),
        "n_cells_qc": n,
        "frac_qc_pass": n / n_raw if n_raw else np.nan,
        "median_genes_qc": float(np.median(obs.loc[qc_pass, "n_genes_by_counts"])) if n else np.nan,
        "median_umi_qc": float(np.median(obs.loc[qc_pass, "total_counts"])) if n else np.nan,
        "median_pct_mt_qc": float(np.median(obs.loc[qc_pass, "pct_counts_mt"])) if n else np.nan,
        "frac_TH_pos_all": float((th_umi >= DA_MIN_UMI).mean()) if n else np.nan,
        "frac_DA_of_neurons": vc.get("DA_neuron", 0) / neurons if neurons else np.nan,
    }
    for ct in ["DA_neuron", "Other_neuron", "Progenitor", "Proliferating", "Ependymal_like", "Unassigned"]:
        metrics[f"frac_{ct}"] = vc.get(ct, 0) / n if n else np.nan
    for g in REPORT_GENES:
        if g in ad.var_names:
            col = np.asarray(counts[:, ad.var_names.get_loc(g)].todense()).ravel()
            metrics[f"det_{g}"] = float((col > 0).mean()) if n else np.nan
        else:
            metrics[f"det_{g}"] = np.nan

    out.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(out / "cell_labels.parquet", index=False)
    pseudobulk.to_parquet(out / "pseudobulk.parquet")
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2, default=float))

    del adata, ad, counts, cells, pseudobulk
    gc.collect()


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--dry-run", action="store_true", help="only count runs found")
    p.add_argument("--force", action="store_true", help="reprocess runs with existing outputs")
    p.add_argument("--limit", type=int, default=None, help="process at most N runs (testing)")
    a = p.parse_args()

    if "data" not in a.out_dir.resolve().parts:
        sys.exit("--out-dir must be inside the git-ignored data/ folder.")

    runs = find_runs(a.in_dir)
    print(f"Found {len(runs)} filtered .h5 runs.")
    if a.dry_run or not runs:
        return

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "run_info.json").write_text(json.dumps({
        "pipeline_version": PIPELINE_VERSION, "git_commit": git_commit(),
        "started": time.strftime("%Y-%m-%d %H:%M:%S"), "qc": QC, "markers": MARKERS,
        "da_gene": DA_GENE, "da_min_umi": DA_MIN_UMI, "score_seed": SCORE_SEED,
        "min_assign_score": MIN_ASSIGN_SCORE, "scanpy": version("scanpy"), "anndata": version("anndata"),
        "python": sys.version.split()[0],
    }, indent=2))

    todo = list(runs.items())[: a.limit] if a.limit else list(runs.items())
    done = skipped = failed = 0
    t0 = time.time()
    for i, (run, h5) in enumerate(todo, 1):
        out = a.out_dir / "runs" / run
        if (out / "metrics.json").exists() and not a.force:
            skipped += 1
            continue
        try:
            process_run(run, h5, out)
            done += 1
        except Exception as e:  # keep going; error type only, no sample ID on console
            failed += 1
            (a.out_dir / "errors.log").open("a").write(f"{run}\t{type(e).__name__}: {e}\n")
        print(f"[{i}/{len(todo)}] processed={done} skipped={skipped} failed={failed} "
              f"elapsed={time.time() - t0:.0f}s", flush=True)

    rows = [json.loads(f.read_text()) for f in sorted((a.out_dir / "runs").glob("*/metrics.json"))]
    pd.DataFrame(rows).to_csv(a.out_dir / "summary.csv", index=False)
    print(f"Done. summary.csv has {len(rows)} runs. Failures: {failed} (see errors.log).")


if __name__ == "__main__":
    main()
