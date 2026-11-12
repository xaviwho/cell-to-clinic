#!/usr/bin/env python
"""
FOUNDIN-PD day-65 scRNA: per-run QC, marker-based annotation and pseudobulk.

Processes Cell Ranger `filtered_feature_bc_matrix.h5` files ONE RUN AT A TIME, so it
fits in modest RAM (no need to load the 42.7 GB Seurat object).

Dependencies: numpy, scipy, pandas, h5py only. scanpy/anndata/pyarrow are NOT used
because the analysis machine's Windows Application Control policy blocks some of their
compiled dependencies (numba, zarr/numcodecs, pyarrow). The 10x reader and the
gene-set scoring are reimplemented here (see read_10x_h5 and score_genes).

For each run it writes (under --out-dir, which must be inside the git-ignored data/ area):
  runs/<run>/metrics.json          QC + phenotype metrics (one row of the summary)
  runs/<run>/cell_labels.csv.gz    barcode -> QC pass, cell type, TH UMI count
  runs/<run>/pseudobulk.csv.gz     raw-count sums, genes x cell type
and finally:
  summary.csv                      all runs' metrics combined
  run_info.json                    parameters, package versions, git commit

Governance: the console only reports counts of runs processed, timing and memory; it
never prints per-run values or sample IDs. All outputs stay local (data/ is git-ignored).
Re-running skips runs whose outputs already exist (use --force to redo).

IMPORTANT: the QC thresholds and marker sets below are analysis decisions. They are
frozen by git tag before running on the data and described in the OSF registration.
Any change after seeing outputs must be logged in ACTIVITY.md.

Usage:
  python -I scripts/scrn_per_run.py --in-dir data/unpacked \
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

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

PIPELINE_VERSION = "0.1.0"

# ---------------------------------------------------------------- frozen parameters
QC = {
    "min_genes": 500,       # genes detected per cell
    "max_genes": 8000,      # crude doublet guard
    "min_counts": 1000,     # UMIs per cell
    "max_pct_mt": 20.0,     # % mitochondrial UMIs
}

# Broad compartments, scored with score_genes() (fixed seed). Highest score wins.
MARKERS = {
    "Neuron": ["STMN2", "SNAP25", "SYT1", "GAP43", "DCX", "ELAVL4", "INA", "RBFOX3"],
    "Progenitor": ["SOX2", "VIM", "NES", "FOXA2", "LMX1A", "CORIN", "HES1", "OTX2"],
    "Proliferating": ["MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C", "BIRC5"],
    "Ependymal_like": ["FOXJ1", "PIFO", "RSPH1", "CCDC153", "TPPP3", "DYNLRB2"],
}
# Neurons are then split into DA vs other neurons by TH expression.
DA_GENE = "TH"
DA_MIN_UMI = 1              # primary rule
DA_MIN_UMI_SENS = 2         # pre-specified sensitivity rule (ambient-RNA robustness)
# Extra genes whose detection rate is reported per run (phenotype candidates).
# Tuples = aliases: first symbol found in the reference is used (GBA1 was "GBA"
# in GENCODE releases used by Cell Ranger refs of the FOUNDIN era).
REPORT_GENES = ["TH", "DDC", "SLC6A3", "KCNJ6", "NR4A2", "SNCA", "LRRK2", ("GBA1", "GBA")]
SCORE_SEED = 0
SCORE_CTRL_SIZE = 50        # scanpy default
SCORE_N_BINS = 25           # scanpy default
MIN_ASSIGN_SCORE = 0.0      # below this for all compartments -> "Unassigned"

# Validation gate (pre-specified, see ACTIVITY.md; computed by scripts/scrn_gate.py):
# Pearson r between per-line 100*frac_TH_pos_all (primary rule; control line averaged
# over runs) and published per-line "scRNA %TH+ cells" (Bressan et al. 2023, Table S1).
# Secondary, reported but not gating: frac_DA_neuron vs Table S3 "totalDopa".
#   r >= 0.80 -> pass: re-annotation adequate, Seurat object not needed
#   0.60 <= r < 0.80 -> partial: proceed, request published cell metadata as cross-check
#   r < 0.60 -> fail: do not use re-annotation for phenotypes; obtain published metadata
GATE_R_PASS = 0.80
GATE_R_PARTIAL = 0.60
# ----------------------------------------------------------------------------------


def find_runs(in_dir: Path) -> dict[str, Path]:
    """Map run name -> filtered h5. Handles <run>/outs/filtered_feature_bc_matrix.h5,
    <run>/filtered_feature_bc_matrix.h5 and flat <run>_filtered_feature_bc_matrix.h5."""
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


def read_10x_h5(path: Path) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Read a Cell Ranger v3+ filtered matrix. Returns (cells x genes CSR counts,
    barcodes, gene symbols), keeping 'Gene Expression' features only, with gene
    symbols made unique like anndata's var_names_make_unique ('-1', '-2', ...)."""
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        n_genes, n_cells = g["shape"][:]
        m = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=(n_genes, n_cells))
        barcodes = g["barcodes"][:].astype(str)
        names = g["features/name"][:].astype(str)
        ftype = g["features/feature_type"][:].astype(str)
    keep = ftype == "Gene Expression"
    X = m[keep].T.tocsr().astype(np.float32)
    names = names[keep]
    seen: dict[str, int] = {}
    uniq = []
    dup = pd.Index(names).duplicated(keep="first")
    for nm, d in zip(names, dup):
        if d:
            seen[nm] = seen.get(nm, 0) + 1
            uniq.append(f"{nm}-{seen[nm]}")
        else:
            uniq.append(nm)
    return X, barcodes, np.array(uniq)


def score_genes(Xn: sp.csr_matrix, gene_idx: np.ndarray, rng: np.random.Generator,
                n_bins: int = SCORE_N_BINS, ctrl_size: int = SCORE_CTRL_SIZE) -> np.ndarray:
    """Seurat/scanpy-style module score (Satija et al. 2015; scanpy tl.score_genes
    algorithm): genes are binned by mean expression; for each bin touched by the gene
    set, ctrl_size random genes from that bin form the control set; score = mean
    expression of the set minus mean of the control set. Not bit-identical to scanpy
    (different RNG), same method."""
    gene_means = np.asarray(Xn.mean(axis=0)).ravel()
    n = len(gene_means)
    n_items = int(np.round(n / (n_bins - 1)))
    ranks = pd.Series(gene_means).rank(method="min").to_numpy()
    cut = (ranks // n_items).astype(int)
    set_mask = np.zeros(n, dtype=bool)
    set_mask[gene_idx] = True
    ctrl = set()
    for b in np.unique(cut[gene_idx]):
        pool = np.flatnonzero((cut == b) & ~set_mask)
        if len(pool):
            ctrl.update(rng.choice(pool, size=min(ctrl_size, len(pool)), replace=False).tolist())
    ctrl_idx = np.array(sorted(ctrl), dtype=int)
    s = np.asarray(Xn[:, gene_idx].mean(axis=1)).ravel()
    c = np.asarray(Xn[:, ctrl_idx].mean(axis=1)).ravel() if len(ctrl_idx) else 0.0
    return s - c


def process_run(run: str, h5: Path, out: Path) -> None:
    X, barcodes, genes = read_10x_h5(h5)
    gidx = {g: i for i, g in enumerate(genes)}
    n_raw = X.shape[0]

    total = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.diff(X.indptr)
    mt = np.char.startswith(genes.astype(str), "MT-")
    mt_counts = np.asarray(X[:, np.flatnonzero(mt)].sum(axis=1)).ravel()
    pct_mt = np.divide(mt_counts * 100, total, out=np.zeros_like(total), where=total > 0)
    qc_pass = (
        (n_genes >= QC["min_genes"]) & (n_genes <= QC["max_genes"])
        & (total >= QC["min_counts"]) & (pct_mt <= QC["max_pct_mt"])
    )

    counts = X[qc_pass]
    n = counts.shape[0]

    def gene_col(sym: str) -> np.ndarray | None:
        i = gidx.get(sym)
        return None if i is None else np.asarray(counts[:, i].todense()).ravel()

    th_umi = gene_col(DA_GENE)
    if th_umi is None:
        th_umi = np.zeros(n)

    # Normalised copy for scoring only (normalize_total 1e4 + log1p)
    tot_q = total[qc_pass]
    scale = np.divide(1e4, tot_q, out=np.zeros_like(tot_q), where=tot_q > 0)
    Xn = sp.diags(scale.astype(np.float32)) @ counts
    Xn = Xn.tocsr()
    Xn.data = np.log1p(Xn.data)

    rng = np.random.default_rng(SCORE_SEED)
    scores = np.full((n, len(MARKERS)), -np.inf)
    for j, (ct, gl) in enumerate(MARKERS.items()):
        present = np.array([gidx[g] for g in gl if g in gidx], dtype=int)
        if len(present) >= 2 and n:
            scores[:, j] = score_genes(Xn, present, rng)
    del Xn

    if n:
        best = scores.argmax(axis=1)
        labels = np.array(list(MARKERS), dtype=object)[best]
        labels[scores.max(axis=1) <= MIN_ASSIGN_SCORE] = "Unassigned"
    else:
        labels = np.array([], dtype=object)
    is_neuron = labels == "Neuron"
    labels[is_neuron & (th_umi >= DA_MIN_UMI)] = "DA_neuron"
    labels[is_neuron & (th_umi < DA_MIN_UMI)] = "Other_neuron"

    # Per-cell table (all cells, QC-failed ones labelled as such)
    full_label = np.full(n_raw, "QC_fail", dtype=object)
    full_label[qc_pass] = labels
    full_th = np.zeros(n_raw, dtype=np.int32)
    full_th[qc_pass] = th_umi.astype(np.int32)
    cells = pd.DataFrame({
        "barcode": barcodes, "qc_pass": qc_pass, "cell_type": full_label,
        "th_umi": full_th, "n_genes": n_genes, "total_counts": total.astype(np.int64),
        "pct_mt": pct_mt.round(3),
    })

    # Pseudobulk: raw counts summed per cell type (genes x cell type)
    pb = {}
    for ct in sorted(set(labels)):
        pb[ct] = np.asarray(counts[labels == ct].sum(axis=0)).ravel()
    pb["All_QC_pass"] = np.asarray(counts.sum(axis=0)).ravel()
    pseudobulk = pd.DataFrame(pb, index=pd.Index(genes, name="gene")).round().astype(np.int64)

    # Metrics
    vc = pd.Series(labels).value_counts()
    neurons = int(vc.get("DA_neuron", 0) + vc.get("Other_neuron", 0))
    med = lambda v: float(np.median(v)) if n else np.nan  # noqa: E731
    metrics = {
        "run": run,
        "pipeline_version": PIPELINE_VERSION,
        "n_cells_raw": int(n_raw),
        "n_cells_qc": int(n),
        "frac_qc_pass": n / n_raw if n_raw else np.nan,
        "median_genes_qc": med(n_genes[qc_pass]),
        "median_umi_qc": med(total[qc_pass]),
        "median_pct_mt_qc": med(pct_mt[qc_pass]),
        "frac_TH_pos_all": float((th_umi >= DA_MIN_UMI).mean()) if n else np.nan,
        "frac_TH_pos_all_sens": float((th_umi >= DA_MIN_UMI_SENS).mean()) if n else np.nan,
        "frac_DA_of_neurons": vc.get("DA_neuron", 0) / neurons if neurons else np.nan,
        "frac_DA_neuron_sens": float((is_neuron & (th_umi >= DA_MIN_UMI_SENS)).mean()) if n else np.nan,
    }
    for ct in ["DA_neuron", "Other_neuron", "Progenitor", "Proliferating", "Ependymal_like", "Unassigned"]:
        metrics[f"frac_{ct}"] = vc.get(ct, 0) / n if n else np.nan
    for g in REPORT_GENES:
        aliases = g if isinstance(g, tuple) else (g,)
        col = next((c for c in (gene_col(s) for s in aliases) if c is not None), None)
        metrics[f"det_{aliases[0]}"] = float((col > 0).mean()) if (col is not None and n) else np.nan

    out.mkdir(parents=True, exist_ok=True)
    cells.to_csv(out / "cell_labels.csv.gz", index=False)
    pseudobulk.to_csv(out / "pseudobulk.csv.gz")
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2, default=float))

    del X, counts, cells, pseudobulk
    gc.collect()


def peak_rss_gb() -> float:
    """Peak working set of this process (Windows), in GB; NaN elsewhere."""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL("kernel32")
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
        if not k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
            return float("nan")
        return pmc.PeakWorkingSetSize / 1024**3
    except Exception:
        return float("nan")


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
        "da_gene": DA_GENE, "da_min_umi": DA_MIN_UMI, "da_min_umi_sens": DA_MIN_UMI_SENS,
        "report_genes": [list(g) if isinstance(g, tuple) else g for g in REPORT_GENES],
        "score_seed": SCORE_SEED, "score_ctrl_size": SCORE_CTRL_SIZE, "score_n_bins": SCORE_N_BINS,
        "min_assign_score": MIN_ASSIGN_SCORE, "gate_r_pass": GATE_R_PASS, "gate_r_partial": GATE_R_PARTIAL,
        "numpy": version("numpy"), "scipy": version("scipy"), "pandas": version("pandas"),
        "h5py": version("h5py"), "python": sys.version.split()[0],
    }, indent=2))

    todo = list(runs.items())[: a.limit] if a.limit else list(runs.items())
    done = skipped = failed = 0
    t0 = time.time()
    for i, (run, h5) in enumerate(todo, 1):
        out = a.out_dir / "runs" / run
        if (out / "metrics.json").exists() and not a.force:
            skipped += 1
            continue
        t1 = time.time()
        try:
            process_run(run, h5, out)
            done += 1
        except Exception as e:  # keep going; error type only, no sample ID on console
            failed += 1
            with (a.out_dir / "errors.log").open("a") as fh:
                fh.write(f"{run}\t{type(e).__name__}: {e}\n")
        print(f"[{i}/{len(todo)}] processed={done} skipped={skipped} failed={failed} "
              f"run_s={time.time() - t1:.0f} elapsed_s={time.time() - t0:.0f} "
              f"peak_rss_gb={peak_rss_gb():.2f}", flush=True)

    rows = [json.loads(f.read_text()) for f in sorted((a.out_dir / "runs").glob("*/metrics.json"))]
    pd.DataFrame(rows).to_csv(a.out_dir / "summary.csv", index=False)
    print(f"Done. summary.csv has {len(rows)} runs. Failures: {failed} (see errors.log).")


if __name__ == "__main__":
    main()
