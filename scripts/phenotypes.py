"""Phenotypes P1-P4 per line, exactly as registered (OSF osf.io/s6w9c, Part 3). Config: config/phenotypes.json.

P1  logit(Table S3 "Dopaminergic Neurons" / 100), clipped to [0.005, 0.995] before the logit.
P2  log2(CPM + 1) of SNCA in the all-cell pseudobulk (column All_QC_pass of scrn-v0.1.0 output).
P3  mean over HALLMARK_OXIDATIVE_PHOSPHORYLATION genes of gene-wise z-scored log2(CPM + 1).
P4  as P3 over KEGG_LYSOSOME.
    Genes with mean CPM < 10 across the 79 donor lines are dropped before scoring; z-scores use the
    79 donor lines' mean and SD (population SD, ddof=0).
P2-P4 (and every P3/P4 member gene, for gate A0-5) are residualised across the 79 donor lines on
    [1, logit P1, logit(sum of Table S3 progenitor clusters)]; the residual is the phenotype.
Standardisation to mean 0, SD 1 over the analysed PD donors happens later, in models.py.

Also computed (registered sensitivity analyses):
  S3a  logit %TH+ (>= 1 UMI), S3b logit %TH+ (>= 2 UMI), S3c logit Table S3 totalDopa   (replace P1)
  S4   P2-P4 residualised on our marker-label composition (frac_DA_neuron; frac_Progenitor +
       frac_Proliferating) instead of Table S3.

Control line (4 runs, needed only for gate A0-2): Table S3 holds 4 rows for it, but which row is
which run is not recorded. P1 within-line variance uses the 4 rows as they are (order does not
matter). For P2-P4 the control runs are residualised with the line's mean Table S3 composition,
which leaves any composition-driven run-to-run variation in the within-line variance; this can only
lower the ICC. The S4 variant uses each run's own marker-label composition.

Outputs (all inside data/, never committed):
  <out>/lines.csv        one row per run: line id, is_control, n_cells_qc, P1..P4, variants
  <out>/genes_resid.csv  residualised member genes of P3/P4 (rows = runs), for gate A0-5
  <out>/phenotypes_info.json  gene counts kept per set, git commit, versions
"""
from __future__ import annotations

import argparse, glob, json, os, re, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUN_RE = re.compile(r"^SCRN_(PPMI\d+)(B\w+?)?_(\d+)_da65_v\d+$")
VERSION = "1.0.0"


def load_config(path=ROOT / "config" / "phenotypes.json"):
    with open(path) as f:
        return json.load(f)


def read_geneset(path):
    with open(ROOT / path if not os.path.isabs(path) else path) as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]


def logit_clip(p, lo, hi):
    p = np.clip(np.asarray(p, float), lo, hi)
    return np.log(p / (1 - p))


def parse_runs(runs):
    rows = []
    for r in runs:
        m = RUN_RE.match(r)
        if not m:
            raise ValueError("unrecognised run name format")
        rows.append({"run": r, "line": m.group(1), "cdi": m.group(3), "is_control": m.group(2) is not None})
    return pd.DataFrame(rows)


def read_table_s3(path, cfg):
    s3 = pd.read_excel(path, header=1)
    c = cfg["table_s3"]
    need = [c["id_column"], c["P1_da_column"], c["totaldopa_column"], *c["progenitor_columns"]]
    miss = [x for x in need if x not in s3.columns]
    if miss:
        raise ValueError(f"Table S3 columns missing: {miss}")
    s3 = s3.dropna(subset=[c["id_column"]]).copy()
    s3["cdi"] = s3[c["id_column"]].astype(str).str.extract(r"(\d+)")[0]
    out = pd.DataFrame({"cdi": s3["cdi"],
                        "s3_da": s3[c["P1_da_column"]].astype(float) / 100,
                        "s3_totaldopa": s3[c["totaldopa_column"]].astype(float) / 100,
                        "s3_prog": s3[c["progenitor_columns"]].astype(float).sum(1) / 100})
    return out


def residualise(Y, Z_fit, Z_apply):
    """Fit Y (lines x k) on Z_fit with intercept; return residuals for Z_apply rows' Y."""
    A = np.c_[np.ones(len(Z_fit)), Z_fit]
    B, *_ = np.linalg.lstsq(A, Y["fit"], rcond=None)
    return Y["apply"] - np.c_[np.ones(len(Z_apply)), Z_apply] @ B


def compute(scrn_dir, s3_path, cfg):
    scrn_dir = Path(scrn_dir)
    summ = pd.read_csv(scrn_dir / "summary.csv")
    runs = parse_runs(summ["run"]).merge(summ, on="run")
    if runs["is_control"].sum() < 2 or runs.loc[~runs.is_control, "line"].duplicated().any():
        raise ValueError("expected >= 2 control runs and one run per donor line")
    lo, hi = cfg["logit_clip"]

    # Table S3: donors by CDI; control rows kept separately (order-free use only)
    s3 = read_table_s3(s3_path, cfg)
    ctrl_cdi = set(runs.loc[runs.is_control, "cdi"])
    s3_ctrl = s3[s3.cdi.isin(ctrl_cdi)]
    s3_don = s3[~s3.cdi.isin(ctrl_cdi)]
    if s3_don.cdi.duplicated().any():
        raise ValueError("duplicate Table S3 rows for a donor line")
    don = runs[~runs.is_control].merge(s3_don, on="cdi", how="left")
    ctl = runs[runs.is_control].copy()
    for col in ["s3_da", "s3_totaldopa", "s3_prog"]:
        ctl[col] = s3_ctrl[col].mean() if len(s3_ctrl) else np.nan
    lines = pd.concat([don, ctl], ignore_index=True)

    # pseudobulk -> CPM -> log2(CPM+1)
    col = cfg["pseudobulk_column"]
    counts = {}
    for r in lines["run"]:
        pb = pd.read_csv(scrn_dir / "runs" / r / "pseudobulk.csv.gz", usecols=["gene", col]).set_index("gene")[col]
        counts[r] = pb
    C = pd.DataFrame(counts).fillna(0.0).T.loc[lines["run"]]  # runs x genes
    cpm = C.div(C.sum(1), axis=0) * 1e6
    isd = (~lines["is_control"]).to_numpy()
    keep = cpm.loc[isd].mean(0) >= cfg["min_mean_cpm"]
    L = np.log2(cpm.loc[:, keep] + 1)
    mu, sd = L.loc[isd].mean(0), L.loc[isd].std(0, ddof=0)
    Zg = (L - mu) / sd.replace(0, np.nan)

    info = {"version": VERSION, "n_genes_total": int(C.shape[1]), "n_genes_kept": int(keep.sum())}
    gs = {k: read_geneset(v) for k, v in cfg["gene_sets"].items()}
    raw = pd.DataFrame(index=lines["run"])
    g2 = cfg["P2_gene"]
    raw["P2"] = L[g2] if g2 in L.columns else np.nan
    info["P2_gene_kept"] = bool(g2 in L.columns)
    members = {}
    for k, genes in gs.items():
        g = [x for x in genes if x in Zg.columns]
        members[k] = g
        raw[k] = Zg[g].mean(1)
        info[f"{k}_genes_in_set"] = len(genes)
        info[f"{k}_genes_used"] = len(g)

    # P1 and composition covariates
    lines["P1"] = logit_clip(lines["s3_da"], lo, hi)
    lines["S3c_totaldopa"] = logit_clip(lines["s3_totaldopa"], lo, hi)
    lines["S3a_thpos"] = logit_clip(lines["frac_TH_pos_all"], lo, hi)
    lines["S3b_thpos2"] = logit_clip(lines["frac_TH_pos_all_sens"], lo, hi)
    comp_s3 = np.c_[lines["P1"], logit_clip(lines["s3_prog"], lo, hi)]
    mk = cfg["marker_label_columns_for_sensitivity_4"]
    comp_mk = np.c_[logit_clip(lines[mk["da"]], lo, hi),
                    logit_clip(lines[mk["progenitor"]].sum(1), lo, hi)]

    gene_cols = sorted(set(members["P3"]) | set(members["P4"]))
    for label, comp in [("", comp_s3), ("S4_", comp_mk)]:
        Y = raw[["P2", "P3", "P4"]].to_numpy()
        ok = np.isfinite(Y).all(0)
        res = np.full_like(Y, np.nan)
        res[:, ok] = residualise({"fit": Y[isd][:, ok], "apply": Y[:, ok]}, comp[isd], comp)
        for j, k in enumerate(["P2", "P3", "P4"]):
            lines[label + k] = res[:, j]
        if label == "":
            Yg = Zg[gene_cols].to_numpy()
            Rg = residualise({"fit": Yg[isd], "apply": Yg}, comp[isd], comp)
            genes_resid = pd.DataFrame(Rg, index=lines["run"], columns=gene_cols)
    keep_cols = ["run", "line", "is_control", "n_cells_qc", "P1", "P2", "P3", "P4",
                 "S3a_thpos", "S3b_thpos2", "S3c_totaldopa", "S4_P2", "S4_P3", "S4_P4"]
    info["members"] = members
    info["n_donor_lines"] = int(isd.sum())
    info["n_control_runs"] = int((~isd).sum())
    info["n_control_s3_rows"] = int(len(s3_ctrl))
    info["n_donor_lines_missing_s3"] = int(lines.loc[isd, "s3_da"].isna().sum())
    return lines[keep_cols], genes_resid, s3_ctrl.reset_index(drop=True), info


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrn-dir", default=str(ROOT / "data" / "derived" / "scrn"))
    ap.add_argument("--s3", required=True, help="path to Table S3 (mmc4.xlsx)")
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "phenotypes"))
    a = ap.parse_args(argv)
    out = Path(a.out).resolve()
    if (ROOT / "data") not in out.parents and not os.environ.get("C2C_SYNTHETIC"):
        sys.exit("output directory must be inside data/ (set C2C_SYNTHETIC=1 for synthetic tests)")
    out.mkdir(parents=True, exist_ok=True)
    lines, genes, s3_ctrl, info = compute(a.scrn_dir, a.s3, load_config())
    lines.to_csv(out / "lines.csv", index=False)
    genes.to_csv(out / "genes_resid.csv")
    s3_ctrl.to_csv(out / "control_s3_rows.csv", index=False)
    info["git_commit"] = git_commit()
    info["numpy"], info["pandas"] = np.__version__, pd.__version__
    with open(out / "phenotypes_info.json", "w") as f:
        json.dump(info, f, indent=1)
    print(f"phenotypes: {info['n_donor_lines']} donor lines, {info['n_control_runs']} control runs, "
          f"{info['n_genes_kept']} genes kept; P3 {info['P3_genes_used']}/{info['P3_genes_in_set']}, "
          f"P4 {info['P4_genes_used']}/{info['P4_genes_in_set']}")


if __name__ == "__main__":
    main()
