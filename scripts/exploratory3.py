"""EXPLORATORY, POST HOC: explore-v0.3.0. Does composition adjustment make bulk SNCA reliable?
Approved by the user on 2026-11-12 after E5 showed bulk SNCA ICC ~ 0 against scRNA
(composition-residualised) 0.94. ONE analysis only; specified and tagged before running.

Data: FOUNDIN RNAB bulk count table; IDs mapped to symbols as in explore-v0.2.1; log2(CPM + 1).
Composition markers, fixed in advance:
  TH        log2(CPM + 1) of TH (DA neurons)
  NEURO     mean of gene-wise z-scored log2(CPM + 1) of STMN2, SNAP25, SYT1 (pan-neuronal);
            z-scores use the donor lines of that day.
Adjusted SNCA = residual of SNCA on [1, TH, NEURO], fitted across the donor lines of that day by
least squares and applied to all samples of that day (donors and control runs), as in the scRNA
residualisation of P2.

Outputs, unadjusted vs adjusted SNCA:
  Day 65  ICC = 1 - S_w^2 / S_d^2 (5 control runs vs donor lines; 79 scRNA donor lines and all
          donor lines), with the F-based 95% CI (as in A0-2 and E5b).
  Day 25  share of control-run variance that is between batches = 1 - var(9 batch-1 runs) /
          var(all 13 control runs), with a 95% percentile bootstrap CI (B = 2000; batch-1 runs and
          the other 4 runs resampled separately; seed fixed), plus ICC_all_runs as in E5c.
No genes are filtered by expression here: the 5 marker genes and SNCA are all highly expressed.
Output: data/derived/exploratory3.json (aggregates only).
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exploratory2 import load_bulk, to_symbols  # noqa: E402
from gates_prelink import icc  # noqa: E402
from models import _clean  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261114
B = 2000
NEURO = ["STMN2", "SNAP25", "SYT1"]
NEED = ["SNCA", "TH", *NEURO]


def day_scores(raw, meta, day):
    m = meta[meta.day == day]
    c = raw[m.col]
    L = np.log2(c.div(c.sum(0), axis=1) * 1e6 + 1)
    miss = [g for g in NEED if g not in L.index]
    if miss:
        raise KeyError(f"marker genes missing: {miss}")
    don = list(m.loc[~m.ctrl, "col"])
    Zn = L.loc[NEURO]
    Zn = Zn.sub(Zn[don].mean(1), axis=0).div(Zn[don].std(1, ddof=0), axis=0)
    df = pd.DataFrame({"SNCA": L.loc["SNCA"], "TH": L.loc["TH"], "NEURO": Zn.mean(0)})
    A = np.c_[np.ones(len(don)), df.loc[don, ["TH", "NEURO"]]]
    beta, *_ = np.linalg.lstsq(A, df.loc[don, "SNCA"].to_numpy(), rcond=None)
    df["SNCA_adj"] = df["SNCA"] - np.c_[np.ones(len(df)), df[["TH", "NEURO"]]] @ beta
    r2 = 1 - np.var(df.loc[don, "SNCA_adj"]) / np.var(df.loc[don, "SNCA"])
    return df, m, {"composition_R2_in_donors": float(r2), "coef_TH": float(beta[1]), "coef_NEURO": float(beta[2])}


def batch_share(vals_b1, vals_other, rng):
    def share(b1, ot):
        a = np.r_[b1, ot]
        return max(0.0, 1 - np.var(b1, ddof=1) / np.var(a, ddof=1))
    pt = share(vals_b1, vals_other)
    bs = [share(rng.choice(vals_b1, len(vals_b1)), rng.choice(vals_other, len(vals_other))) for _ in range(B)]
    return {"share": float(pt), "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}


def run(bulk_path, h5_path, scrn_lines):
    rng = np.random.default_rng(SEED)
    raw, meta = load_bulk(bulk_path)
    raw, mapinfo = to_symbols(raw, h5_path)
    res = {"label": "EXPLORATORY, POST HOC", "version": "0.3.0", "seed": SEED, "markers": {"TH": "TH", "NEURO": NEURO},
           "bulk_id_mapping": mapinfo}
    df, m, fit = day_scores(raw, meta, 65)
    ctrl, don = m[m.ctrl], m[~m.ctrl]
    res["day65"] = {"fit": fit, "n_control_runs": int(len(ctrl))}
    for label, dset in [("scrn_79", don[don.line.isin(set(scrn_lines))]), ("all", don)]:
        res["day65"][label] = {k: icc(df.loc[dset.col, k], df.loc[ctrl.col, k]) for k in ("SNCA", "SNCA_adj")}
    df, m, fit = day_scores(raw, meta, 25)
    ctrl, don = m[m.ctrl], m[~m.ctrl]
    b1 = ctrl[ctrl.batch == "B1"].col; ot = ctrl[ctrl.batch != "B1"].col
    res["day25"] = {"fit": fit, "n_batch1_runs": int(len(b1)), "n_other_runs": int(len(ot))}
    for k in ("SNCA", "SNCA_adj"):
        bs = batch_share(df.loc[b1, k].to_numpy(), df.loc[ot, k].to_numpy(), rng)
        bs["icc_all_runs"] = icc(df.loc[don.col, k], df.loc[ctrl.col, k])
        res["day25"][k] = bs
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", default=str(ROOT / "data" / "unpacked_rnab" / "FOUNDIN-PD_150.9_RNAB" / "aggregated_expression" / "countTable.tsv"))
    ap.add_argument("--phen-dir", default=str(ROOT / "data" / "derived" / "phenotypes"))
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "exploratory3.json"))
    a = ap.parse_args(argv)
    lines = pd.read_csv(Path(a.phen_dir) / "lines.csv")
    scrn_lines = lines.loc[~lines.is_control, "line"].str.extract(r"(\d+)")[0]
    h5 = next((ROOT / "data" / "unpacked").rglob("filtered_feature_bc_matrix.h5"))
    res = run(a.bulk, h5, scrn_lines)
    with open(a.out, "w") as f:
        json.dump(_clean(res), f, indent=1)
    print("written exploratory3.json")


if __name__ == "__main__":
    main()
