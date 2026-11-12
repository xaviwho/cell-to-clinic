"""Gates A0-2 (technical stability) and A0-3 (cell sufficiency), OSF osf.io/s6w9c Part 5. No clinical data.

A0-2  For each phenotype: ICC = sigma_b^2 / (sigma_b^2 + sigma_w^2), where
        S_d^2 = variance (ddof=1) of the phenotype over the 79 donor lines (one run each, so this
                is between-line plus within-line variance),
        S_w^2 = variance (ddof=1) over the control line's 4 runs (3 df),
        ICC   = 1 - S_w^2 / S_d^2 (reported unclipped and clipped at 0).
      95% CI from (1 - ICC) * S_d^2/S_w^2 ~ F(n_d - 1, n_w - 1).
      Pass: ICC >= 0.5 (point estimate). A failing phenotype is dropped and reported.
      P1 uses the control line's 4 Table S3 rows; P2-P4 use its 4 runs (see phenotypes.py).
      Registered phenotypes: P1-P4. Also reported, not gating: S3a/S3b/S3c, S4_P2..S4_P4.
A0-3  Each line needs >= 1,000 QC-passing cells in its run; others are excluded. Reported as counts
      over all donor lines here; the count among analysed PD donors is reported by models.py.

Output: <phen-dir>/gates_prelink.json (aggregates only, safe to summarise in ACTIVITY.md).
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenotypes import load_config, logit_clip  # noqa: E402

REGISTERED = ["P1", "P2", "P3", "P4"]
REPORT_ONLY = ["S3a_thpos", "S3b_thpos2", "S3c_totaldopa", "S4_P2", "S4_P3", "S4_P4"]


def icc(donor_vals, ctrl_vals):
    d = np.asarray(donor_vals, float); d = d[np.isfinite(d)]
    w = np.asarray(ctrl_vals, float); w = w[np.isfinite(w)]
    if len(d) < 3 or len(w) < 2:
        return {"icc": None, "n_donor": int(len(d)), "n_within": int(len(w))}
    sd2, sw2 = d.var(ddof=1), w.var(ddof=1)
    R = sd2 / sw2 if sw2 > 0 else np.inf
    est = 1 - 1 / R
    f_lo, f_hi = stats.f.ppf([0.025, 0.975], len(d) - 1, len(w) - 1)
    return {"icc": float(est), "icc_clipped": float(max(est, 0.0)),
            "ci95": [float(1 - f_hi / R), float(1 - f_lo / R)],
            "n_donor": int(len(d)), "n_within": int(len(w)), "df_within": int(len(w) - 1)}


def run(phen_dir, cfg):
    phen_dir = Path(phen_dir)
    lines = pd.read_csv(phen_dir / "lines.csv")
    s3c = pd.read_csv(phen_dir / "control_s3_rows.csv")
    lo, hi = cfg["logit_clip"]
    don, ctl = lines[~lines.is_control], lines[lines.is_control]
    within = {k: ctl[k].to_numpy() for k in REGISTERED + REPORT_ONLY if k in ctl}
    within["P1"] = logit_clip(s3c["s3_da"], lo, hi)
    within["S3c_totaldopa"] = logit_clip(s3c["s3_totaldopa"], lo, hi)
    thr = cfg["a0_2_min_icc"]
    a02 = {}
    for k in REGISTERED + REPORT_ONLY:
        r = icc(don[k], within[k])
        r["gating"] = k in REGISTERED
        if r["icc"] is not None:
            r["pass"] = bool(r["icc"] >= thr)
        a02[k] = r
    n_min = cfg["a0_3_min_qc_cells"]
    a03 = {"rule": f">= {n_min} QC-passing cells", "min_qc_cells": int(n_min),
           "n_donor_lines": int(len(don)),
           "n_donor_lines_excluded": int((don.n_cells_qc < n_min).sum()),
           "n_control_runs_below": int((ctl.n_cells_qc < n_min).sum())}
    return {"A0_2": a02, "A0_2_threshold": thr,
            "A0_2_note": "within-line variance has 3 df (control line's 4 runs); CIs are wide",
            "A0_2_phenotypes_retained": [k for k in REGISTERED if a02[k].get("pass")],
            "A0_3": a03}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phen-dir", default=str(Path(__file__).resolve().parents[1] / "data" / "derived" / "phenotypes"))
    a = ap.parse_args(argv)
    res = run(a.phen_dir, load_config())
    with open(Path(a.phen_dir) / "gates_prelink.json", "w") as f:
        json.dump(res, f, indent=1)
    for k, r in res["A0_2"].items():
        if r["icc"] is not None:
            print(f"A0-2 {k:14s} ICC {r['icc']:+.2f} (95% CI {r['ci95'][0]:+.2f} to {r['ci95'][1]:+.2f})"
                  f" {'PASS' if r['pass'] else 'FAIL'}{'' if r['gating'] else ' [report only]'}")
    print("A0-3:", res["A0_3"])


if __name__ == "__main__":
    main()
