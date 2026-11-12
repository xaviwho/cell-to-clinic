#!/usr/bin/env python
"""
Aggregate summary + validation gate for scrn_per_run.py outputs.

Prints AGGREGATES ONLY (counts, medians/IQRs, correlations, n). No IDs, no per-line values.

Gate (pre-specified, thresholds imported from scrn_per_run.py):
  primary   Pearson r: 100*frac_TH_pos_all  vs  Table S1 "scRNA %TH+ cells"  (matched by PATNO)
  secondary Pearson r: 100*frac_DA_neuron   vs  Table S3 "totalDopa"          (matched by FOUNDIN_ID
            digits = CDI suffix in run name; reported, not gating)
Lines with several runs (the control line) are averaged on both sides before correlating.

Usage:
  python -I scripts/scrn_gate.py --summary data/derived/scrn/summary.csv \
         --table-s1 <mmc2.xlsx> --table-s3 <mmc4.xlsx> [--out data/derived/scrn/gate.json]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrn_per_run import GATE_R_PARTIAL, GATE_R_PASS  # noqa: E402

RUN_RE = re.compile(r"^SCRN_(PPMI\d+)(B\w+?)?_(\d+)_da65_v\d+$")


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.replace(",", ".", regex=False).str.strip(), errors="coerce")


def read_table(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    hdr = next(i for i in range(10) if raw.iloc[i].notna().sum() >= 3)
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(c).strip() for c in raw.iloc[hdr]]
    return df.dropna(how="all")


def iqr(v: pd.Series) -> str:
    v = v.dropna()
    return f"median {v.median():.3g} (IQR {v.quantile(.25):.3g}-{v.quantile(.75):.3g}), n={len(v)}"


def corr(a: pd.Series, b: pd.Series) -> dict:
    d = pd.concat([a, b], axis=1).dropna()
    if len(d) < 3:
        return {"n": int(len(d)), "pearson_r": None, "spearman_rho": None}
    return {"n": int(len(d)), "pearson_r": round(float(d.iloc[:, 0].corr(d.iloc[:, 1])), 3),
            "spearman_rho": round(float(d.iloc[:, 0].corr(d.iloc[:, 1], method="spearman")), 3)}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--table-s1", type=Path, required=True)
    p.add_argument("--table-s3", type=Path, required=True)
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args()

    s = pd.read_csv(a.summary)
    parsed = s["run"].str.extract(RUN_RE)
    s["patno"], s["cdi"] = parsed[0], parsed[2]
    unparsed = int(s["patno"].isna().sum())

    res: dict = {"runs": int(len(s)), "runs_unparsed": unparsed,
                 "distinct_lines": int(s["patno"].nunique()),
                 "lines_with_multiple_runs": int((s.groupby("patno").size() > 1).sum())}
    agg = {
        "cells_raw": iqr(s["n_cells_raw"]), "cells_qc": iqr(s["n_cells_qc"]),
        "frac_qc_pass": iqr(s["frac_qc_pass"]),
        "pct_TH_pos (>=1 UMI)": iqr(100 * s["frac_TH_pos_all"]),
        "pct_TH_pos (>=2 UMI)": iqr(100 * s["frac_TH_pos_all_sens"]),
        "pct_DA_neuron": iqr(100 * s["frac_DA_neuron"]),
        "pct_Other_neuron": iqr(100 * s["frac_Other_neuron"]),
        "pct_Progenitor": iqr(100 * s["frac_Progenitor"]),
        "pct_Proliferating": iqr(100 * s["frac_Proliferating"]),
        "pct_Ependymal_like": iqr(100 * s["frac_Ependymal_like"]),
        "pct_Unassigned": iqr(100 * s["frac_Unassigned"]),
        "lines_lt10pct_TH (>=1 UMI)": int((s.groupby("patno")["frac_TH_pos_all"].mean() < .10).sum()),
        "lines_lt5pct_TH (>=1 UMI)": int((s.groupby("patno")["frac_TH_pos_all"].mean() < .05).sum()),
    }
    res["aggregates"] = agg

    line = s.groupby("patno")[["frac_TH_pos_all", "frac_TH_pos_all_sens", "frac_DA_neuron"]].mean() * 100

    s1 = read_table(a.table_s1)
    s1["patno"] = s1["PATNO"].astype(str).str.extract(r"(PPMI\d+)")[0]
    s1["th"] = num(s1["scRNA %TH+ cells"])
    s1l = s1.groupby("patno")["th"].mean()
    res["table_s1"] = {"lines_with_scRNA_TH": int(s1l.notna().sum()),
                       "our_lines_matched": int(line.index.isin(s1l.dropna().index).sum()),
                       "published_lines_not_in_ours": int((~s1l.dropna().index.isin(line.index)).sum())}
    prim = corr(line["frac_TH_pos_all"], s1l)
    res["gate_primary"] = prim
    res["gate_sensitivity_2umi"] = corr(line["frac_TH_pos_all_sens"], s1l)

    s3 = read_table(a.table_s3)
    s3["cdi"] = s3["FOUNDIN_ID"].astype(str).str.extract(r"(\d+)")[0]
    s3["dopa"] = num(s3["totalDopa"])
    s3c = s3.groupby("cdi")["dopa"].mean()
    cdi_line = s.dropna(subset=["cdi"]).groupby("cdi")["frac_DA_neuron"].mean() * 100
    res["table_s3"] = {"rows": int(len(s3)), "distinct_ids": int(s3c.size),
                       "our_lines_matched": int(cdi_line.index.isin(s3c.index).sum())}
    res["secondary_DA_vs_totalDopa"] = corr(cdi_line, s3c)

    r = prim["pearson_r"]
    res["gate_thresholds"] = {"pass": GATE_R_PASS, "partial": GATE_R_PARTIAL}
    res["gate_verdict"] = ("undetermined" if r is None else "pass" if r >= GATE_R_PASS
                           else "partial" if r >= GATE_R_PARTIAL else "fail")

    txt = json.dumps(res, indent=2)
    print(txt)
    if a.out:
        a.out.write_text(txt)


if __name__ == "__main__":
    main()
