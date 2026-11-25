"""EXPLORATORY, POST HOC (explore-II-v0.1.0). Approved by the user 2026-11-25 after the Cell-to-Clinic II
primary (H1-DA, D2 = SNCA in DA neurons) weakened from beta -0.57 (p 0.040) to -0.37 (p 0.20) when the
media covariate was replaced by five-level batch dummies (registered sensitivity 4). ONE analysis,
specified here and tagged before running. Aggregates only.

Question: is batch a confounder (batches differ in both D2 and progression), or does the batch
covariate mainly cost precision / absorb donor signal?
Sample: the 47 PD donors of the primary analysis; slopes need >= 3 OFF visits (the A0-1 rule).
  (a) D2 by batch: one-way ANOVA of D2 (z over the 47) across the 5 FOUNDIN batches: R^2, F, p.
  (b) Progression by batch: per-donor OLS slope of OFF-state MDS-UPDRS III on years since diagnosis
      (donors with >= 3 OFF visits); one-way ANOVA across batches: R^2, F, p.
  (c) Where the D2-slope association lives:
        overall       Pearson r(D2, slope)
        within-batch  r of batch-centred D2 and batch-centred slope (pooled; df reduced by 4)
        between-batch r of the 5 batch means (descriptive only, n = 5)
  (d) Media vs batch: share of the batch-mean variance of D2 and of slope that the media split
      (E6 = batches 4-5) already captures: R^2(media) / R^2(batch).
  Also reported: donors per batch (counts only).
Reading, fixed in advance:
  - confounding is plausible if BOTH (a) and (b) have p < 0.10 and the between-batch r has the same sign
    as the overall r;
  - otherwise, if the within-batch r is at least 75% of the overall r, the D2 association is mainly
    within batches and the attenuation reflects lost precision / adjustment noise;
  - anything else: "unresolved".
These are simple per-donor summaries, not the registered LMM, and are for interpretation only.
Output: data/derived/explore_batch.json
"""
from __future__ import annotations

import argparse, json, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from models import donor_table, A01_MIN_VISITS, _clean  # noqa: E402
from models_da import assemble, MIN_DA  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.1.0"


def anova(y, g):
    y, g = np.asarray(y, float), np.asarray(g)
    groups = [y[g == k] for k in np.unique(g)]
    f, p = stats.f_oneway(*groups)
    ss_b = sum(len(x) * (x.mean() - y.mean()) ** 2 for x in groups)
    return {"R2": float(ss_b / ((y - y.mean()) ** 2).sum()), "F": float(f), "p": float(p), "k": len(groups), "n": int(len(y))}


def r2_of(y, g):
    y = np.asarray(y, float); g = np.asarray(g)
    m = pd.Series(y).groupby(g).transform("mean").to_numpy()
    return float(((m - y.mean()) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def slopes(off, ids):
    out = {}
    for pid, v in off[off.PATNO.isin(ids)].dropna(subset=["total", "t_dx"]).groupby("PATNO"):
        if len(v) >= A01_MIN_VISITS and v.t_dx.nunique() > 1:
            out[pid] = np.polyfit(v.t_dx, v.total, 1)[0]
    return pd.Series(out)


def analyse(d, off):
    res = {"version": VERSION, "label": "EXPLORATORY, POST HOC"}
    d = d.copy()
    d["z"] = (d.D2 - d.D2.mean()) / d.D2.std(ddof=1)
    res["donors_per_batch"] = {str(int(k)): int(v) for k, v in d.batch.value_counts().sort_index().items()}
    res["a_D2_by_batch"] = anova(d.z, d.batch)
    s = slopes(off, d.index)
    dd = d.loc[s.index].assign(slope=s)
    res["b_slope_by_batch"] = anova(dd.slope, dd.batch)
    x, y, b = dd.z.to_numpy(), dd.slope.to_numpy(), dd.batch.to_numpy()
    xc = x - pd.Series(x).groupby(b).transform("mean").to_numpy()
    yc = y - pd.Series(y).groupby(b).transform("mean").to_numpy()
    r_all = stats.pearsonr(x, y)
    r_w = float(np.corrcoef(xc, yc)[0, 1]); df_w = len(x) - 2 - (len(np.unique(b)) - 1)
    t_w = r_w * np.sqrt(df_w / (1 - r_w ** 2))
    mb = dd.groupby("batch")[["z", "slope"]].mean()
    res["c_association"] = {"n": int(len(x)), "overall_r": float(r_all[0]), "overall_p": float(r_all[1]),
                            "within_batch_r": r_w, "within_batch_p": float(2 * stats.t.sf(abs(t_w), df_w)), "within_df": int(df_w),
                            "between_batch_r_n5": float(np.corrcoef(mb.z, mb.slope)[0, 1])}
    med = dd.media_e6.to_numpy()
    res["d_media_share"] = {"D2": r2_of(x, med) / r2_of(x, b), "slope": r2_of(y, med) / r2_of(y, b)}
    a, bb, c = res["a_D2_by_batch"], res["b_slope_by_batch"], res["c_association"]
    if a["p"] < 0.10 and bb["p"] < 0.10 and np.sign(c["between_batch_r_n5"]) == np.sign(c["overall_r"]):
        res["reading"] = "confounding plausible"
    elif abs(c["within_batch_r"]) >= 0.75 * abs(c["overall_r"]) and np.sign(c["within_batch_r"]) == np.sign(c["overall_r"]):
        res["reading"] = "association mainly within batches"
    else:
        res["reading"] = "unresolved"
    return res


def main(argv=None):
    from ppmi_load import build
    from scrn_gate import read_table
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", required=True)
    ap.add_argument("--ppmi-dir", required=True)
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "explore_batch.json"))
    a = ap.parse_args(argv)
    lines, _ = assemble(ROOT / "data/derived/phenotypes", ROOT / "data/derived/phenotypes_da")
    gates = json.load(open(ROOT / "data/derived/phenotypes/gates_prelink.json"))
    s1 = read_table(Path(a.s1))
    d0 = donor_table(lines, s1, pd.DataFrame(columns=["PATNO", "sex_female", "age_onset", "dur_baseline", "saa_pos", "nsd_iss"]))
    clin = build(a.ppmi_dir, list(d0.index))
    d = donor_table(lines, s1, clin["donors"])
    d = d[(d.n_cells_qc >= gates["A0_3"]["min_qc_cells"]) & (d.n_da_v2 >= MIN_DA)]
    res = analyse(d, clin["updrs3"][clin["updrs3"].state == "OFF"])
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump(_clean(res), open(a.out, "w"), indent=1)
    print(json.dumps(_clean(res), indent=1))


if __name__ == "__main__":
    main()
