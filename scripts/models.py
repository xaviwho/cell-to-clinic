"""Gate A0-1, hypotheses H1-H3, gate A0-5 and sensitivity analyses 1-9 (OSF osf.io/s6w9c, Part 5).

Run only after the registration is approved and the code release is tagged and archived.
Inputs: phenotypes (scripts/phenotypes.py), pre-link gates (scripts/gates_prelink.py), public
Table S1, and tidy clinical tables (scripts/ppmi_load.py). Output: one JSON of aggregate results.

Sample: PD donors (Table S1 Disease_Status "Idiopathic PD" or "Genetic PD") with a day-65 scRNA run,
minus A0-3 exclusions. A donor missing a covariate is excluded from that model and counted.

Covariates C (main effects; in H1/H2 also x time): age at onset, sex (female = 1), genetic PD
(= 1), PRS PD_SCORE_Nalls2019_no_GBA_LRRK2, media (E6 = batches 4-5 = 1), disease duration at PPMI
baseline. Continuous covariates (age at onset, PRS, duration) and phenotypes are z-scored over the
donors in each model. Time is not centred. Outcomes are untransformed.

H1  OFF-state MDS-UPDRS III ~ 1 + t + P + P:t + C + C:t + (1 + t | donor), REML, Satterthwaite df
    on P:t; t = years since diagnosis. Donors need >= 3 OFF visits (A0-1 rule). One model per
    phenotype retained by A0-2; A0-5 can drop P3/P4; BH over the remaining tests.
H2  Same model, Y = mean putamen SBR, all donors with >= 1 scan; separate BH family.
H3  OLS: P ~ SAA + C, and P ~ NSD-ISS + C; separate BH families. SAA is descriptive only if
    < 5 donors are SAA-negative.
Every test: beta, SE, df, p, q, partial r = t / sqrt(t^2 + df), TOST with 90% CI inside +/-0.40,
class supported / equivalent / inconclusive.
A0-5 For P3/P4: refit H1 with each member gene (residualised as the composite) in place of P;
    the composite stays in H1 if >= 80% of genes have the composite's sign of beta3.
Sensitivity analyses (registered numbering):
    1 idiopathic only (genetic-PD covariate dropped as constant)   [H1, H2, H3]
    2 five-level batch dummies instead of media                     [H1, H2, H3]
    3 P1 replaced by logit %TH+ (>=1 UMI, >=2 UMI) and totalDopa    [H1, H2, H3]
    4 P2-P4 residualised on our marker labels                       [H1, H2, H3]
    5 leave-one-donor-out influence on beta3                        [H1]
    6 ON-state UPDRS III with LEDD (z-scored) as time-varying cov.  [H1]
    7 visits within 5 years of diagnosis                            [H1, H2]
    8 MDS-UPDRS Part II total as outcome (all visits)               [H1]
    9 time from PPMI baseline instead of diagnosis                  [H1, H2]
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lmm import lmm_fit_fallback, satterthwaite, ols, tost_partial_r, classify, bh  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PHENOS = ["P1", "P2", "P3", "P4"]
COVS = ["age_onset", "sex_female", "genetic_pd", "prs", "media_e6", "dur_baseline"]
CONT = {"age_onset", "prs", "dur_baseline"}
A01_MIN_DONORS, A01_MIN_VISITS, A01_MIN_SPAN = 30, 3, 2.0
A05_MIN_SHARE = 0.80
SAA_MIN_NEG = 5
EQ_BOUND = 0.40


def num(x):
    return pd.to_numeric(pd.Series(x).astype(str).str.replace(",", ".", regex=False), errors="coerce").to_numpy()


def donor_table(lines, s1, clin_donors):
    s1 = s1.copy()
    s1["PATNO"] = s1["PATNO"].astype(str).str.extract(r"(\d+)")[0]
    s1 = s1.dropna(subset=["PATNO"]).groupby("PATNO").first()
    pd_rows = s1[s1["Disease_Status"].isin(["Idiopathic PD", "Genetic PD"])]
    d = lines[~lines.is_control].copy()
    d["PATNO"] = d["line"].str.extract(r"(\d+)")[0]
    d = d[d.PATNO.isin(pd_rows.index)].set_index("PATNO")
    d["genetic_pd"] = (pd_rows.loc[d.index, "Disease_Status"] == "Genetic PD").astype(float)
    d["prs"] = num(pd_rows.loc[d.index, "PD_SCORE_Nalls2019_no_GBA_LRRK2"])
    d["batch"] = num(pd_rows.loc[d.index, "BATCH"])
    d["media_e6"] = d["batch"].isin([4, 5]).astype(float).where(d["batch"].notna())
    c = clin_donors.copy()
    c["PATNO"] = c["PATNO"].astype(str)
    d = d.join(c.set_index("PATNO")[["sex_female", "age_onset", "dur_baseline", "saa_pos", "nsd_iss"]], how="left")
    return d


def z(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)


def design_donor(d, pcol, covs):
    """Per-donor columns, z-scored over these donors."""
    out = pd.DataFrame(index=d.index)
    out["P"] = z(d[pcol])
    for c in covs:
        out[c] = z(d[c]) if c in CONT else d[c].astype(float)
    return out


def batch_dummies(d):
    b = d["batch"]
    levels = sorted(b.dropna().unique())
    return pd.DataFrame({f"batch{int(l)}": (b == l).astype(float) for l in levels[1:]}, index=d.index)


def covariate_frame(d, covs, use_batch=False):
    base = [c for c in covs if not (use_batch and c == "media_e6")]
    f = d[base].copy()
    if use_batch:
        f = f.join(batch_dummies(d))
    return f


def fit_long(d, visits, pcol, covs, tcol="t_dx", min_visits=1, use_batch=False, visit_cov=None, keep_fit=False):
    cf = covariate_frame(d, covs, use_batch)
    ok = d[pcol].notna() & cf.notna().all(1)
    v = visits[visits.PATNO.isin(d.index[ok])].dropna(subset=["total", tcol])
    if visit_cov:
        v = v.dropna(subset=[visit_cov])
    cnt = v.groupby("PATNO").size()
    keep = cnt.index[cnt >= min_visits]
    v = v[v.PATNO.isin(keep)]
    dd = d.loc[keep]
    if len(dd) < 5:
        return {"n_donors": int(len(dd)), "error": "too few donors"}
    cfd = covariate_frame(dd, covs, use_batch)
    cfd = cfd.loc[:, cfd.nunique() > 1]  # drop constant columns (e.g. genetic_pd in idiopathic-only)
    P = z(dd[pcol])
    C = np.column_stack([z(cfd[c]) if c in CONT else cfd[c].to_numpy(float) for c in cfd.columns]) if cfd.shape[1] else np.zeros((len(dd), 0))
    idx = pd.Series(np.arange(len(dd)), index=dd.index)
    gi = idx.loc[v.PATNO].to_numpy()
    t = v[tcol].to_numpy(float)
    Pv, Cv = P[gi], C[gi]
    cols = [np.ones(len(v)), t, Pv, Pv * t, Cv, Cv * t[:, None]]
    if visit_cov:
        cols.append(z(v[visit_cov]))
    X = np.column_stack(cols)
    fit = lmm_fit_fallback(v["total"].to_numpy(float), X, gi, t)
    L = np.zeros(X.shape[1]); L[3] = 1
    s = satterthwaite(fit, L)
    s.update(tost_partial_r(s["t"], s["df"], EQ_BOUND))
    s.update({"n_donors": int(len(dd)), "n_obs": int(len(v)), "re": fit["re"],
              "fallback_log": fit["fallback_log"], "covariates": list(cfd.columns)})
    if keep_fit:
        s["_ids"] = list(dd.index)
    return s


def fit_cross(d, pcol, predictor, covs, use_batch=False):
    cf = covariate_frame(d, covs, use_batch)
    ok = d[pcol].notna() & d[predictor].notna() & cf.notna().all(1)
    dd, cf = d[ok], cf[ok]
    cf = cf.loc[:, cf.nunique() > 1]
    if len(dd) < len(cf.columns) + 5:
        return {"n_donors": int(len(dd)), "error": "too few donors"}
    C = np.column_stack([z(cf[c]) if c in CONT else cf[c].to_numpy(float) for c in cf.columns])
    pred = dd[predictor].to_numpy(float)
    X = np.column_stack([np.ones(len(dd)), pred, C])
    o = ols(z(dd[pcol]), X)
    s = {"estimate": float(o["beta"][1]), "se": float(o["se"][1]), "t": float(o["t"][1]),
         "df": o["df"], "p": float(o["p"][1]), "n_donors": int(o["n"])}
    s.update(tost_partial_r(s["t"], s["df"], EQ_BOUND))
    if predictor == "saa_pos":
        s["n_saa_neg"] = int((pred == 0).sum()); s["n_saa_pos"] = int((pred == 1).sum())
        if s["n_saa_neg"] < SAA_MIN_NEG:
            ph = z(dd[pcol])
            s["descriptive_only"] = True
            s["mean_by_group"] = {"neg": float(np.mean(ph[pred == 0])) if s["n_saa_neg"] else None,
                                  "pos": float(np.mean(ph[pred == 1]))}
    return s


def family(results, drop=()):
    keys = [k for k, r in results.items() if "p" in r and k not in drop and not r.get("descriptive_only")]
    q = bh([results[k]["p"] for k in keys]) if keys else []
    for k, qq in zip(keys, q):
        results[k]["q"] = float(qq)
        results[k]["class"] = classify(qq, results[k])
    for k in results:
        if k in drop:
            results[k]["class"] = "dropped_by_A0-5"
    return results


def a01(d, upd3):
    off = upd3[(upd3.state == "OFF") & upd3.PATNO.isin(d.index)].dropna(subset=["t_dx"])
    g = off.groupby("PATNO")["t_dx"].agg(["size", "min", "max"])
    ok = (g["size"] >= A01_MIN_VISITS) & ((g["max"] - g["min"]) >= A01_MIN_SPAN)
    n = int(ok.sum())
    return {"n_pd_donors": int(len(d)), "n_meeting_rule": n, "pass": n >= A01_MIN_DONORS,
            "rule": f">= {A01_MIN_DONORS} donors with >= {A01_MIN_VISITS} OFF visits spanning >= {A01_MIN_SPAN} y"}


def a05(d, off, genes, members, comp_sign):
    out = {}
    for k in ["P3", "P4"]:
        if k not in comp_sign:
            continue
        signs = []
        dg = d.join(genes[members[k]], how="left")
        for g in members[k]:
            r = fit_long(dg, off, g, COVS, min_visits=A01_MIN_VISITS)
            if "estimate" in r:
                signs.append(np.sign(r["estimate"]))
        share = float(np.mean(np.array(signs) == comp_sign[k])) if signs else np.nan
        out[k] = {"n_genes": len(signs), "share_same_sign": share, "pass": bool(share >= A05_MIN_SHARE)}
    return out


def run(lines, genes, members, gates, s1, clin, progress=print):
    res = {"version": "1.0.0"}
    d = donor_table(lines, s1, clin["donors"])
    res["sample"] = {"n_pd_with_scrna": int(len(d))}
    d = d[d.n_cells_qc >= gates["A0_3"]["min_qc_cells"]]
    res["sample"]["n_after_A0_3"] = int(len(d))
    res["sample"]["n_missing_covariate"] = int(d[COVS].isna().any(axis=1).sum())
    phenos = [p for p in PHENOS if p in gates["A0_2_phenotypes_retained"]]
    res["phenotypes_after_A0_2"] = phenos

    res["A0_1"] = a01(d, clin["updrs3"])
    if not res["A0_1"]["pass"]:
        res["stopped"] = "A0-1 failed: data audit only"
        return res
    off = clin["updrs3"][clin["updrs3"].state == "OFF"]
    on = clin["updrs3"][clin["updrs3"].state == "ON"]

    progress("H1")
    h1 = {p: fit_long(d, off, p, COVS, min_visits=A01_MIN_VISITS) for p in phenos}
    comp_sign = {k: np.sign(h1[k]["estimate"]) for k in ["P3", "P4"] if k in h1 and "estimate" in h1[k]}
    progress("A0-5")
    res["A0_5"] = a05(d, off, genes, members, comp_sign) if genes is not None else {}
    drop = [k for k, r in res["A0_5"].items() if not r["pass"]]
    res["H1"] = family(h1, drop)
    progress("H2")
    res["H2"] = family({p: fit_long(d, clin["datscan"], p, COVS) for p in phenos})
    progress("H3")
    res["H3_SAA"] = family({p: fit_cross(d, p, "saa_pos", COVS) for p in phenos})
    res["H3_NSDISS"] = family({p: fit_cross(d, p, "nsd_iss", COVS) for p in phenos})

    sens = {}
    progress("sensitivity 1-4")
    idio = d[d.genetic_pd == 0]
    variants = {
        "S1_idiopathic": dict(dd=idio, ph={p: p for p in phenos}),
        "S2_batch": dict(dd=d, ph={p: p for p in phenos}, use_batch=True),
        "S3_P1_alternatives": dict(dd=d, ph={k: k for k in ["S3a_thpos", "S3b_thpos2", "S3c_totaldopa"]}),
        "S4_marker_residualised": dict(dd=d, ph={f"S4_{p}": f"S4_{p}" for p in ["P2", "P3", "P4"] if p in phenos}),
    }
    for name, v in variants.items():
        ub = v.get("use_batch", False)
        sens[name] = {
            "H1": family({k: fit_long(v["dd"], off, c, COVS, min_visits=A01_MIN_VISITS, use_batch=ub) for k, c in v["ph"].items()}),
            "H2": family({k: fit_long(v["dd"], clin["datscan"], c, COVS, use_batch=ub) for k, c in v["ph"].items()}),
            "H3_SAA": family({k: fit_cross(v["dd"], c, "saa_pos", COVS, use_batch=ub) for k, c in v["ph"].items()}),
            "H3_NSDISS": family({k: fit_cross(v["dd"], c, "nsd_iss", COVS, use_batch=ub) for k, c in v["ph"].items()}),
        }
    progress("sensitivity 5 (leave-one-donor-out)")
    s5 = {}
    for p in phenos:
        full = h1[p]
        if "estimate" not in full:
            continue
        ids = fit_long(d, off, p, COVS, min_visits=A01_MIN_VISITS, keep_fit=True)["_ids"]
        est = []
        for i in ids:
            r = fit_long(d.drop(index=i), off, p, COVS, min_visits=A01_MIN_VISITS)
            if "estimate" in r:
                est.append(r["estimate"])
        est = np.array(est)
        s5[p] = {"beta3_full": full["estimate"], "loo_min": float(est.min()), "loo_max": float(est.max()),
                 "n_sign_changes": int((np.sign(est) != np.sign(full["estimate"])).sum()),
                 "max_abs_change_in_se_units": float(np.max(np.abs(est - full["estimate"])) / full["se"])}
    sens["S5_leave_one_out"] = s5
    progress("sensitivity 6-9")
    sens["S6_on_state_ledd"] = {"H1": family({p: fit_long(d, on, p, COVS, min_visits=A01_MIN_VISITS, visit_cov="ledd") for p in phenos})}
    sens["S7_first_5y"] = {"H1": family({p: fit_long(d, off[off.t_dx <= 5], p, COVS, min_visits=A01_MIN_VISITS) for p in phenos}),
                           "H2": family({p: fit_long(d, clin["datscan"][clin["datscan"].t_dx <= 5], p, COVS) for p in phenos})}
    sens["S8_part2"] = {"H1": family({p: fit_long(d, clin["updrs2"], p, COVS, min_visits=A01_MIN_VISITS) for p in phenos})}
    sens["S9_time_from_baseline"] = {"H1": family({p: fit_long(d, off, p, COVS, tcol="t_bl", min_visits=A01_MIN_VISITS) for p in phenos}),
                                     "H2": family({p: fit_long(d, clin["datscan"], p, COVS, tcol="t_bl") for p in phenos})}
    res["sensitivity"] = sens
    return res


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def main(argv=None):
    from phenotypes import load_config
    from ppmi_load import build
    ap = argparse.ArgumentParser()
    ap.add_argument("--phen-dir", default=str(ROOT / "data" / "derived" / "phenotypes"))
    ap.add_argument("--s1", required=True)
    ap.add_argument("--ppmi-dir", required=True)
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "models.json"))
    a = ap.parse_args(argv)
    from scrn_gate import read_table
    pdir = Path(a.phen_dir)
    lines = pd.read_csv(pdir / "lines.csv")
    genes = pd.read_csv(pdir / "genes_resid.csv", index_col=0)
    info = json.load(open(pdir / "phenotypes_info.json"))
    gates = json.load(open(pdir / "gates_prelink.json"))
    # genes_resid is indexed by run; re-key to PATNO for donors
    run2pat = lines.set_index("run")["line"].str.extract(r"(\d+)")[0]
    genes = genes.loc[lines.loc[~lines.is_control, "run"]]
    genes.index = run2pat.loc[genes.index].to_numpy()
    s1 = read_table(Path(a.s1))
    d0 = donor_table(lines, s1, pd.DataFrame({"PATNO": [], "sex_female": [], "age_onset": [], "dur_baseline": [], "saa_pos": [], "nsd_iss": []}))
    clin = build(a.ppmi_dir, list(d0.index))
    res = run(lines, genes, info["members"], gates, s1, clin)
    with open(a.out, "w") as f:
        json.dump(_clean(res), f, indent=1)
    print("written", a.out)


if __name__ == "__main__":
    main()
