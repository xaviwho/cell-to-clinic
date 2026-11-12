"""Exploratory analyses after the confirmatory run (OSF osf.io/s6w9c, Part 5 "Exploratory").
Every result here is EXPLORATORY and must be labelled so. Specified and frozen (git tag
explore-v0.1.0) before running on real data. Inputs as for models.py.

E1  Component genes of P3 (OXPHOS) and P4 (Lysosome), sent here by gate A0-5. The H1 model with each
    residualised gene in place of the phenotype; BH within each gene set. Reports the number of genes
    with q < 0.05, the sign share, the distribution of partial r and the 10 strongest genes.
E2  Time from diagnosis to clinical milestones (registered as exploratory):
      M1 Hoehn & Yahr >= 3 at an OFF-state Part III exam (NHY; values > 5 = not rated -> missing)
      M2 freezing of gait: MDS-UPDRS Part II item NP2FREZ >= 1
      M3 cognitive impairment: investigator diagnosis cogstate >= 2 (MCI or dementia), Curated Cut
      (falls: no variable in the downloaded data -> not analysed)
    Cox proportional hazards with delayed entry: entry = first assessment of that milestone (years
    since diagnosis); event = first assessment meeting the criterion; otherwise censored at the last
    assessment. Donors meeting the criterion at their first assessment (prevalent) are excluded.
    Model: phenotype (z) + age at onset (z) + genetic PD + media E6. Breslow ties. Milestones with
    < 10 incident events are descriptive only. BH over all phenotype x milestone tests. HR per SD.
E3  Manifesting vs non-manifesting carriers (registered as exploratory): Table S1 "Genetic PD"
    vs "Genetic Uffected" (sic) lines. OLS: phenotype (z within these lines) ~ manifesting + LRRK2
    carrier + media E6 + sex female. BH over phenotypes; TOST +/-0.40 on partial r.
E4  Design analysis: what a cell-to-clinic study needs. No phenotype-outcome association is used.
    - Phenotype reliability R1 = A0-2 ICC (point estimate); mean of k independent differentiation runs
      per line: R_k = k R1 / (1 + (k - 1) R1) (Spearman-Brown).
    - Outcome-slope reliability: LMM of the outcome on 1 + t + C + C:t + (1 + t | donor) with no
      phenotype; lambda_i = tau1^2 / (tau1^2 + s_i^2), where s_i^2 = sigma^2 / sum_j (t_ij - tbar_i)^2;
      lambda = mean over donors. Computed for OFF-state UPDRS III and putamen SBR.
    - Observable correlation r_obs = r_true * sqrt(R_k * lambda). Donors needed for 80% power
      (two-sided alpha 0.05): N = ((z_.975 + z_.80) / atanh(r_obs))^2 + 3 + 6 covariates.
    - Grid: r_true in {0.2, 0.3, 0.4, 0.5}, k in {1, 2, 3, 4}, per phenotype (P1-P4) and outcome.
    - Also: the true correlation detectable at 80% power with this study's n and k = 1.
Output: data/derived/exploratory.json (aggregates only).
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lmm import lmm_fit, ols, tost_partial_r, bh  # noqa: E402
from models import donor_table, fit_long, z, COVS, CONT, A01_MIN_VISITS, _clean  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MIN_EVENTS = 10
R_TRUE = [0.2, 0.3, 0.4, 0.5]
K_RUNS = [1, 2, 3, 4]
COX_COVS = ["age_onset", "genetic_pd", "media_e6"]


# ---------------- Cox PH with delayed entry ----------------
def cox_fit(entry, time, event, X, max_iter=50):
    entry, time, event, X = map(np.asarray, (entry, time, event, X))
    X = X.astype(float)
    n, p = X.shape
    beta = np.zeros(p)
    ev_times = np.unique(time[event == 1])

    def ll_grad_hess(b):
        eta = X @ b
        w = np.exp(eta - eta.max())
        ll, g, H = 0.0, np.zeros(p), np.zeros((p, p))
        for t in ev_times:
            risk = (entry < t) & (time >= t)
            d = (time == t) & (event == 1)
            nd = d.sum()
            wr = w[risk]
            S0 = wr.sum(); S1 = wr @ X[risk]; S2 = (X[risk] * wr[:, None]).T @ X[risk]
            ll += (eta[d]).sum() - nd * (np.log(S0) + eta.max())
            g += X[d].sum(0) - nd * S1 / S0
            H -= nd * (S2 / S0 - np.outer(S1, S1) / S0 ** 2)
        return ll, g, H

    ll, g, H = ll_grad_hess(beta)
    for _ in range(max_iter):
        step = np.linalg.solve(-H, g)
        s = 1.0
        while True:
            nb = beta + s * step
            nll, ng, nH = ll_grad_hess(nb)
            if nll >= ll - 1e-12 or s < 1e-4:
                break
            s /= 2
        beta, conv = nb, abs(nll - ll) < 1e-10
        ll, g, H = nll, ng, nH
        if conv:
            break
    cov = np.linalg.inv(-H)
    se = np.sqrt(np.diag(cov))
    zv = beta / se
    return {"beta": beta, "se": se, "p": 2 * stats.norm.sf(np.abs(zv)), "loglik": float(ll)}


def milestone_table(assess, crit_col, threshold):
    """assess: PATNO, t_dx, value. Returns per-donor entry, time, event; prevalent removed."""
    a = assess.dropna(subset=["t_dx", crit_col]).sort_values(["PATNO", "t_dx"])
    rows = []
    n_prev = 0
    for p, g in a.groupby("PATNO"):
        hit = g[crit_col].to_numpy() >= threshold
        if hit[0]:
            n_prev += 1
            continue
        t = g["t_dx"].to_numpy()
        if hit.any():
            rows.append((p, t[0], t[np.argmax(hit)], 1))
        elif len(t) > 1:
            rows.append((p, t[0], t[-1], 0))
    return pd.DataFrame(rows, columns=["PATNO", "entry", "time", "event"]), n_prev


# ---------------- helpers to read milestone variables (same OFF rule as ppmi_load) ----------------
def read_milestones(ppmi_dir, patnos, donors):
    from ppmi_load import parse_date, load_colmap
    cm = load_colmap()
    pdir = Path(ppmi_dir)
    pats = set(map(str, patnos))
    dx = donors.set_index("PATNO")["dx_date"]

    def dates(s):
        return s.map(lambda v: parse_date(v, cm["date_formats"], cm["date_day_if_month_only"]))

    def tdx(df, datecol):
        d = dates(df[datecol])
        return (d - df["PATNO"].map(dx)).dt.days / 365.25

    u3 = pd.read_csv(next(pdir.glob("MDS-UPDRS_Part_III*.csv")), low_memory=False)
    u3 = u3[u3.PATNO.astype(str).isin(pats)].copy(); u3["PATNO"] = u3.PATNO.astype(str)
    tr = pd.to_numeric(u3.PDTRTMNT, errors="coerce"); st = u3.PDSTATE.astype(str).str.upper().str.strip()
    u3 = u3[(tr == 0) | (st == "OFF")].copy()
    u3["nhy"] = pd.to_numeric(u3.NHY, errors="coerce").where(lambda v: v <= 5)
    u3["t_dx"] = tdx(u3, "INFODT")
    u2 = pd.read_csv(next(pdir.glob("MDS_UPDRS_Part_II__Patient_Questionnaire*.csv")), low_memory=False)
    u2 = u2[u2.PATNO.astype(str).isin(pats)].copy(); u2["PATNO"] = u2.PATNO.astype(str)
    u2["frz"] = pd.to_numeric(u2.NP2FREZ, errors="coerce"); u2["t_dx"] = tdx(u2, "INFODT")
    cur = next(pdir.glob("PPMI_Curated_Data_Cut_*.xlsx"))
    xl = pd.ExcelFile(cur)
    cg = pd.read_excel(xl, xl.sheet_names[0], usecols=["PATNO", "visit_date", "cogstate"], dtype={"PATNO": str, "visit_date": str})
    cg = cg[cg.PATNO.isin(pats)].copy()
    cg["cog"] = pd.to_numeric(cg.cogstate, errors="coerce"); cg["t_dx"] = tdx(cg, "visit_date")
    return {"M1_HY3": (u3[["PATNO", "t_dx", "nhy"]], "nhy", 3),
            "M2_freezing": (u2[["PATNO", "t_dx", "frz"]], "frz", 1),
            "M3_cognitive": (cg[["PATNO", "t_dx", "cog"]], "cog", 2)}


# ---------------- E4 ----------------
def slope_reliability(d, visits, tcol="t_dx", min_visits=3):
    cf = d[COVS]
    v = visits[visits.PATNO.isin(d.index)].dropna(subset=["total", tcol])
    cnt = v.groupby("PATNO").size(); keep = cnt.index[cnt >= min_visits]
    v = v[v.PATNO.isin(keep)]; dd = d.loc[keep]
    C = np.column_stack([z(dd[c]) if c in CONT else dd[c].to_numpy(float) for c in COVS])
    idx = pd.Series(np.arange(len(dd)), index=dd.index); gi = idx.loc[v.PATNO].to_numpy()
    t = v[tcol].to_numpy(float)
    X = np.column_stack([np.ones(len(v)), t, C[gi], C[gi] * t[:, None]])
    f = lmm_fit(v.total.to_numpy(float), X, gi, t, "slope_corr")
    tau1 = f["re_cov"][1, 1]; s2 = f["sigma"] ** 2
    ssd = v.assign(t=t).groupby("PATNO")["t"].apply(lambda x: ((x - x.mean()) ** 2).sum())
    lam = tau1 / (tau1 + s2 / ssd)
    return {"n_donors": int(len(dd)), "slope_sd": float(np.sqrt(tau1)), "resid_sd": float(np.sqrt(s2)),
            "lambda_mean": float(lam.mean()), "lambda_median": float(lam.median()), "converged": f["converged"]}


def n_needed(r_obs, n_cov=6, alpha=0.05, power=0.8):
    if r_obs <= 0:
        return None
    zz = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    return float((zz / np.arctanh(r_obs)) ** 2 + 3 + n_cov)


def detectable_r(n, rel, n_cov=6, alpha=0.05, power=0.8):
    zz = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    r_obs = np.tanh(zz / np.sqrt(max(n - 3 - n_cov, 1)))
    return float(r_obs / np.sqrt(rel)) if rel > 0 else None


def design(icc, rel_out):
    out = {}
    for ph, R1 in icc.items():
        for oname, lam in rel_out.items():
            cells = []
            for k in K_RUNS:
                Rk = k * R1 / (1 + (k - 1) * R1) if R1 > 0 else 0.0
                for r in R_TRUE:
                    ro = r * np.sqrt(Rk * lam)
                    cells.append({"k_runs": k, "r_true": r, "reliability_product": float(Rk * lam),
                                  "r_observable": float(ro), "n_donors_80pct": n_needed(ro)})
            out[f"{ph}|{oname}"] = cells
    return out


# ---------------- main ----------------
def run(lines, genes, members, gates, s1, clin, ppmi_dir, progress=print):
    res = {"label": "EXPLORATORY", "version": "0.1.0"}
    d_all = donor_table(lines, s1, clin["donors"])
    d = d_all[d_all.n_cells_qc >= gates["A0_3"]["min_qc_cells"]]
    off = clin["updrs3"][clin["updrs3"].state == "OFF"]
    phenos = [p for p in ["P2", "P3", "P4"] if p in gates["A0_2_phenotypes_retained"]]

    progress("E1 gene-level")
    e1 = {}
    dg = d.join(genes, how="left")
    for k in ["P3", "P4"]:
        rows = []
        for g in members[k]:
            r = fit_long(dg, off, g, COVS, min_visits=A01_MIN_VISITS)
            if "estimate" in r:
                rows.append((g, r["estimate"], r["p"], r["r"]))
        t = pd.DataFrame(rows, columns=["gene", "beta3", "p", "r"])
        t["q"] = bh(t.p.to_numpy())
        comp = np.sign(fit_long(d, off, k, COVS, min_visits=A01_MIN_VISITS)["estimate"])
        e1[k] = {"n_genes": int(len(t)), "n_q05": int((t.q < 0.05).sum()), "n_p05": int((t.p < 0.05).sum()),
                 "expected_p05_under_null": float(0.05 * len(t)),
                 "share_negative_beta3": float((t.beta3 < 0).mean()), "composite_sign": float(comp),
                 "r_quantiles": {q: float(t.r.quantile(q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
                 "binomial_sign_test_p": float(stats.binomtest(int((t.beta3 < 0).sum()), len(t)).pvalue),
                 "top10": t.nsmallest(10, "p")[["gene", "beta3", "r", "p", "q"]].round(4).to_dict("records")}
    res["E1_gene_level"] = e1

    progress("E2 milestones")
    ms = read_milestones(ppmi_dir, d.index, clin["donors"].assign(PATNO=lambda x: x.PATNO.astype(str)))
    e2 = {"falls": "not available in the downloaded data"}
    tests = []
    for mname, (assess, col, thr) in ms.items():
        tab, n_prev = milestone_table(assess, col, thr)
        tab = tab.set_index("PATNO").join(d[phenos + COX_COVS], how="inner").dropna()
        info = {"n_donors": int(len(tab)), "n_prevalent_excluded": int(n_prev), "n_events": int(tab.event.sum()),
                "median_followup_y": float((tab.time - tab.entry).median()) if len(tab) else None}
        if info["n_events"] < MIN_EVENTS:
            info["descriptive_only"] = True
            e2[mname] = info
            continue
        for ph in phenos:
            X = np.column_stack([z(tab[ph]), z(tab["age_onset"]), tab["genetic_pd"], tab["media_e6"]])
            f = cox_fit(tab.entry, tab.time, tab.event, X)
            info[ph] = {"HR_per_SD": float(np.exp(f["beta"][0])),
                        "CI95": [float(np.exp(f["beta"][0] - 1.96 * f["se"][0])), float(np.exp(f["beta"][0] + 1.96 * f["se"][0]))],
                        "p": float(f["p"][0])}
            tests.append((mname, ph))
        e2[mname] = info
    if tests:
        q = bh([e2[m][p]["p"] for m, p in tests])
        for (m, p), qq in zip(tests, q):
            e2[m][p]["q"] = float(qq)
    res["E2_milestones"] = e2

    progress("E3 carriers")
    s1c = s1.copy(); s1c["PATNO"] = s1c.PATNO.astype(str).str.extract(r"(\d+)")[0]
    s1c = s1c.dropna(subset=["PATNO"]).groupby("PATNO").first()
    L = lines[~lines.is_control].copy(); L["PATNO"] = L.line.str.extract(r"(\d+)")[0]
    L = L[L.n_cells_qc >= gates["A0_3"]["min_qc_cells"]].set_index("PATNO")
    car = L.join(s1c[["Disease_Status", "Genetic_Status", "BATCH"]], how="inner")
    car = car[car.Disease_Status.isin(["Genetic PD", "Genetic Uffected"])].copy()
    car["manifest"] = (car.Disease_Status == "Genetic PD").astype(float)
    car["lrrk2"] = car.Genetic_Status.astype(str).str.contains("LRRK2").astype(float)
    car["media_e6"] = pd.to_numeric(car.BATCH, errors="coerce").isin([4, 5]).astype(float)
    sex = clin["all_sex"]
    car["sex_female"] = car.index.map(sex)
    e3 = {"n_manifesting": int(car.manifest.sum()), "n_non_manifesting": int((1 - car.manifest).sum()),
          "n_missing_sex": int(car.sex_female.isna().sum())}
    ps = []
    for ph in phenos:
        c = car.dropna(subset=[ph, "sex_female"])
        X = np.column_stack([np.ones(len(c)), c.manifest, c.lrrk2, c.media_e6, c.sex_female])
        o = ols(z(c[ph]), X)
        r = {"beta_SD": float(o["beta"][1]), "se": float(o["se"][1]), "p": float(o["p"][1]), "df": o["df"], "n": int(o["n"])}
        r.update(tost_partial_r(float(o["t"][1]), o["df"]))
        e3[ph] = r; ps.append(r["p"])
    for ph, qq in zip(phenos, bh(ps)):
        e3[ph]["q"] = float(qq)
    res["E3_carriers"] = e3

    progress("E4 design")
    icc = {k: v["icc"] for k, v in gates["A0_2"].items() if k in ("P1", "P2", "P3", "P4") and v["icc"] is not None}
    rel = {"UPDRS3_OFF": slope_reliability(d, off), "putamen_SBR": slope_reliability(d, clin["datscan"], min_visits=2)}
    res["E4_design"] = {
        "phenotype_reliability_icc": icc,
        "phenotype_icc_ci95": {k: gates["A0_2"][k]["ci95"] for k in icc},
        "outcome_slope_reliability": rel,
        "grid": design(icc, {k: v["lambda_mean"] for k, v in rel.items()}),
        "this_study_detectable_true_r": {
            f"{ph}|{o}": detectable_r(rel[o]["n_donors"], max(icc[ph], 0) * rel[o]["lambda_mean"])
            for ph in icc for o in rel},
    }
    return res


def main(argv=None):
    from ppmi_load import build, parse_date, load_colmap
    from scrn_gate import read_table
    ap = argparse.ArgumentParser()
    ap.add_argument("--phen-dir", default=str(ROOT / "data" / "derived" / "phenotypes"))
    ap.add_argument("--s1", required=True)
    ap.add_argument("--ppmi-dir", required=True)
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "exploratory.json"))
    a = ap.parse_args(argv)
    pdir = Path(a.phen_dir)
    lines = pd.read_csv(pdir / "lines.csv")
    genes = pd.read_csv(pdir / "genes_resid.csv", index_col=0)
    info = json.load(open(pdir / "phenotypes_info.json")); gates = json.load(open(pdir / "gates_prelink.json"))
    run2pat = lines.set_index("run")["line"].str.extract(r"(\d+)")[0]
    genes = genes.loc[lines.loc[~lines.is_control, "run"]]; genes.index = run2pat.loc[genes.index].to_numpy()
    s1 = read_table(Path(a.s1))
    d0 = donor_table(lines, s1, pd.DataFrame(columns=["PATNO", "sex_female", "age_onset", "dur_baseline", "saa_pos", "nsd_iss"]))
    clin = build(a.ppmi_dir, list(d0.index))
    # sex for all FOUNDIN donor lines (needed for E3 non-manifesting carriers)
    cm = load_colmap()["tables"]["demographics"]
    dem = pd.read_csv(next(Path(a.ppmi_dir).glob(cm["file"])), low_memory=False)
    dem["PATNO"] = dem[cm["PATNO"]].astype(str)
    clin["all_sex"] = dem.groupby("PATNO")[cm["sex"]].first().eq(cm["sex_female_value"]).astype(float).where(
        dem.groupby("PATNO")[cm["sex"]].first().notna())
    res = run(lines, genes, info["members"], gates, s1, clin, a.ppmi_dir)
    with open(a.out, "w") as f:
        json.dump(_clean(res), f, indent=1)
    print("written", a.out)


if __name__ == "__main__":
    main()
