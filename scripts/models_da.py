"""Cell-to-Clinic II (OSF osf.io/4pzxa; text = tag prereg-II-v1.0.0): linkage of DA-neuron phenotypes.

Run ONLY after (1) the registration is filed and (2) gate G-eQTL has run (data/derived/eqtl/g_eqtl.json
must exist). Model code is NOT re-implemented: donor_table, fit_long, fit_cross and family are imported
unchanged from scripts/models.py (prereg-v1.0.1): same REML LMM, Satterthwaite df, covariates, OFF-state
definition, time = years since diagnosis, fallbacks, z-scoring over the donors in each model, TOST
(+/-0.40, 90% CI) and the supported / equivalent / inconclusive classes.

Sample: registration-1 PD donors with a day-65 scRNA run, minus A0-3 (>= 1000 QC cells), minus G-cells
(>= 100 new-rule DA cells). A0-1 as registration 1.
Phenotypes (scripts/phenotypes_da.py, da-phen-v1.0.0): D2 (SNCA, DA pseudobulk), D4 (lysosome, DA pseudobulk).
  D1 and D3 failed G-icc: "not tested"; their H1 estimates are reported descriptively.
H1-DA (primary, single test, alpha 0.05): D2 x time on OFF-state MDS-UPDRS III (>= 3 OFF visits).
Secondary family, BH over 7: H1 D4 (only if A0-5 passes), H2 D2/D4 (putamen SBR), H3 D2/D4 by SAA,
  H3 D2/D4 by NSD-ISS. SAA tests are descriptive only if < 5 donors are SAA-negative (registration-1 rule).
A0-5 (D4 only): refit H1 with each KEGG lysosome member gene (gene-z in the DA pseudobulk; D4 is not
  residualised, so neither are its genes); D4 stays if >= 80% share the sign of D4's beta3.
H4 (single test, alpha 0.05): HALLMARK OXPHOS mean gene-z in day-65 bulk RNA (FOUNDIN RNAB), x time on
  OFF-state UPDRS III. "Composition-adjusted as in registration 1" is implemented literally as P3 was:
  genes with mean CPM >= 10 over the scRNA donor lines, z by those lines' mean and population SD, mean
  over set genes, then residualised across those lines on [1, logit Table S3 DA, logit Table S3
  progenitor sum]. One bulk sample per line (checked; the script stops otherwise).
G-eQTL: read from g_eqtl.json. If it failed, D2/D4 are still tested but labelled "marker-defined neurons".
Sensitivity (registered list): 1 old DA label (D*_old); 2 composition-adjusted (D*_adj); 3 idiopathic PD
  only; 4 five-level batch covariate; 5 leave-one-donor-out (H1-DA); 6 first 5 years after diagnosis
  (H1, H2). 1-4 rerun the primary test and the secondary family.
Output: data/derived/models_da.json (aggregates only).
"""
from __future__ import annotations

import argparse, json, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lmm import classify  # noqa: E402
from models import (COVS, A01_MIN_VISITS, A05_MIN_SHARE, donor_table, fit_long, fit_cross,  # noqa: E402
                    family, a01, _clean)

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.0.0"
MIN_DA = 100
TESTED, NOT_TESTED = ["D2", "D4"], ["D1", "D3"]
REGISTRATION = "osf.io/4pzxa"


def single(r):
    """A single pre-specified test at alpha 0.05 (no multiplicity correction)."""
    if "p" in r:
        r["class"] = classify(r["p"], r)
    return r


def bulk_oxphos(bulk_path, h5_path, s3_path, lines, cfg):
    """H4 phenotype per scRNA donor run (controls NaN); see module docstring."""
    from exploratory2 import load_bulk, to_symbols, bulk_scores
    from phenotypes import read_table_s3, parse_runs, logit_clip, residualise
    raw, meta = load_bulk(bulk_path)
    raw, _ = to_symbols(raw, h5_path)
    don = lines[~lines.is_control].copy()
    don["pat"] = don["line"].str.extract(r"(\d+)")[0]
    m = meta[(meta.day == 65) & (~meta.ctrl) & meta.line.isin(set(don.pat))]
    if m.line.duplicated().any() or len(m) != len(don):
        raise SystemExit("H4 needs exactly one day-65 bulk sample per scRNA donor line")
    sc, nkeep = bulk_scores(raw[m.col], list(m.col), cfg, extra_genes=())
    ox = sc["B-OXPHOS"].rename(index=dict(zip(m.col, m.line)))
    s3 = read_table_s3(s3_path, cfg).drop_duplicates("cdi").set_index("cdi")
    cdi = parse_runs(don["run"]).set_index("run")["cdi"]
    lo, hi = cfg["logit_clip"]
    comp = np.c_[logit_clip(s3.loc[cdi.loc[don.run], "s3_da"], lo, hi), logit_clip(s3.loc[cdi.loc[don.run], "s3_prog"], lo, hi)]
    y = ox.loc[don.pat].to_numpy(float)[:, None]
    adj = residualise({"fit": y, "apply": y}, comp, comp)[:, 0]
    out = pd.Series(np.nan, index=lines["run"])
    out.loc[don.run] = adj
    return out, {"n_lines": int(len(m)), "n_genes_kept": nkeep}


def secondary(d, off, dat, ph, use_batch=False, drop=()):
    s = {"H1_D4": fit_long(d, off, ph["D4"], COVS, min_visits=A01_MIN_VISITS, use_batch=use_batch)}
    for k in TESTED:
        s[f"H2_{k}"] = fit_long(d, dat, ph[k], COVS, use_batch=use_batch)
        s[f"H3_SAA_{k}"] = fit_cross(d, ph[k], "saa_pos", COVS, use_batch=use_batch)
        s[f"H3_NSDISS_{k}"] = fit_cross(d, ph[k], "nsd_iss", COVS, use_batch=use_batch)
    return family(s, drop)


def a05_d4(d, off, members, d4_sign):
    signs = []
    for g in members:
        r = fit_long(d, off, g, COVS, min_visits=A01_MIN_VISITS)
        if "estimate" in r:
            signs.append(np.sign(r["estimate"]))
    share = float(np.mean(np.array(signs) == d4_sign)) if signs else np.nan
    return {"n_genes": len(signs), "share_same_sign": share, "pass": bool(share >= A05_MIN_SHARE),
            "rule": f">= {A05_MIN_SHARE:.0%} of member genes share the sign of D4's beta3"}


def run(lines, s1, clin, geqtl, members_d4, min_qc_cells, progress=print):
    """lines: registration-1 lines.csv joined (on run) with lines_da.csv, the D4 member gene columns
    (prefixed 'g_') and H4 ('H4_bulk_oxphos')."""
    res = {"version": VERSION, "registration": REGISTRATION}
    d = donor_table(lines, s1, clin["donors"])
    res["sample"] = {"n_pd_with_scrna": int(len(d))}
    d = d[d.n_cells_qc >= min_qc_cells]
    res["sample"]["n_after_A0_3"] = int(len(d))
    d = d[d.n_da_v2 >= MIN_DA]
    res["sample"]["n_after_G_cells"] = int(len(d))
    res["sample"]["n_missing_covariate"] = int(d[COVS].isna().any(axis=1).sum())
    res["G_eQTL"] = {"pass": geqtl.get("G_eQTL_pass"),
                     "label": "DA neurons" if geqtl.get("G_eQTL_pass") else "marker-defined neurons"}
    res["A0_1"] = a01(d, clin["updrs3"])
    if not res["A0_1"]["pass"]:
        res["stopped"] = "A0-1 failed: data audit only"
        return res
    off = clin["updrs3"][clin["updrs3"].state == "OFF"]
    dat = clin["datscan"]

    progress("H1-DA primary")
    res["H1_DA_primary"] = single(fit_long(d, off, "D2", COVS, min_visits=A01_MIN_VISITS))
    progress("A0-5 (D4)")
    h1d4 = fit_long(d, off, "D4", COVS, min_visits=A01_MIN_VISITS)
    res["A0_5_D4"] = a05_d4(d, off, members_d4, np.sign(h1d4["estimate"])) if "estimate" in h1d4 else {"pass": False, "error": "H1 D4 not estimable"}
    drop = () if res["A0_5_D4"]["pass"] else ("H1_D4",)
    progress("secondary family")
    ph = {k: k for k in TESTED}
    res["secondary"] = secondary(d, off, dat, ph, drop=drop)
    res["secondary_n_in_BH"] = sum("q" in r for r in res["secondary"].values())
    progress("H4 bulk OXPHOS")
    res["H4_bulk_oxphos"] = single(fit_long(d, off, "H4_bulk_oxphos", COVS, min_visits=A01_MIN_VISITS))
    res["not_tested_descriptive_H1"] = {k: fit_long(d, off, k, COVS, min_visits=A01_MIN_VISITS) for k in NOT_TESTED}

    progress("sensitivity 1-4")
    sens = {}
    variants = {"S1_old_label": dict(dd=d, ph={k: f"{k}_old" for k in TESTED}),
                "S2_composition_adjusted": dict(dd=d, ph={k: f"{k}_adj" for k in TESTED}),
                "S3_idiopathic": dict(dd=d[d.genetic_pd == 0], ph=ph),
                "S4_batch": dict(dd=d, ph=ph, use_batch=True)}
    for name, v in variants.items():
        ub = v.get("use_batch", False)
        sens[name] = {"H1_DA_primary": single(fit_long(v["dd"], off, v["ph"]["D2"], COVS, min_visits=A01_MIN_VISITS, use_batch=ub)),
                      "secondary": secondary(v["dd"], off, dat, v["ph"], use_batch=ub, drop=drop)}
    progress("sensitivity 5 (leave-one-donor-out)")
    full = res["H1_DA_primary"]
    if "estimate" in full:
        ids = fit_long(d, off, "D2", COVS, min_visits=A01_MIN_VISITS, keep_fit=True)["_ids"]
        est = np.array([r["estimate"] for i in ids
                        if "estimate" in (r := fit_long(d.drop(index=i), off, "D2", COVS, min_visits=A01_MIN_VISITS))])
        sens["S5_leave_one_out"] = {"beta3_full": full["estimate"], "loo_min": float(est.min()), "loo_max": float(est.max()),
                                    "n_sign_changes": int((np.sign(est) != np.sign(full["estimate"])).sum()),
                                    "max_abs_change_in_se_units": float(np.max(np.abs(est - full["estimate"])) / full["se"])}
    progress("sensitivity 6 (first 5 years)")
    off5, dat5 = off[off.t_dx <= 5], dat[dat.t_dx <= 5]
    sens["S6_first_5y"] = {"H1_DA_primary": single(fit_long(d, off5, "D2", COVS, min_visits=A01_MIN_VISITS)),
                           "H1_D4": fit_long(d, off5, "D4", COVS, min_visits=A01_MIN_VISITS),
                           "H2": {k: fit_long(d, dat5, k, COVS) for k in TESTED}}
    res["sensitivity"] = sens
    return res


def assemble(phen_dir, phen_da_dir, h4=None):
    lines = pd.read_csv(Path(phen_dir) / "lines.csv")
    da = pd.read_csv(Path(phen_da_dir) / "lines_da.csv").drop(columns=["line", "is_control"])
    genes = pd.read_csv(Path(phen_da_dir) / "genes_da.csv", index_col=0)
    from phenotypes import read_geneset, load_config
    members = [g for g in read_geneset(load_config()["gene_sets"]["P4"]) if g in genes.columns]
    g = genes[members].add_prefix("g_")
    out = lines.merge(da, on="run", how="left").merge(g, left_on="run", right_index=True, how="left")
    if h4 is not None:
        out["H4_bulk_oxphos"] = out["run"].map(h4)
    return out, ["g_" + m for m in members]


def main(argv=None):
    from phenotypes import load_config
    from ppmi_load import build
    from scrn_gate import read_table
    ap = argparse.ArgumentParser()
    ap.add_argument("--phen-dir", default=str(ROOT / "data" / "derived" / "phenotypes"))
    ap.add_argument("--phen-da-dir", default=str(ROOT / "data" / "derived" / "phenotypes_da"))
    ap.add_argument("--geqtl", default=str(ROOT / "data" / "derived" / "eqtl" / "g_eqtl.json"))
    ap.add_argument("--s1", required=True)
    ap.add_argument("--s3", required=True)
    ap.add_argument("--ppmi-dir", required=True)
    ap.add_argument("--bulk", default=str(ROOT / "data" / "unpacked_rnab" / "FOUNDIN-PD_150.9_RNAB" / "aggregated_expression" / "countTable.tsv"))
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "models_da.json"))
    a = ap.parse_args(argv)
    if not Path(a.geqtl).exists():
        raise SystemExit("linkage is held until gate G-eQTL has run (g_eqtl.json missing)")
    geqtl = json.load(open(a.geqtl))
    cfg = load_config()
    lines0 = pd.read_csv(Path(a.phen_dir) / "lines.csv")
    h5 = next((ROOT / "data" / "unpacked").rglob("filtered_feature_bc_matrix.h5"))
    h4, h4info = bulk_oxphos(a.bulk, h5, a.s3, lines0, cfg)
    lines, members = assemble(a.phen_dir, a.phen_da_dir, h4)
    gates = json.load(open(Path(a.phen_dir) / "gates_prelink.json"))
    s1 = read_table(Path(a.s1))
    d0 = donor_table(lines, s1, pd.DataFrame(columns=["PATNO", "sex_female", "age_onset", "dur_baseline", "saa_pos", "nsd_iss"]))
    clin = build(a.ppmi_dir, list(d0.index))
    res = run(lines, s1, clin, geqtl, members, gates["A0_3"]["min_qc_cells"])
    res["H4_phenotype"] = h4info
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    with open(a.out, "w") as f:
        json.dump(_clean(res), f, indent=1)
    print("written", a.out)


if __name__ == "__main__":
    main()
