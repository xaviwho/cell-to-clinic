"""End-to-end synthetic test: phenotypes -> gates_prelink -> ppmi_load -> models. No real data.

Builds fake scrn-v0.1.0 outputs, fake Tables S1/S3 and fake PPMI-format CSVs in a temp directory,
plants a P1 x time effect on OFF-state UPDRS III (and nothing for P2-P4), and checks that the
pipeline recovers it and that the gates behave as specified.
Run: C2C_SYNTHETIC=1 python -I tests/test_pipeline_synthetic.py <empty_tmp_dir> [--fast]
"""
import gzip, json, os, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd

import phenotypes, gates_prelink, ppmi_load, models

rng = np.random.default_rng(42)
TMP = Path(sys.argv[1]); FAST = "--fast" in sys.argv
assert TMP.exists() and not any(TMP.iterdir()), "give an empty directory"
cfg = phenotypes.load_config()

# ---------------- donors ----------------
n_pd_i, n_pd_g, n_other = 28, 19, 32
n_don = n_pd_i + n_pd_g + n_other
pat = [str(90000 + i) for i in range(n_don)]
status = ["Idiopathic PD"] * n_pd_i + ["Genetic PD"] * n_pd_g + ["Healthy Control"] * 20 + ["Prodromal"] * 12
batch = rng.integers(1, 6, n_don)
lat_da = rng.normal(size=n_don)          # drives DA fraction (P1)
lat_ox = rng.normal(size=n_don)          # drives OXPHOS genes (P3)
ctrl_pat, ctrl_cdi = "99999", "0001"

# ---------------- Table S3 ----------------
def comp(l):
    a = np.exp(np.r_[rng.normal(0, .3, 5) + [0, 1.5, 0, 3, 3], -0.5 + 0, 2.5 + 0.8 * l])
    return 100 * a / a.sum()
cols = cfg["table_s3"]["progenitor_columns"][:2] + ["Ependymal-like Cells"] + cfg["table_s3"]["progenitor_columns"][2:] + \
       ["Immature Dopaminergic Neurons", "Dopaminergic Neurons"]
s3_rows = []
cdis = [f"{1000 + i:04d}" for i in range(n_don)]
for i in range(n_don):
    c = comp(lat_da[i]); s3_rows.append([f"FOUNDIN_{cdis[i]}", *c, c[5] + c[6]])
ctrl_lat = 0.3
for b in range(4):
    c = comp(ctrl_lat + rng.normal(0, 0.1)); s3_rows.append([f"FOUNDIN_{ctrl_cdi}", *c, c[5] + c[6]])
s3 = pd.DataFrame(s3_rows, columns=["FOUNDIN_ID", *cols, "totalDopa"])
s3_path = TMP / "s3.xlsx"
with pd.ExcelWriter(s3_path) as w:
    pd.DataFrame([["Supplemental Table 3. synthetic"]]).to_excel(w, header=False, index=False, startrow=0)
    s3.to_excel(w, index=False, startrow=1)

# ---------------- scRNA outputs ----------------
p3 = phenotypes.read_geneset(cfg["gene_sets"]["P3"]); p4 = phenotypes.read_geneset(cfg["gene_sets"]["P4"])
genes = sorted(set(p3) | set(p4) | {"SNCA"} | {f"FILL{i}" for i in range(300)})
base = rng.normal(np.log(30), 1.2, len(genes))
ox_idx = np.array([g in p3 for g in genes])
scrn = TMP / "scrn"; (scrn / "runs").mkdir(parents=True)
summary = []
def write_run(run, lat_o, lat_d, ncell):
    mu = np.exp(base + 0.3 * lat_o * ox_idx + rng.normal(0, 0.05, len(genes)) + 0.1 * lat_d * (np.array(genes) == "SNCA"))
    cnt = rng.poisson(mu * ncell / 2)
    pb = pd.DataFrame({"gene": genes, "DA_neuron": cnt // 3, "Ependymal_like": 0, "Other_neuron": 0, "Progenitor": cnt - cnt // 3,
                       "Proliferating": 0, "Unassigned": 0, "All_QC_pass": cnt})
    (scrn / "runs" / run).mkdir()
    pb.to_csv(scrn / "runs" / run / "pseudobulk.csv.gz", index=False)
    da = 1 / (1 + np.exp(-(-1 + 0.8 * lat_d)))
    summary.append({"run": run, "n_cells_qc": ncell, "frac_TH_pos_all": da * 0.9, "frac_TH_pos_all_sens": da * 0.7,
                    "frac_DA_neuron": da, "frac_Progenitor": 0.4, "frac_Proliferating": 0.05})
ncells = rng.integers(1500, 8000, n_don); ncells[3] = 700   # one PD donor fails A0-3
for i in range(n_don):
    write_run(f"SCRN_PPMI{pat[i]}_{cdis[i]}_da65_v1", lat_ox[i], lat_da[i], int(ncells[i]))
for b in ["B1", "B2", "B3", "B5E6"]:
    write_run(f"SCRN_PPMI{ctrl_pat}{b}_{ctrl_cdi}_da65_v1", 0.2 + rng.normal(0, 0.1), ctrl_lat, 4000)
pd.DataFrame(summary).to_csv(scrn / "summary.csv", index=False)

# ---------------- Table S1 ----------------
s1 = pd.DataFrame({"PATNO": pat, "BATCH": batch, "Disease_Status": status,
                   "PD_SCORE_Nalls2019_no_GBA_LRRK2": [f"{v:.4f}".replace(".", ",") for v in rng.normal(size=n_don)]})

# ---------------- phenotypes + gates ----------------
t0 = time.time()
lines, genes_resid, s3_ctrl, info = phenotypes.compute(scrn, s3_path, cfg)
pdir = TMP / "phen"; pdir.mkdir()
lines.to_csv(pdir / "lines.csv", index=False); s3_ctrl.to_csv(pdir / "control_s3_rows.csv", index=False)
don = lines[~lines.is_control].set_index("line")
exp_p1 = phenotypes.logit_clip(s3.iloc[:n_don]["Dopaminergic Neurons"].to_numpy() / 100, .005, .995)
assert np.allclose(don.loc[[f"PPMI{p}" for p in pat], "P1"], exp_p1), "P1 mismatch"
assert info["P3_genes_used"] > 150 and info["P4_genes_used"] > 80 and info["P2_gene_kept"]
r = np.corrcoef(don.loc[[f"PPMI{p}" for p in pat], "P3"], lat_ox)[0, 1]
print(f"phenotypes ok: P1 exact; corr(P3, planted OXPHOS factor) = {r:.2f}; residual means ~0: "
      f"{np.abs(don[['P2','P3','P4']].mean()).max():.1e}")
assert r > 0.8
g = gates_prelink.run(pdir, cfg)
print("A0-2:", {k: round(v["icc"], 2) for k, v in g["A0_2"].items() if v["icc"] is not None and v["gating"]}, "A0-3:", g["A0_3"])
assert g["A0_3"]["n_donor_lines_excluded"] == 1 and g["A0_2"]["P1"]["icc"] > 0.5 and g["A0_2"]["P3"]["icc"] > 0.5

# ---------------- fake PPMI files ----------------
ppmi = TMP / "ppmi"; ppmi.mkdir()
pdi = [i for i in range(n_don) if "PD" in status[i] and status[i] != "Prodromal"]
pd_pats = [pat[i] for i in pdi]
def md(y, m): return f"{m:02d}/{y}"
dem, dxh, u3, u2, dsc, led, saa, nsd = [], [], [], [], [], [], [], []
z1 = (lat_da[pdi] - lat_da[pdi].mean()) / lat_da[pdi].std(ddof=1)
truth = {}
for k, i in enumerate(pdi):
    p = pat[i]; by = int(rng.integers(1940, 1965)); dxy = int(rng.integers(2005, 2012)); bly = dxy + int(rng.integers(0, 3))
    dem.append({"PATNO": p, "SEX": int(rng.integers(0, 2)), "BIRTHDT": md(by, 6)})
    dxh.append({"PATNO": p, "PDDXDT": md(dxy, 3), "SXDT": md(dxy - 1, 3)})
    slope = 2.0 + 0.9 * z1[k] + rng.normal(0, 0.8); icpt = 20 + rng.normal(0, 4)
    nvis = int(rng.integers(3, 10)) if k != 0 else 2   # donor 0 fails the >=3 OFF visit rule
    truth[p] = nvis
    for v in range(nvis):
        y, ev = bly + v, ("BL" if v == 0 else f"V{v:02d}")
        tdx = (y - dxy) + (3 - 3) / 12
        val = icpt + slope * tdx + rng.normal(0, 3)
        treated = int(v >= 2)
        u3.append({"PATNO": p, "EVENT_ID": ev, "INFODT": md(y, 3), "NP3TOT": round(val), "PDTRTMNT": treated, "PDSTATE": "OFF" if treated else ""})
        if treated:
            u3.append({"PATNO": p, "EVENT_ID": ev, "INFODT": md(y, 3), "NP3TOT": round(val - 8), "PDTRTMNT": 1, "PDSTATE": "ON"})
            led.append({"PATNO": p, "EVENT_ID": ev, "INFODT": md(y, 3), "LEDD": 300 + 50 * v})
        u2.append({"PATNO": p, "EVENT_ID": ev, "INFODT": md(y, 3), "NP2PTOT": round(5 + 1.0 * tdx + rng.normal(0, 2))})
        if v % 2 == 0:
            sbr = 1.2 - 0.05 * tdx + rng.normal(0, 0.1)
            dsc.append({"PATNO": p, "INFODT": md(y, 4), "DATSCAN_PUTAMEN_L": sbr, "DATSCAN_PUTAMEN_R": sbr + 0.02})
    saa.append({"PATNO": p, "RUNDATE": md(bly + 1, 1), "SAA_Status": "Negative" if k % 6 == 0 else "Positive"})
    saa.append({"PATNO": p, "RUNDATE": md(bly + 3, 1), "SAA_Status": "Positive"})
    nsd.append({"PATNO": p, "INFODT": md(bly, 5), "NSD_ISS_STAGE": "2B"})
    nsd.append({"PATNO": p, "INFODT": md(bly + 2, 5), "NSD_ISS_STAGE": "3" if k % 2 else "2B"})
# a non-FOUNDIN participant that must be dropped on read
dem.append({"PATNO": "12345", "SEX": 1, "BIRTHDT": md(1950, 1)})
for name, rows in [("Demographics_2026.csv", dem), ("PD_Diagnosis_History_2026.csv", dxh), ("MDS-UPDRS_Part_III_2026.csv", u3),
                   ("MDS_UPDRS_Part_II__Patient_Questionnaire_2026.csv", u2), ("DaTScan_Analysis_2026.csv", dsc),
                   ("LEDD_2026.csv", led), ("SAA_2026.csv", saa), ("PPMI_NSD_ISS_2026.csv", nsd)]:
    pd.DataFrame(rows).to_csv(ppmi / name, index=False)

clin = ppmi_load.build(ppmi, pd_pats)
dn = clin["donors"].set_index("PATNO")
assert "12345" not in dn.index and len(dn) == len(pd_pats)
off = clin["updrs3"][clin["updrs3"].state == "OFF"]
assert all(off.groupby("PATNO").size().loc[p] == truth[p] for p in pd_pats), "OFF visit count"
assert (clin["updrs3"].state == "ON").sum() == sum(max(0, n - 2) for n in truth.values())
assert dn["nsd_iss"].isin([2.0, 3.0]).all() and dn["saa_pos"].sum() == sum(1 for k in range(len(pdi)) if k % 6)
assert np.allclose(dn["dur_baseline"].dropna().min(), 0, atol=0.01)
assert clin["updrs3"]["ledd"].notna().sum() == (clin["updrs3"].state == "ON").sum() + (clin["updrs3"].state == "OFF")[clin["updrs3"].event.isin([f"V{v:02d}" for v in range(2, 10)])].sum()
print("ppmi_load ok: OFF/ON split, earliest SAA, latest NSD-ISS, LEDD by EVENT_ID, non-FOUNDIN rows dropped")

# ---------------- models ----------------
genes_d = genes_resid.loc[lines.loc[~lines.is_control, "run"]]
genes_d.index = lines.loc[~lines.is_control, "line"].str.extract(r"(\d+)")[0].to_numpy()
members = info["members"]
if FAST:
    members = {k: v[:20] for k, v in members.items()}
res = models.run(lines, genes_d, members, g, s1, clin, progress=lambda s: print("  ..", s, f"{time.time()-t0:.0f}s"))
res = models._clean(res)
json.dump(res, open(TMP / "models_synthetic.json", "w"), indent=1)
print("sample:", res["sample"], "A0-1:", res["A0_1"])
for fam in ["H1", "H2", "H3_SAA", "H3_NSDISS"]:
    print(fam, {k: (round(v.get("estimate", np.nan), 2), round(v.get("q", np.nan), 3), v.get("class"), v.get("re", "")) for k, v in res[fam].items()})
print("A0-5:", res["A0_5"])
print("S5:", {k: v["n_sign_changes"] for k, v in res["sensitivity"]["S5_leave_one_out"].items()})
for s in ["S1_idiopathic", "S2_batch", "S3_P1_alternatives", "S6_on_state_ledd", "S7_first_5y", "S8_part2", "S9_time_from_baseline"]:
    print(s, {k: (round(v.get("estimate", np.nan), 2), v.get("class")) for k, v in res["sensitivity"][s]["H1"].items()})
assert res["sample"]["n_pd_with_scrna"] == 47 and res["sample"]["n_after_A0_3"] == 46
assert res["H1"]["P1"]["class"] == "supported" and res["H1"]["P1"]["estimate"] > 0
assert res["H1"]["P1"]["n_donors"] == 45            # A0-3 drop + one donor with 2 OFF visits
assert res["H3_SAA"]["P1"].get("n_saa_neg", 0) >= 5
assert res["sensitivity"]["S3_P1_alternatives"]["H1"]["S3a_thpos"]["class"] == "supported"
print(f"ALL PIPELINE TESTS PASSED ({time.time()-t0:.0f}s)")
