"""Synthetic test for scripts/models_da.py (Cell-to-Clinic II linkage). No real data.

Plants a D2 x time effect on OFF-state UPDRS III and nothing for D4/H4, and checks: sample rules
(A0-3, G-cells), the primary single test, the secondary BH family size, A0-5 dropping, the G-eQTL
label, and that sensitivity analyses run.
Run: python -I tests/test_models_da.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd

import models_da as md


def synth(seed=1, n_pd=50, effect=1.2, d4_genes_coherent=True):
    rng = np.random.default_rng(seed)
    n = n_pd + 6
    pat = [str(70000 + i) for i in range(n)]
    status = ["Idiopathic PD"] * 30 + ["Genetic PD"] * (n_pd - 30) + ["Healthy Control"] * 6
    s1 = pd.DataFrame({"PATNO": pat, "Disease_Status": status, "PD_SCORE_Nalls2019_no_GBA_LRRK2": rng.normal(size=n),
                       "BATCH": rng.integers(1, 6, n)})
    D = rng.normal(size=(n, 4))
    lines = pd.DataFrame({"run": [f"R{i}" for i in range(n)], "line": [f"PPMI{p}" for p in pat], "is_control": False,
                          "n_cells_qc": 5000, "n_da_v2": 400, "D1": D[:, 0], "D2": D[:, 1], "D3": D[:, 2], "D4": D[:, 3]})
    lines.loc[0, "n_cells_qc"] = 500          # removed by A0-3
    lines.loc[1, "n_da_v2"] = 50              # removed by G-cells
    for k in ("D2", "D4"):
        lines[f"{k}_old"] = lines[k] + rng.normal(0, 0.3, n)
        lines[f"{k}_adj"] = lines[k] + rng.normal(0, 0.3, n)
    for j in range(10):
        sgn = 1 if (d4_genes_coherent or j < 5) else -1
        lines[f"g_G{j}"] = sgn * lines["D4"] + rng.normal(0, 0.2, n)
    lines["H4_bulk_oxphos"] = rng.normal(size=n)
    lines = pd.concat([lines, pd.DataFrame({"run": ["C1", "C2"], "line": ["PPMI99999"] * 2, "is_control": True,
                                            "n_cells_qc": 5000, "n_da_v2": 400})], ignore_index=True)
    donors = pd.DataFrame({"PATNO": pat, "sex_female": rng.integers(0, 2, n).astype(float), "age_onset": rng.normal(60, 8, n),
                           "dur_baseline": rng.gamma(2, 1, n), "saa_pos": (rng.random(n) < 0.85).astype(float),
                           "nsd_iss": rng.integers(1, 4, n).astype(float)})
    rows, scans = [], []
    for i, p in enumerate(pat):
        t0 = rng.uniform(0, 3); slope = 2 + effect * D[i, 1] + rng.normal(0, 0.7); icpt = rng.normal(20, 6)
        for v in range(6):
            t = t0 + v
            rows.append({"PATNO": p, "state": "OFF", "total": icpt + slope * t + rng.normal(0, 3), "t_dx": t, "t_bl": v})
            scans.append({"PATNO": p, "total": 1.5 - 0.08 * t + rng.normal(0, 0.1), "t_dx": t, "t_bl": v})
    clin = {"donors": donors, "updrs3": pd.DataFrame(rows), "datscan": pd.DataFrame(scans)}
    return lines, s1, clin, [f"g_G{j}" for j in range(10)]


def test_planted_d2_effect():
    lines, s1, clin, mem = synth()
    r = md.run(lines, s1, clin, {"G_eQTL_pass": True}, mem, 1000, progress=lambda *_: None)
    assert r["sample"] == {"n_pd_with_scrna": 50, "n_after_A0_3": 49, "n_after_G_cells": 48, "n_missing_covariate": 0}, r["sample"]
    assert r["A0_1"]["pass"] and r["G_eQTL"]["label"] == "DA neurons"
    p = r["H1_DA_primary"]
    assert p["class"] == "supported" and p["estimate"] > 0 and "q" not in p
    assert r["A0_5_D4"]["pass"] and r["secondary_n_in_BH"] == 7
    assert r["secondary"]["H1_D4"]["class"] != "supported"
    assert r["H4_bulk_oxphos"]["class"] in ("equivalent", "inconclusive")
    assert set(r["not_tested_descriptive_H1"]) == {"D1", "D3"} and "class" not in r["not_tested_descriptive_H1"]["D1"]
    s = r["sensitivity"]
    assert set(s) == {"S1_old_label", "S2_composition_adjusted", "S3_idiopathic", "S4_batch", "S5_leave_one_out", "S6_first_5y"}
    assert s["S3_idiopathic"]["H1_DA_primary"]["n_donors"] < p["n_donors"]
    assert s["S5_leave_one_out"]["n_sign_changes"] == 0


def test_a05_fail_drops_d4_and_geqtl_fail_label():
    lines, s1, clin, mem = synth(seed=2, effect=0.0, d4_genes_coherent=False)
    # make D4 carry a real slope so its sign is defined, while half its genes are flipped
    r = md.run(lines, s1, clin, {"G_eQTL_pass": False}, mem, 1000, progress=lambda *_: None)
    assert r["G_eQTL"]["label"] == "marker-defined neurons"
    assert not r["A0_5_D4"]["pass"] and r["A0_5_D4"]["share_same_sign"] <= 0.6
    assert r["secondary"]["H1_D4"]["class"] == "dropped_by_A0-5" and r["secondary_n_in_BH"] == 6


def test_bulk_oxphos_h4_phenotype():
    """Fake RNAB count table, Cell Ranger feature list and Table S3: H4 recovers a planted OXPHOS
    latent and is orthogonal to the Table S3 composition it is residualised on."""
    import tempfile, h5py
    from phenotypes import load_config, read_geneset
    cfg = load_config()
    rng = np.random.default_rng(5)
    ox = read_geneset(cfg["gene_sets"]["P3"])
    genes = ox + [f"FILL{i}" for i in range(200)]
    ens = [f"ENSG{i:011d}.{1 + i % 3}" for i in range(len(genes))]
    n = 40
    pats = [str(80000 + i) for i in range(n)]; cdis = [f"{3000 + i:04d}" for i in range(n)]
    lat, da, prog = rng.normal(size=n), rng.uniform(0.1, 0.5, n), rng.uniform(0.1, 0.4, n)
    base = rng.normal(np.log(80), 1, len(genes)); isox = np.array([g in ox for g in genes])
    cols, data = [], []
    for i in range(n):
        mu = np.exp(base + isox * (0.3 * lat[i] + 0.5 * np.log(da[i] / (1 - da[i]))))
        cols.append(f"RNAB_PPMI{pats[i]}_{cdis[i]}_da65_v1"); data.append(rng.poisson(mu * 50))
    cols.append(f"RNAB_PPMI{pats[0]}_{cdis[0]}_da25_v1"); data.append(rng.poisson(np.exp(base) * 50))  # other day: ignored
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pd.DataFrame(np.array(data).T, index=ens, columns=cols).to_csv(tmp / "counts.tsv", sep="\t")
        with h5py.File(tmp / "f.h5", "w") as f:
            f["matrix/features/id"] = np.array([e.split(".")[0] for e in ens], dtype="S")
            f["matrix/features/name"] = np.array(genes, dtype="S")
        c = cfg["table_s3"]
        s3 = pd.DataFrame({c["id_column"]: [f"FOUNDIN_{x}" for x in cdis], c["P1_da_column"]: 100 * da,
                           c["totaldopa_column"]: 100 * da})
        for j, pc in enumerate(c["progenitor_columns"]):
            s3[pc] = 100 * prog / len(c["progenitor_columns"])
        with pd.ExcelWriter(tmp / "s3.xlsx") as w:
            pd.DataFrame([["synthetic"]]).to_excel(w, header=False, index=False)
            s3.to_excel(w, index=False, startrow=1)
        lines = pd.DataFrame({"run": [f"SCRN_PPMI{p}_{c_}_da65_v1" for p, c_ in zip(pats, cdis)],
                              "line": [f"PPMI{p}" for p in pats], "is_control": False})
        lines = pd.concat([lines, pd.DataFrame({"run": ["SCRN_PPMI99999B1_0001_da65_v1"], "line": ["PPMI99999"], "is_control": [True]})],
                          ignore_index=True)
        h4, info = md.bulk_oxphos(tmp / "counts.tsv", tmp / "f.h5", tmp / "s3.xlsx", lines, cfg)
    v = h4.loc[lines.run[:n]].to_numpy()
    assert info["n_lines"] == n and np.isnan(h4.iloc[-1])
    assert abs(np.corrcoef(v, np.log(da / (1 - da)))[0, 1]) < 1e-6 and abs(np.corrcoef(v, np.log(prog / (1 - prog)))[0, 1]) < 1e-6
    assert np.corrcoef(v, lat)[0, 1] > 0.9


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f(); print("ok", k, flush=True)
    print("ALL PASSED")
