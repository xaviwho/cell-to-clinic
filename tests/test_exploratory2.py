"""Synthetic checks for scripts/exploratory2.py. Run: python -I tests/test_exploratory2.py <empty_tmp_dir>"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import h5py
import numpy as np
import pandas as pd
import exploratory2 as ex
from phenotypes import load_config, read_geneset

TMP = pathlib.Path(sys.argv[1]); assert TMP.exists() and not any(TMP.iterdir())
rng = np.random.default_rng(3)
cfg = load_config()
genes = sorted(set(read_geneset(cfg["gene_sets"]["P3"])) | set(read_geneset(cfg["gene_sets"]["P4"])) | {"SNCA", "TH"} | {f"F{i}" for i in range(200)})
ids = [f"ENSG{i:011d}.{i % 7 + 1}" for i in range(len(genes))]  # versioned, as in Cell Ranger
with h5py.File(TMP / "f.h5", "w") as f:
    g = f.create_group("matrix/features")
    g["id"] = np.array([x.encode() for x in ids]); g["name"] = np.array([x.encode() for x in genes])
ox = np.array([x in read_geneset(cfg["gene_sets"]["P3"]) for x in genes])
base = rng.normal(np.log(50), 1, len(genes))
cols, data = [], []
def sample(name, line_eff, run_eff):
    mu = np.exp(base + ox * (line_eff + run_eff) + (np.array(genes) == "SNCA") * line_eff)
    cols.append(name); data.append(rng.poisson(mu * 20))
line_eff = rng.normal(0, 0.5, 60)
for day in (65, 25):
    for i in range(60):
        sample(f"RNAB_PPMI{3000+i}_{1000+i}_da{day}_v1", line_eff[i], rng.normal(0, 0.1))
    # control: batch-1 replicates share a batch effect; other batches add their own
    b1 = rng.normal(0, 0.3)
    reps = ["B1"] * (9 if day == 25 else 1) + ["B2", "B3", "B5E6", "B5E8"]
    for j, b in enumerate(reps):
        be = b1 if b == "B1" else rng.normal(0, 0.3)
        sample(f"RNAB_PPMI9999{b}_0001_da{day}_v{j+1}", 0.2, be + rng.normal(0, 0.05))
ct = pd.DataFrame(np.array(data).T, index=ids, columns=cols)
ct.index.name = "Geneid"
ct.to_csv(TMP / "counts.tsv", sep="\t")
res = ex.e5b_e5c(TMP / "counts.tsv", [str(3000 + i) for i in range(50)], cfg, TMP / "f.h5")
print("mapping", res["bulk_id_mapping"])
b = res["E5b_day65"]
print("E5b ICC scrn subset:", {k: round(v["icc"], 2) for k, v in b["scrn_79"].items()}, "n_donor", b["scrn_79"]["B-SNCA"]["n_donor"], "n_ctrl", b["n_control_runs"])
assert b["scrn_79"]["B-SNCA"]["n_donor"] == 50 and b["all"]["B-SNCA"]["n_donor"] == 60 and b["n_control_runs"] == 5
assert b["all"]["B-SNCA"]["icc"] > 0.8          # SNCA: line effect, no run effect
c = res["E5c_day25"]["B-OXPHOS"]
print("E5c OXPHOS:", {k: round(v, 3) if isinstance(v, float) else v for k, v in c.items()})
assert c["n_batch1_runs"] == 9 and c["between_batch_share"] > 0.5 and c["icc_within_batch"] > c["icc_all_runs"]

# E5a with small B on synthetic clinical data
ex.B, ex.B_LAMBDA = 200, 15
assert res["bulk_id_mapping"]["n_mapped"] == res["bulk_id_mapping"]["n_ids"]
n = 40
d = pd.DataFrame({"age_onset": rng.normal(60, 8, n), "sex_female": rng.integers(0, 2, n).astype(float),
                  "genetic_pd": rng.integers(0, 2, n).astype(float), "prs": rng.normal(size=n),
                  "media_e6": rng.integers(0, 2, n).astype(float), "dur_baseline": rng.uniform(0, 2, n)},
                 index=[f"p{i}" for i in range(n)])
rows = []
for i, p in enumerate(d.index):
    s = 2 + rng.normal(0, 1.5)
    for t in np.arange(rng.integers(3, 9)) + rng.uniform(0, 2):
        rows.append((p, t, 20 + s * t + rng.normal(0, 5), "OFF"))
v = pd.DataFrame(rows, columns=["PATNO", "t_dx", "total", "state"])
clin = {"updrs3": v, "datscan": v.assign(total=v.total / 100)}
gates = {"A0_2": {p: {"icc": ic, "n_donor": 79, "n_within": 4} for p, ic in
                  [("P1", 0.33), ("P2", 0.94), ("P3", 0.55), ("P4", 0.64)]}}
r = ex.e5a(gates, d, clin, rng)
print("E5a lambda:", r["lambda_boot"]); print("E5a P2 UPDRS k1 r0.3:", r["P2"]["UPDRS3_OFF|k1|r0.3"], "P1 P(icc<=0):", r["P1"]["P_icc_le_0"])
assert 0 < r["lambda_boot"]["UPDRS3_OFF"]["median"] < 1 and r["P1"]["P_icc_le_0"] > r["P2"]["P_icc_le_0"]
assert r["P2"]["UPDRS3_OFF|k1|r0.3"]["N_p10"] <= r["P2"]["UPDRS3_OFF|k1|r0.3"]["N_median"]
print("ALL E5 TESTS PASSED")
