"""Synthetic check for scripts/exploratory3.py. Run: python -I tests/test_exploratory3.py <empty_tmp_dir>"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import h5py
import numpy as np
import pandas as pd
import exploratory3 as ex

TMP = pathlib.Path(sys.argv[1]); assert TMP.exists() and not any(TMP.iterdir())
rng = np.random.default_rng(11)
genes = ["SNCA", "TH", "STMN2", "SNAP25", "SYT1"] + [f"F{i}" for i in range(300)]
ids = [f"ENSG{i:011d}.{i % 5 + 1}" for i in range(len(genes))]
with h5py.File(TMP / "f.h5", "w") as f:
    g = f.create_group("matrix/features")
    g["id"] = np.array([x.encode() for x in ids]); g["name"] = np.array([x.encode() for x in genes])
base = rng.normal(np.log(40), 1, len(genes)); base[:5] = np.log(400)
neuro_load = np.zeros(len(genes)); neuro_load[:5] = 1.0          # SNCA and markers scale with neuron fraction
cols, data = [], []
def sample(name, line_eff, neuro):
    mu = np.exp(base + neuro_load * neuro + (np.arange(len(genes)) == 0) * line_eff)
    cols.append(name); data.append(rng.poisson(mu * 30))
line_eff = rng.normal(0, 0.3, 60)
for day in (65, 25):
    for i in range(60):
        sample(f"RNAB_PPMI{3000+i}_{1000+i}_da{day}_v1", line_eff[i], rng.normal(0, 0.6))
    batch_neuro = {b: rng.normal(0, 0.6) for b in ["B1", "B2", "B3", "B5E6", "B5E8"]}
    reps = ["B1"] * (9 if day == 25 else 1) + ["B2", "B3", "B5E6", "B5E8"]
    for j, b in enumerate(reps):
        sample(f"RNAB_PPMI9999{b}_0001_da{day}_v{j+1}", 0.1, batch_neuro[b] + rng.normal(0, 0.02))
pd.DataFrame(np.array(data).T, index=pd.Index(ids, name="Geneid"), columns=cols).to_csv(TMP / "c.tsv", sep="\t")
r = ex.run(TMP / "c.tsv", TMP / "f.h5", [str(3000 + i) for i in range(50)])
d65 = r["day65"]["all"]
print("day65 ICC unadjusted", round(d65["SNCA"]["icc"], 2), "adjusted", round(d65["SNCA_adj"]["icc"], 2), "fit", r["day65"]["fit"])
print("day25 batch share unadj", r["day25"]["SNCA"]["share"], "adj", r["day25"]["SNCA_adj"]["share"], r["day25"]["SNCA_adj"]["ci95"])
assert d65["SNCA_adj"]["icc"] > d65["SNCA"]["icc"] and d65["SNCA_adj"]["icc"] > 0.9
assert r["day25"]["SNCA"]["share"] > 0.8 and r["day25"]["SNCA_adj"]["share"] < r["day25"]["SNCA"]["share"]
assert r["day65"]["scrn_79"]["SNCA"]["n_donor"] == 50 and r["day25"]["n_batch1_runs"] == 9
print("ALL E6 TESTS PASSED")
