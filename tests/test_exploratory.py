"""Synthetic checks for scripts/exploratory.py (Cox with delayed entry, milestone table, design maths)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd
from exploratory import cox_fit, milestone_table, n_needed, design

rng = np.random.default_rng(7)
# left-truncated exponential survival with known log-HRs (0.5, -0.3)
n = 4000
X = np.c_[rng.normal(size=n), rng.integers(0, 2, n)]
T = rng.exponential(1 / (0.1 * np.exp(X @ [0.5, -0.3])))
E = rng.uniform(0, 5, n)
C = E + rng.uniform(2, 12, n)
keep = T > E                      # only those event-free at entry are observed
time = np.minimum(T, C)[keep]; ev = (T <= C)[keep].astype(int)
f = cox_fit(E[keep], time, ev, X[keep])
print("Cox delayed entry: beta", f["beta"].round(3), "se", f["se"].round(3), "(truth 0.5, -0.3)")
assert np.all(np.abs(f["beta"] - [0.5, -0.3]) < 3 * f["se"])

a = pd.DataFrame({"PATNO": ["a"] * 3 + ["b"] * 3 + ["c"] * 2 + ["d"],
                  "t_dx": [1, 2, 3, 1, 2, 3, 0.5, 4, 2], "v": [1, 2, 3, 3, 3, 3, 0, 1, 0]})
tab, prev = milestone_table(a, "v", 3)
assert prev == 1 and list(tab.PATNO) == ["a", "c"] and tab.set_index("PATNO").loc["a"].tolist() == [1, 3, 1]
assert tab.set_index("PATNO").loc["c"].tolist() == [0.5, 4, 0]   # censored at last; 'd' has 1 visit -> dropped

assert abs(n_needed(0.3) - (((1.959964 + 0.841621) / np.arctanh(0.3)) ** 2 + 9)) < 1e-3
g = design({"P": 1.0}, {"O": 1.0})["P|O"]
assert all(abs(c["r_observable"] - c["r_true"]) < 1e-12 for c in g)
g = design({"P": 0.5}, {"O": 0.5})["P|O"]
c = [x for x in g if x["k_runs"] == 1 and x["r_true"] == 0.4][0]
assert abs(c["r_observable"] - 0.2) < 1e-12
c4 = [x for x in g if x["k_runs"] == 4 and x["r_true"] == 0.4][0]
assert abs(c4["reliability_product"] - 0.8 * 0.5) < 1e-12
print("ALL EXPLORATORY TESTS PASSED")
