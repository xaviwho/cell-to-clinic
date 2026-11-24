"""Synthetic checks for scripts/phenotypes_da.py (no real data)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd
import phenotypes_da as pdz


def test_scores_exclusion_and_z():
    rng = np.random.default_rng(0)
    genes = ["SNCA"] + [f"G{i}" for i in range(40)]
    C = pd.DataFrame(rng.poisson(200, (12, len(genes))), columns=genes, index=[f"r{i}" for i in range(12)]).astype(float)
    isd = pd.Series([True] * 10 + [False] * 2, C.index)
    ok = pd.Series([True] * 9 + [False] + [True] * 2, C.index)
    cfg = {"min_mean_cpm": 10, "P2_gene": "SNCA"}
    gs = {"P3": [f"G{i}" for i in range(20)], "P4": [f"G{i}" for i in range(20, 40)] + ["MISSING"]}
    out, Z, n, kept = pdz.scores(C, isd, ok, gs, cfg)
    assert out.loc["r9"].isna().all()                      # excluded run
    assert n == {"D3": 20, "D4": 20} and kept == len(genes)
    assert np.allclose(Z.loc[isd & ok].mean(0), 0, atol=1e-9)
    L = np.log2(C.div(C.sum(1), axis=0) * 1e6 + 1)
    assert np.isclose(out.loc["r0", "D2"], L.loc["r0", "SNCA"])


if __name__ == "__main__":
    test_scores_exclusion_and_z(); print("ALL PASSED")
