"""Synthetic checks for scripts/explore_batch.py (no real data). Run: python -I tests/test_explore_batch.py"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd
import explore_batch as eb


def synth(confound, seed=0, n=60):
    rng = np.random.default_rng(seed)
    batch = rng.integers(1, 6, n)
    shift = np.array([0, -1.5, -0.5, 0.5, 1.5])[batch - 1] if confound else np.zeros(n)
    d2 = rng.normal(size=n) + shift
    slope = 2 + (-1.0 * shift if confound else -0.8 * d2) + rng.normal(0, 0.5, n)
    ids = [str(i) for i in range(n)]
    d = pd.DataFrame({"D2": d2, "batch": batch.astype(float), "media_e6": np.isin(batch, [4, 5]).astype(float)}, index=ids)
    rows = [{"PATNO": i, "t_dx": t, "total": 20 + s * t + rng.normal(0, 1)} for i, s in zip(ids, slope) for t in range(5)]
    return d, pd.DataFrame(rows)


def test_confounded():
    r = eb.analyse(*synth(True))
    assert r["reading"] == "confounding plausible", r
    assert r["a_D2_by_batch"]["p"] < 0.01 and r["b_slope_by_batch"]["p"] < 0.01


def test_within_batch():
    r = eb.analyse(*synth(False))
    assert r["reading"] == "association mainly within batches", r
    c = r["c_association"]
    assert c["within_batch_r"] < -0.5 and c["within_df"] == c["n"] - 6
    assert sum(r["donors_per_batch"].values()) == 60


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f(); print("ok", k)
    print("ALL PASSED")
