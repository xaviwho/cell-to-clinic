"""Synthetic check of the A3 variance decomposition (no real data)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts")); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import jerber_a3 as a3
from test_jerber_reliability import _design


def test_decompose_recovers_components():
    rng = np.random.default_rng(3)
    out = []
    for _ in range(150):
        lines, pools = _design(rng, n_lines=200, n_rep=40)
        n = rng.integers(50, 800, len(lines)).astype(float)
        c = 20.0                                                   # per-cell noise constant
        u = {l: rng.normal(0, 1.0) for l in np.unique(lines)}      # line var 1
        b = {p: rng.normal(0, 0.5) for p in np.unique(pools)}      # batch var 0.25
        y = np.array([u[l] + b[p] for l, p in zip(lines, pools)]) + rng.normal(0, np.sqrt(c / n))
        out.append(a3.decompose(c / n, n, y, lines, pools))
    assert abs(np.median([o["sigma_line2"] for o in out]) - 1.0) < 0.15
    assert abs(np.median([o["sigma_batch2"] for o in out]) - 0.25) < 0.1
    assert abs(np.median([o["icc_by_cells"]["inf"] for o in out]) - 0.8) < 0.07


if __name__ == "__main__":
    test_decompose_recovers_components(); print("ALL PASSED")
