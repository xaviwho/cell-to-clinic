"""Synthetic checks for scripts/jerber_reliability.py (no real data)."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import numpy as np
import jerber_reliability as jr
from cell2clinic.gates import stability_icc


def _design(rng, n_lines=150, n_rep=30, n_pools=12):
    lines, pools = [], []
    for i in range(n_lines):
        ps = rng.choice(n_pools, 2 if i < n_rep else 1, replace=False)
        for p in ps:
            lines.append(f"L{i:03d}"); pools.append(f"P{p:02d}")
    return np.array(lines), np.array(pools)


def test_vectorised_icc_matches_scalar():
    rng = np.random.default_rng(0)
    lines, pools = _design(rng)
    Y = rng.normal(size=(len(lines), 5))
    v = jr.icc_genes(Y, lines, pools)[0]
    for j in range(5):
        dv, reps = jr.split(lines, pools, Y[:, j])
        assert np.isclose(v[j], stability_icc(dv, reps)["icc"])


def test_icc_recovers_truth_with_pool_effect():
    rng = np.random.default_rng(1)
    est = []
    for _ in range(200):
        lines, pools = _design(rng)
        u = {l: rng.normal() for l in np.unique(lines)}           # line effect, var 1
        b = {p: rng.normal(0, 0.5) for p in np.unique(pools)}     # pool (batch) effect, var 0.25
        y = np.array([u[l] + b[p] + rng.normal(0, 0.5) for l, p in zip(lines, pools)])
        dv, reps = jr.split(lines, pools, y)
        est.append(stability_icc(dv, reps)["icc"])
    assert abs(np.median(est) - 1 / 1.5) < 0.05     # across-batch truth: 1 / (1 + 0.25 + 0.25)


def test_resid_matrix():
    rng = np.random.default_rng(2)
    Z = rng.normal(size=(50, 2)); Y = Z @ [[1, 2], [3, 4]] + 0.1 * rng.normal(size=(50, 2))
    R = jr.resid(Y, Z)
    assert np.allclose(np.c_[np.ones(50), Z].T @ R, 0, atol=1e-8)


if __name__ == "__main__":
    for n, fn in list(globals().items()):
        if n.startswith("test_"):
            fn(); print("ok", n)
    print("ALL PASSED")
