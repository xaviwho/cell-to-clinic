"""Simulation check of the registered power statements (OSF osf.io/s6w9c, Part 5). Synthetic data only.

Uses the registered H1 analysis exactly: REML LMM with random intercept + slope (registered
fallbacks), Satterthwaite t test on P:t, BH over 4 phenotypes, TOST on partial r with bound 0.40.

Data-generating assumptions (set before any clinical data is seen; the real visit structure is
unknown until linkage):
  - donors n in {47, 30}; each has 3-9 annual OFF-state UPDRS III visits starting 0-3 years after
    diagnosis (uniform), visit times jittered by +/- 0.1 y;
  - intercept 20 (SD 6), mean slope 2.0 points/y (SD 1.5), residual SD 5;
  - 6 covariates (3 continuous, 3 binary), each with main effect 0.5 and slope effect 0.1;
  - 4 independent standard-normal phenotypes; phenotype 1 has P:t effect beta3 (grid), the others 0.
Reported per cell: mean partial r of phenotype 1, power at p < 0.05 (single test), power at BH
q < 0.05 among the 4 tests, rate of "equivalent" (90% CI of partial r inside +/-0.40), and how
often a fallback random-effects structure was used.

Run: python -I scripts/power_sim.py --out results/power_sim.json [--nsim 300] [--workers 14]
"""
from __future__ import annotations

import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lmm import lmm_fit_fallback, satterthwaite, tost_partial_r, bh  # noqa: E402

SEED = 20261112
BETAS = [0.0, 0.4, 0.6, 0.8, 1.0, 1.2]
NS = [47, 30]


def one_rep(args):
    n, beta3, seed = args
    rng = np.random.default_rng(seed)
    P = rng.normal(size=(n, 4))
    C = np.c_[rng.normal(size=(n, 3)), rng.integers(0, 2, (n, 3))].astype(float)
    a0 = 20 + rng.normal(0, 6, n); a1 = 2.0 + rng.normal(0, 1.5, n)
    g, t, y = [], [], []
    for i in range(n):
        m = int(rng.integers(3, 10)); start = rng.uniform(0, 3)
        tt = start + np.arange(m) + rng.uniform(-0.1, 0.1, m)
        mu = a0[i] + a1[i] * tt + beta3 * P[i, 0] * tt + 0.5 * C[i].sum() + 0.1 * C[i].sum() * tt
        g += [i] * m; t += list(tt); y += list(mu + rng.normal(0, 5, m))
    g, t, y = np.array(g), np.array(t), np.array(y)
    Cz = C.copy(); Cz[:, :3] = (Cz[:, :3] - Cz[:, :3].mean(0)) / Cz[:, :3].std(0, ddof=1)
    out = []
    for k in range(4):
        pk = (P[:, k] - P[:, k].mean()) / P[:, k].std(ddof=1)
        X = np.column_stack([np.ones(len(y)), t, pk[g], pk[g] * t, Cz[g], Cz[g] * t[:, None]])
        f = lmm_fit_fallback(y, X, g, t)
        L = np.zeros(X.shape[1]); L[3] = 1
        s = satterthwaite(f, L)
        e = tost_partial_r(s["t"], s["df"])
        out.append((s["p"], e["r"], e["equivalent"], f["re"] != "slope_corr", s["df"]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--nsim", type=int, default=300)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    a = ap.parse_args(argv)
    ss = np.random.SeedSequence(SEED)
    cells = [(n, b) for n in NS for b in BETAS]
    seeds = ss.spawn(len(cells))
    t0 = time.time()
    res = {"seed": SEED, "nsim": a.nsim, "assumptions": __doc__.split("Data-generating assumptions")[1].split("Reported")[0].strip(),
           "cells": []}
    with ProcessPoolExecutor(a.workers) as ex:
        for (n, b), s in zip(cells, seeds):
            reps = list(ex.map(one_rep, [(n, b, int(x.generate_state(1)[0])) for x in s.spawn(a.nsim)], chunksize=4))
            p1 = np.array([r[0][0] for r in reps]); r1 = np.array([r[0][1] for r in reps])
            q1 = np.array([bh([x[0] for x in r])[0] for r in reps])
            eq1 = np.array([r[0][2] for r in reps]); fb = np.mean([any(x[3] for x in r) for r in reps])
            df1 = np.array([r[0][4] for r in reps])
            null_p = np.array([x[0] for r in reps for x in r[1:]])
            cell = {"n": n, "beta3": b, "mean_partial_r": float(r1.mean()), "power_p05": float((p1 < 0.05).mean()),
                    "power_bh4": float((q1 < 0.05).mean()), "equivalent_rate": float(eq1.mean()),
                    "type1_null_phenotypes": float((null_p < 0.05).mean()), "fallback_rate": float(fb),
                    "median_df": float(np.median(df1))}
            res["cells"].append(cell)
            print(f"n={n} beta3={b:.1f}: r={cell['mean_partial_r']:.2f} power={cell['power_p05']:.2f} "
                  f"BH={cell['power_bh4']:.2f} equiv={cell['equivalent_rate']:.2f} typeI={cell['type1_null_phenotypes']:.3f} "
                  f"fallback={fb:.2f} df={cell['median_df']:.1f} [{time.time()-t0:.0f}s]", flush=True)
    for n in NS:
        cs = [c for c in res["cells"] if c["n"] == n]
        r = np.array([c["mean_partial_r"] for c in cs])
        for key in ["power_p05", "power_bh4"]:
            pw = np.maximum.accumulate(np.array([c[key] for c in cs]))  # monotone for interpolation
            res[f"r_at_80pct_{key}_n{n}"] = float(np.interp(0.8, pw, r)) if pw.max() >= 0.8 else None
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(res, f, indent=1)
    print({k: v for k, v in res.items() if k.startswith("r_at")})


if __name__ == "__main__":
    main()
