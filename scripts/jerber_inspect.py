"""A1: inspect Jerber et al. 2021 day-52 h5 (public, CC-BY) - structure and counts only.

Step 1 (--tree): list groups/datasets with shapes and dtypes (no values).
Step 2 (--obs): summarise categorical cell metadata: per column, number of levels; for columns that look
like donor / pool / treatment, the counts needed for the go/no-go rule:
  >= 10 lines with cells in >= 2 independent pools, untreated condition.
"""
import argparse, json, sys
import h5py
import numpy as np


def tree(f):
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"D {name} {obj.shape} {obj.dtype}")
        else:
            print(f"G {name} attrs={list(obj.attrs.keys())[:6]}")
    f.visititems(visit)


def _decode(a):
    a = np.asarray(a)
    if a.dtype.kind in "SO":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return a


def read_col(g, name):
    """AnnData-style column: categorical group (codes+categories) or plain dataset, or older
    __categories layout."""
    obj = g[name]
    if isinstance(obj, h5py.Group) and "codes" in obj:
        cats = _decode(obj["categories"][()]); codes = obj["codes"][()]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    v = _decode(obj[()])
    if "__categories" in g and name in g["__categories"]:
        cats = _decode(g["__categories"][name][()])
        return np.where(v.astype(int) >= 0, cats[np.clip(v.astype(int), 0, None)], "NA")
    return v


def obs_summary(f, obs_path):
    g = f[obs_path]
    cols = [k for k in g.keys() if not k.startswith("__") and k not in ("_index", "index")]
    out = {}
    for c in cols:
        try:
            v = read_col(g, c)
        except Exception as e:
            out[c] = {"error": type(e).__name__}; continue
        if v.ndim != 1:
            continue
        lv = np.unique(v)
        out[c] = {"n_levels": int(len(lv)), "dtype": str(v.dtype)}
        if len(lv) <= 12:            # small categorical: show levels (study design labels, not people)
            out[c]["levels"] = {str(l): int((v == l).sum()) for l in lv}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("h5"); ap.add_argument("--tree", action="store_true")
    ap.add_argument("--obs"); ap.add_argument("--donor"); ap.add_argument("--pool"); ap.add_argument("--treat")
    ap.add_argument("--untreated"); ap.add_argument("--celltype"); ap.add_argument("--out")
    a = ap.parse_args()
    with h5py.File(a.h5, "r") as f:
        if a.tree:
            tree(f); sys.exit()
        if a.obs and not a.donor:
            print(json.dumps(obs_summary(f, a.obs), indent=1)); sys.exit()
        g = f[a.obs]
        d, p = read_col(g, a.donor), read_col(g, a.pool)
        keep = np.ones(len(d), bool)
        if a.treat:
            keep &= read_col(g, a.treat) == a.untreated
        res = {"n_cells": int(len(d)), "n_cells_untreated": int(keep.sum())}
        pairs = {}
        for di, pi in zip(d[keep], p[keep]):
            pairs.setdefault(di, {}).setdefault(pi, 0); pairs[di][pi] += 1
        res["n_lines_untreated"] = len(pairs)
        res["n_pools_untreated"] = int(len(np.unique(p[keep])))
        for thr in (1, 50, 100):
            res[f"lines_in_>=2_pools_min{thr}cells"] = int(sum(sum(n >= thr for n in v.values()) >= 2 for v in pairs.values()))
        res["pools_per_line_hist"] = {int(k): int(n) for k, n in zip(*np.unique([len(v) for v in pairs.values()], return_counts=True))}
        if a.celltype:
            ct = read_col(g, a.celltype)[keep]
            res["celltypes_untreated"] = {str(k): int(n) for k, n in zip(*np.unique(ct, return_counts=True))}
        print(json.dumps(res, indent=1))
        if a.out:
            json.dump(res, open(a.out, "w"), indent=1)
