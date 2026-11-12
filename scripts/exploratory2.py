"""EXPLORATORY E5: how uncertain is the design guide (E4), and what do more replicates say?
Specified and frozen (git tag explore-v0.2.0) before running on real data. Aggregates only.

E5a  Uncertainty in E4. For each of B = 1000 draws:
       - phenotype ICC from its sampling distribution: (1 - ICC) * R_obs ~ F(n_d - 1, n_w - 1), where
         R_obs = S_d^2 / S_w^2 (the pivot behind the A0-2 CI), so ICC* = 1 - F* / R_obs;
       - outcome-slope reliability lambda* from a donor bootstrap (donors resampled with
         replacement, LMM refitted; B_lambda = 200 per outcome, recycled across the 1000 draws).
     Reported: N donors for 80% power at r_true = 0.3 and 0.4, k = 1 and 4 runs per line, as the
     median and 10th-90th percentile, and P(ICC* <= 0), where the phenotype carries no line signal.
E5b  Replicate reliability in bulk RNA-seq (FOUNDIN RNAB, day 65; an independent assay of separate
     cultures). Phenotypes on log2(CPM + 1) of the bulk count table, genes with mean CPM >= 10 over
     the donor lines: B-SNCA; B-OXPHOS and B-LYSO (mean gene-wise z over the same MSigDB sets);
     B-TH (DA-neuron marker, a proxy for P1). Not residualised, since bulk has no composition.
     Bulk IDs are versioned Ensembl IDs, mapped to symbols with the Cell Ranger feature list
     (version stripped; IDs sharing a symbol summed; unmapped dropped).
     ICC as in A0-2: donor lines (one run each) vs the control line's 5 day-65 runs (4 df), with the
     F-based 95% CI. Two donor sets: the 79 scRNA donor lines, and all bulk day-65 donor lines.
E5c  Where does run-to-run noise come from? Day 25, control line: 9 runs within batch 1 and 4 runs in
     batches 2, 3, 5E6 and 5E8. For each E5b phenotype (genes filtered and z-scored over day-25
     donor lines):
       within-batch variance  = variance over the 9 batch-1 runs (8 df);
       all-run variance       = variance over all 13 control runs;
       between-batch share    = 1 - within / all-run (clipped at 0);
       ICC_within_batch       = 1 - within / S_d^2 (day-25 donor lines).
     Answers whether replicate differentiations need to span batches to average out noise.
Output: data/derived/exploratory2.json
"""
from __future__ import annotations

import argparse, json, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exploratory import slope_reliability, n_needed, K_RUNS  # noqa: E402
from gates_prelink import icc  # noqa: E402
from models import donor_table, _clean  # noqa: E402
from phenotypes import load_config, read_geneset  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261113
B = 1000
B_LAMBDA = 200
RNAB_RE = re.compile(r"^RNAB_PPMI(\d+)(B\w+?)?_(\d+)_da(\d+)_v(\d+)$")


def boot_lambda(d, visits, min_visits, rng, n_boot):
    out = []
    ids = np.array(d.index)
    for _ in range(n_boot):
        pick = rng.choice(ids, len(ids), replace=True)
        newd, newv = [], []
        for j, p in enumerate(pick):
            lab = f"b{j}"
            newd.append(d.loc[[p]].rename(index={p: lab}))
            vv = visits[visits.PATNO == p].copy(); vv["PATNO"] = lab; newv.append(vv)
        dd = pd.concat(newd); vv = pd.concat(newv)
        try:
            out.append(slope_reliability(dd, vv, min_visits=min_visits)["lambda_mean"])
        except Exception:
            out.append(np.nan)
    return np.array(out)


def e5a(gates, d, clin, rng):
    off = clin["updrs3"][clin["updrs3"].state == "OFF"]
    lam = {"UPDRS3_OFF": boot_lambda(d, off, 3, rng, B_LAMBDA),
           "putamen_SBR": boot_lambda(d, clin["datscan"], 2, rng, B_LAMBDA)}
    res = {"lambda_boot": {k: {"median": float(np.nanmedian(v)), "p10": float(np.nanpercentile(v, 10)),
                               "p90": float(np.nanpercentile(v, 90)), "n_failed": int(np.isnan(v).sum())}
                           for k, v in lam.items()}}
    zz = stats.norm.ppf(0.975) + stats.norm.ppf(0.8)
    for ph in ["P1", "P2", "P3", "P4"]:
        g = gates["A0_2"][ph]
        n_d, n_w = g["n_donor"], g["n_within"]
        R = 1 / (1 - g["icc"])
        F = stats.f.rvs(n_d - 1, n_w - 1, size=B, random_state=rng)
        icc_s = 1 - F / R
        res[ph] = {"icc_point": g["icc"], "icc_p10": float(np.percentile(icc_s, 10)),
                   "icc_p90": float(np.percentile(icc_s, 90)), "P_icc_le_0": float((icc_s <= 0).mean())}
        for o, lv in lam.items():
            lv = lv[np.isfinite(lv)]
            ld = rng.choice(lv, B, replace=True)
            for k in (1, 4):
                rk = np.where(icc_s > 0, k * icc_s / (1 + (k - 1) * np.clip(icc_s, 0, None)), 0.0)
                for r in (0.3, 0.4):
                    ro = r * np.sqrt(rk * ld)
                    N = np.where(ro > 0, (zz / np.arctanh(np.clip(ro, 1e-9, 0.999))) ** 2 + 9, np.inf)
                    res[ph][f"{o}|k{k}|r{r}"] = {"N_median": float(np.median(N)), "N_p10": float(np.percentile(N, 10)),
                                                 "N_p90": float(np.percentile(N, 90)) if np.isfinite(np.percentile(N, 90)) else None}
    return res


def bulk_scores(counts, cols_fit, cfg, extra_genes=("SNCA", "TH")):
    cpm = counts.div(counts.sum(0), axis=1) * 1e6
    keep = cpm[cols_fit].mean(1) >= cfg["min_mean_cpm"]
    L = np.log2(cpm[keep] + 1)
    mu, sd = L[cols_fit].mean(1), L[cols_fit].std(1, ddof=0)
    Z = L.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)
    out = pd.DataFrame(index=counts.columns)
    for g in extra_genes:
        out[f"B-{g}"] = L.loc[g] if g in L.index else np.nan
    for name, k in [("B-OXPHOS", "P3"), ("B-LYSO", "P4")]:
        gs = [x for x in read_geneset(cfg["gene_sets"][k]) if x in Z.index]
        out[name] = Z.loc[gs].mean(0)
    return out, int(keep.sum())


def load_bulk(path):
    raw = pd.read_csv(path, sep="\t", index_col=0)
    meta = []
    for c in raw.columns:
        m = RNAB_RE.match(c)
        meta.append({"col": c, "line": m.group(1), "ctrl": m.group(2) is not None, "batch": m.group(2),
                     "day": int(m.group(4))})
    return raw, pd.DataFrame(meta)


def to_symbols(raw, h5_path):
    """Map versioned Ensembl IDs to symbols using the Cell Ranger feature list (id -> name);
    unmapped IDs are dropped; counts of IDs sharing a symbol are summed."""
    import h5py
    with h5py.File(h5_path, "r") as f:
        ft = f["matrix/features"]
        ids = [x.decode() for x in ft["id"][:]]; names = [x.decode() for x in ft["name"][:]]
    mp = dict(zip([re.sub(r"\.\d+$", "", x) for x in ids], names))  # v0.2.1: strip versions on both sides
    sym = raw.index.astype(str).str.replace(r"\.\d+$", "", regex=True).map(mp)
    out = raw[sym.notna()].copy()
    out.index = sym[sym.notna()]
    return out.groupby(level=0).sum(), {"n_ids": int(len(raw)), "n_mapped": int(sym.notna().sum()),
                                         "n_symbols": int(out.index.nunique())}


def e5b_e5c(bulk_path, scrn_lines, cfg, h5_path):
    raw, meta = load_bulk(bulk_path)
    raw, mapinfo = to_symbols(raw, h5_path)
    res = {"bulk_id_mapping": mapinfo}
    scrn_pat = set(scrn_lines)
    for day in (65, 25):
        m = meta[meta.day == day]
        don = m[~m.ctrl]
        counts = raw[m.col]
        sc, nkeep = bulk_scores(counts, list(don.col), cfg)
        ctrl = m[m.ctrl]
        out = {"n_genes_kept": nkeep, "n_donor_lines_all": int(len(don)), "n_control_runs": int(len(ctrl))}
        if day == 65:
            for label, dset in [("scrn_79", don[don.line.isin(scrn_pat)]), ("all", don)]:
                out[label] = {ph: icc(sc.loc[dset.col, ph], sc.loc[ctrl.col, ph]) for ph in sc.columns}
            res["E5b_day65"] = out
        else:
            b1 = ctrl[ctrl.batch == "B1"]
            for ph in sc.columns:
                w = sc.loc[b1.col, ph].var(ddof=1); a = sc.loc[ctrl.col, ph].var(ddof=1)
                sd2 = sc.loc[don.col, ph].var(ddof=1)
                out[ph] = {"within_batch_var": float(w), "all_run_var": float(a),
                           "between_batch_share": float(max(0.0, 1 - w / a)),
                           "icc_within_batch": float(1 - w / sd2), "icc_all_runs": float(1 - a / sd2),
                           "n_batch1_runs": int(len(b1))}
            res["E5c_day25"] = out
    return res


def main(argv=None):
    from ppmi_load import build
    from scrn_gate import read_table
    ap = argparse.ArgumentParser()
    ap.add_argument("--phen-dir", default=str(ROOT / "data" / "derived" / "phenotypes"))
    ap.add_argument("--s1", required=True)
    ap.add_argument("--ppmi-dir", required=True)
    ap.add_argument("--bulk", default=str(ROOT / "data" / "unpacked_rnab" / "FOUNDIN-PD_150.9_RNAB" / "aggregated_expression" / "countTable.tsv"))
    ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "exploratory2.json"))
    a = ap.parse_args(argv)
    rng = np.random.default_rng(SEED)
    cfg = load_config()
    pdir = Path(a.phen_dir)
    lines = pd.read_csv(pdir / "lines.csv"); gates = json.load(open(pdir / "gates_prelink.json"))
    s1 = read_table(Path(a.s1))
    d0 = donor_table(lines, s1, pd.DataFrame(columns=["PATNO", "sex_female", "age_onset", "dur_baseline", "saa_pos", "nsd_iss"]))
    clin = build(a.ppmi_dir, list(d0.index))
    d = donor_table(lines, s1, clin["donors"])
    d = d[d.n_cells_qc >= gates["A0_3"]["min_qc_cells"]]
    res = {"label": "EXPLORATORY", "version": "0.2.0", "seed": SEED}
    scrn_lines = lines.loc[~lines.is_control, "line"].str.extract(r"(\d+)")[0]
    h5 = next((ROOT / "data" / "unpacked").rglob("filtered_feature_bc_matrix.h5"))
    res.update(e5b_e5c(a.bulk, scrn_lines, cfg, h5))
    print("E5b/E5c done", flush=True)
    res["E5a"] = e5a(gates, d, clin, rng)
    with open(a.out, "w") as f:
        json.dump(_clean(res), f, indent=1)
    print("written", a.out)


if __name__ == "__main__":
    main()
