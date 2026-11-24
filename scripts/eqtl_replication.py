"""C2 / gate G-eQTL (Cell-to-Clinic II draft, section 4): do Jerber DA-neuron lead eQTLs replicate in
FOUNDIN DA-neuron pseudobulk? NO clinical data. Freeze with a tag before running on FOUNDIN data.

Step `leads` (public data only). DEVIATION from osf.io/4pzxa, decided by the user 2026-11-24 before any
  FOUNDIN data was used: the released qtl_results_all files carry empirical_feature_p_value = -1 on every
  row, so the gene-level p cannot be taken from them. Leads therefore come from Jerber et al. 2021
  Supplementary Table 7 (publisher file 41588_2021_801_MOESM3_ESM.xlsx: "all lead eQTL SNP-gene pairs"),
  context label DA_D52 (day 52, untreated), and the registered rule is applied to its columns:
  BH_fdr (BH over the gene-level empirical p) < 0.05. Cross-check: every lead is looked up in the D52
  untreated DA qtl_results_all file (same SNP-gene pair; beta and p agreement reported).
  The original all-SNP route (`--source allsnp`) is kept for the record; it cannot produce leads.
  -> data/derived/eqtl/jerber_leads.tsv, jerber_leads.json
Step `replicate` (FOUNDIN, aggregates only):
  expression  log2(CPM+1) of the new-rule DA pseudobulk (da-labels-v1.0.0), 79 donor lines;
              genes with mean CPM >= 10; symbols mapped to Ensembl via the Cell Ranger features.
  genotype    data/derived/geno dosage (count of a1), matched on chr:pos (GRCh37) and alleles;
              strand-ambiguous A/T and C/G SNPs dropped; aligned to Jerber's assessed allele.
  model       OLS: expression ~ dosage + 3 genotype PCs + 5 expression PCs (top-2,000-variance genes).
  statistics  sign concordance with Jerber beta; pi1 = 1 - pi0 (Storey, lambda = 0.5).
  null        200 sets: each lead gene paired with a random FOUNDIN SNP on another chromosome in the
              same MAF bin (0.05-0.1-0.2-0.3-0.5), compared with the lead's Jerber sign.
  PASS        concordance >= 0.70 AND pi1 >= 0.30 AND both above the null's 95th percentile.
  sensitivity the same with the old-label DA pseudobulk and the all-cell pseudobulk (report only).
Output: data/derived/eqtl/g_eqtl.json (aggregates only).
"""
from __future__ import annotations

import argparse, gzip, json, re, subprocess, sys, tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "derived" / "eqtl"
MAF_BINS = [0.05, 0.1, 0.2, 0.3, 0.5]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
N_NULL, SEED = 200, 0


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def pi1(p, lam=0.5):
    p = np.asarray(p, float); p = p[np.isfinite(p)]
    return float(1 - min(1.0, np.mean(p > lam) / (1 - lam))) if len(p) else np.nan


def align(jref, jalt, jassessed, a1, a2):
    """+1 if FOUNDIN dosage (count of a1) counts Jerber's assessed allele, -1 if it counts the other
    allele, 0 if the alleles do not match or the SNP is strand-ambiguous."""
    if {jref, jalt} in ({"A", "T"}, {"C", "G"}):
        return 0
    for flip in (False, True):
        b1, b2 = (COMP.get(a1, "?"), COMP.get(a2, "?")) if flip else (a1, a2)
        if {b1, b2} == {jref, jalt}:
            return 1 if b1 == jassessed else -1
    return 0


def ols_t(y, x, C):
    X = np.c_[np.ones(len(y)), x, C]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; df = len(y) - X.shape[1]
    s2 = r @ r / df
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    t = b[1] / se
    return b[1], 2 * stats.t.sf(abs(t), df)


# archive members are <Day>.<celltype>.<untreated|ROT_treated>.qtl_results_all...; the registration uses untreated
MEMBER = r"/D52\.DA\.untreated\.qtl_results_all\.sorted\.txt\.gz$"

COLS = ["feature_id", "snp_id", "p_value", "beta", "empirical_feature_p_value",
        "snp_chromosome", "snp_position", "assessed_allele", "maf"]


def lead_table(chunks):
    """per gene: the min-p SNP row, with the gene's empirical p (the value that is not -1)."""
    tops, emps = [], []
    for ch in chunks:
        tops.append(ch.sort_values("p_value").drop_duplicates("feature_id"))
        emps.append(ch.groupby("feature_id").empirical_feature_p_value.max())
    T = pd.concat(tops).sort_values("p_value").drop_duplicates("feature_id").set_index("feature_id")
    T["empirical_feature_p_value"] = pd.concat(emps).groupby(level=0).max()
    T = T[T.empirical_feature_p_value >= 0].reset_index()
    T["q"] = bh(T.empirical_feature_p_value)
    return T


def leads(archive, pattern):
    """Original registered route (all-SNP file). Kept for the record: the file has no gene-level p."""
    OUT.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r|gz") as t:
        for m in t:
            if m.isfile() and re.search(pattern, m.name):
                src = gzip.GzipFile(fileobj=t.extractfile(m))
                L = lead_table(pd.read_csv(src, sep="\t", chunksize=2_000_000, usecols=COLS))
                info = {"member": m.name}
                break
        else:
            raise SystemExit(f"no archive member matches {pattern}")
    if L.empty:
        raise SystemExit("no gene-level empirical p in this file (all -1): use --source table7")
    lead = L[L.q < 0.05].copy()
    al = lead.snp_id.str.split("_", expand=True)
    lead["ref"], lead["alt"] = al[2], al[3]
    lead.to_csv(OUT / "jerber_leads.tsv", sep="\t", index=False)
    info.update({"n_genes_tested": int(len(L)), "n_lead_fdr05": int(len(lead))})
    json.dump(info, open(OUT / "jerber_leads.json", "w"), indent=1)
    print(info)


T7_SHEET, T7_LABEL = "SupplementaryTable7", "DA_D52"


def table7_leads(xlsx, label=T7_LABEL):
    T = pd.read_excel(xlsx, T7_SHEET, header=3)
    ctx = T[T.label == label]
    if ctx.ensembl_gene_id.duplicated().any():
        raise SystemExit("Table 7: more than one lead per gene in this context")
    lead = ctx[ctx.BH_fdr < 0.05].rename(columns={"ensembl_gene_id": "feature_id"}).copy()
    al = lead.snp_id.str.split("_", expand=True)
    lead["ref"], lead["alt"] = al[2], al[3]
    info = {"source": "Jerber 2021 Supplementary Table 7", "label": label, "n_in_context": int(len(ctx)),
            "n_q_lt_05": int((ctx.q_value < 0.05).sum()), "n_lead_bh_fdr05": int(len(lead))}
    return lead[["feature_id", "snp_id", "p_value", "beta", "empirical_feature_p_value", "BH_fdr", "q_value",
                 "snp_chromosome", "snp_position", "assessed_allele", "maf", "ref", "alt", "hgnc_symbol"]], info


def crosscheck(lead, archive, pattern):
    """Look the Table 7 pairs up in the all-SNP file of the same context."""
    want = set(zip(lead.feature_id, lead.snp_id))
    hits = []
    with tarfile.open(archive, "r|gz") as t:
        for m in t:
            if m.isfile() and re.search(pattern, m.name):
                for ch in pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(m)), sep="\t", chunksize=2_000_000,
                                      usecols=["feature_id", "snp_id", "p_value", "beta"]):
                    k = pd.Series(list(zip(ch.feature_id, ch.snp_id)), index=ch.index).isin(want)
                    hits.append(ch[k])
                break
    H = pd.concat(hits).merge(lead[["feature_id", "snp_id", "beta", "p_value"]], on=["feature_id", "snp_id"], suffixes=("_all", "_t7"))
    return {"n_pairs_found_in_allsnp_file": int(len(H)), "share_found": float(len(H) / len(lead)),
            "beta_max_abs_diff": float((H.beta_all - H.beta_t7).abs().max()),
            "beta_corr": float(np.corrcoef(H.beta_all, H.beta_t7)[0, 1]) if len(H) > 2 else None}


def leads_t7(xlsx, archive, pattern):
    OUT.mkdir(parents=True, exist_ok=True)
    lead, info = table7_leads(xlsx)
    info["crosscheck_member"] = pattern
    info.update(crosscheck(lead, archive, pattern))
    lead.to_csv(OUT / "jerber_leads.tsv", sep="\t", index=False)
    json.dump(info, open(OUT / "jerber_leads.json", "w"), indent=1)
    print(info)


def expression_pcs(L, k=5, top=2000):
    v = L.var(0).sort_values(ascending=False).index[:top]
    Z = (L[v] - L[v].mean()) / L[v].std().replace(0, 1)
    U, s, _ = np.linalg.svd(Z.to_numpy(), full_matrices=False)
    return U[:, :k] * s[:k]


def replicate_one(lead, expr, G, snps, gpcs, rng):
    """expr: lines x Ensembl log2CPM (donor lines aligned to G rows)."""
    C = np.c_[gpcs[:, :3], expression_pcs(expr)]
    key = {(c, p): j for j, (c, p) in enumerate(zip(snps.chr, snps.pos))}
    G = G.astype(float); G[G < 0] = np.nan
    af = np.nanmean(G, 0) / 2; maf = np.minimum(af, 1 - af)
    mbin = np.digitize(maf, MAF_BINS[1:-1])
    obs, tested = [], []
    for r in lead.itertuples(index=False):
        if r.feature_id not in expr.columns:
            continue
        j = key.get((int(r.snp_chromosome), int(r.snp_position)))
        if j is None:
            continue
        s = align(r.ref, r.alt, r.assessed_allele, snps.a1[j], snps.a2[j])
        if s == 0:
            continue
        x = G[:, j] if s == 1 else 2 - G[:, j]
        ok = np.isfinite(x)
        if ok.sum() < 30 or np.nanstd(x[ok]) == 0:
            continue
        b, p = ols_t(expr[r.feature_id].to_numpy()[ok], x[ok], C[ok])
        obs.append((np.sign(b) == np.sign(r.beta), p)); tested.append((r.feature_id, int(r.snp_chromosome), mbin[j], np.sign(r.beta)))
    if not obs:
        return {"n_tested": 0}
    conc = float(np.mean([o[0] for o in obs])); p1 = pi1([o[1] for o in obs])
    null_c, null_p = [], []
    chrom = snps.chr.to_numpy()
    for _ in range(N_NULL):
        cc, pp = [], []
        for gene, ch, b_, sg in tested:
            cand = np.flatnonzero((chrom != ch) & (mbin == b_))
            j = rng.choice(cand); x = G[:, j]; ok = np.isfinite(x)
            b, p = ols_t(expr[gene].to_numpy()[ok], x[ok], C[ok])
            cc.append(np.sign(b) == sg); pp.append(p)
        null_c.append(np.mean(cc)); null_p.append(pi1(pp))
    c95, p95 = float(np.percentile(null_c, 95)), float(np.percentile(null_p, 95))
    return {"n_tested": len(obs), "sign_concordance": conc, "pi1": p1,
            "null95_concordance": c95, "null95_pi1": p95,
            "pass": bool(conc >= 0.70 and p1 >= 0.30 and conc > c95 and p1 > p95)}


def symbol_to_ensembl():
    """Map the unique gene symbols used by scrn-v0.1.0 / da-labels-v1.0.0 to Ensembl IDs, from the
    Cell Ranger features of one run (same reference for all runs)."""
    import h5py
    sys.path.insert(0, str(ROOT / "scripts"))
    from scrn_per_run import find_runs
    h5 = sorted(find_runs(ROOT / "data" / "unpacked").values())[0]
    with h5py.File(h5, "r") as f:
        ids = f["matrix/features/id"][:].astype(str); names = f["matrix/features/name"][:].astype(str)
        ft = f["matrix/features/feature_type"][:].astype(str)
    ids, names = ids[ft == "Gene Expression"], names[ft == "Gene Expression"]
    seen, uniq = {}, []
    for nm, d in zip(names, pd.Index(names).duplicated(keep="first")):
        if d:
            seen[nm] = seen.get(nm, 0) + 1; uniq.append(f"{nm}-{seen[nm]}")
        else:
            uniq.append(nm)
    return dict(zip(uniq, [i.split(".")[0] for i in ids]))


def foundin_inputs(kind, sym2ens):
    """lines x Ensembl log2(CPM+1) for donor lines, with matching dosage rows and genotype PCs."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import phenotypes_da as pdz
    d = pd.read_csv(ROOT / "data/derived/phenotypes_da/lines_da.csv")
    d = d[~d.is_control]
    runs = d.run.tolist()
    if kind in ("new", "old"):
        C = pdz.load_counts(runs, kind)
    else:
        C = pd.DataFrame({r: pd.read_csv(pdz.SCRN / "runs" / r / "pseudobulk.csv.gz", usecols=["gene", "All_QC_pass"])
                          .set_index("gene")["All_QC_pass"] for r in runs}).fillna(0.0).T.loc[runs]
    cpm = C.div(C.sum(1), axis=0) * 1e6
    L = np.log2(cpm.loc[:, cpm.mean(0) >= 10] + 1)
    L = L[[g for g in L.columns if g in sym2ens]]
    L.columns = [sym2ens[g] for g in L.columns]
    L = L.loc[:, ~L.columns.duplicated()].reset_index(drop=True)
    glines = [x.strip() for x in open(ROOT / "data/derived/geno/lines.txt") if x.strip()]
    gi = {l: i for i, l in enumerate(glines)}
    G = np.load(ROOT / "data/derived/geno/dosage.npz")["G"][[gi[l] for l in d.line]]
    pcs = pd.read_csv(ROOT / "data/derived/geno/pcs.tsv", sep="\t", index_col=0)
    pcs.index = pcs.index.astype(str)
    return L, G, pcs.loc[d.line.astype(str)].to_numpy()


def replicate(doi):
    lead = pd.read_csv(OUT / "jerber_leads.tsv", sep="\t")
    snps = pd.read_csv(ROOT / "data/derived/geno/snps.tsv", sep="\t")
    sym2ens = symbol_to_ensembl()
    res = {"registration_doi": doi, "n_jerber_leads": int(len(lead)),
           "rule": "PASS if concordance >= 0.70 and pi1 >= 0.30 and both above the null 95th percentile (primary: new label)"}
    for kind in ("new", "old", "all"):
        expr, G, gpcs = foundin_inputs(kind, sym2ens)
        res[kind] = replicate_one(lead, expr, G, snps, gpcs, np.random.default_rng(SEED))
        print(kind, res[kind], flush=True)
    res["G_eQTL_pass"] = res["new"].get("pass")
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump(res, open(OUT / "g_eqtl.json", "w"), indent=1, default=float)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["leads", "replicate"])
    ap.add_argument("--archive", default="D:/ST-DT4/data/external/jerber/eqtl/eqtl_summary_stats.tar.gz")
    ap.add_argument("--member", default=MEMBER)
    ap.add_argument("--source", choices=["table7", "allsnp"], default="table7")
    ap.add_argument("--table7", default="D:/ST-DT4/data/external/jerber/supp/esm3.xlsx")
    ap.add_argument("--registration-doi", help="required for `replicate`: the filed Cell-to-Clinic II OSF DOI")
    a = ap.parse_args(argv)
    if a.step == "leads":
        if a.source == "table7":
            leads_t7(a.table7, a.archive, a.member)
        else:
            leads(a.archive, a.member)
        return
    if not a.registration_doi:
        raise SystemExit("replicate is held until the Cell-to-Clinic II registration is filed: pass --registration-doi")
    replicate(a.registration_doi)


if __name__ == "__main__":
    main()
