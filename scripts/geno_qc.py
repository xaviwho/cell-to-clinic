"""Phase 2, Workstream C1: genotype QC for the FOUNDIN-PD NeuroChip array (PLINK .bed/.bim/.fam).
No expression or clinical data is used here. Aggregates only are printed.

Steps (fixed in advance):
 1. Read PLINK (SNP-major .bed) with numpy; 0/1/2 = count of A1 alleles, missing = NaN.
 2. Sample IDs are <LINE>_day<D>[_<batch>]. Each line is genotyped at two culture days; identity check =
    genotype concordance between a line's samples over SNPs called in both (report min / median).
    A line is flagged if any pair < 0.98.
 3. Per line, one genotype vector: the mean over its samples where called (identical samples agree),
    rounded to 0/1/2 (discordant calls -> missing).
 4. SNP QC on all genotyped FOUNDIN lines (95; includes the 80 scRNA lines): autosomes 1-22; call rate >= 0.95; MAF >= 0.05;
    HWE exact-test p >= 1e-6 (Wigginton et al. 2005).
 5. Line QC: call rate >= 0.95 after SNP QC.
 6. Genotype PCs: top 10 from LD-thinned SNPs (every SNP >= 50 kb from the previous kept one on the
    same chromosome), z-scored dosages, missing set to the mean.
Output (data/derived/geno/): dosage.npz (lines x SNPs, int8, -1 = missing), snps.tsv (chr, id, pos,
A1, A2), lines.txt, pcs.tsv, geno_qc.json (counts only).
"""
from __future__ import annotations

import json, re, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GDIR = ROOT / "data" / "unpacked_geno"
OUT = ROOT / "data" / "derived" / "geno"
ID_RE = re.compile(r"^(PPMI\d+)_day(\d+)(?:_(B\w+))?$")


def read_plink(prefix):
    fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, names=["fid", "iid", "p", "m", "sex", "ph"])
    bim = pd.read_csv(f"{prefix}.bim", sep="\t", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
    n, m = len(fam), len(bim)
    raw = np.fromfile(f"{prefix}.bed", dtype=np.uint8)
    assert raw[0] == 0x6C and raw[1] == 0x1B and raw[2] == 1, "not SNP-major PLINK bed"
    nb = (n + 3) // 4
    b = raw[3:].reshape(m, nb)
    codes = np.stack([(b >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(m, nb * 4)[:, :n]
    # PLINK: 00 hom A1 (2 copies of A1), 01 missing, 10 het, 11 hom A2
    lut = np.array([2, np.nan, 1, 0], dtype=float)
    return fam, bim, lut[codes].T  # samples x SNPs


def hwe_p(n_aa, n_ab, n_bb):
    """Exact SNP-HWE test (Wigginton, Cutler & Abecasis 2005)."""
    n_ab, n_aa, n_bb = int(n_ab), int(n_aa), int(n_bb)
    n = n_aa + n_ab + n_bb
    rare = min(2 * n_aa + n_ab, 2 * n_bb + n_ab)
    if n == 0:
        return 1.0
    probs = np.zeros(rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if (rare - mid) % 2:
        mid += 1
    probs[mid] = 1.0
    hom_r = (rare - mid) // 2; hom_c = n - mid - hom_r
    for het in range(mid, 1, -2):
        probs[het - 2] = probs[het] * het * (het - 1) / (4.0 * (hom_r + 1) * (hom_c + 1)); hom_r += 1; hom_c += 1
    hom_r = (rare - mid) // 2; hom_c = n - mid - hom_r
    for het in range(mid, rare - 1, 2):
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2) * (het + 1)); hom_r -= 1; hom_c -= 1
    probs /= probs.sum()
    obs = n_ab
    return float(min(1.0, probs[probs <= probs[obs] + 1e-12].sum()))


def main():
    prefix = str(next(GDIR.rglob("*.bed"))).removesuffix(".bed")
    fam, bim, G = read_plink(prefix)
    ids = fam["iid"].astype(str)
    m = ids.str.extract(ID_RE)
    if m[0].isna().any():
        sys.exit("unrecognised sample ID format")
    fam["line"] = m[0]
    info = {"n_samples": int(len(fam)), "n_variants": int(len(bim)), "n_lines": int(fam.line.nunique())}

    # identity concordance within line
    conc = []
    for line, idx in fam.groupby("line").groups.items():
        idx = list(idx)
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                x, y = G[idx[a]], G[idx[b]]
                ok = ~np.isnan(x) & ~np.isnan(y)
                conc.append((line, float((x[ok] == y[ok]).mean())))
    c = np.array([v for _, v in conc])
    info["identity_pairs"] = int(len(c)); info["identity_concordance_min"] = float(c.min())
    info["identity_concordance_median"] = float(np.median(c)); info["lines_flagged_lt_0.98"] = int(len({l for l, v in conc if v < 0.98}))

    # one vector per line
    lines = sorted(fam.line.unique())
    L = np.full((len(lines), G.shape[1]), np.nan)
    for i, line in enumerate(lines):
        S = G[fam.index[fam.line == line]]
        mean = np.nanmean(np.where(np.isnan(S), np.nan, S), axis=0) if len(S) > 1 else S[0]
        agree = np.nanmax(S, 0) == np.nanmin(S, 0) if len(S) > 1 else ~np.isnan(S[0])
        L[i] = np.where(agree, np.round(mean), np.nan)

    # SNP QC
    auto = bim.chr.isin(range(1, 23)).to_numpy()
    call = (~np.isnan(L)).mean(0)
    af = np.nanmean(L, 0) / 2
    maf = np.minimum(af, 1 - af)
    keep = auto & (call >= 0.95) & (maf >= 0.05)
    info["snps_autosomal"] = int(auto.sum()); info["snps_after_call_maf"] = int(keep.sum())
    hw = np.ones(L.shape[1])
    for j in np.where(keep)[0]:
        g = L[:, j]; g = g[~np.isnan(g)]
        hw[j] = hwe_p((g == 2).sum(), (g == 1).sum(), (g == 0).sum())
    keep &= hw >= 1e-6
    info["snps_final"] = int(keep.sum())
    Lk = L[:, keep]; bk = bim[keep].reset_index(drop=True)
    lcall = (~np.isnan(Lk)).mean(1)
    info["lines_call_lt_0.95"] = int((lcall < 0.95).sum())

    # PCs on thinned SNPs
    sel, last = [], {}
    for j, (ch, pos) in enumerate(zip(bk.chr, bk.pos)):
        if ch not in last or pos - last[ch] >= 50_000:
            sel.append(j); last[ch] = pos
    X = Lk[:, sel]
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0); sd[sd == 0] = 1
    Z = np.where(np.isnan(X), 0, (X - mu) / sd)
    U, s, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    pcs = U[:, :10] * s[:10]
    info["pc_snps"] = int(len(sel)); info["pc_var_explained_top3"] = [float(v) for v in (s[:3] ** 2 / (s ** 2).sum())]

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "dosage.npz", G=np.where(np.isnan(Lk), -1, Lk).astype(np.int8))
    bk[["chr", "id", "pos", "a1", "a2"]].to_csv(OUT / "snps.tsv", sep="\t", index=False)
    (OUT / "lines.txt").write_text("\n".join(lines))
    pd.DataFrame(pcs, index=lines, columns=[f"PC{k+1}" for k in range(10)]).to_csv(OUT / "pcs.tsv", sep="\t")
    json.dump(info, open(OUT / "geno_qc.json", "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
