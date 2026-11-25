"""DESCRIPTIVE, POST HOC (added 2026-11-25 after G-eQTL ran): where the Jerber leads drop out before
testing in the new-label DA pseudobulk. Uses the g-eqtl-v1.1.0 functions unchanged. Counts only.
Output: results/c2c2_eqtl_attrition.json"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import json
import numpy as np, pandas as pd
import eqtl_replication as er

lead = pd.read_csv(er.OUT / "jerber_leads.tsv", sep="\t")
snps = pd.read_csv(er.ROOT / "data/derived/geno/snps.tsv", sep="\t")
expr, G, gpcs = er.foundin_inputs("new", er.symbol_to_ensembl())
key = {(c, p): j for j, (c, p) in enumerate(zip(snps.chr, snps.pos))}
n = {"leads": len(lead)}
a = lead[lead.feature_id.isin(expr.columns)]; n["gene_expressed_cpm10"] = len(a)
pos = [key.get((int(c), int(p))) for c, p in zip(a.snp_chromosome, a.snp_position)]
b = a[[j is not None for j in pos]]; jj = [j for j in pos if j is not None]; n["snp_on_array_same_pos"] = len(b)
al = [er.align(r.ref, r.alt, r.assessed_allele, snps.a1[j], snps.a2[j]) for r, j in zip(b.itertuples(), jj)]
n["alleles_aligned_not_ambiguous"] = int(sum(s != 0 for s in al))
# how many leads would be on the array at all, regardless of expression
n["snp_on_array_any_gene"] = int(sum(key.get((int(c), int(p))) is not None for c, p in zip(lead.snp_chromosome, lead.snp_position)))
n["array_snps_total"] = len(snps)
json.dump({k: int(v) for k, v in n.items()}, open(er.ROOT / "results" / "c2c2_eqtl_attrition.json", "w"), indent=1)
print(n)
