"""Synthetic checks for scripts/eqtl_replication.py (no real data)."""
import sys, pathlib
from pathlib import Path
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pandas as pd
import eqtl_replication as er


def test_align():
    assert er.align("G", "A", "A", "A", "G") == 1      # a1 is the assessed allele
    assert er.align("G", "A", "A", "G", "A") == -1     # a1 is the other allele
    assert er.align("G", "A", "A", "T", "C") == 1      # opposite strand, a1 = comp(A)
    assert er.align("A", "T", "T", "T", "A") == 0      # ambiguous
    assert er.align("G", "A", "A", "C", "T") == -1     # opposite strand, a1 = comp(G)


def test_pi1_and_bh():
    rng = np.random.default_rng(0)
    p = np.r_[rng.uniform(size=5000), rng.beta(0.05, 1, 5000)]
    assert abs(er.pi1(p) - 0.5) < 0.05
    assert abs(er.pi1(rng.uniform(size=20000))) < 0.03
    q = er.bh([0.01, 0.04, 0.03, 0.5]); assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.5], atol=1e-6)


def test_lead_table():
    ch1 = pd.DataFrame({"feature_id": ["g1", "g1", "g2"], "snp_id": ["1_1_A_G", "1_2_A_G", "1_3_C_T"],
                        "p_value": [1e-3, 1e-8, 0.2], "beta": [0.1, 0.5, 0.1], "empirical_feature_p_value": [-1, 1e-6, -1],
                        "snp_chromosome": 1, "snp_position": [1, 2, 3], "assessed_allele": "G", "maf": 0.2})
    ch2 = pd.DataFrame({"feature_id": ["g2"], "snp_id": ["1_4_C_T"], "p_value": [0.3], "beta": [0.0],
                        "empirical_feature_p_value": [0.9], "snp_chromosome": 1, "snp_position": [4],
                        "assessed_allele": "T", "maf": 0.3})
    T = er.lead_table([ch1, ch2]).set_index("feature_id")
    assert T.loc["g1", "snp_id"] == "1_2_A_G" and T.loc["g1", "empirical_feature_p_value"] == 1e-6
    assert T.loc["g2", "snp_id"] == "1_3_C_T" and T.loc["g2", "empirical_feature_p_value"] == 0.9


def test_replicate_detects_true_eqtls_and_null_does_not():
    rng = np.random.default_rng(1)
    n, ngene = 79, 60
    snps = pd.DataFrame({"chr": np.repeat(np.arange(1, 23), 40), "pos": np.tile(np.arange(40) * 1000 + 5, 22),
                         "a1": "A", "a2": "G"})
    maf = rng.uniform(0.1, 0.5, len(snps))
    G = rng.binomial(2, maf, (n, len(snps))).astype(np.int8)
    gpcs = rng.normal(size=(n, 10))
    expr = pd.DataFrame(rng.normal(size=(n, 2500)), columns=[f"E{i}" for i in range(2500)])
    rows = []
    for g in range(ngene):
        j = g * 7; beta = rng.choice([-1, 1]) * 0.8
        expr[f"E{g}"] += beta * (G[:, j] - G[:, j].mean())          # dosage counts A = assessed
        rows.append({"feature_id": f"E{g}", "snp_chromosome": snps.chr[j], "snp_position": snps.pos[j],
                     "ref": "G", "alt": "A", "assessed_allele": "A", "beta": beta})
    er.N_NULL = 20
    res = er.replicate_one(pd.DataFrame(rows), expr, G, snps, gpcs, rng)
    assert res["n_tested"] == ngene and res["sign_concordance"] > 0.9 and res["pi1"] > 0.6 and res["pass"]
    assert res["null95_concordance"] < 0.7 and res["null95_pi1"] < 0.3


def test_member_pattern():
    import re
    from eqtl_replication import MEMBER
    d = "eqtl_summary_stats_renamed/"
    assert re.search(MEMBER, d + "D52.DA.untreated.qtl_results_all.sorted.txt.gz")
    for bad in ("D52.DA.ROT_treated.qtl_results_all.sorted.txt.gz", "D52.DA.untreated.qtl_results_all.sorted.txt.gz.tbi",
                "D52.Astro.untreated.qtl_results_all.sorted.txt.gz", "D11.DA.untreated.qtl_results_all.sorted.txt.gz"):
        assert not re.search(MEMBER, d + bad)


def test_table7_leads():
    import tempfile
    from eqtl_replication import table7_leads, T7_SHEET
    cols = ["snp_id", "p_value", "beta", "beta_se", "empirical_feature_p_value", "feature_chromosome", "feature_start",
            "feature_end", "n_samples", "snp_chromosome", "snp_position", "assessed_allele", "maf", "hwe_p", "BH_fdr",
            "q_value", "label", "ensembl_gene_id", "hgnc_symbol"]
    rows = [["1_100_C_T", 1e-6, 0.3, .05, .001, 1, 1, 2, 170, 1, 100, "T", .2, .5, .01, .01, "DA_D52", "ENSG1", "A"],
            ["1_200_G_A", 1e-5, -.2, .05, .02, 1, 1, 2, 170, 1, 200, "A", .2, .5, .08, .04, "DA_D52", "ENSG2", "B"],
            ["2_300_A_G", 1e-7, 0.4, .05, .001, 2, 1, 2, 170, 2, 300, "G", .3, .5, .01, .01, "DA_D52_ROT", "ENSG3", "C"],
            ["2_400_A_C", 1e-7, 0.4, .05, .001, 2, 1, 2, 170, 2, 400, "C", .3, .5, .01, .01, "DA_D30", "ENSG4", "D"]]
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "t7.xlsx"
        with pd.ExcelWriter(f) as w:
            pd.DataFrame([["List of eQTL discoveries"], ["All lead eQTL SNP-gene pairs are provided."]]).to_excel(
                w, sheet_name=T7_SHEET, header=False, index=False)
            pd.DataFrame(rows, columns=cols).to_excel(w, sheet_name=T7_SHEET, index=False, startrow=3)
        lead, info = table7_leads(f)
    assert list(lead.feature_id) == ["ENSG1"] and info["n_in_context"] == 2 and info["n_q_lt_05"] == 2
    assert lead.ref.iloc[0] == "C" and lead.alt.iloc[0] == "T" and lead.assessed_allele.iloc[0] == "T"


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f(); print("ok", k)
    print("ALL PASSED")

