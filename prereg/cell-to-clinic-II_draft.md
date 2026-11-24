# Cell-to-Clinic II: DA-neuron phenotypes and clinical progression (DRAFT for OSF)

Status: **final for filing**, reviewed and approved by the user 2026-11-24 (git tag `prereg-II-v1.0.0`). To be filed on OSF as a new registration linked to osf.io/s6w9c; the OSF text must match this file. **No clinical data will be linked to the phenotypes below until this is filed.**

## 1. Aims

Registration 1 (osf.io/s6w9c) tested all-cell day-65 phenotypes and found nothing supported (0 of 10 tests; 6 equivalent, 4 inconclusive). Since then, work without clinical data showed two things:
- the earlier DA label (TH ≥ 1 UMI) was not DA-specific;
- across-batch reliability is much lower than the within-batch control runs suggested (Jerber et al. 2021, 199 lines).

This registration tests phenotypes measured **inside DA neurons**, selected on reliability alone. It also adds a genetic positive-control gate for the new DA labels.

Aims:
1. **Confirmatory:** does SNCA expression in iPSC-derived DA neurons relate to motor progression in the same donors?
2. **Secondary:** lysosome score in DA neurons, DaTscan decline, SAA status and NSD-ISS stage.
3. **Framework:** report every gate and the design consequence, whatever the clinical result.

**Expected power (stated before any linkage).** With n = 47 donors, the slope reliability λ = 0.64 measured in E4, and phenotype reliability 0.62–0.83, the smallest true correlation detectable at 80% power is **about 0.58–0.68** for D2. Plausible effects (r 0.2–0.4) are out of reach. The tests are therefore mainly **effect-size bounding**: equivalence at ±0.40 has limited power (58% at n = 47 in the registration-1 simulation). We register them because a pre-specified, gated test is still informative, and because the gates and design numbers are an outcome in their own right.

## 2. Data and prior exposure (complete list)

- **Data:** FOUNDIN-PD Project 150 processed data only (day-65 scRNA, 83 runs: 79 donor lines plus 4 control-line runs; day-0/25/65 bulk RNA; NeuroChip genotypes). PPMI Study Data 2026-11-12 and Curated Cut v.2026-08-10. Jerber et al. 2021 (Zenodo 4651413, public).
- **Already seen with clinical data (registration 1):**
  - all-cell SNCA (P2), OXPHOS (P3) and lysosome (P4) against UPDRS III slope, SBR slope, SAA and NSD-ISS;
  - the exploratory analyses E1–E5.
  - In particular, all-cell P2 vs UPDRS III slope gave r −0.18 (inconclusive). The DA-only phenotypes below are **correlated with these**, so this is not a blind test of SNCA biology in general. It tests whether measuring inside DA neurons changes the answer.
- **Seen without clinical data, before this draft:**
  - new DA labels (B0, `da-labels-v1.0.0`);
  - DA-only phenotypes and their FOUNDIN control-run ICCs (B1, `da-phen-v1.0.0`);
  - Jerber across-batch ICCs (A1/A2, `jerber-a2-v0.1.0`);
  - genotype QC (C1).
  - **Phenotype selection below uses these reliability values, which were known when this was written.** No DA-only phenotype has been linked to any clinical, imaging or SAA variable.

## 3. Phenotypes (frozen in `scripts/phenotypes_da.py`, tag `da-phen-v1.0.0`)

- **DA neurons:** neurons (scrn-v0.1.0 compartment) with a 9-marker DA score above t = −0.10. The threshold was calibrated on public Table S3 (B0).
- **Run requirement:** ≥ 100 DA cells. No donor line is excluded by this rule.

| ID | Definition | FOUNDIN ICC (4 control runs, 3 df) | Jerber across-batch ICC (18 df) | Role |
|---|---|---|---|---|
| **D2** | SNCA log2(CPM+1), DA pseudobulk | 0.83 | 0.62 (95% CI 0.14–0.79) | **Primary** |
| D4 | KEGG lysosome mean gene-z, DA pseudobulk | 0.88 | 0.43 | Secondary |
| D3 | HALLMARK OXPHOS mean gene-z, DA pseudobulk | 0.46 | 0.18 | Exploratory only (fails both) |
| D1 | logit DA fraction (new rule) | −0.24 | 0.38 (P1 analogue) | Not tested (fails) |

**Retention rule:** FOUNDIN ICC ≥ 0.5 **and** Jerber across-batch ICC ≥ 0.40, both point estimates. The second condition is new. It is justified by E5c and A2, which showed that within-batch replicates overstate reliability.

**Variants (sensitivity, not confirmatory):**
- **D*_old:** the same phenotypes from the old TH ≥ 1 DA label.
- **D*_adj:** residualised on new-rule DA fraction and marker progenitor fraction.

**Standardisation:** to mean 0 and SD 1 over the analysed PD donors, as in registration 1.

## 4. Gates (pass rules fixed here)

| Gate | Question | Rule | Status at drafting |
|---|---|---|---|
| G-cells | Enough DA cells per run? | ≥ 100 new-rule DA cells | Pass (0 lines excluded) |
| G-icc | Reliable across batches? | As in §3 | D2 and D4 pass; D1 and D3 fail |
| **G-eQTL** (new) | Do the new DA-neuron pseudobulks carry known DA-neuron genetic regulation? | Lead eQTLs from Jerber DA neurons (day 52, untreated; per gene the min-p SNP, gene-level empirical p, BH FDR < 5%), tested in FOUNDIN DA pseudobulk (79 donor lines; log2 CPM ~ dosage + 3 genotype PCs + 5 expression PCs; SNPs matched on GRCh37 position and alleles, strand-ambiguous SNPs dropped). Pass: sign concordance ≥ 70% **and** π1 ≥ 0.3 (Storey, λ = 0.5), both above the 95th percentile of a null of 200 sets pairing the same lead genes with random SNPs on other chromosomes in the same MAF bin. Code: `scripts/eqtl_replication.py`, tagged before the run. | **Decided 2026-11-24: run after filing and before any linkage.** The FOUNDIN step refuses to run without the filed registration DOI. |
| A0-1 | Linkage | As registration 1 | Pass (44/47) |
| A0-5 | Coherence (D4 only) | ≥ 80% of member genes share the composite's sign in the H1 model | Evaluated at linkage |

**If G-eQTL fails:** D2 and D4 are still tested, but the paper reports them as phenotypes of "marker-defined neurons" and does not claim DA identity. The old-label results are given equal prominence.

## 5. Hypotheses and models

- **Model code:** frozen from registration 1 (`scripts/models.py`, tag `prereg-v1.0.1`), the same REML LMM, Satterthwaite df, covariates, OFF-state definition, time axis and fallbacks, with D2/D4 in place of P2–P4. Any change needed only to read the new phenotype file goes into a new tag before linkage.
- **H1-DA (primary, single test, α = 0.05 two-sided):** the D2 × time interaction on OFF-state MDS-UPDRS III.
- **Secondary family (BH across 7 tests):**
  - H1-DA on D4 (if A0-5 passes);
  - H2-DA, D2 and D4 × time on putamen SBR;
  - H3-DA, D2 and D4 by SAA status, and D2 and D4 by NSD-ISS stage.
- **Classes:** supported (q < 0.05), equivalent (90% CI of partial r within ±0.40), inconclusive. TOST uses the same mapping as registration 1.

## 6. Additional confirmatory test (D1, bulk OXPHOS)

The exploratory finding from registration 1 was all-cell OXPHOS vs UPDRS III slope, r −0.33. It was seen in data, failed A0-5, and was not tested confirmatorily.
- **H4 (single test, α = 0.05, same model):** HALLMARK OXPHOS mean gene-z in **day-65 bulk RNA** (an independent assay of the same lines), composition-adjusted as in registration 1, × time on UPDRS III.
- **Not an independent-donor replication.** It tests whether the signal survives a change of assay. An independent-cohort replication (D2 of the Phase 2 plan) is not available.

## 7. Sensitivity analyses
1. Old DA label (D*_old).
2. Composition-adjusted (D*_adj).
3. Idiopathic PD only.
4. Five-level batch covariate.
5. Leave one donor out (H1-DA).
6. First 5 years after diagnosis.

## 8. Reporting commitments
- Every gate result is reported, including failures.
- Tests a gate removes are reported as "not tested", with their estimates given descriptively.
- The design consequence is reported with the result: the donors needed for r = 0.3 at the measured reliabilities, from the A2 and A4 design guide.
- Code tags are frozen before each step. Deviations are listed with reasons.
- Only aggregate results are published; there are no participant-level values.
