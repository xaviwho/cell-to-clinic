# Activity log: FOUNDIN-PD cell-to-clinic

Progress record for the FOUNDIN-PD (PPMI Project 150) cell-to-clinic analysis. The newest entry is first.

> **Public repository.** This log holds no participant-level information: no PATNOs, cell-line (CDI) numbers, per-line values, credentials or data. Controlled-access data stays on the local machine under `data/`, which is git-ignored.

> **Two sessions share this log.** The **local session** (Claude in VS Code, on the machine with the data) is the **main environment**: from 2026-11-12 it does all the work, including scripts, analyses, commits and tags. The **planning session** (Claude on claude.ai, no data access) reads this log and can raise requests or reviews. Each entry is tagged with which session wrote it. Requests go under **"For the local session"** / **"For the planning session"** headings. Pull before starting, commit and push when done. Only aggregate, non-identifying results go here.

---

## 2026-11-12 (local session): figures started: Fig. 3 and Fig. 1 done

Built from aggregate outputs only (`models.json`, gate JSONs). Scripts in `scripts/figures/`; PNG (300 dpi) and PDF in `paper/figures/`.

- **Fig. 3** (`fig3_forest`): equivalence forest plot. All 10 confirmatory tests (partial r, 90% CI) in four hypothesis blocks, with the ±0.40 equivalence region shaded. The two H1 rows dropped by A0-5 are shown as hollow grey markers, labelled "not tested". Each row has n and its class/q in a separate text column.
  - Colours: equivalent = blue, inconclusive = orange (reference palette slots 1–2; validated: CVD ΔE 24.7, normal-vision ΔE 33.6).
  - Class is also shown by marker shape and text, so the figure reads in greyscale and for colour-blind readers.
- **Fig. 1** (`fig1_gates`):
  - (A) framework schematic: pipeline-validity gates, and the phenotype flow 4 → A0-2 → 3 → H2/H3, with A0-5 → P2 only in H1; the result is 10 tests, 0 supported, 6 equivalent, 4 inconclusive.
  - (B) every gate's measured value against its threshold (A0-6, A0-4 ×2 with the partial threshold, A0-2 ×4, A0-5 ×2, A0-1).
  - Pass / partial / fail use the reserved status colours, always with a marker shape and a text label.
- **Next:** Fig. 2 (A0-6 scatter + A0-4 ROC; this needs per-fold aggregate ROC points, so a small export step will be added), Fig. 4 (variance decomposition, E5b/E5c + `explore-v0.3.0`), Fig. 5 (design curves with E5a uncertainty bands).

---

## 2026-11-12 (local session): **EXPLORATORY, POST HOC**: bulk SNCA composition check (`explore-v0.3.0`)

Approved by the user. One analysis, specified in the header of `scripts/exploratory3.py` and frozen as tag `explore-v0.3.0` (commit `f6e3267`) before running. A synthetic test plants a batch → composition → SNCA effect, and adjustment recovers it.

**Design:**
- Bulk log2(CPM + 1) SNCA, residualised across donor lines on TH and a pan-neuronal score (mean z of STMN2, SNAP25, SYT1). Markers fixed in advance.
- Day 65: ICC with the control line's 5 runs; F-based 95% CI.
- Day 25: between-batch share of control-run variance (9 batch-1 runs + 4 others); stratified bootstrap 95% CI.

| | Unadjusted bulk SNCA | Composition-adjusted bulk SNCA | scRNA SNCA (P2, single-cell composition) |
|---|---|---|---|
| Composition R² across donor lines | n/a | **0.82** at day 65; 0.48 at day 25 | n/a |
| Day-65 ICC, 79 scRNA lines (95% CI) | −0.15 (−8.6 to 0.61) | **0.36** (−4.4 to 0.78) | 0.94 |
| Day-65 ICC, all 91 lines | 0.02 | **0.39** (−4.1 to 0.79) | n/a |
| Day-25 between-batch share of run variance | 0.97 (0.72–0.99) | **0.71** (0.26–0.93) | n/a |
| Day-25 ICC, all 13 control runs | 0.54 (−0.28 to 0.78) | **0.87** (0.62 to 0.94) | n/a |

### Interpretation (post hoc)
- **Composition explains most of the variation in bulk SNCA:** 82% of between-line variance at day 65.
- **Adjusting for it moves reliability in the predicted direction at both time points:**
  - day 65: ≈ 0 → 0.36–0.39;
  - day 25: 0.54 → 0.87, with the batch share falling 0.97 → 0.71.
- **Bulk marker adjustment recovers only part of the reliability.** The day-65 ICC stays below 0.5 and far below the single-cell-adjusted 0.94; the CIs are very wide (4 df). The likely reason is that three pan-neuronal genes and TH approximate composition far more crudely than single-cell cluster proportions.
- **Wording for the paper:** "composition is the main route by which batch affects SNCA, and **single-cell composition** adjustment is needed to make it a reliable donor trait; marker-based bulk adjustment helps but is not sufficient." This is post hoc, from one control line.
- No variants were run, and none will be.

### Also done
- **E4 wording fix applied** in the E4 entry and the status summary: "given the measured reliabilities, the smallest detectable true correlation was 0.55 (SNCA) to 0.92 (DA yield)"; separately, n = 47 alone limits detection to about r ≥ 0.43.
- **Data handling from now on:** no `data/` listings or file names in chat; counts or hashed IDs only.

### Next
Figures, from aggregate outputs only, in the planned order: Fig. 3 (equivalence forest plot), then Fig. 1, 2, 4, 5. Saved to `paper/figures/`; scripts in `scripts/figures/`.

---

## 2026-11-12 (planning session): user decisions: bulk SNCA check approved; paper direction confirmed

### Decided by the user
1. **Approved:** one post-hoc exploratory check, composition-adjusted bulk SNCA (`explore-v0.3.0`).
2. **Confirmed:** the paper direction and outline in the planning review below (TCBB; two contributions: the gated framework and the empirical design guide).

### For the local session, in order
1. **`explore-v0.3.0`: bulk SNCA composition check.** Write the specification into the script header and tag it before running.
   - **Question:** does adjusting bulk SNCA for neuronal/DA composition raise its run-to-run reliability (bulk ICC ≈ 0, against scRNA composition-adjusted 0.94)?
   - **Design:** residualise bulk log-CPM SNCA on bulk composition markers. Fix the marker set in advance (e.g. TH plus a pan-neuronal marker score such as STMN2/SNAP25/SYT1). Then recompute the ICC as in E5b (day 65; 5 control runs) and the batch share as in E5c (day 25).
   - **Report:** ICC and batch share, unadjusted vs adjusted, with CIs. One analysis only; no variants after seeing results.
   - **Label it post hoc.** Log the result here, whichever way it comes out.
2. **Apply the E4 wording fix** (planning review, point 3) wherever the summary is reused.
3. **Start the figures**, from aggregate outputs only, in this order: Fig. 3 (equivalence forest plot), Fig. 1 (framework and gate flow), then Figs 2, 4, 5. Save to `paper/figures/` with the scripts in `scripts/figures/`.
4. **Data handling:** no directory listings or file names from `data/` in chat. Use counts or hashed IDs.

### For the user
- Check whether the PPMI DUA requires reporting the participant identifier printed into the local chat (disclosure 4).
- Aspera access ends about 2026-12-11. Nothing more is planned to download, but the published per-cell labels request (LONI/DZNE) is still open if you want it.

---

## 2026-11-12 (planning session): review of the confirmatory and exploratory results; paper outline

Careful, well-documented work. Every analysis was frozen before it ran, the bug fix is logged, and failures are kept. I agree with the summary and the recommendation (experimental phase complete), with the points below.

### Review notes
1. **Results and classes are correct as registered:** 10 tests, 0 supported, 6 equivalent, 4 inconclusive; gates applied mechanically.
2. **OXPHOS in H1 (r −0.33, 90% CI −0.55 to −0.04).** Its CI excludes zero, but it is not a confirmatory result: A0-5 dropped the composite before testing. In the paper:
   - show it in the "reported only" row, with the gate decision beside it;
   - describe the gene-level direction (E1) as "consistent but not significant";
   - name it as the one candidate signal for replication.

   Don't put it in the abstract as a finding.
3. **Wording fix (E4, point 1).** "Even a perfectly measured study of 47 donors could only have detected r ≥ 0.55" mixes two ideas. 0.55 already includes attenuation by the measured reliabilities. Say instead: "given the measured reliabilities, the smallest detectable true correlation was 0.55 (SNCA) to 0.92 (DA yield)". Separately, n = 47 alone limits detection to about r ≥ 0.43.
4. **E5 rests on one control line** (4–13 runs). The batch-share and design numbers describe that line. Present them as an empirical example with wide uncertainty (E5a gives it), not as general constants.
5. **Disclosure 4 (a PATNO and CDI printed into the local chat) is more than cosmetic.** A participant identifier went to a third-party service.
   - Record it in the prior-exposure and data-handling section.
   - Check whether your DUA requires reporting it.
   - From now on: no directory listings of `data/` in chat; use counts or hashed IDs.
6. **The optional bulk SNCA composition check (summary 6.2).** It's worth doing only because it tests the paper's central design claim directly ("composition adjustment is what makes SNCA reliable"). Freeze it first, keep it to one analysis, and stop there. User's call.
7. **Release:** agree on `paper-v1.0` + Zenodo, with aggregate outputs only. Also cite the OSF registration and both A0-4 versions.

### Paper outline (TCBB): *Cell-to-Clinic: A Gated Framework Linking iPSC Neuron Phenotypes to Parkinson's Progression*
1. **Introduction.** The promise of patient iPSC "avatars"; the untested assumption that dish phenotypes track the patient's clinical course; why a pre-registered, gated test is needed.
2. **The framework (main methodological contribution).** Gates A0-1 to A0-6, each with a pass rule and a positive control; how they order the analysis; equivalence testing; pre-registration plus frozen code. Written to apply beyond FOUNDIN-PD.
3. **Data.** FOUNDIN-PD Project 150 (80 lines, day-65 scRNA, day-0/25/65 bulk) linked to PPMI; access, timeline and firewall (download dates, OSF, Zenodo).
4. **Results I: the gates.**
   - A0-6 reproduction;
   - A0-4 two versions (fail, then partial);
   - A0-2 drops DA yield;
   - A0-5 drops both composites.

   4 phenotypes in, 1 tested in H1.
5. **Results II: confirmatory tests.** 10 tests: 0 supported, 6 equivalent (DaTscan decline; 2 NSD-ISS), 4 inconclusive; the 9 sensitivity analyses.
6. **Results III: why (exploratory, labelled).**
   - E4/E5 measurement properties: run-to-run noise is batch-driven, composition adjustment is the key, lysosome is robust;
   - design curves;
   - E1–E3 briefly.
7. **Discussion.**
   - What the null means and doesn't mean;
   - the design guide for future studies (5 points);
   - limitations: one control line, n = 47, power discrepancy, iPSC aging reset, prior exposure and disclosures.
8. **Reproducibility statement.** OSF, Zenodo DOIs, tags, ACTIVITY.md.

### Figures
1. **Framework schematic + gate flow:** 4 phenotypes → A0-2 → A0-5 → tests, with each gate's measured value against its threshold.
2. **Reproduction and positive control:** A0-6 scatter (r = 0.999); A0-4 ROC for v0.1.0 and v0.2.0 against the published value.
3. **Equivalence forest plot:** all 10 confirmatory tests (partial r, 90% CI, ±0.40 bounds, class colour); the reported-only rows in grey.
4. **Variance decomposition:** within-batch vs between-batch run variance per phenotype (E5c), and scRNA vs bulk ICC (E5b).
5. **Design curves:** donors needed against true r, per phenotype, for 1 vs 4 runs per line, with uncertainty bands (E4/E5a).
6. **Supplementary:** sensitivity-analysis grid; E1 gene-level volcano; E2 milestone hazard ratios; E3 carriers; power simulation vs analytic.

Tables: 1, cohort; 2, gates; 3, confirmatory results; S-tables, sensitivity analyses and exploratory results.

### For the local session
- Apply the E4 wording fix (point 3) wherever the summary is reused.
- If the user approves point 6, freeze and run the bulk SNCA composition check (`explore-v0.3.0`) as one analysis only.
- Then start the figures in the order above (Fig. 3 and Fig. 1 first), from aggregate outputs only.

---

## 2026-11-12 (local session): **PROJECT STATUS SUMMARY AND RECOMMENDATION**: all registered and exploratory analyses complete

A single overview of everything so far; the detailed entries are below. Aggregates only.

### 1. Timeline and frozen versions
| Tag | What it froze | Run on real data |
|---|---|---|
| `scrn-v0.1.0` | scRNA QC, %TH+, cell types, pseudobulk (83 runs, 80 lines) | yes: gate A0-6 |
| `day0-gate-v0.1.0` | A0-4 v1 (our %TH+ outcome) | yes: **fail** (CV AUC 0.57) |
| `day0-gate-v0.2.0` | A0-4 v2 (exact replication of the published classifier) | yes: **partial** (CV AUC 0.72 vs 0.85; 5/5 published genes replicate) |
| `prereg-v1.0.0` | All registered code (P1–P4, gates, H1–H3, sensitivity 1–9); synthetic-tested; **Zenodo 10.5281/zenodo.24156791** | yes, after OSF approval |
| `prereg-v1.0.1` | Name-only PPMI config fixes and an xlsx→csv helper; analysis code byte-identical | yes |
| `explore-v0.1.0` | Exploratory E1–E4 | yes |
| `explore-v0.2.0` / `v0.2.1` | Exploratory E5 (`v0.2.1` = ID-mapping bug fix) | yes |

**Registration:** OSF [osf.io/s6w9c](https://osf.io/s6w9c), filed and approved 2026-11-12, embargoed to 2027-11-12.

**Clinical data:** PPMI Study Data downloaded 2026-11-12, plus Curated Data Cut v.2026-08-10. Downloaded only after OSF approval and the Zenodo archive, as registered.

### 2. Gates
| Gate | Rule | Result |
|---|---|---|
| A0-6 Reproduction | r ≥ 0.80 with Table S1 %TH+ | **pass**, r = 0.999 (n = 79) |
| A0-4 Pipeline validity | CV AUC ≥ 0.75 (partial 0.65–0.75, p < 0.05) | v0.1.0 **fail**; v0.2.0 **partial**, 0.72 (p = 0.02); proceed, reported as a limitation |
| A0-3 Cell sufficiency | ≥ 1,000 QC cells | 1 of 79 lines excluded (not a PD donor) |
| A0-2 Technical stability | ICC ≥ 0.5 | **P1 DA fraction fails (0.33), dropped**; P2 0.94, P3 0.55, P4 0.64 pass |
| A0-1 Linkage | ≥ 30 PD donors with ≥ 3 OFF visits over ≥ 2 years | **pass**, 44 / 47 |
| A0-5 Composite coherence | ≥ 80% of genes share the composite's sign | **P3 fails (75%), P4 fails (62%)**: both dropped from H1 |

### 3. Confirmatory results (n = 47 PD donors; 10 tests)
| Hypothesis | Tests | Supported | Equivalent | Inconclusive |
|---|---|---|---|---|
| H1 motor progression (OFF-state UPDRS III slope) | P2 only | 0 | 0 | 1 (r −0.18) |
| H2 DaTscan putamen decline (n = 38) | P2, P3, P4 | 0 | **3** | 0 |
| H3 SAA status (7 negative / 38 positive) | P2, P3, P4 | 0 | 1 | 2 |
| H3 NSD-ISS stage (n = 39) | P2, P3, P4 | 0 | 2 | 1 |
| **Total** | 10 | **0** | **6** | **4** |

- **Sensitivity analyses 1–9:** nothing reaches q < 0.05; leaving out any one donor never changes a sign.
- **Power:** the registered equivalence power was optimistic (58% by simulation vs 83% registered at n = 47). Reported as a discrepancy.

### 4. Exploratory results (labelled as such)
- **E1, genes:** no single OXPHOS or lysosome gene is significant; the nominal hits match chance. 75% of OXPHOS genes point toward "higher OXPHOS, slower progression", led by ATP-synthase and complex III subunits. This is one weak, coherent signal, not a finding.
- **E2, milestones** (H&Y ≥ 3, freezing, cognitive impairment): no associations; HR CIs about 0.55–1.7 per SD. Falls were not available.
- **E3, carriers:** manifesting (19) vs non-manifesting (21) carriers don't differ in P2, P3 or P4 (all **equivalent**).
- **E4, design:**
  - Reliability of each patient's progression rate is λ = 0.64 (motor) and 0.54 (DaTscan).
  - Given the measured reliabilities, the smallest detectable true correlation was 0.55 (SNCA) to 0.92 (DA yield). *(Wording corrected 2026-11-12.)* Separately, n = 47 alone limits detection to about r ≥ 0.43.
  - Plausible effects (r 0.2–0.4) were out of reach.
- **E5a, design uncertainty:** SNCA needs about **164 donors (141–234)** for r = 0.3. For DA yield, OXPHOS and lysosome the upper bound is unbounded, because 4 control runs can't exclude zero line signal (P = 43%, 26%, 23%).
- **E5b, bulk RNA replicates (5 control runs):** lysosome is reliable (ICC **0.89**). Unadjusted SNCA, TH and OXPHOS are not (ICC ≈ 0–0.2).
- **E5c, noise source (day 25; 9 within-batch runs):** **85–97% of run-to-run variance is between batches** for SNCA, TH and OXPHOS, against 21% for lysosome. Within one batch, replicates are near-identical (ICC 0.87–0.99).

### 5. What the study shows
1. **The gated framework worked as intended.** Of the 4 planned phenotypes, the gates removed 3 before or during the primary test, for documented, pre-specified reasons:
   - DA yield is unstable across runs;
   - the two pathway composites were incoherent at gene level.

   Without the gates, a DA-yield × progression result would have been uninterpretable.
2. **No evidence that day-65 iPSC neuron phenotypes track PD progression in 47 donors.** For DaTscan decline, NSD-ISS stage and carrier status, meaningful associations (|r| ≥ 0.40) are **ruled out** for several phenotypes. For motor progression the result is inconclusive.
3. **The null results are what the measurement properties predict.** They are not evidence against cell-to-clinic links in general. Effects of a plausible size need roughly 3–5 times more donors.
4. **The main limiter is batch, acting through cell composition.** Run noise is batch-driven. SNCA becomes reliable only after single-cell composition adjustment. Lysosomal expression is the most robust line-intrinsic trait.

### 6. Recommendation

**For this project:**
1. **Treat the experimental phase as complete.** All registered analyses have run exactly as tagged, and the planned exploratory set (E1–E5) is done. Don't add more analyses to look for a positive result.
2. **One optional, clearly post-hoc check worth doing:** composition-adjust bulk SNCA using bulk neuronal/TH markers, to test directly whether composition explains its unreliability (section 5, point 4). It is small, it is mechanistic, and it would firm up the paper's main design message. Freeze it as `explore-v0.3.0` before running.
3. **Write the paper around two contributions:**
   - (a) a pre-registered, gated cell-to-clinic analysis, with an honest null and equivalence result;
   - (b) an empirical design guide for iPSC cell-to-clinic studies, from E4 and E5.

   Report every failed version (A0-4 v0.1.0, explore-v0.2.0), the power discrepancy, the 9 donors without DaTscan, and the prior-exposure disclosures.
4. **Release the analysis code and aggregate results** with the paper: tag `paper-v1.0` and archive it on Zenodo. No participant-level data.

**For future iPSC cell-to-clinic studies** (the practical take-home):
1. Measure cell-type composition per run (single-cell, or a validated deconvolution), and **adjust phenotypes for composition**.
2. **Spread replicate differentiations across batches.** Same-batch replicates give a false sense of reliability. Balance batches across groups.
3. **Estimate run-to-run reliability before sizing the donor cohort,** using several lines with ≥ 3 cross-batch runs each, not one control line.
4. **Prefer composition-robust molecular phenotypes** (lysosomal programme; composition-adjusted SNCA) over DA-neuron yield as donor traits.
5. **Plan for about 150–300 donors** to detect correlations of r ≈ 0.3 with clinical progression, with regular clinical follow-up of ≥ 5 years.

### 7. Open items
- **User:**
  - decide on the optional bulk SNCA composition check (section 6.2);
  - confirm the paper's direction.
- **Planning session:**
  - review this summary;
  - draft the paper outline and the figure list. Candidate figures:
    - gate flow (4 phenotypes → 1 tested);
    - equivalence forest plot;
    - design curves: N against runs per line, per phenotype;
    - the batch-variance decomposition.
- **Still open from earlier, now lower priority:**
  - asking LONI/DZNE for the published per-cell labels (would enable cell-type-specific pseudobulk, exploratory);
  - the LONI follow-up email.
- **The Aspera account expires about 2026-12-11.** Any further processed-data download must happen before then.

### 8. Disclosures (cumulative, all in the local chat only, none in this log)
1. CDI numbers printed from a README (early).
2. Public Table S1/S3 rows with IDs printed, on two occasions.
3. One aggregate A0-3 count seen before the tag.
4. One scRNA run directory name (PATNO and CDI) printed during E5.

---

## 2026-11-12 (local session): **EXPLORATORY** E5: design-guide uncertainty, bulk-RNA replicates, batch vs within-batch noise

Code: `scripts/exploratory2.py`, frozen as `explore-v0.2.0` (commit `1a52f41`) before running.
- **Bug found on the first real-data run:** bulk Ensembl IDs were mapped to symbols with the version stripped on one side only, so E5b and E5c returned all-NaN.
- **Fix:** commit `ef33d61`, tag **`explore-v0.2.1`**; strips versions on both sides. The test now uses versioned IDs. No analysis choice changed.
- **E5a** does not use the mapping and has a fixed seed; its results are identical in both runs.

### E5a: uncertainty in the design guide (E4)
- ICC drawn from its F-based sampling distribution (4 control runs); slope reliability from a donor bootstrap (200 refits per outcome).
- **Slope reliability λ**, median (10th–90th percentile): OFF-state UPDRS III **0.60 (0.49–0.69)**; putamen SBR **0.45 (0.14–0.60)**.
- **Donors needed for 80% power, true r = 0.3, OFF-state UPDRS III slope**, median (10th percentile); the 90th percentile is unbounded for all except P2:

| Phenotype | P(no line signal, ICC ≤ 0) | k = 1 run per line | k = 4 runs per line |
|---|---|---|---|
| P2 SNCA | 2% | 164 (141); 90th: 234 | 153 (134); 90th: 198 |
| P4 Lysosome | 23% | 283 (175) | 192 (146) |
| P3 OXPHOS | 26% | 352 (186) | 208 (152) |
| P1 DA yield | 43% | 1,053 (214) | 392 (163) |

### E5b: replicate reliability in bulk RNA-seq, day 65
Independent library preparation; the control line has **5 runs (4 df)**; composition is **not** adjusted.

| Bulk phenotype | ICC, 79 scRNA donor lines (95% CI) | ICC, all 91 lines | scRNA ICC (A0-2, composition-residualised) |
|---|---|---|---|
| SNCA | **−0.15** (−8.6 to 0.61) | 0.02 | **0.94** |
| TH (DA proxy) | 0.02 (−7.2 to 0.67) | 0.23 | P1: 0.33 |
| OXPHOS | 0.21 (−5.6 to 0.73) | 0.19 | 0.55 |
| **Lysosome** | **0.89 (0.12 to 0.96)** | **0.89** | 0.64 |

### E5c: where run-to-run noise comes from (bulk, day 25)
The control line has 9 runs in batch 1 and 4 in other batches.

| Phenotype | Share of run variance between batches | ICC using within-batch-1 variance | ICC using all 13 runs |
|---|---|---|---|
| SNCA | **97%** | 0.99 | 0.54 |
| TH | **94%** | 0.98 | 0.69 |
| OXPHOS | **85%** | 0.92 | 0.48 |
| Lysosome | 21% | 0.87 | 0.83 |

### What E5 adds (exploratory)
1. **Run-to-run noise is almost entirely batch-to-batch.** Within one batch, replicate differentiations are nearly identical (within-batch ICC 0.87–0.99). Across batches, 85–97% of the run variance appears for SNCA, TH and OXPHOS. Two consequences:
   - replicates within a single batch overstate reliability and don't help;
   - to average out the noise, replicate differentiations must **span batches**, or batch must be designed out (balanced across groups) and modelled.
2. **Composition adjustment is what makes SNCA reliable.** Unadjusted bulk SNCA has no line signal (ICC ≈ 0). The scRNA SNCA residualised on cell-type composition has ICC 0.94. The batch effect on SNCA runs mostly through how many neurons a run makes. Without single-cell composition, SNCA would not be a usable donor trait.
3. **The lysosomal programme is the most robust line-intrinsic signal.** It is reliable in both assays (bulk 0.89; scRNA 0.64), and the least batch-driven (21% between-batch). It is the best candidate phenotype for a future, larger cell-to-clinic study. In this study it was dropped from H1 only by A0-5, the composite-coherence gate.
4. **DA yield (TH, P1) is unreliable in both assays.** This supports the A0-2 decision.

Caveats:
- Within-line variance comes from one control line throughout.
- Day 25 is a progenitor stage.
- E5b is not composition-adjusted, by design.
- All results are exploratory.

### Disclosure
A file listing printed one scRNA run directory name (PATNO and CDI) into the local chat. Nothing went to this log.

### For the planning session
- E4 and E5 together give a concrete, data-based design recommendation for iPSC cell-to-clinic studies:
  1. single-cell composition adjustment;
  2. replicate differentiations across batches, not within;
  3. lysosomal and SNCA phenotypes in preference to DA yield;
  4. about 150–300 donors for effects of r ≈ 0.3.

  This would make a strong second contribution alongside the gated-framework result.
- A possible follow-up (exploratory, post hoc; only on request): composition-adjusting bulk SNCA using the bulk TH/neuronal markers, to test point 2 directly in the bulk data.

---

## 2026-11-12 (local session): **EXPLORATORY** analyses E1–E4 (frozen as `explore-v0.1.0` before running)

Code: `scripts/exploratory.py`, commit `2db113d`, tag `explore-v0.1.0`; synthetic tests in `tests/test_exploratory.py` (Cox with delayed entry recovers planted hazard ratios). Same data and sample as the confirmatory run (47 PD donors). **Everything below is exploratory and hypothesis-generating.**

### E1: component genes of P3 and P4 (H1 model per gene; BH within each set)
| Set | Genes | q < 0.05 | p < 0.05 (null expectation) | Share with β₃ < 0 | Median r |
|---|---|---|---|---|---|
| OXPHOS | 189 | **0** | 8 (9.5) | 75% | −0.12 |
| Lysosome | 85 | **0** | 3 (4.3) | 62% | −0.05 |

- **No single gene is associated with motor progression after correction.** The number of nominal hits matches chance.
- **OXPHOS direction is consistent:** three quarters of the genes point the same way as the composite ("higher OXPHOS, slower OFF-state progression"). The nominally strongest are core ATP-synthase and complex III subunits (ATP5PB, ATP5F1A, ATP5MC3, UQCRC1) and PDP1, all at r ≈ −0.36 to −0.40.
- **Caveat:** the genes are co-expressed, so the sign share is not 189 independent votes. A naive sign test (p ≈ 3e-12) greatly overstates the evidence. Read it as one weak, directionally consistent signal, not a finding.

### E2: time from diagnosis to milestones (Cox with delayed entry; phenotype + age at onset + genetic PD + media)
| Milestone | Donors | Events | Median follow-up | HR per SD (95% CI): P2 / P3 / P4 | Min q |
|---|---|---|---|---|---|
| H&Y ≥ 3 (OFF) | 46 | 15 | 6.7 y | 1.01 (0.61–1.68) / 0.98 (0.60–1.59) / 1.17 (0.71–1.92) | 0.95 |
| Freezing (NP2FREZ ≥ 1) | 47 | 37 | 4.7 y | 1.06 (0.75–1.51) / 1.19 (0.85–1.65) / 1.05 (0.78–1.41) | 0.95 |
| Cognitive impairment (MCI or dementia) | 46 | 26 | 7.1 y | 0.85 (0.56–1.30) / 0.78 (0.54–1.12) / 1.11 (0.78–1.58) | 0.95 |

- Falls are not available in the downloaded data.
- **No association with any milestone.** The CIs are wide: hazard ratios of about 0.6–1.7 per SD are not excluded.

### E3: manifesting (n = 19) vs non-manifesting (n = 21) carriers (OLS; adjusted for LRRK2, media, sex)
| | Difference (SD) | q | partial r [90% CI] | |
|---|---|---|---|---|
| P2 SNCA | +0.18 | 0.70 | +0.09 [−0.19, 0.35] | **equivalent** |
| P3 OXPHOS | +0.14 | 0.70 | +0.07 [−0.21, 0.33] | **equivalent** |
| P4 Lysosome | −0.14 | 0.70 | −0.07 [−0.34, 0.21] | **equivalent** |

- **Day-65 neuron phenotypes don't distinguish carriers who developed PD from those who haven't.** All three differences fall inside the ±0.40 equivalence bound.
- One donor line was dropped by A0-3, giving 21 non-manifesting carriers rather than 22.

### E4: what a cell-to-clinic study needs (uses reliabilities only, no phenotype–outcome association)
**Measured reliabilities:**
- **Phenotypes** (A0-2 ICC, one run per line): P1 0.33, P2 0.94, P3 0.55, P4 0.64. These rest on 3 df, so they are very uncertain.
- **Outcome slopes** (per donor, from the LMM without phenotype):
  - OFF-state UPDRS III: λ = 0.64 (slope SD 1.45 points/year; residual SD 7.2);
  - putamen SBR: λ = 0.54.

**Smallest true correlation this study could detect** (80% power, single test):

| | OFF-state UPDRS III | Putamen SBR |
|---|---|---|
| P1 | 0.92 | > 1 |
| P2 | 0.55 | 0.67 |
| P3 | 0.71 | 0.88 |
| P4 | 0.66 | 0.81 |

**Donors needed for 80% power** (OFF-state UPDRS III slope; k = differentiation runs per line):

| True r | P2 (k = 1 → 4) | P3 (k = 1 → 4) | P1, DA yield (k = 1 → 4) |
|---|---|---|---|
| 0.3 | 149 → 142 | 249 → 167 | 416 → 209 |
| 0.4 | 85 → 81 | 142 → 96 | 235 → 119 |

**Interpretation:**
1. **Given the measured reliabilities, the smallest detectable true correlation was 0.55 (SNCA) to 0.92 (DA yield).** *(Wording corrected 2026-11-12: 0.55 already includes attenuation. Separately, n = 47 alone limits detection to about r ≥ 0.43.)* Plausible biological effects (r 0.2–0.4) were out of reach. The null results are expected, not evidence against cell-to-clinic links.
2. **For noisy phenotypes (DA yield, OXPHOS), repeating differentiations helps more than recruiting donors.** Four runs per line roughly halves the donors needed for DA yield (416 → 209 at r = 0.3).
3. **For reliable phenotypes (SNCA), the limit is the clinical slope.** Here λ = 0.64, so more visits per donor or longer follow-up help, and extra runs do not.
4. **A realistic target is about 150–250 donors with 2–4 runs per line, for effects of r ≈ 0.3.** That is 3–5 times FOUNDIN-PD's PD arm.

**Caveats:**
- The ICCs come from one control line with 4 runs.
- λ assumes PPMI-like visit schedules.
- The Fisher-z sample sizes assume about 6 covariates.

### For the planning session
- E4 is the most useful new result: an empirical design guide for iPSC-to-clinic studies, built on measured reliabilities. It would make a strong figure: required N against runs per line, one curve per phenotype.
- Possible next steps, all exploratory; I need a request before running any:
  1. A bootstrap CI for E4, resampling the donor lines and the 4 control runs, to show how uncertain the ICC is.
  2. Using the replicate structure that FOUNDIN does have across batches for other modalities. The bulk RNA day-0/25/65 control runs may give more within-line runs, and so a better ICC.
  3. An OXPHOS replication check, if an independent iPSC cohort with clinical follow-up is accessible.

---

## 2026-11-12 (local session): **CONFIRMATORY RESULTS**: H1–H3 on P2–P4. No hypothesis supported; A0-5 drops P3 and P4 from H1

Run at `prereg-v1.0.1` (analysis code identical to `prereg-v1.0.0`), PPMI Study Data 2026-11-12 plus Curated Cut v.2026-08-10. Output `data/derived/models.json` (local only). No model needed a random-effects fallback in the primary analyses; all df are Satterthwaite. Classes: **supported** (q < 0.05), **equivalent** (90% CI of partial r inside ±0.40), **inconclusive**. Phenotypes are per SD.

### Gate A0-5: do the composites move with their genes? (H1 model, ≥ 80% of genes share the composite's sign)
| Composite | Genes | Same sign as composite | Result |
|---|---|---|---|
| P3 OXPHOS | 189 | 75.1% | **FAIL**: dropped from H1; components exploratory |
| P4 Lysosome | 85 | 62.4% | **FAIL**: dropped from H1; components exploratory |

### H1 (primary): phenotype × time on OFF-state MDS-UPDRS III (n = 47 donors, 527 visits)
| | β₃ (points/year per SD) | SE | df | p | q | partial r [90% CI] | Class |
|---|---|---|---|---|---|---|---|
| **P2 SNCA** | −0.27 | 0.28 | 30.7 | 0.33 | 0.33 | −0.18 [−0.44, 0.13] | **inconclusive** |
| P3 OXPHOS (dropped by A0-5; reported only) | −0.49 | 0.26 | 29.3 | 0.069 | n/a | −0.33 [−0.55, −0.04] | not tested |
| P4 Lysosome (dropped by A0-5; reported only) | −0.18 | 0.26 | 27.3 | 0.51 | n/a | −0.13 [−0.41, 0.19] | not tested |

After A0-2 (P1 dropped) and A0-5 (P3, P4 dropped), the confirmatory H1 family is the single test on P2.

### H2: phenotype × time on putamen SBR (n = 38, 135 scans)
| | β₃ (SBR/year per SD) | p | q | partial r [90% CI] | Class |
|---|---|---|---|---|---|
| P2 | +0.002 | 0.65 | 0.88 | +0.08 [−0.21, 0.34] | **equivalent** |
| P3 | +0.003 | 0.50 | 0.88 | +0.12 [−0.17, 0.38] | **equivalent** |
| P4 | +0.001 | 0.88 | 0.88 | +0.03 [−0.29, 0.34] | **equivalent** |

### H3: phenotype by SAA status (n = 45; 7 negative) and by NSD-ISS stage (n = 39)
| | SAA: β (SD), q, r [90% CI] | Class | NSD-ISS: β (SD per stage), q, r [90% CI] | Class |
|---|---|---|---|---|
| P2 | +0.69, 0.36, +0.25 [−0.02, 0.47] | inconclusive | −0.06, 0.90, −0.03 [−0.32, 0.26] | **equivalent** |
| P3 | −0.45, 0.36, −0.16 [−0.40, 0.11] | inconclusive | −0.04, 0.90, −0.02 [−0.31, 0.27] | **equivalent** |
| P4 | −0.40, 0.36, −0.15 [−0.40, 0.12] | **equivalent** | +0.23, 0.90, +0.14 [−0.17, 0.40] | inconclusive |

### Sensitivity analyses (registered); nothing reaches q < 0.05 anywhere
- **1 Idiopathic only:**
  - H1 n = 28; H2 n = 22. All classes inconclusive or equivalent.
  - SAA is descriptive only, because idiopathic donors have < 5 SAA-negatives.
- **2 Five-level batch:** the same pattern as the primary analysis (H1 P2 r −0.15; P3 r −0.32, q 0.25).
- **3 Alternatives to P1** (all failed A0-2, so exploratory):
  - H1: %TH+ r −0.28 and −0.32 (q 0.21); totalDopa r +0.02 (*equivalent*).
  - Nothing supported in H2 or H3.
- **4 Marker-label residualisation:** the same as the primary analysis (H1 P2 r −0.22, P3 r −0.31), so composition residualisation doesn't drive the results.
- **5 Leave one donor out, H1:** no change of sign for P2, P3 or P4. The largest shift is 1.16 SE, for P3.
- **6 ON-state with LEDD** (n = 42): all inconclusive (|r| ≤ 0.18).
- **7 First 5 years after diagnosis:**
  - H1 n = 44: |r| ≤ 0.16.
  - H2 n = 35: inconclusive.
- **8 Part II slope** (n = 47, 913 visits): r −0.09 / −0.14 / +0.21; P2 *equivalent*.
- **9 Time from PPMI baseline:** the same as the primary analysis.

### Plain-language summary
- **No pre-registered hypothesis is supported.**
- **The day-65 phenotypes don't track the rate of DaTscan decline:** all three are *equivalent* to a negligible effect.
- **NSD-ISS stage:** SNCA and OXPHOS are *equivalent* (negligible association).
- **Motor progression (H1), SAA status and the remaining tests:** the data can neither show an effect nor rule one out (*inconclusive*), as the power simulation predicted.
- **The only consistent pattern is exploratory:** higher OXPHOS score goes with slower OFF-state motor progression (r ≈ −0.3) in the primary, batch, marker-residualised and baseline-time models. It is not significant, and the composite failed A0-5.

### Not yet done (exploratory, to be labelled as such)
- Gene-level P3/P4 components (A0-5 sends them there).
- Milestone time-to-event analyses.
- Any batch-adjusted P1.
- The archived DaTscan file for the 9 donors without Xing scans.

None of these will be run without a request here.

### For the planning session
1. Please review how to present these results: the gates removed 3 of 4 phenotypes before or during H1, and the 10 confirmatory tests (H1 1, H2 3, H3-SAA 3, H3-NSD-ISS 3) gave 0 supported, 6 equivalent and 4 inconclusive.
2. Decide which exploratory analyses (if any) to run, and in what order.
3. Paper outline: the gated-framework results (A0-2, A0-4, A0-5 decisions), the equivalence results for DaTscan, and the power discrepancy.

---

## 2026-11-12 (local session): PPMI data checked; name-only fixes (`prereg-v1.0.1`); **A0-1 PASSES**; donor counts per analysis

### Data
- Unzipped to `data/ppmi/2026-11-12/`: 8 study-data CSVs (dated 12Nov2026), the Curated Data Cut **v.2026-08-10**, the data dictionary and code list, and 2 method PDFs. Only file names, column headers, code lists and aggregate value formats were inspected.
- **Matched the config unchanged:**
  - Demographics (SEX 0 = female, per the code list);
  - PD Diagnosis History;
  - MDS-UPDRS Part III (PDTRTMNT 0/1; PDSTATE OFF/ON);
  - Part II (NP2PTOT);
  - SAA (all rows are CSF; "Inconclusive" is treated as missing, as coded).
  - All dates are MM/YYYY (set to day 15), except SAA RUNDATE (YYYY-MM-DD).
- **Name-only fixes** (`config/ppmi_columns.json`, commit `a97b8f7`, tag **`prereg-v1.0.1`**):
  - DaTscan: the file is "Xing Core Lab - Quant SBR"; date `DATSCAN_DATE`; putamen `PUTAMEN_L_REF_CWM` / `PUTAMEN_R_REF_CWM` (reference: cerebral white matter). 31 TRODAT scans exist in the file; **none** belong to FOUNDIN-PD donors, so no filter is needed.
  - LEDD: the study-data file is a per-drug medication log, not per visit. As the checklist anticipated, **per-visit LEDD is taken from the Curated Data Cut** (`LEDD`, matched by EVENT_ID).
  - NSD-ISS: from the Curated Data Cut (`NSD_STAGE`; "Not NSD" = not stageable).
  - New helper `scripts/prep_curated.py` converts the 5 needed curated columns from xlsx to CSV (format only), because the frozen loader reads CSV.
  - **`ppmi_load.py`, `models.py`, `lmm.py`, `phenotypes.py` and `gates_prelink.py` are byte-identical to `prereg-v1.0.0`** (checked with git diff).

### A0-1 and the sample (computed before any model was fitted)
| | n |
|---|---|
| PD donors with day-65 scRNA | 47 (28 idiopathic, 19 genetic) |
| After A0-3 (the excluded line is not a PD donor) | 47 |
| Missing any covariate | 0 |
| **A0-1:** ≥ 3 OFF visits spanning ≥ 2 years | **44 → PASS** (threshold 30) |
| H1: ≥ 3 OFF visits, complete covariates | 47 (median 12 OFF exams per donor; median span 8 years) |
| Sensitivity 6: donors with ON-state exams | 46 |
| **H2:** ≥ 1 DaTscan / ≥ 2 DaTscans | **38 / 38** |
| H3 SAA: negative / positive / missing | **7 / 38 / 2** (≥ 5 negatives, so it is tested, not descriptive) |
| H3 NSD-ISS: stageable | 39 (stages 3: 11, 4: 22, 5: 6) |
| Sensitivity 8: ≥ 3 Part II visits | 47 |

### Flag for the planning session: DaTscan coverage
- 9 of 47 PD donors have **no scan** in the Xing Core Lab file. The archived "DaTScan SBR Analysis", an older processing, was not downloaded. **H2 runs as registered on the 38**; every one of them has ≥ 2 scans, so the effective n for the slope is also 38.
- Any use of the archived file would be a post-hoc, exploratory addition, because it mixes processing pipelines. Not done.

### Next
H1–H3, A0-5 and sensitivity analyses 1–9 on P2–P4 (P1 was dropped by A0-2) are running now at `prereg-v1.0.1`.

---

## 2026-11-12 (planning session): PPMI clinical data cut downloaded

The user logged in to LONI IDA; the planning session selected and downloaded the files via the user's browser. Download date: **2026-11-12** (after OSF approval and the Zenodo archive, as registered).

- **Location:** `D:\ST-DT4\data\ppmi\` (inside the git-ignored `data/` folder; the original zip is kept there). Note: a single folder, not the dated subfolder suggested in the checklist. Record the date from this entry.
- **Source:** LONI IDA → PPMI → Download → Study Data (one combined download).
- **Data-cut version:** Study Data files current as of 2026-11-12; **Curated data cut v.2026-08-10**. The latest quarterly freeze listed on IDA is 2026Q2 (2026-07-01); it was not downloaded.

| # | Checklist item | File on IDA |
|---|---|---|
| 1 | Demographics | Demographics |
| 2 | Diagnosis date | PD Diagnosis History |
| 3 | H1 outcome | MDS-UPDRS Part III Treatment Determination and Part III: Motor Examination |
| 4 | Sensitivity 8 | MDS-UPDRS Part II Patient Questionnaire: Motor Aspects of Experiences of Daily Living (M-EDL) |
| 5 | H2 outcome | Xing Core Lab - Quant SBR, plus DaTSCAN SPECT Quantitative Analysis Results README v1.0 |
| 6 | Sensitivity 6 | LEDD Concomitant Medication Log, plus Method for Calculating Levodopa Equivalent Dose in PPMI Data |
| 7 | H3 outcome | SAA Biospecimen Analysis Results |
| 8 | H3 outcome (NSD-ISS) | Curated data cut (v.2026-08-10) |
| 9 | Name checks | Data Dictionary; Code List |

- **Not downloaded:** the archived "DaTScan SBR Analysis" (older processing). If Xing Core Lab Quant SBR lacks early visits, report it here before deciding anything.
- Nothing was opened or viewed during the download.

### For the local session
1. Unzip if not already done, list file names and **column headers only**, compare with `config/ppmi_columns.json`, and log any name-only fixes here.
2. Report anything missing or structured differently **before** running models (e.g. DaTscan file layout, NSD-ISS location in the curated cut).
3. Then A0-1, the donor counts per analysis (H1; H2 with ≥ 2 scans; H3 SAA−/SAA+), and H1–H3 exactly as tagged in `prereg-v1.0.0`.

---

## 2026-11-12 (planning session): review of `prereg-v1.0.0`, power check and the A0-2 result

Good work: frozen and archived before data, gates applied exactly as registered, and the deviations and disclosures logged openly.

### Fixed
- Claude Doc: the "Coding and model details" bullet now says time = years since PD diagnosis (baseline = sensitivity 9). The OSF text already said diagnosis (checked on the submitted Review page). The filing checklist is ticked and the version DOI is added.

### A0-2: P1 fails. How to frame it in the paper
1. **Present it as the framework working, not as a setback.** The gate stopped the planned primary phenotype from going to clinical testing because its run-to-run variation (across batches and media) was as large as the variation between donors. Without the gate, a P1 × progression result would have been uninterpretable. This is the central demonstration for the "gated framework" contribution, so give it a figure.
2. **Wording: "not shown to be stable", not "shown to be unstable".** ICC 0.33 with 95% CI −8.4 to 0.80 rests on 3 df. The registered rule uses the point estimate, so the decision stands. But state the CI and the 3-df limit in the same sentence.
3. **Connect the evidence.** Run-to-run instability (A0-2), media and batch sensitivity of %TH+, and day-0 predictability of the published fraction (A0-4) together point to the same picture: DA yield has a line-intrinsic component, but in this protocol it is strongly modulated by differentiation run and media. That's a useful message for iPSC disease-modelling studies that treat yield as a donor trait.
4. **No rescue in the confirmatory analysis.** Don't add a batch-adjusted P1 to H1–H3. If you want one, it is a post-hoc exploratory analysis, labelled as such and reported after the registered results.
5. **Knock-on effect to note:** P2–P4 are still residualised on logit P1 as registered, but P1 is noisy, so some composition signal may remain in P2–P4. Sensitivity 4 (marker-label residualisation; ICCs 0.62/0.51/0.80) is the check; discuss it alongside the main results.

### Power discrepancy
- **Agree it's real; report it, don't fix it.** The analytic figures used n − 3 − 6 df and ignored the extra slope parameters and slope-estimation noise. The simulation (equivalence 58% at n = 47, 16% at n = 30) is the honest figure.
- In the paper: put both side by side, explain the cause in one sentence, and say before results are known that most nulls will be classed inconclusive. Note that the filed analytic figures were computed before the simulation existed.

### The 6 decisions made in code
1–3, 5: fine; standard choices, frozen before data. Two comments:
- **(4) H2 with donors having ≥ 1 DaTscan.** Acceptable in a mixed model (single-scan donors inform the intercept, not the slope), but report how many H2 donors have ≥ 2 scans. That number is the effective sample for the slope.
- **(5) Age at onset from symptom-onset date.** Fine and conventional. Name it explicitly in the methods ("age at symptom onset").
- (6) Milestones not coded: fine, exploratory.

### Disclosures
- Both are minor (public Table S3 rows; one aggregate count). Keep them in the prior-exposure section of the paper as logged.

### For the local session
- Proceed with the clinical-data steps in your checklist once the user downloads the data cut. Before running H1–H3, log the A0-1 result and the number of PD donors per analysis (H1, H2 with ≥ 2 scans, H3 SAA−/SAA+).

---

## 2026-11-12 (local session): PPMI clinical data download checklist (for the user)

The OSF registration is approved and the code is archived (10.5281/zenodo.24156791), so the clinical data cut **may now be downloaded**. LONI IDA needs the user's login, so the user does the download; the local session does everything after that.

### Where to save
- On the analysis machine only: **`D:\ST-DT4\data\ppmi6-11-12\`**, or whatever the actual download date is.
- `data/` is git-ignored. Never put these files in cloud storage, email or the repo.
- Unzip into that folder. Keep the original zip files there too, for provenance.

### What to download
From LONI IDA → PPMI → **Download → Study Data**. Download the **whole file** for each item; the loader keeps only the FOUNDIN-PD PD donors and drops everyone else when it reads each file.

| # | Needed for | PPMI study-data file (usual name) | Must contain |
|---|---|---|---|
| 1 | Sex, age at onset | **Demographics** | PATNO, SEX, BIRTHDT |
| 2 | Time since diagnosis, onset, duration | **PD Diagnosis History** | PATNO, PDDXDT (diagnosis date), SXDT (symptom onset date) |
| 3 | **H1 primary outcome**, A0-1, sensitivity 6 | **MDS-UPDRS Part III** (Motor Examination) | PATNO, EVENT_ID, INFODT, NP3TOT, PDTRTMNT, PDSTATE |
| 4 | Sensitivity 8 | **MDS-UPDRS Part II: Patient Questionnaire** | PATNO, EVENT_ID, INFODT, NP2PTOT |
| 5 | **H2 outcome** | **DaTScan Analysis** (SBR results) | PATNO, INFODT or scan date, DATSCAN_PUTAMEN_L, DATSCAN_PUTAMEN_R |
| 6 | Sensitivity 6 (LEDD) | **LEDD Concomitant Medication Log**, or a per-visit LEDD table | PATNO, EVENT_ID, LEDD |
| 7 | **H3 outcome** | **α-synuclein SAA results** (CSF; e.g. "SAA Biospecimen Analysis Results") | PATNO, assay or collection date, result (Positive/Negative) |
| 8 | **H3 outcome** | **NSD-ISS staging**, often in the **Curated Data Cut** (`PPMI_Curated_Data_Cut_*.xlsx`) | PATNO, visit date, NSD-ISS stage |
| 9 | Name checks | **Data dictionaries / codebooks** for the files above, and for the Curated Data Cut if used | |

**Also useful:** the **PPMI Curated Data Cut** spreadsheet and its dictionary. It often has per-visit LEDD and NSD-ISS already derived, and is the easiest source for items 6 and 8.

**Do not download** imaging (raw DaTscan/MRI), genetics, biospecimen inventories beyond the SAA result, or any other cohort's data. They are not needed and not in scope.

### What to record (reply here or tell the local session)
- The **download date**.
- The **data-cut version or label** shown on IDA (e.g. the date stamp in the file names, or the Curated Data Cut version).
- Anything that wasn't available.

### What the local session does next
1. Lists the file names and **column headers only**, and compares them with `config/ppmi_columns.json`. File and column-name fixes are logged here as name-only changes. **The analysis logic in `scripts/ppmi_load.py` and `scripts/models.py` stays exactly as tagged in `prereg-v1.0.0`.**
2. If a needed item is missing or structured differently (e.g. LEDD only as a medication log), it reports here **before** running any model, and asks how to proceed.
3. Runs A0-1, the linkage gate (≥ 30 PD donors with ≥ 3 OFF-state UPDRS III visits spanning ≥ 2 years). If it fails, it stops with a data audit only.
4. Runs H1–H3 on P2–P4 (P1 was dropped by A0-2), A0-5, and sensitivity analyses 1–9. Then it logs aggregate results here: no PATNOs, no per-donor values.

---

## 2026-11-12 (local session): OSF approved; pre-linkage gates run. **A0-2 FAILS for P1**; P2–P4 pass

- **OSF registration osf.io/s6w9c approved**, as reported by the user on 2026-11-12. The code archive (DOI 10.5281/zenodo.24156791) was finished before this.
- Ran on real data at tag `prereg-v1.0.0`; the working tree was confirmed identical to the tag for `scripts/`, `config/` and `tests/`. Commands: `phenotypes.py --s3 <Table S3>`, then `gates_prelink.py`. **No clinical data yet.**

### Phenotype construction
- 79 donor lines and 4 control runs.
- 9,340 genes pass mean CPM ≥ 10.
- Genes used: P3 OXPHOS 189 of 200; P4 Lysosome 85 of 121.
- SNCA (P2) passes the filter.

### Gate A0-2: technical stability (ICC ≥ 0.5; within-line variance from the control line, 3 df)
| Phenotype | ICC | 95% CI | Result |
|---|---|---|---|
| **P1** DA-neuron fraction (Table S3) | **0.33** | −8.4 to 0.80 | **FAIL → dropped** |
| P2 SNCA | 0.94 | 0.12 to 0.98 | pass |
| P3 OXPHOS | 0.55 | −5.2 to 0.86 | pass |
| P4 Lysosome | 0.64 | −4.0 to 0.89 | pass |
| *Report only:* %TH+ ≥1 UMI / ≥2 UMI / totalDopa | 0.30 / 0.19 / 0.41 | all lower bounds < −7 | below 0.5 |
| *Report only:* P2/P3/P4 residualised on our marker labels | 0.62 / 0.51 / 0.80 | | ≥ 0.5 |

### Gate A0-3: cell sufficiency (≥ 1,000 QC cells)
1 of 79 donor lines is excluded. Whether it is a PD donor is counted when the models run.

### What this means, per the registration
- **P1 is dropped from H1–H3** ("Drop that phenotype; report it"). The H1, H2 and H3 families each become **3 tests (P2, P3, P4)**, with BH over 3.
- `models.py` already reads `A0_2_phenotypes_retained` from the gate output, so nothing in the code changes.
- Sensitivity 3 (alternatives to P1) is coded to run regardless. Its P1 replacements also fall below 0.5 here, so those results will be reported as **exploratory**, with this caveat.
- **Interpretation, for the paper; nothing here was decided on data:**
  - All the DA-composition measures (P1, %TH+, totalDopa) vary about as much between differentiation runs of the *same* line as between lines.
  - The control line's 4 runs span batches 1, 2, 3 and 5, so E8 and E6 media are both included. Its within-line variance therefore includes the batch and media effects that A0-4 and the batch analysis had already flagged for DA yield.
  - The pathway phenotypes (P2–P4) are stable.
  - With 3 df the CIs are very wide, as the registration anticipated. The rule uses the point estimate, so the decision stands.
- This is a substantive negative result in its own right: DA-neuron yield in this protocol is, at this sample size, not reliable enough as a donor-level trait.

### Next: PPMI clinical data cut
The user downloads it from LONI IDA (an authenticated portal; the local session has no access) into `data/ppmi/<download-date>/`. Then the local session records the download date and data-cut version here, checks the file and column names (name-only fixes logged), and runs A0-1 → H1–H3.

### For the planning session
- Please review the P1 failure. The registered rule applies mechanically; the question is how to frame it in the paper. P1 was the primary phenotype, and A0-4 and A0-6 supported its construction, but not its run-to-run stability.

---

## 2026-11-12 (local session): Zenodo archive of `prereg-v1.0.0`, DOI minted

- The user switched on the Zenodo GitHub integration and published the GitHub release [`prereg-v1.0.0`](https://github.com/xaviwho/cell-to-clinic/releases/tag/prereg-v1.0.0) (not a pre-release), on 2026-11-12 at 09:18 UTC.
- Zenodo archive: **version DOI [10.5281/zenodo.24156791](https://doi.org/10.5281/zenodo.24156791)** (this exact code). Concept DOI covering all versions: [10.5281/zenodo.24156790](https://doi.org/10.5281/zenodo.24156790).
- This completes the registration's commitment to a "tagged, Zenodo-archived release of the code before it is applied to real data". Cite the version DOI in the paper.

### Still to happen before any real-data run
1. OSF shows osf.io/s6w9c as **approved**. The user approves via the OSF email; please record the approval date here.
2. Then the local session downloads the PPMI clinical data cut, records the download date and data-cut version here, and runs, at `prereg-v1.0.0`: phenotypes → A0-2/A0-3 → A0-1 → H1–H3.

### For the planning session
- The Claude Doc's checklist items "Make a Zenodo archive…" and "Write the phenotype and gate scripts…" can be ticked, and the version DOI added to Part 5 of the doc. The OSF record itself doesn't need it.

---

## 2026-11-12 (local session): pre-registration code release `prereg-v1.0.0`, synthetic-tested only

Done the planning session's tasks 1–2. **Nothing has been run on real FOUNDIN-PD or clinical data.** Code frozen as tag [`prereg-v1.0.0`](https://github.com/xaviwho/cell-to-clinic/releases/tag/prereg-v1.0.0) (commit `c61fb9b`).

### What is in the release
| File | Contents |
|---|---|
| `config/genesets/` | HALLMARK_OXIDATIVE_PHOSPHORYLATION (200 genes) and KEGG_LYSOSOME (121 genes), MSigDB v2023.2.Hs (KEGG set from the `c2.cp.kegg_legacy` file) |
| `config/phenotypes.json` | Table S3 columns: P1 = "Dopaminergic Neurons"; progenitor sum = Neuroepithelial-like Cells + Proliferating Floor Plate Progenitors + Early neuron Progenitor + Late neuron Progenitor (ependymal and both DA clusters excluded); thresholds |
| `scripts/phenotypes.py` | P1–P4 as in Part 3, plus the registered variants (sensitivity 3: %TH+ at ≥1 and ≥2 UMI, totalDopa; sensitivity 4: marker-label residualisation) |
| `scripts/gates_prelink.py` | A0-2 ICC = 1 − S²(control runs)/S²(79 donor lines), with F-based 95% CI; A0-3 ≥ 1,000 QC cells |
| `scripts/lmm.py` | REML mixed model with Satterthwaite df (lmerTest method), registered fallbacks, OLS, partial r, TOST, BH; numpy/scipy only |
| `scripts/ppmi_load.py` + `config/ppmi_columns.json` | Tidy clinical tables. Logic fixed in code (OFF = untreated or PDSTATE OFF; time = years since diagnosis; earliest SAA; latest NSD-ISS; mean L/R putamen). File and column names live in the config and are checked on download |
| `scripts/models.py` | A0-1, H1, H2, H3 (SAA, NSD-ISS), A0-5, sensitivity analyses 1–9 |
| `scripts/power_sim.py`, `results/power_sim.json` | Simulation check of the registered power figures (synthetic) |
| `tests/` | `test_lmm.py`, `test_pipeline_synthetic.py` |

### Validation (synthetic)
- **LMM:**
  - fixed effects and variance components match statsmodels MixedLM REML to ≤ 5e-4; SEs agree within 0.7%, the expected difference between the two SE definitions;
  - Satterthwaite df on a balanced design is 39.000, exactly the theoretical value (47 − 8);
  - type I error 0.047 over 300 null simulations;
  - with no simulated slope variance, the fallback rule switches as specified.
- **End-to-end** (fake scRNA outputs, Tables S1/S3, and PPMI-format CSVs; P1 × time effect planted):
  - P1 is recovered as *supported* in H1 and in sensitivity analyses 1, 2, 3, 6, 7 and 9;
  - P2–P4 come out *equivalent* or *inconclusive*;
  - A0-3 drops the planted low-cell line, and the ≥ 3-OFF-visit rule drops the planted 2-visit donor;
  - A0-5 runs over all 200 OXPHOS and 121 lysosome genes;
  - the loader keeps OFF/ON, SAA, NSD-ISS and LEDD-by-visit correct and drops non-FOUNDIN participants.

### Power simulation (registered H1 analysis, 300 reps per cell)
Assumptions (made before linkage): 3–9 annual OFF visits starting 0–3 years after diagnosis; slope 2 ± 1.5 points/year; residual SD 5.

| | n = 47 | n = 30 |
|---|---|---|
| \|partial r\| at 80% power, single test | **0.45** (registered: 0.43) | > 0.52, power 74% at r = 0.52 (registered: ≈ 0.54) |
| same, BH over 4, one true effect | **0.53** (registered: ≈ 0.49) | > 0.52, power 53% at r = 0.52 |
| P(equivalent) when the true effect is 0 | **58%** (registered: ≈ 83%) | **16%** (registered: ≈ 62%) |
| Type I error, null phenotypes | 0.04–0.06 | 0.04–0.07 |
| Median Satterthwaite df | ≈ 33 | ≈ 19 |

**The registered equivalence-power figures are optimistic.** The analytic figures appear to assume about n − 3 df. The registered model spends 8 df on slope terms (time, phenotype × time, 6 covariates × time), and REML/Satterthwaite lowers the effective df further. Detection power is close to the registered figures. The registration is filed, so this goes into the paper as a reported discrepancy, not a change of plan. Expect most null results at n = 47 to be classed *inconclusive* rather than *equivalent*, and nearly all at n = 30.

### Decisions made in code (not spelled out in the registration; frozen in the tag before any data)
1. **"Fails to converge"** = optimiser failure, OR a singular fit (random-effect SD < 1e-4 × residual SD, or |correlation| > 0.999, as in lme4 `isSingular`), OR a non-positive-definite deviance Hessian.
2. **TOST:** the 90% CI of β₃ (t distribution, Satterthwaite df) is mapped to partial r by r = t/√(t² + df).
3. **A0-2 ICC:**
   - P1 within-line variance uses the control line's 4 Table S3 rows. The note "the control line has no usable Table S3 linkage" holds only for the mapping of row to run, and a variance doesn't depend on that mapping.
   - For P2–P4, the control runs are residualised with the line's mean Table S3 composition. This can only lower the ICC.
4. **Sample rules by analysis:**
   - H2 uses all donors with ≥ 1 DaTscan.
   - Sensitivity 8 (Part II) uses the same ≥ 3-visit rule as H1.
   - In sensitivity 1 (idiopathic only), the genetic-PD covariate is dropped because it is constant.
   - In sensitivity 6, LEDD is attached to the UPDRS III exam with the same EVENT_ID.
5. **Coding:**
   - NSD-ISS numeric = integer stage, with subtype letters dropped (2B → 2).
   - Age at onset uses the symptom-onset date.
   - Month/year dates are set to day 15.
6. **Milestone time-to-event analyses** (exploratory) are not coded yet. They stay labelled exploratory.

### Disclosures (this session, before the tag)
- A structure check printed **2 rows of public Table S3** (IDs and proportions) into the local chat. It is public data, but it breaks our own masking rule. Nothing went to this log.
- A column check displayed **one aggregate count:** 1 of the 83 runs has fewer than 1,000 QC cells (minimum 517). That is effectively the overall A0-3 count. The run was not identified and not related to diagnosis or any donor characteristic.

### For the planning session
- **Inconsistency in the Claude Doc:** the "Coding and model details" bullet still says *"time in years since PPMI baseline, not centred"*, while the filed change makes time = years since diagnosis (baseline is sensitivity 9). The code follows the filed version. Please check that the OSF text matches, and fix the doc.
- Please review the power discrepancy and the 6 coded decisions above.

### For the user
1. **Zenodo:** log in to zenodo.org with GitHub and switch the toggle **on** for `xaviwho/cell-to-clinic`, under Account → GitHub. Then create a GitHub release from tag `prereg-v1.0.0`; Zenodo archives it and mints the DOI. I can create the release with `gh` once the toggle is on. A release made before that gets no DOI.
2. Once OSF shows the registration as **approved** and the DOI is logged, the local session will download the PPMI data cut and run: phenotypes → A0-2/A0-3 → A0-1 → H1–H3, exactly as tagged.

---

## 2026-11-12 (planning session): OSF pre-registration SUBMITTED

- **Registered 2026-11-12** on OSF: [osf.io/s6w9c](https://osf.io/s6w9c), template "OSF Preregistration", **embargoed until 2027-11-12**. Status at submission: pending admin approval (the user approves via the OSF email). The OSF record is authoritative; the Claude Doc mirrors it.
- Foreknowledge level certified: *"Authors have observed the data, but have not performed the proposed analyses."*
- **Changes made at filing (differ from the drafted entry below):**
  - **Time axis = years since PD diagnosis** at each visit (was years since PPMI baseline). Disease duration at PPMI baseline stays as a covariate.
  - New sensitivity analyses: (7) H1/H2 restricted to the first 5 years after diagnosis; (8) MDS-UPDRS Part II slope as outcome; (9) time from PPMI baseline instead of diagnosis.
  - Time to clinical milestones (falls, freezing, cognitive impairment, H&Y 3) added as exploratory.
  - Coding fixed: sex female = 1; genetic PD = 1; media E6 = 1; SAA positive = 1; NSD-ISS numeric; continuous covariates z-scored; outcomes untransformed. Convergence fallback: drop random-effect correlation, then random intercept only. Missing visits handled by likelihood (MAR); no imputation.
  - The Zenodo DOI is **not** in the registration. It commits to a tagged, Zenodo-archived release of the code **before it is applied to real data**, with the DOI reported in the paper.

### For the local session
1. The pre-registration tasks below still apply (phenotypes P1–P4, gates A0-2/A0-3/A0-5, H1–H3 models), now with **time = years since diagnosis** and sensitivity analyses 7–9. Test on synthetic data only.
2. Tag the release and archive it on Zenodo **before** running anything on real data; log the tag and DOI here.
3. **The PPMI clinical data cut may be downloaded once OSF shows the registration as approved.** Record the download date and data-cut version here. Then run A0-2 and A0-3 (no clinical data), then A0-1, then H1–H3 exactly as registered.

---

## 2026-11-12 (planning session): pre-registration drafted

The full registration text is in the Claude Doc [Cell-to-Clinic: OSF Pre-registration](https://claude.ai/artifact/ED1Q5tAFzQ3H3WTupWX3G6), laid out to match the OSF Secondary Data Preregistration template. It supersedes the earlier working draft.

### Decided by the user
- The Tier 3 commitment is to register **before linking to clinical data**. Processing so far is disclosed (Part 4), not a deviation.
- Primary group: all PD donors with day-65 scRNA, n = 47; idiopathic-only as sensitivity.
- Phenotypes: P1 Table S3 "Dopaminergic Neurons" fraction (logit); P2 SNCA; P3 OXPHOS score (MSigDB Hallmark); P4 lysosomal score (MSigDB KEGG Lysosome). P2–P4 from all-cell pseudobulk, residualised on composition across the 79 donor lines.
- Carrier comparison and non-scRNA modalities: exploratory.

### Proposed, awaiting the user
- **A0-3 changed** to ≥ 1,000 QC-passing cells per line (was ≥ 200 DA cells). P1 is now a composition measure, so excluding low-DA lines would truncate it.
- Six covariates: age at onset, sex, genetic PD, PRS without GBA/LRRK2, media (E6/E8), disease duration at baseline.

### For the local session (do not run on real data until the registration is filed)
1. Save gene lists for P3/P4 (MSigDB v2023.2.Hs) and the Table S3 progenitor-cluster column names under `config/`.
2. Write `scripts/phenotypes.py` (P1–P4 exactly as in Part 3 of the doc), `scripts/gates_prelink.py` (A0-2 ICC from the control line's 4 runs vs 79 donor lines; A0-3) and the H1–H3 model code with A0-5. Test on synthetic data only.
3. Tag the release (e.g. `prereg-v1.0.0`) and archive it on Zenodo; report the DOI here.
4. The PPMI clinical data cut is downloaded only after the user records the OSF registration date here.

---

## 2026-11-12 (planning session): review of A0-4 v0.1.0 / v0.2.0 and the two decisions

Paper title noted (`paper/main.tex`): *Cell-to-Clinic: A Gated Framework Linking iPSC Neuron Phenotypes to Parkinson's Progression.* Handling the failed gate, choosing route (a), freezing before fitting and reporting both versions is exactly how this should go.

### Decision 1: partial A0-4 → **agree: treat as met for pipeline validity**, with conditions
- A0-4's purpose is "does the pipeline recover a known line-intrinsic signal?" The evidence that it does is the **gene-level replication**: 5/5 published genes selected, all in the published direction at FDR < 1%. That, not apparent AUC (uninformative, as your null run showed), is the basis. Say so in the registration.
- Be explicit that this is a **second attempt chosen after a failure**, and that the "partial = met" reading is being **decided now, before any clinical linkage**. Record it as a rule change with its date, and keep v0.1.0's fail in the gates table.
- Report CV AUC 0.72 vs published 0.85 plainly. Control-line leakage is a plausible explanation, but **don't test it** to rescue the number; mention it only as an untested possibility.
- Fix the gates table: add the "partial" definition (AUC ≥ 0.65, p < 0.05) and what it allows. Proposed wording: *"Partial: proceed, with the shortfall reported as a limitation."*

### Decision 2: primary DA phenotype → **agree: switch to Table S3 DA-neuron fraction**
- The evidence is good and was gathered without clinical data: the published DA fraction is predictable from day 0 and not batch-driven, whereas %TH+ is batch-sensitive and not predictable from day 0. A phenotype dominated by batch would be a poor cell-to-clinic predictor.
- **Pick one column as primary now:** Table S3 **"Dopaminergic Neurons"**, the same outcome used in A0-4, so the gate and the phenotype match. `totalDopa` (immature + mature) becomes a sensitivity analysis, and %TH+ a secondary phenotype and the A0-6 reproduction gate.
- Declare the reason for the switch in the registration (batch sensitivity and day-0 predictability, both measured without clinical data).
- Knock-on effects to update:
  - **A0-3 (≥ 200 DA cells per line)** needs a DA-cell count. Table S3 gives proportions; multiply by our per-run QC cell count, and say so.
  - **DA-neuron fraction as a covariate** in the outcome model now means the Table S3 value.
  - The control line has no usable Table S3 linkage. That matters little, since it's not a PD donor.

### Registration wording to fix (for the user)
- The draft says the registration will be filed *"before any FOUNDIN-PD data are downloaded."* That is no longer true. Check the exact commitment in the approved Tier 3 form:
  - if it promised registration before download/examination, declare this as a **deviation** with dates and a summary of what was done (this log is the record), and consider a short note to MJFF;
  - if it only covers linking to clinical outcomes, just correct the draft's wording.
- Either way: **no new phenotype-level analyses beyond the gates and the power simulation until the registration is filed.**

### For the local session
1. Apply decisions 1 and 2 to the registration text above (gates table partial rule, phenotype source, A0-3 cell count method).
2. Next: power simulation with n = 47 and the Table S3 DA fraction as primary, then OSF filing.

---

## 2026-11-12 (local session): A0-4 v0.2.0 (exact replication of the published classifier): PARTIAL; published analysis reproduced

The user chose route (a). The published method was taken **verbatim** from the STAR Methods (Europe PMC full-text XML; PMC itself served a CAPTCHA), implemented, tested on synthetic data and frozen as tag [`day0-gate-v0.2.0`](https://github.com/xaviwho/cell-to-clinic/releases/tag/day0-gate-v0.2.0) (commit `ee99751`) **before** fitting. Script: `scripts/day0_gate_v2.py`. The failed v0.1.0 result stands and is reported alongside.

### What differed from v0.1.0 (and is now replicated)
| | v0.1.0 (ours) | Published, used in v0.2.0 |
|---|---|---|
| Outcome | our %TH+, median split | Table S3 **"Dopaminergic Neurons"** cluster %, **≥ 15.7% (Q25)**: high/low |
| Features | VST, top 1,000 by variance | DESeq2 **normalised counts**, rows with mean ≥ 50 (16,415 kept of 94,946) |
| Model | L2 logistic, C tuned | binomial **elastic net, alpha 0.5, lambda 0.22** (glmnet objective, standardised x) |
| CV | 20 × 5-fold, pooled AUC | **100 × 5-fold**, fold-averaged AUC (as cv.glmnet) |

- **Check that the outcome is right:** Q25 of all 83 Table S3 rows is 15.66%, which reproduces the published 15.7% cut and the **62 high / 21 low** split exactly.
- **Deviation:** 79 donor lines rather than 83 rows (57 high / 22 low). The control line is excluded because Table S3 can't link its 4 rows to its batch-specific day-0 samples.
- scikit-learn is blocked, so elastic net was implemented in numpy: a working-set FISTA solver, checked against the KKT optimality conditions (max violation 2e-7 on the real fit).

### Result
| | Ours (v0.2.0) | Published |
|---|---|---|
| Apparent (in-sample) AUC | **0.935** | 0.93 |
| Non-zero genes | 14 | 10 |
| Published 5 genes among selected | **5 / 5** | n/a |
| Repeated 5-fold CV AUC | **0.72** (SD 0.05); permutation p = 0.02 (null mean 0.51); pooled out-of-fold 0.67 | 0.85 (SD 0.03) |

The five published genes (Spearman ρ of day-0 normalised counts vs DA %, n = 79) all replicate in the published direction with FDR < 1%:

| Gene | ρ | FDR | Published sign |
|---|---|---|---|
| HNRNPH3 | +0.41 | 4e-4 | + |
| SRSF5 | +0.35 | 2e-3 | + |
| HSD17B6 | +0.40 | 4e-4 | + |
| ZSWIM8 | −0.44 | 2e-4 | − |
| ARSA | −0.34 | 2e-3 | − |

**Verdict under the frozen rule: partial** (AUC ≥ 0.65 and p < 0.05, but < 0.75).

**Reading:** the pipeline **reproduces the published analysis**: same in-sample AUC, the same key genes and the same gene-level correlations. The shortfall is only in cross-validated AUC (0.72 vs 0.85). The synthetic null run gave an in-sample AUC of 0.99 with no signal, so the published 0.93 is not informative by itself; the CV figure is the real benchmark. Possible (untested) reasons for the CV gap:
- the paper's 83 rows include the control line's 4 near-identical runs, which can leak across folds;
- cv.glmnet's AUC details;
- n of 79 vs 83.

### Secondary results (pre-specified, not gating)
- **Batch-only baseline** for the published outcome: CV AUC **0.53**. Unlike %TH+ (v0.1.0: batch-only 0.62), DA-cluster fraction is not batch-driven. This is consistent with the paper's "not related to batch".
- **Same model, outcome = our %TH+ split at Q25: CV AUC 0.35.** Day-0 expression does **not** predict %TH+. AUC < 0.5 in small-n CV with no signal is a known artefact, so read it as "no signal", not as a reversed effect.
- **Implication:** %TH+ (TH ≥1 UMI, any cell) and published DA-neuron cluster fraction are **different constructs**. %TH+ is batch-sensitive and not predictable from day 0. The DA-cluster fraction is predictable and not batch-driven. This supports the planning decision to use **Table S3 proportions as primary cell-type phenotypes**, and argues for treating %TH+ as secondary.

### For the user / planning session: decisions
1. **What a partial A0-4 means.** The draft gates table defines only pass (AUC ≥ 0.75) and fail ("fix pipeline before any clinical test"); partial is undefined.
   - Proposal: treat A0-4 as **met for pipeline validity**, because the published analysis is reproduced: apparent AUC, all 5 genes, gene correlations, p = 0.02.
   - Report CV AUC 0.72 vs 0.85 honestly, with v0.1.0's failure alongside.
   - Fix this interpretation in the registration before any clinical step.
2. **Primary DA phenotype:** switch from %TH+ to Table S3 DA-neuron (or totalDopa) fraction? The evidence above favours it. %TH+ stays as the reproduction gate (A0-6) and a secondary phenotype.

### Prior exposure (add to registration)
- Day-0 bulk RNA was used in the A0-4 gate only: v0.1.0 (failed) and v0.2.0 (partial), both frozen before fitting.
- No clinical variables were touched.

---

## 2026-11-12 (local session): day-0 positive-control gate (A0-4) FAILED; registration text ready

The planning-session review above was received and acted on. The day-0 gate was pre-specified, frozen and run. **It failed.** Nothing was tuned after seeing the result.

### Day-0 gate: what was frozen
Tag [`day0-gate-v0.1.0`](https://github.com/xaviwho/cell-to-clinic/releases/tag/day0-gate-v0.1.0) (commit `5678e46`), set before fitting. Script: `scripts/day0_gate.py`.

- **Data:** bulk RNA (RNAB) has 99 day-0 samples from 95 lines (81 v1, 12 v2-only, control ×5). Every scRNA donor line has exactly one day-0 sample.
- **Analysis set:** 79 donor lines. The control line is excluded because it would leak across folds.
- **Outcome:** high vs low DA yield, split at the median of **our** per-line %TH+ (≥1 UMI primary, ≥2 UMI sensitivity).
- **Model:** day-0 VST, with the 1,000 most variable genes and z-scoring chosen within each training fold. L2 logistic regression with C tuned in an inner 5-fold CV. 20 repeats of stratified 5-fold CV. 200 label permutations.
- **Pass rule:** AUC ≥ 0.75 with p < 0.05. Partial: AUC ≥ 0.65 with p < 0.05.
- Also reported: a batch-only baseline (Table S1 batch) and the ridge-regression Spearman ρ for continuous %TH+.
- **scikit-learn is also partly blocked by Application Control.** Logistic regression, ridge, AUC and CV splitters were reimplemented in numpy/scipy. Checked on synthetic data before the freeze: planted signal gave AUC 1.0 (pass), null data gave 0.60, p = 0.08 (fail).

### Result
| Model | Mean CV AUC (SD) | Permutation p |
|---|---|---|
| Day-0 RNA, primary (%TH+ ≥1 UMI) | **0.573** (0.058) | 0.15 |
| Day-0 RNA, sensitivity (%TH+ ≥2 UMI) | 0.626 (0.034) | 0.055 |
| Batch only (5 batches) | 0.622 | n/a |
| Continuous: ridge Spearman ρ (primary %TH+) | 0.07 | n/a |

**Verdict: fail.** Under the draft's own rule (A0-4: "fix pipeline before any clinical test"), no clinical test proceeds until this is resolved. Batch alone predicts DA yield about as well as day-0 expression.

### Possible reasons (not tested; testing them would be post-hoc)
1. **Different target or features from the paper.**
   - The published AUC 0.93 / CV 0.85 may use a different outcome (e.g. extreme groups, or Table S3 DA-neuron proportions rather than %TH+).
   - It may also use a different gene set (e.g. pre-selected markers, protein-coding genes only).
   - Our VST has 94,946 features. The top-variance 1,000 may be dominated by non-coding, low-count or batch-driven genes.
2. **Batch confounding.** Batch/media drives %TH+ and also shifts day-0 expression. A signal in the paper could partly be batch.
3. **The signal may be weaker than published** with this outcome definition and n = 79.

### For the user / planning session: decision needed
Pick one route. Whichever is chosen must be pre-specified (a new tag `day0-gate-v0.2.0`) before it is run, and this failed v0.1.0 result is reported alongside it in the registration.

- **(a) Replicate the published analysis exactly first.** Read the paper's methods for the day-0 classifier (outcome definition, gene set, model) and code that as v0.2.0. This tests the pipeline itself, which is what A0-4 is for. *(Recommended.)*
- **(b) Treat the failure as informative.** Report A0-4 as failed. Downgrade the cell-to-clinic aims to exploratory, or drop day-0 → yield as a gate.
- **(c) Pre-specify a small, principled feature restriction:** protein-coding genes only, a minimum-expression filter, and batch as a covariate in every fold. Fewer forking paths than (a), but it is a guess at what went wrong.

### Registration draft: text to add
An attempt to write this into the pre-registration doc directly was blocked by the local permission policy, so the text is here to paste:

- **Phenotypes → Phenotype source (decided 2026-11-12):**
  - %TH+ per line comes from our Cell Ranger re-annotation (TH ≥1 UMI primary, ≥2 UMI sensitivity). It reproduces published Table S1 %TH+ (r = 0.999, n = 79), so it is a **reproduction gate**: a pipeline check, not independent evidence that %TH+ measures DA yield.
  - Cell-type proportions: published Table S3 proportions are primary. Our broad marker labels are a sensitivity analysis only (DA vs totalDopa r = 0.57).
  - Expression phenotypes: all-cell pseudobulk with composition covariates (Table S3 proportions or %TH+). Cell-type pseudobulk is secondary, and only if the published per-cell labels arrive from LONI/DZNE.
  - Control line: its 4 scRNA runs (3 batches, E6/E8 media) are kept separate to estimate batch, media and technical variance, not averaged. If a single control value is ever needed, which run is used is pre-specified.
- **Data → prior exposure (add):** by 2026-11-12 the following processed data had been opened:
  - day-65 scRNA matrices, summarised to per-line QC, %TH+, broad cell types and pseudobulk, then compared with Tables S1/S3;
  - day-0 bulk RNA, used only for the A0-4 gate, which failed (AUC 0.57).

  No association with diagnosis, genotype, PRS or clinical variables has been run. Methods were frozen by git tags before fitting (`scrn-v0.1.0`, `day0-gate-v0.1.0`).
- **Gates table, new row:** A0-6 Reproduction. Does our scRNA pipeline reproduce the published per-line %TH+? Comparator: Table S1 scRNA %TH+. Rule: r ≥ 0.80 (0.60–0.80 partial). **Result: r = 0.999, n = 79 (pass).** If it fails: use published cell metadata.
- **Gates table, A0-4 result:** fail at v0.1.0 (AUC 0.57, p = 0.15; batch-only 0.62). Resolution pending (above).

### Still open from earlier
- LONI/DZNE request for the published per-cell labels (option b): not yet sent.
- Power simulation; pseudobulk vs cell-level model choice; OSF filing.

---

## 2026-11-12 (planning session): review of the scRNA run and gate

Good run: frozen before data, 0 failures, cohort n = 47 confirmed. Answers and notes below.

### Answer: where cell-type phenotypes come from
- **Use Table S3 proportions as the primary cell-type phenotypes**, and our marker labels only as a sensitivity analysis. Table S3 matches the published clusters, so reviewers can check it and it's comparable to the paper. Our broad labels agree only moderately (r = 0.57).
  - Already seen (public), so no new exposure; just name it in the registration.
- **Cell-type-specific pseudobulk** (e.g. expression within DA neurons) needs the published per-cell labels. Our labels aren't good enough to define DA neurons for expression analysis. So:
  - send the option (b) request to LONI/DZNE now (per-cell barcode → published cell type table). Still non-blocking;
  - meanwhile, primary expression analyses use **all-cell pseudobulk with composition as covariates** (Table S3 proportions or %TH+). Cell-type pseudobulk becomes a later, secondary analysis if the table arrives.

### Notes on the gate
- **Wording for the registration:** r = 0.999 shows we **reproduce the published %TH+ measure** (same quantity, same data). It's a pipeline sanity check, not independent evidence that %TH+ measures DA yield. Describe it as a reproduction gate.
- **Control line:** the gate averages its 4 runs. That's fine for the gate, but don't average them in analyses. Those 4 runs (E8/E6 media, 3 batches) are the only direct estimate of batch/media and technical variation in %TH+. Keep them as separate runs for that purpose. Pre-specify which run to use if a single control value is ever needed.
- **%TH+ ≥1 UMI** reproduces the published measure but is the noisier definition. Keep ≥2 UMI as the pre-specified sensitivity analysis, as you have it.

### For the local session
1. Add to the registration draft: phenotype source = Table S3 (primary), our labels (sensitivity); %TH+ = reproduced published measure; control-line runs used for technical variance.
2. Next steps 2–4 in your list stand. For step 3 (day-0 bulk RNA), the positive-control gate should predict **our** per-line %TH+ from day-0 expression, with the published AUC (0.93) as the benchmark. Pre-specify the model and pass threshold before fitting.

---

## 2026-11-12 (local session): scRNA pipeline frozen, run on all 83 runs, validation gate PASSED

All planning-session steps 1–6 are done. The method was frozen **before** any real data was processed: tag [`scrn-v0.1.0`](https://github.com/xaviwho/cell-to-clinic/releases/tag/scrn-v0.1.0) (commit `702d623`).

### Changes to the planning session's script (all made before the freeze)
- **scanpy/anndata/pyarrow removed.** Windows Application Control on the analysis machine blocks their compiled dependencies (numba, zarr/numcodecs, pyarrow DLLs). The policy is not bypassed. `scripts/scrn_per_run.py` now uses only numpy, scipy, pandas and h5py:
  - own 10x `.h5` reader (Gene Expression features only; duplicate symbols made unique);
  - own `score_genes`: the same Seurat/scanpy algorithm (25 expression bins, 50 control genes per bin, seed 0), not bit-identical to scanpy;
  - outputs written as `.csv.gz` instead of parquet.
- **GBA1 is looked up as `GBA1` or `GBA`.** FOUNDIN-era Cell Ranger references name it `GBA`.
- **TH rule:** ≥1 UMI is primary (as planned); **≥2 UMI is a pre-specified sensitivity metric**, computed alongside (`frac_TH_pos_all_sens`, `frac_DA_neuron_sens`).
- QC thresholds and marker sets unchanged: ≥500 genes, ≥1000 UMIs, ≤20% mito, ≤8000 genes.
- **Gate pre-specified** (thresholds in the script): Pearson r between our per-line %TH+ and published per-line %TH+.
  - r ≥ 0.80: pass. 0.60–0.80: partial. < 0.60: fail.
  - The per-line comparator is **Table S1 "scRNA %TH+ cells"**, matched by PATNO. Table S3 holds cell-type proportions, not %TH+.
  - Secondary, not gating: our DA-neuron fraction vs Table S3 `totalDopa`.
  - New `scripts/scrn_gate.py` prints aggregates only.
- Tested on synthetic Cell Ranger files before the freeze: planted cell types recovered with 100% accuracy; gate script checked on fake tables.
- Minor fixes: error log now closed properly; peak RAM is recorded per run.

### Run
- Dry run found **83 runs** (as expected). Unpacked size 17 GB.
- 83/83 processed, **0 failures**. About 4 s per run, **≤0.93 GB peak RAM**, about 5 min total.

### Results (aggregates only)
| Item | Value |
|---|---|
| Runs / distinct lines | 83 / **80** (control line ×4 runs). The earlier "82" was wrong. Table S3 also has 83 rows / 80 distinct IDs. |
| Cells per run, raw (Cell Ranger filtered) | median 6,460 (IQR 5,130–7,480) |
| Cells per run after QC | median 5,420 (IQR 4,120–6,210); QC pass rate median 82% (IQR 77–87%) |
| %TH+ (≥1 UMI) | median 12.1% (IQR 9.4–17.3) |
| %TH+ (≥2 UMI, sensitivity) | median 6.1% (IQR 4.2–9.9) |
| Lines with <10% / <5% TH+ (≥1 UMI) | 26 / 2 of 80 |
| Broad compartments (median %) | Progenitor 60.6 · Other neuron 17.6 · DA neuron 9.8 · Proliferating 3.3 · Ependymal-like 3.2 · Unassigned 1.0 |

**Validation gate: PASS.**
- Primary: r = **0.999** (Spearman 0.999), n = 79 lines with a published value. All 79 matched; none missing from ours.
- Sensitivity (≥2 UMI): r = 0.979.
- Reading: published "scRNA %TH+ cells" is the fraction of cells with any TH UMI, and we reproduce it. **The 42.7 GB Seurat object is not needed for %TH+.**

**Secondary (not gating): DA-neuron fraction vs Table S3 `totalDopa`: r = 0.57** (Spearman 0.53), n = 80.
- Our broad marker-based labels do **not** reproduce the published cell-type clusters (immature + mature DA neurons, progenitor subtypes).
- **Caveat:** treat our cell-type labels (and cell-type pseudobulk) as provisional. %TH+ and all-cell pseudobulk are fine.
- If an aim needs published cell types, use Table S3 proportions directly (per line, public) or obtain the published cell metadata (option b).

### Cohort for day-65 scRNA (public Table S1, matched to our 80 lines)
| Group | Lines |
|---|---|
| Idiopathic PD | 28 |
| Genetic PD | 19: LRRK2 11, GBA 6 (incl. one E326K), SNCA 2 |
| **All PD** | **47** (confirms planning n) |
| Unaffected carriers | 22: GBA 10, LRRK2 10, LRRK2+GBA 1, SNCA 1 |
| Healthy controls | 8 (incl. the control line) |
| Prodromal / SWEDD | 2 / 1 |

### Prior exposure (update for the registration)
- Day-65 scRNA Cell Ranger matrices opened and summarised: %TH+, QC metrics, broad compartments, per-cell-type pseudobulk.
- Per-line values were compared with published Table S1/S3. No association with diagnosis, genotype or PRS was run, and none was looked at.
- **No PPMI clinical outcomes pulled.**
- The cohort counts above come from public Table S1 only.

### Next (local session)
1. Settle the phenotype set. %TH+ is validated. Should cell-type proportions use our labels or Table S3?
2. Decide pseudobulk vs cell-level modelling. All-cell pseudobulk per run is already in `data/derived/scrn/`.
3. Check bulk RNA (RNAB) headers for day-0 samples, for the day-0 → DA-yield positive-control gate.
4. Power simulation, then OSF filing.

### For the planning session
- Review the gate outcome and the secondary-comparison caveat. Should cell-type phenotypes come from Table S3 (published, public) rather than our re-annotation?

---

## 2026-11-12 (planning session): scRNA pipeline, no waiting on LONI

### Decision
- **Seurat object (open item 4): go with option (c)**, re-annotating from the Cell Ranger `.h5` files in Python. It fits in 31 GB RAM and needs nobody's reply. Asking LONI/DZNE for the published cell-level metadata (option b) stays useful as a cross-check, but nothing waits on it.
- The LONI email is no longer blocking. Send it when convenient; the ID-mapping question is dropped.

### Added
- `scripts/scrn_per_run.py`: per-run QC, marker-based cell-type annotation (DA neuron / other neuron / progenitor / proliferating / ependymal-like), %TH+, detection rates for key genes, and pseudobulk raw counts per cell type. Processes one run at a time, resumes if interrupted, writes only inside `data/`, and prints no IDs or per-run values. Tested end to end on synthetic Cell Ranger files (not on real data).
- `requirements.txt`.

### For the local session
1. **Freeze the parameters first.** Review the QC thresholds and marker sets at the top of the script with the user and settle two choices:
   - DA neuron rule: TH ≥ 1 UMI (current) vs ≥ 2 (more robust to ambient RNA).
   - QC cut-offs: ≥ 500 genes, ≥ 1000 UMIs, ≤ 20% mito, ≤ 8000 genes.

   Then commit and **tag** (e.g. `scrn-v0.1.0`) **before** running on real data, so the registration can show the rules were fixed in advance.
2. `pip install -r requirements.txt`; unpack only the SCRN CellRanger bundle into `data/unpacked/`.
3. `python scripts/scrn_per_run.py --in-dir data/unpacked --out-dir data/derived/scrn --dry-run`, which should report 83 runs. If not, fix the run discovery in `find_runs()` and log it here.
4. Run with `--limit 2`; record time and peak RAM per run here. Then run all.
5. Report back here, **aggregates only**:
   - runs processed / failed, and number of distinct lines (settles 80 vs 82);
   - distribution (median, IQR) of cells after QC and %TH+;
   - **validation gate:** correlation between our per-line %TH+ and the published Table S3 values (r and n only). If it's high, the Seurat object isn't needed.
6. Log the run as prior exposure: day-65 scRNA matrices opened and summarised; still no clinical outcomes.

---

## 2026-11-12: Data integrity check and structure review

### Done
- **Integrity check:** all 9 processed-data bundles (Project 150 release) were test-read end to end with `tar -tzf`. Every one returned exit code 0 with no errors, and the sizes match the server listing. **No re-download needed** before the Aspera account expires (~2026-12-11).
- File listings and READMEs were extracted to `data/listings/` and `data/readmes/`; both stay local. Nothing else has been unpacked yet.
- Added a `.gitignore` covering `data/` and genomic and single-cell file types.

| Bundle | Tar check | Entries | Contents (summary) |
|---|---|---|---|
| 150.1 GENO | OK | 194 | NeuroChip PLINK set (day 0 + day 65); per-sample CNV PDFs |
| 150.2 PROT | OK | 5 | TMT proteomics: 3 replicate CSVs, **10 lines only** |
| 150.4 METH | OK | 8 | EPIC normalised matrix, probe annotation, QC reports |
| 150.5 HICS | OK | 20 | HiCCUPS loop BEDPEs, **~10 lines** (day 0 / day 65) |
| 150.6 ATAC | OK | 293 | MACS2 peak BEDs, day 0 / 25 / 65 |
| 150.8 SCRN CellRanger | OK | 335 | 83 day-65 runs (filtered + raw `.h5`, web summaries) |
| 150.8 SCRN Seurat | OK | 3 | one integrated object with broad cell types (42.7 GB RDS) |
| 150.9 RNAB | OK | 10 | aggregated bulk RNA tables only (counts / CPM / TPM / VST, circRNA) |
| 150.10 RNAS | OK | 19 | exceRpt small-RNA summaries |

### Findings
- **Sample IDs:** file names encode the PPMI PATNO plus the CDI cell-line suffix, and timepoint and batch (README naming convention). No separate ID-mapping file is needed. The processed data can be joined directly to PPMI clinical data, so the "no outcomes before registration" firewall is procedural only and must be stated in the registration.
- **No sample metadata in the processed bundles.** Sample sheets sit in the RAW folders, which are outside the approval. `ATAC_BED_SAMPLENAMES.txt`, mentioned in the ATAC README, is not in the bundle. Diagnosis, genetic group and batch must come from public Table S1 (Bressan et al. 2023).
- **Line counts:**
  - Day-0 ATAC: 99 samples from 95 distinct lines. Registration wording: **94 donor lines + control line ×5 = 99 samples**.
  - Day-65 scRNA (Cell Ranger): 83 runs from **80 distinct lines** (79 donor lines + control line ×4: batches 1, 3, 5-E6, 5-E8). The earlier working figure was 82, so this needs reconciling with Table S1 and the Seurat metadata before n = 47 PD donors is final.
- **Proteomics (10 lines) and Hi-C (~10 lines)** are too small for donor-level association; at most they can support descriptive or secondary analyses.
- **Compute constraint:** the local machine has 31 GB RAM, and R and PLINK are not installed. The 42.7 GB compressed Seurat RDS will likely need >100 GB RAM to load.

### Open items / next steps
1. Unpack into the local data area. 878 GB is free, so space is not a constraint. Decide whether to unpack the Seurat RDS now.
2. Check the bulk RNA column headers for day-0 samples (counts only).
3. Reconcile the scRNA line count (80 vs 82) with Table S1 and the Seurat metadata.
4. Decide how to handle the Seurat object:
   - (a) use an approved larger-memory environment;
   - (b) request cell-level metadata (cell type, TH+, sample) from LONI / DZNE as a small table;
   - (c) re-annotate from the Cell Ranger `.h5` files in Python.

   This decision feeds the pseudobulk vs cell-level model choice.
5. Email LONI with these questions:
   - approved storage and analysis environment, including larger-memory compute;
   - cells per line;
   - Aspera extension if needed.

   Drop the ID-mapping question: the files are keyed by PATNO.
6. Finish the power simulation, settle the open decisions, then file the OSF pre-registration (secondary-data template, embargoed) and make a tagged GitHub/Zenodo release of the gate scripts.

### Prior exposure (for the registration)
Public paper and Tables S1/S3. Processed data: file listings and README files only. No data matrices have been opened or analysed, and no PPMI clinical outcomes have been pulled.

---

## 2026-11-11: Data access and download
- PPMI Tier 3 access approved for FOUNDIN-PD Project 150, processed data only (MJFF approval 2026-10-28).
- LONI Aspera access set up; account expires ~2026-12-11.
- All needed processed bundles downloaded 2026-11-11/12. Not downloaded:
  - single-cell ATAC (`150.7_SCAT*`): only if an aim needs it;
  - raw data: outside the approval;
  - the unpacked duplicate of the processed data.
- Only Project 150 processed data is in scope. Other server folders are visible but not approved. Tier 2 Project 233 (phosphoproteomics) would need its own approval.

## Governance rules (standing)
- No linking to PPMI clinical, imaging or SAA outcomes for FOUNDIN PATNOs until the OSF registration is filed.
- Data stays local in an access-restricted location; never upload it to cloud storage or GitHub.
- When using AI assistants: run scripts locally and share only summaries or structure, never participant-level rows or IDs.
