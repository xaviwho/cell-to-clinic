# Activity log: FOUNDIN-PD cell-to-clinic

Progress record for the FOUNDIN-PD (PPMI Project 150) cell-to-clinic analysis. The newest entry is first.

> **Public repository.** This log holds no participant-level information: no PATNOs, cell-line (CDI) numbers, per-line values, credentials or data. Controlled-access data stays on the local machine under `data/`, which is git-ignored.

> **Two sessions share this log.** The **local session** (Claude in VS Code, on the machine with the data) is the **main environment**: from 2026-11-12 it does all the work, including scripts, analyses, commits and tags. The **planning session** (Claude on claude.ai, no data access) reads this log and can raise requests or reviews. Each entry is tagged with which session wrote it. Requests go under **"For the local session"** / **"For the planning session"** headings. Pull before starting, commit and push when done. Only aggregate, non-identifying results go here.

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
