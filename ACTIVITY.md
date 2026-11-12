# Activity log: FOUNDIN-PD cell-to-clinic

Progress record for the FOUNDIN-PD (PPMI Project 150) cell-to-clinic analysis. The newest entry is first.

> **Public repository.** This log holds no participant-level information: no PATNOs, cell-line (CDI) numbers, per-line values, credentials or data. Controlled-access data stays on the local machine under `data/`, which is git-ignored.

> **Two sessions share this log.** The **local session** (Claude in VS Code, on the machine with the data) is the **main environment**: from 2026-11-12 it does all the work, including scripts, analyses, commits and tags. The **planning session** (Claude on claude.ai, no data access) reads this log and can raise requests or reviews. Each entry is tagged with which session wrote it. Requests go under **"For the local session"** / **"For the planning session"** headings. Pull before starting, commit and push when done. Only aggregate, non-identifying results go here.

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
