# Hyperparathyroidism scRNA-seq Atlas — PhD Thesis Code Reproduction Repository

**Thesis chapters:** Chapter 2 (Single-Cell Transcriptome Atlas), Chapter 3 (SHPT Functional Experiments), Chapter 4 (PHPT Functional Experiments)
**Institution:** Xiangya School of Medicine, Central South University
**Author:** HUI Ouyang — ouyanghui950319@gmail.com
**Status:** Thesis in preparation (2026) — code provided for reviewer reproduction and reader reference.

---

## 1. What this repository is

This repository contains the **complete analytical source code** used to produce every figure in Chapters 2–4 of my PhD dissertation. It is a *code-reproduction* repository, not a data repository — raw sequencing data and large intermediate objects are hosted separately (see §7) and are downloaded on demand by the scripts below.

**Chapter 2** builds the first integrated single-cell transcriptome reference of human parathyroid tissue across three disease states (primary hyperparathyroidism, PHPT; secondary hyperparathyroidism, SHPT; and histologically normal parathyroid, Normal). Starting from 15 fresh surgical specimens we generated **170,852 high-quality cells** and annotated **16 cell types** spanning endocrine, stromal, vascular, and immune compartments. On top of this reference atlas we developed and benchmarked **scMMR**, an in-house multi-task deep neural network framework that unifies cell-type annotation, low-dimensional embedding, gene/pathway/transcription-factor importance scoring, proportion-perturbation testing, expression-perturbation testing, and bulk RNA-seq deconvolution. Key biological outputs of the chapter include: (i) identification of disease-specific parathyroid subpopulations (PHPT Cluster 3 and SHPT Cluster 2) whose abundance correlates with both tumour size and serum PTH and is independently validated on an external bulk cohort by BayesPrism; (ii) three parathyroid differentiation lineages reconstructed by Monocle3 / RNA velocity / CytoTRACE2 that jointly recapitulate the normal → SHPT and normal → PHPT trajectories; and (iii) the CASR↓ / PTH↑ transcriptional program that accompanies pseudotime progression in both disease lineages.

**Chapter 3** focuses on SHPT (secondary hyperparathyroidism) functional experiments: differential gene screening between SHPT and Normal parathyroid cells, gene validation via expression-level analysis, GSEA enrichment of proliferation-related pathways, and pathway activity scoring.

**Chapter 4** focuses on PHPT (primary hyperparathyroidism) functional experiments: parallel gene screening and validation for PHPT vs Normal, GSEA and pathway analysis, and CellChat-based cell-cell communication analysis revealing disease-specific ligand-receptor interactions.

## 2. Developed tools released alongside this thesis

Two R/Python packages were written specifically for this thesis and are released under separate repositories:

- **scMMR** v0.3.0 — multi-task DNN framework for scRNA-seq; six functional modules (`DNN_annotation`, `PlotMAP`, `IntegratedGradients`, `RankPercent`, `RankPerturbation`, `DNN_deconvolution`). Repository: https://github.com/HUI950319/scMMR
- **UtilsR** v0.5.0 — lightweight visualization / utility helper used by every figure script. Repository: https://github.com/HUI950319/UtilsR
All figure scripts in `Ch2_Fig02_20_Visualization/` import scMMR and UtilsR. Installing these two packages is the main dependency for running the visualization layer.

## 3. Repository layout

```
PhDThesis-Scripts/
├── 00_scRNA_upstream/            # Cell Ranger v8.0.1 pipeline (fastq → filtered matrices)
├── 00_bulk_upstream/             # fastp + STAR v2.7.10a + featureCounts v2.0.3 (bulk RNA-seq)
├── 00_scVelo/                    # Velocyto + scVelo v0.3.4 RNA velocity
├── 00_scenic_pipe/               # pySCENIC v0.12.1 gene regulatory network inference
│   └── cistarget/DOWNLOAD.md     # download links for the 4× feather + 2× motif TBL databases (~1.25 GB total, not tracked)
├── 00_GeneList/                  # Proliferation-related pathway gene list construction (msigdbr)
├── Ch2_Fig02_20_Visualization/   # R scripts that draw every figure in Chapter 2 (Fig 2-2 through Fig 2-20)
├── Ch3_Fig21_24_SHPT/            # R scripts for Chapter 3 figures (SHPT gene screening, GSEA, pathway analysis)
├── Ch4_Fig25_29_PHPT/            # R scripts for Chapter 4 figures (PHPT gene screening, GSEA, pathway, CellChat)
├── .gitignore                    # ignores *.rds / *.h5ad / *.bam / backup_* / cisTarget references / IDE junk
└── README.md                     # this file
```

Most figures follow a two-file `*_Data.R` → `*_Plot.R` convention:
- `*_Data.R` — reads the integrated Seurat / AnnData object, runs the statistical computation, saves a compact intermediate `.rds` under the object's `extdata/` slot.
- `*_Plot.R` — loads that intermediate and draws the final publication figure with `UtilsR` + `ggplot2` helpers.

Three layout variants exist in practice and are reflected one-to-one in the §5 table below:
- **Single-file figures** (`Fig02_QualityControl.R`, `Fig03_CellType.R`) — data prep and plotting are trivial enough to live in one script.
- **One-to-one pairs** (`Fig04_Benchmark.R` + `Fig04_Benchmark_Plot.R`, `Fig05_EvalEmbedding_*`, `Fig06_ROC_Calibration_*`, `Fig19_ClinicalCorrelation_*`, `Fig20_BayesPrism_*`) — one Data script feeds exactly one Plot script. Note that `Fig04_Benchmark.R` is the Data half despite lacking the `_Data` suffix (kept for historical consistency with the benchmark harness).
- **Grouped pairs shared across tightly related panels** (`Fig07_08_Data.R` feeds both `Fig07_PlotMAP_Plot.R` and `Fig08_QC_Scatter_Plot.R`; `Fig09_11_Interpretation_*` covers Fig 2-9/10/11; `Fig12_15_Perturbation_*` covers Fig 2-12/13/14/15; `Fig16_18_Trajectory_*` covers Fig 2-16/17/18) — a single intermediate `.rds` is reused to draw several panels, which is why the Data script name spans a figure range.

This split means a reviewer can regenerate any panel in seconds once the upstream object is on disk, without re-running Harmony integration or the DNN.

Ch3 and Ch4 follow the same two-file `*_Data.R` → `*_Plot.R` convention as Ch2, with a shared data-preparation script (`Fig21_22_GeneScreening_Data.R` for Ch3, `Fig25_26_GeneScreening_Data.R` for Ch4) feeding multiple downstream plot scripts.

## 4. Datasets used in Chapter 2

| # | Dataset | Assay | Species | Samples | Role in chapter | Source |
|---|---------|-------|---------|---------|-----------------|--------|
| 1 | **Self-Para** (this study) | 10× 3′ v3.1 scRNA-seq | Human | 15 (PHPT 6, SHPT 6, Normal 3) | Reference atlas, Benchmark, perturbation, trajectory, deconvolution reference | GSA-Human (accession on publication) |
| 2 | **GSE190773** | 10× 3′ scRNA-seq | Human | 5 (PHPT 5) | Benchmark, PlotMAP projection, ROC/calibration, deconvolution | GEO GSE190773 |
| 3 | **GSE233962** | 10× 3′ scRNA-seq | Human / M. mulatta | 8 (Human 4, Macaque 4) | Benchmark, PlotMAP projection, ROC/calibration | GEO GSE233962 (Venkat et al. 2024) |
| 4 | **Baron** | inDrop scRNA-seq | Human | 4 donors, 8,569 cells (pancreas, non-diabetic) | Cross-tissue benchmark | scRNAseq R pkg (Baron et al. 2016) |
| 5 | **Muraro** | CEL-Seq2 scRNA-seq | Human | 4 donors, 3,072 cells (pancreas, cadaver donors) | Cross-tissue benchmark | scRNAseq R pkg (Muraro et al. 2016) |
| 6 | **Segerstolpe** | Smart-Seq2 scRNA-seq | Human | 10 donors, 3,514 cells (pancreas, 6 healthy + 4 T2D) | Cross-tissue benchmark | scRNAseq R pkg (Segerstolpe et al. 2016) |
| 7 | **Xin** | Fluidigm C1 scRNA-seq | Human | 18 donors, 1,600 cells (pancreas, 12 non-diabetic + 6 T2D) | Cross-tissue benchmark | scRNAseq R pkg (Xin et al. 2016) |
| 8 | **Self-PHPT-Bulk** | Bulk RNA-seq (NovaSeq PE150) | Human | 12 (PHPT 12, with tumour size and PTH) | Deconvolution and clinical correlation, BayesPrism validation | GSA (accession on publication) |
| 9 | **Self-SHPT-Bulk** | Bulk RNA-seq (NovaSeq PE150) | Human | 12 (SHPT 12, with tumour size and PTH) | Deconvolution and clinical correlation, BayesPrism validation | GSA (accession on publication) |
| 10 | **PRJNA516535** | Bulk RNA-seq (HiSeq 2500) | Human | 10 (PHPT adenoma 10, with clinical info) | External deconvolution validation, BayesPrism validation | NCBI SRA / GSE125433 (He et al. 2019) |

**Sample provenance.** All parathyroid samples were collected at Xiangya Hospital of Central South University between 2023 and 2025. PHPT samples were solitary adenomas removed during curative parathyroidectomy; SHPT samples came from end-stage renal-disease patients undergoing total parathyroidectomy with forearm autotransplantation; Normal samples were histologically normal parathyroid glands incidentally resected during thyroid lobectomy. Libraries were prepared with the 10× Chromium Next GEM Single Cell 3′ Kit v3.1 (PN-1000268) and sequenced on Illumina NovaSeq 6000 PE150 by OE Biotech (Shanghai). Bulk libraries were prepared and sequenced on the same platform.

## 5. Chapter 2 figure index

The table below links every thesis figure to the script that produces it and the biological message conveyed by that figure. Script paths are relative to `Ch2_Fig02_20_Visualization/`.

| Figure | Script(s) | What it shows |
|--------|-----------|----------------|
| Fig 2-2 | `Fig02_QualityControl.R` | Per-sample QC: nFeature_RNA, nCount_RNA, percent.mt violin plots; DoubletFinder doublet calls; decontX ambient-RNA contamination score. Cells passing all filters (200–10,000 features, 1,000–50,000 UMIs, <40% mt, non-doublet, decontX <0.5) advance to integration. |
| Fig 2-3 | `Fig03_CellType.R` | UMAP of the full 170,852-cell integrated reference coloured by the 16 annotated cell types, accompanied by a marker-gene dot/heatmap panel. |
| Fig 2-4 | `Fig04_Benchmark.R` + `Fig04_Benchmark_Plot.R` | Accuracy / Macro-F1 benchmark of 12 annotation methods (XGBoost, RF, SVM, Elastic-Net, LDA, Naive Bayes, KNN, NNET, SingleR, Seurat Label Transfer, CellTypist, **scMMR**) on 7 datasets under 10 × 10-fold stratified CV. scMMR ranks first on both metrics (>0.95 / >0.90). |
| Fig 2-5 | `Fig05_EvalEmbedding_Data.R` + `Fig05_EvalEmbedding_Plot.R` | Silhouette score and standardized pairwise distance for scMMR 512-d embeddings vs. PCA embeddings — scMMR yields tighter and better-separated clusters. |
| Fig 2-6 | `Fig06_ROC_Calibration_Data.R` + `Fig06_ROC_Calibration_Plot.R` | Macro-averaged ROC and per-class AUC on the in-house reference (AUC 0.9961, Acc 0.9925), GSE190773 (Acc 0.9774) and GSE233962 (Acc 0.9655). |
| Fig 2-7 | `Fig07_08_Data.R` + `Fig07_PlotMAP_Plot.R` | Projection of GSE190773 and GSE233962 cells onto the reference UMAP using scMMR's `PlotMAP` module. |
| Fig 2-8 | `Fig07_08_Data.R` + `Fig08_QC_Scatter_Plot.R` | Pearson correlation between scMMR prediction-confidence and DecontX contamination (r = −0.577 and −0.443) — low-confidence cells are the ambient-contaminated ones. |
| Fig 2-9 | `Fig09_11_Interpretation_Data.R` + `Fig09_11_Interpretation_Plot.R` | Cell-type-specific gene importance via Integrated Gradients: PTH and CASR dominate parathyroid cells; CD3D / CD3E dominate T cells; consistent with canonical markers. |
| Fig 2-10 | `Fig09_11_Interpretation_Data.R` + `Fig09_11_Interpretation_Plot.R` | KEGG pathway importance aggregated from IG scores — parathyroid hormone synthesis / calcium signaling are the top parathyroid pathways. |
| Fig 2-11 | `Fig09_11_Interpretation_Data.R` + `Fig09_11_Interpretation_Plot.R` | Transcription-factor importance (IG on pySCENIC regulon activity): GCM2 and MAFB are ranked highest for parathyroid cells; TCF7 / LEF1 for T cells. |
| Fig 2-12 | `Fig12_15_Perturbation_Data.R` + `Fig12_15_Perturbation_Plot.R` | Alluvial plot of cell-type composition across Normal / PHPT / SHPT, with parathyroid epithelial fractions 57.1 % / 71.4 % / 73.8 %. |
| Fig 2-13 | `Fig12_15_Perturbation_Data.R` + `Fig12_15_Perturbation_Plot.R` | miloR differential-abundance neighbourhoods correlated with scMMR `RankPercent` proportion perturbations (Spearman ρ = 0.859) — two independent methods agree. |
| Fig 2-14 | `Fig12_15_Perturbation_Data.R` + `Fig12_15_Perturbation_Plot.R` | scMMR `RankPerturbation` expression-perturbation scores vs. Augur cell-type prioritization — parathyroid cells are the top-perturbed population in both diseases. |
| Fig 2-15 | `Fig12_15_Perturbation_Data.R` + `Fig12_15_Perturbation_Plot.R` | PHPT vs. SHPT perturbation correlation (rs = 0.724 for proportion, 0.657 for expression) — shared and distinct disease programs. |
| Fig 2-16 | `Fig16_18_Trajectory_Data.R` + `Fig16_18_Trajectory_Plot.R` | Monocle3 pseudotime tree overlaid with scVelo RNA-velocity arrows; three lineages emerge: Lineage 1 Normal-dominant, Lineage 2 SHPT-dominant, Lineage 3 PHPT-dominant. |
| Fig 2-17 | `Fig16_18_Trajectory_Data.R` + `Fig16_18_Trajectory_Plot.R` | CytoTRACE2 developmental potential jointly cross-validated against Monocle3 pseudotime and velocity latent-time — all three orderings agree along each lineage. |
| Fig 2-18 | `Fig16_18_Trajectory_Data.R` + `Fig16_18_Trajectory_Plot.R` | CASR / VDR expression decreases and PTH expression increases along SHPT and PHPT pseudotime lineages — the transcriptional signature of the clinical disease phenotype. |
| Fig 2-19 | `Fig19_ClinicalCorrelation_Data.R` + `Fig19_ClinicalCorrelation_Plot.R` | PHPT Cluster 3 fraction correlates with tumour size (R = 0.667, P = 2.28 × 10⁻⁵) and serum PTH (R = 0.694, P = 7.42 × 10⁻⁶); SHPT Cluster 2 fraction correlates with tumour size (R = 0.738, P = 4.67 × 10⁻⁴) and PTH (R = 0.539, P = 0.021). |
| Fig 2-20 | `Fig20_BayesPrism_Data.R` + `Fig20_BayesPrism_Validation_Plot.R` | Independent BayesPrism deconvolution of 12 self + 10 public PHPT bulks and 12 self SHPT bulks confirms scMMR's cluster-fraction estimates (PHPT C3 R = 0.901, P = 1.1 × 10⁻⁸; SHPT C2 R = 0.928, P = 1.33 × 10⁻⁵). |

## 5b. Chapter 3 figure index (SHPT functional experiments)

Scripts are in `Ch3_Fig21_24_SHPT/`. The shared gene-list dependency is in `00_GeneList/Fig23_ProlifGeneList_Data.R`.

| Figure | Script(s) | What it shows |
|--------|-----------|----------------|
| Fig 3-21 | `Fig21_22_GeneScreening_Data.R` + `Fig21_VolcanoVenn_Plot.R` | Differential gene screening between SHPT and Normal parathyroid cells: volcano plot of DEGs and Venn diagram of overlap with proliferation-related genes. |
| Fig 3-22 | `Fig21_22_GeneScreening_Data.R` + `Fig22_GeneValidation_Plot.R` | Gene validation: expression-level box/violin plots of selected proliferation-related DEGs across SHPT vs Normal, with statistical comparisons. |
| Fig 3-23 | `Fig23_GSEA_Data.R` + `Fig23_GSEA_Plot.R` | GSEA enrichment analysis of 34 proliferation-related pathways (Hallmark + KEGG_LEGACY + Reactome/WikiPathways supplements) in SHPT vs Normal. |
| Fig 3-24 | `Fig21_22_GeneScreening_Data.R` + `Fig24_PathwayAnalysis_Plot.R` | Pathway activity scoring: AUCell-based pathway activity across parathyroid cell clusters, showing enrichment of proliferation pathways in SHPT-specific subpopulations. |

## 5c. Chapter 4 figure index (PHPT functional experiments)

Scripts are in `Ch4_Fig25_29_PHPT/`. The shared gene-list dependency is in `00_GeneList/Fig23_ProlifGeneList_Data.R`.

| Figure | Script(s) | What it shows |
|--------|-----------|----------------|
| Fig 4-25 | `Fig25_26_GeneScreening_Data.R` + `Fig25_VolcanoVenn_Plot.R` | Differential gene screening between PHPT and Normal parathyroid cells: volcano plot and Venn diagram with proliferation gene overlap. |
| Fig 4-26 | `Fig25_26_GeneScreening_Data.R` + `Fig26_GeneValidation_Plot.R` | Gene validation: expression-level plots of selected proliferation-related DEGs across PHPT vs Normal. |
| Fig 4-27 | `Fig27_GSEA_Data.R` + `Fig27_GSEA_Plot.R` | GSEA enrichment of 34 proliferation-related pathways in PHPT vs Normal parathyroid cells. |
| Fig 4-28 | `Fig25_26_GeneScreening_Data.R` + `Fig28_PathwayAnalysis_Plot.R` | Pathway activity scoring via AUCell across parathyroid cell clusters in PHPT. |
| Fig 4-29 | `Fig29_CellChat_Data.R` + `Fig29_CellChat_Plot.R` | CellChat cell-cell communication analysis: ligand-receptor interaction networks across cell types in PHPT, including signaling pathway comparison between PHPT and Normal, and visualization of disease-specific communication patterns. |

## 5d. Shared upstream: proliferation gene list

`00_GeneList/Fig23_ProlifGeneList_Data.R` constructs the master proliferation-related pathway gene list used by both Ch3 and Ch4. It queries the MSigDB database (via the `msigdbr` R package, v7.5.1) to assemble 34 pathways: 17 Hallmark, 12 KEGG_LEGACY, 4 Reactome supplements, and 1 WikiPathways supplement. The output `.rds` file is consumed by the GSEA and pathway-analysis scripts in both chapters.

## 6. Software environment (pinned versions)

All versions below are the exact releases used in the thesis. Upstream alignment and SCENIC were run on an internal cloud server (西柚云); all downstream R / Python analyses were run on a local workstation with WSL Ubuntu 22.04.

**Upstream (cloud).** Cell Ranger v8.0.1 with reference GRCh38-2020-A; STAR v2.7.10a with GRCh38 primary + GENCODE v44; fastp (adapter and low-quality trimming); featureCounts v2.0.3 (subread); Velocyto (spliced / unspliced counting for RNA velocity); pySCENIC v0.12.1 (with the four hg38 cisTarget feather databases + v10 motif annotation — see `00_scenic_pipe/cistarget/DOWNLOAD.md`).

**Downstream (local WSL).** R v4.4.3 on Ubuntu 22.04; Python v3.11.14; PyTorch v2.9.1 with CUDA for scMMR training; scanpy v1.11.5; Seurat v5.1.0; Harmony v1.2.1; mlr3verse v0.3.0 (the 12-method benchmark harness); **scMMR v0.3.0**, **UtilsR v0.5.0**, and BayesPrism v2.0; Monocle3 v1.4.26; scVelo v0.3.4; CytoTRACE2 v1.1.0; miloR v2.5.1; Augur v1.0.3; edgeR v4.2.1; CellChat v2.1.0 (Ch4 cell-cell communication); msigdbr (proliferation gene list construction for Ch3/Ch4); AUCell (pathway activity scoring for Ch3/Ch4).

Install the two lab packages first — figure scripts depend on them:

```r
# R ≥ 4.4
if (!require("remotes")) install.packages("remotes")
remotes::install_github("HUI950319/UtilsR")
remotes::install_github("HUI950319/scMMR")
```

## 7. How to reproduce

1. **Obtain the raw data.** Download FASTQs for the in-house ParaDataset (GSA-Human accession on publication), the in-house bulk cohorts, and the public datasets listed in §4. Place them under a user-defined `data/` root.
2. **Run upstream alignment.** Use `00_scRNA_upstream/` for 10× libraries (Cell Ranger v8.0.1, reference GRCh38-2020-A) and `00_bulk_upstream/` for bulk libraries (fastp → STAR → featureCounts). Outputs: filtered count matrices per sample.
3. **Build the integrated reference.** QC + DoubletFinder + decontX filters → merge → SCTransform/NormalizeData → Harmony integration (50 PCs) → UMAP (30 dims) → Louvain clustering (resolution 0.2) → manual annotation of the 16 cell types. This produces the master Seurat object consumed by every downstream script.
4. **Train and apply scMMR.** Binarize the gene expression matrix, train the multi-task DNN (InputBlock → 4× ResNetBlock → ResNetLastBlock → 512-d shared embedding → dual heads for classification + UMAP regression) with PyTorch 2.9.1, then run the six downstream modules.
5. **RNA velocity and SCENIC.** Run Velocyto + scVelo (`00_scVelo/`) and pySCENIC (`00_scenic_pipe/`) using the outputs of Step 3.
6. **Regenerate Chapter 2 figures.** For each figure *N*, execute the matching Data script once — see the §5 table for the exact filename, since a single Data script such as `Fig07_08_Data.R`, `Fig09_11_Interpretation_Data.R`, `Fig12_15_Perturbation_Data.R` or `Fig16_18_Trajectory_Data.R` covers a group of related panels — then run the corresponding Plot script as many times as desired to iterate on cosmetics.
7. **Regenerate Chapter 3/4 figures.** First run `00_GeneList/Fig23_ProlifGeneList_Data.R` to build the proliferation gene list. Then for Ch3, run `Fig21_22_GeneScreening_Data.R` followed by each plot script; for Ch4, run `Fig25_26_GeneScreening_Data.R` followed by each plot script. `Fig29_CellChat_Data.R` (Ch4) requires CellChat v2.1.0 and produces the cell-cell communication intermediate consumed by `Fig29_CellChat_Plot.R`.

Expected wall time: ~2–3 days end-to-end on a single workstation with 256 GB RAM and a single consumer GPU, excluding FASTQ download and cloud upstream.

## 8. Ethics and patient consent

All human parathyroid and bulk samples were collected under a protocol approved by the **Medical Ethics Committee of Xiangya Hospital, Central South University (approval no. 2023030424)**. Written informed consent was obtained from every participant prior to sample collection. Patient identifiers were removed before any computational analysis; only de-identified metadata (disease group, age bracket, serum chemistry, tumour size) is used in the scripts distributed here.

## 9. Data availability

Because the cisTarget motif reference databases (~1.25 GB total) exceed the GitHub per-file 100 MB limit they are **not** tracked in this repository; `00_scenic_pipe/cistarget/DOWNLOAD.md` provides the Aerts-Lab download URLs and the expected filenames.

Raw sequencing data and processed Seurat / AnnData objects will be deposited at the Genome Sequence Archive for Human (GSA-Human) and released upon thesis publication. Until then, processed objects can be obtained directly from the author for academic reviewers.

## 10. What is **not** included in this repository

- Raw FASTQ files and Cell Ranger / STAR binary outputs (too large; hosted separately).
- Intermediate `.rds` / `.h5ad` / `.loom` objects (excluded via `.gitignore`; regeneratable from the scripts).
- cisTarget feather + motif TBL databases (~1.25 GB; see §9).
- Backup / scratch directories (`backup_prefix_*/`, `backup_fixbugs_*/`, `*.bak`) are all excluded.

## 11. How to cite

If you reuse any component of this repository, please cite the thesis (publication pending, 2026) **and** the two lab packages:

> Ouyang H. *Construction of a single-cell transcriptome atlas of hyperparathyroidism.* PhD thesis, Central South University, 2026.

> scMMR: https://github.com/HUI950319/scMMR
> UtilsR: https://github.com/HUI950319/UtilsR

## 12. Contact

Questions, bug reports, and reviewer reproduction issues are welcome.
Author: **HUI Ouyang**
Email: **ouyanghui950319@gmail.com**
Affiliation: Xiangya School of Medicine, Central South University, Changsha, Hunan, China.

---

*This repository is released for academic non-commercial reproduction. No patient identifying information is contained in any file.*
