# Fhl2_manuscript_code

Analysis code for the manuscript "Four and a half LIM domains protein 2 drives cardiomyocyte reprogramming during exercise-induced cardiac recovery".

| Folder | Content |
|---|---|
| `01_snRNAseq_mouse/` | Mouse snRNA-seq (sham, MI + sedentary, MI + exercise): cell-type annotation, Augur cell-type prioritisation, per-cell-type pseudobulk DESeq2 rescue analysis and Hallmark GSEA, cardiomyocyte and epicardial subclustering, scFates pseudotime and gene modules |
| `02_cellchat/` | CellChat cell–cell communication analysis, exercise versus sedentary |
| `03_snRNAseq_Fhl2OE/` | snRNA-seq of cardiomyocyte-specific Fhl2 overexpression: annotation, cardiomyocyte analysis, gene scores, differential expression and enrichment |
| `04_snATACseq/` | snATAC-seq (ArchR): quality control, peak calling, motif enrichment and chromVAR motif activity |
| `05_human_LVRR_snRNAseq/` | Human LVAD cohort: cardiomyocyte selection and per-donor pseudobulk DESeq2 test of FHL2 by LVRR status |
| `06_human_public_cohorts/` | FHL2 in public human cohorts (MAGNet, Kuppe et al. 2022, Reichart et al. 2022): pseudobulk DESeq2 and box plots |

Notebooks are provided without image outputs. Package versions for the human analyses are listed in the `requirements_*.txt` files.

Sequencing data: GEO GSE346720.
