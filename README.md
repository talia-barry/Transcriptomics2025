# Transcriptomic response of the eastern oyster (*Crassostrea virginica*) to blue crab predators

Analysis code for:

> Barry, T. *et al.* [Title – add when published]. *[Journal]*, [year]. [DOI]

## Experimental design

Nested factorial design with 16 oysters (n = 4 per group):

| Factor | Levels |
|---|---|
| Long-term exposure (6 weeks) | **T** = predator-exposed, **U** = unexposed |
| Short-term re-exposure (24 h) | **E** = exposed, **C** = control |

Groups: **TE**, **TC**, **UE**, **UC**. Comparisons are always named treated group first (e.g. TE vs UE).

## Sequencing, mapping and differential expression (Novogene)

RNA extraction, library preparation, sequencing (Illumina, paired-end), quality control,
read mapping to the NCBI *C. virginica* genome (C_virginica-3.0, 2017), gene counting,
DESeq2 differential expression and GSEA for the five pairwise comparisons were performed
by Novogene Co., Ltd. These steps are not included as code here. The interaction-term
DESeq2 model and its GSEA (preranked fgsea) were run by the author (`02a`, `03b`).

All analyses below use Novogene's gene counts and results on the 2017 NCBI genome.
DEGs are defined as padj < 0.05 and |log2 fold change| > 1.5.

## Repository layout

Scripts are numbered in the order they should be run.

| Folder | Script | What it does |
|---|---|---|
| `01_QC_PCA` | `01a_QC_table.R` | Sequencing QC summary table |
| | `01b_PCA.R` | PCA of samples |
| `02_DEGs` | `02a_interaction_term_DESeq2.R` | DESeq2 interaction model (long-term × short-term) |
| | `02b_DEG_summary_table.R` | DEG counts for all comparisons |
| | `02c_volcano_plots.R` | Volcano plots |
| | `02d_venn_diagram.R` | DEG overlap Venn diagram |
| `03_GSEA` | `03a_rrvgo_pairwise_GSEA.R` | rrvgo clustering of GSEA GO terms (pairwise comparisons) |
| | `03b_fGSEA_rrvgo_interaction.R` | fgsea + rrvgo for the interaction term |
| `04_WGCNA` | `04a_soft_power_check.R` | Soft-threshold power selection |
| | `04b_WGCNA_signed_network.R` + `04b_run_WGCNA_amarel.sh` | Signed WGCNA network (run on an HPC cluster) |
| | `04c_WGCNA_figures.R` | WGCNA figures and module summary |
| | `04d_WGCNA_GO_enrichment.R` | GO enrichment of significant modules |
| `05_Guo2022_comparison` | `05a_compare_Guo2022.R` | Comparison with Zeng & Guo (2022) DEGs |

Outputs are written to `results/` (created automatically).

## How to run

1. Clone this repository and open `Transcriptomics2025.Rproj` in RStudio.
   All paths are relative to the repository root, so no `setwd()` is needed.
2. Place the input files in `data/` as described in [`data/README.md`](data/README.md).
3. Run the scripts in numbered order.
   `04b` was run on the Rutgers Amarel cluster with SLURM (`sbatch 04_WGCNA/04b_run_WGCNA_amarel.sh`);
   copy its outputs to `results/04_WGCNA/WGCNA3/` before running `04c` and `04d`.

## Software

R 4.5 with Bioconductor packages: DESeq2, WGCNA, fgsea, rrvgo, clusterProfiler,
AnnotationHub, AnnotationDbi, GO.db, plus dplyr, stringr, ggplot2, ggpubr, ggrepel,
ggVennDiagram, patchwork, treemap.

## Data availability

Raw reads: NCBI SRA BioProject [ACCESSION – add when available].

## Contact

Talia Barry – Rutgers University
