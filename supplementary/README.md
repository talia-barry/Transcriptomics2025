# Supplementary materials

Supplementary tables and figures for:

> Barry, T. *et al.* [Title – add when published].

All files were produced by the scripts in this repository from the final analysis
(NCBI *C. virginica* genome C_virginica-3.0; Novogene counts; final signed WGCNA network).
Comparisons are named treated group first; a positive log2 fold change means higher
expression in the first-named group.

## Tables

| File | Description | Made by |
|---|---|---|
| `tables/Table_S1_sequencing_QC.csv` | **Table S1.** Sequencing and mapping statistics for each library: raw and clean reads, Q20, Q30, GC content, total and uniquely mapped reads and rates. | `01_QC_PCA/01a_QC_table.R` |
| `tables/Table_S2_DEG_results_all_genes.csv` | **Table S2.** Differential expression results for all 39,529 genes: log2 fold change, p-value and adjusted p-value for the five pairwise comparisons (UE vs UC, TC vs UC, TE vs TC, TE vs UC, TE vs UE; DESeq2 by Novogene) and the long-term × short-term interaction term (DESeq2, this study). DEGs: padj < 0.05 and \|log2FC\| > 1.5. | Novogene master file + `02_DEGs/02a_interaction_term_DESeq2.R` |
| `tables/Table_S3_interaction_GSEA_significant_terms.csv` | **Table S3.** All 130 GO gene sets significant (FDR < 0.25) in pre-ranked GSEA (fgsea) of the interaction term, with rrvgo parent term, NES, gene-set size (genes in the ranked list) and leading-edge genes. | `03_GSEA/03b_fGSEA_rrvgo_interaction.R` |
| `tables/Table_S4a_Zeng_Guo_2022_comparison_all_genes.csv` | **Table S4a.** The 104 DEGs from Zeng & Guo (2022; large vs small oysters) compared with each comparison in this study: log2FC, padj, significance and direction concordance. | `05_Guo2022_comparison/05a_compare_Guo2022.R` |
| `tables/Table_S4b_Zeng_Guo_2022_shell_related_genes.csv` | **Table S4b.** Subset of Table S4a flagged as potentially shell-related by Zeng & Guo (2022). | `05_Guo2022_comparison/05a_compare_Guo2022.R` |
| `tables/Table_S5_WGCNA_module_assignments.csv` | **Table S5.** Module assignment (number and colour) for all 26,172 genes in the signed WGCNA network. | `04_WGCNA/04b_WGCNA_signed_network.R` |

## Figures

| File | Description |
|---|---|
| `figures/Figure_S1_WGCNA_soft_threshold.png` | **Figure S1.** Scale-free topology fit (R²) and mean connectivity across soft-thresholding powers 1–30; power 11 was selected (R² = 0.861). |
| `figures/Figure_S2_WGCNA_sample_clustering.png` | **Figure S2.** Hierarchical clustering of the 16 samples used to check for outliers. |
| `figures/Figure_S3_WGCNA_module_dendrogram.png` | **Figure S3.** Gene dendrogram and module colour assignment. |
| `figures/Figure_S4_WGCNA_module_sizes.pdf` | **Figure S4.** Number of genes in each module. |
| `figures/Figure_S5a–e_…pdf` | **Figure S5a–e.** Module membership vs gene significance for the five significant modules: (a) ME12/tan, (b) ME13/salmon, (c) ME19/lightyellow, (d) ME27/white, (e) ME29/saddlebrown. |

Made by `04_WGCNA/04b_WGCNA_signed_network.R` and `04_WGCNA/04c_WGCNA_figures.R`.
