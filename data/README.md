# Input data

Novogene files below are included, trimmed to the 16 samples in this study
(TC, TE, UC, UE). Columns for unrelated samples sequenced in the same batch
were removed; all values for the 16 samples are unchanged.

The Zeng & Guo (2022) DEG table is **not** included, because it belongs to that
publication. To run `05a_compare_Guo2022.R`, download it from the paper's
electronic supplementary material:

> Zeng, D. & Guo, X. (2022). Mantle transcriptome provides insights into
> biomineralization and growth regulation in the eastern oyster
> (*Crassostrea virginica*). *Marine Biotechnology* 24, 82–96.
> https://doi.org/10.1007/s10126-021-10088-7

Save their **Table S2** ("Differentially expressed genes between large and small eastern oysters") as tab-delimited text at `data/Guo2022/DEGs_Zeng_Guo_2022.txt`.
The script skips the first (title) line and expects these columns:
`Gene_ID`, `log2FoldChange`, `pvalue`, `padj`, `significant(LvsS)`,
`NCBI annotation`, `Blast swiss prot`, `Potentially shell related`.

`phenotype/shell_strength_5month_2.5mm.csv` holds the shell measurements for the
180 oysters (5 months old, 2.5 mm start size; 90 control, 90 predator-treated)
from which the transcriptomics oysters were drawn. It is the subset of the
Chapter 3 master data sheet used for Table 1 and Figure 4.

Layout:

```
data/
├── phenotype/
│   └── shell_strength_5month_2.5mm.csv      # shell strength, size and dry weight (Table 1, Figure 4)
├── novogene/
│   ├── gene_count_novogene_2025.txt          # raw read counts per gene (Novogene)
│   ├── MasterFile_novogene_DEGs_2025.xls     # counts, FPKM and DESeq2 results for all comparisons (Novogene)
│   ├── QC/
│   │   ├── FastQC_report_trans2025_novogene.xls
│   │   └── MultiQC_report_trans2025_novogene.xls
│   ├── GO_enrichment/                         # Novogene GO over-representation tests on their DEG lists
│   │   └── TCvsUC / TEvsTC / TEvsUC / UEvsUC .all_GOenrich.xls  (no file was supplied for TEvsUE)
│   └── GSEA/
│       ├── TCvsUC/  gsea_report_for_TCvsUC_TC.tsv, gsea_report_for_TCvsUC_UC.tsv
│       ├── TEvsTC/  gsea_report_for_TEvsTC_TE.tsv, gsea_report_for_TEvsTC_TC.tsv
│       ├── TEvsUC/  gsea_report_for_TEvsUC_TE.tsv, gsea_report_for_TEvsUC_UC.tsv
│       ├── TEvsUE/  gsea_report_for_TEvsUE_TE.tsv, gsea_report_for_TEvsUE_UE.tsv
│       └── UEvsUC/  gsea_report_for_UEvsUC_UE.tsv, gsea_report_for_UEvsUC_UC.tsv
└── Guo2022/
    └── DEGs_Zeng_Guo_2022.txt                # NOT included - see above
```

Raw sequencing reads: NCBI SRA BioProject [ACCESSION – add when available].
