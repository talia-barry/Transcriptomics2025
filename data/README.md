# Input data

Novogene files below are included, trimmed to the 16 samples in this study
(TC, TE, UC, UE). Columns for unrelated samples sequenced in the same batch
were removed; all values for the 16 samples are unchanged.

The Zeng & Guo (2022) DEG table is **not** included. Download it from the
supplementary material of that paper [add citation/link] and save it as
`data/Guo2022/DEGs_Zeng_Guo_2022.txt`.

Layout:

```
data/
├── novogene/
│   ├── gene_count_novogene_2025.txt          # raw read counts per gene (Novogene)
│   ├── MasterFile_novogene_DEGs_2025.xls     # counts, FPKM and DESeq2 results for all comparisons (Novogene)
│   ├── QC/
│   │   ├── FastQC_report_trans2025_novogene.xls
│   │   └── MultiQC_report_trans2025_novogene.xls
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
