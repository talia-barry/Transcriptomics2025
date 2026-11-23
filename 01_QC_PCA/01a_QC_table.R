# ===================================================================
# 01a_QC_table.R
# Purpose: Builds the sequencing QC summary table (Table 1) from Novogene's FastQC/MultiQC reports
# Input:   data/novogene/QC/
# Output:  results/01_QC_PCA/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/01_QC_PCA", recursive = TRUE, showWarnings = FALSE)

####################
#create QC table for results in paper
####novogene data qc report from Box drive
###################
library(dplyr)
library(stringr)
library(tidyr)

# --- Read the QC files (tab-delimited despite .xls extension) ---
fastqc  <- read.delim("data/novogene/QC/FastQC_report_trans2025_novogene.xls",
                      stringsAsFactors = FALSE, check.names = FALSE)
multiqc <- read.delim("data/novogene/QC/MultiQC_report_trans2025_novogene.xls",
                      stringsAsFactors = FALSE, check.names = FALSE)

# --- Helpers to pull numbers out of "12345678(77.56%)" strings ---
extract_count <- function(x) as.numeric(str_extract(x, "^[0-9]+"))
extract_pct   <- function(x) as.numeric(str_extract(x, "(?<=\\()[0-9.]+(?=%\\))"))

multiqc_clean <- multiqc %>%
  transmute(
    sample                     = sample,
    `Total mapped reads`       = extract_count(total_map),
    `Uniquely mapped reads`    = extract_count(unique_map),
    `Uniquely mapping rate (%)`= extract_pct(unique_map)
  )

qc_summary <- fastqc %>%
  select(sample, `Raw reads` = raw_reads, `Clean reads` = clean_reads,
         `Q20 (%)` = Q20, `GC content (%)` = GC_pct) %>%
  left_join(multiqc_clean, by = "sample")

#--- remove the A-D samples
qc_summary <- qc_summary %>%
  filter(str_detect(sample, "^(TE|UE|TC|UC)"))

# --- Transpose to match Table 1 layout (samples as columns, metrics as rows) ---
metric_order <- c("Raw reads", "Clean reads", "Q20 (%)", "GC content (%)",
                  "Total mapped reads", "Uniquely mapped reads", "Uniquely mapping rate (%)")

table1 <- qc_summary %>%
  pivot_longer(-sample, names_to = "Metric", values_to = "value") %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  mutate(Metric = factor(Metric, levels = metric_order)) %>%
  arrange(Metric)

write.csv(table1, "results/01_QC_PCA/Table1_QC_summary.csv", row.names = FALSE)
print(table1, width = Inf)

table1

#########################
#Versions of packages for methods
##############################

packageVersion("dplyr")
packageVersion("DESeq2")
packageVersion("WGCNA")
packageVersion("fgsea")
packageVersion("rrvgo")
packageVersion("ggplot2")
packageVersion("pheatmap")
