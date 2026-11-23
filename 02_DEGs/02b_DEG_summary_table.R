# ===================================================================
# 02b_DEG_summary_table.R
# Purpose: DEG summary table for the 5 pairwise comparisons (Novogene DESeq2) + interaction term; DEG = padj < 0.05 & |log2FC| > 1.5
# Input:   Novogene master file + 02a output
# Output:  results/02_DEGs/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/02_DEGs", recursive = TRUE, showWarnings = FALSE)

##############################################################
# Summary DEG Table - all comparisons including interaction term
# Trans25 - Eastern oyster predator exposure RNA-seq
# VERSION 2: uses the single master wide-format Novogene file
# (one row per gene, all comparisons as columns) instead of
# separate pre-filtered "deg_all" files.
##############################################################

library(dplyr)

# ---- SET WORKING DIRECTORY ----
# (no setwd(): run from the RStudio Project at the repository root)

# ---- SIGNIFICANCE THRESHOLDS ----
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.5

##############################################################
# STEP 1: Load the master file
##############################################################
# This file has gene_id, count columns, fpkm columns, and then
# per-comparison blocks named like:
#   TEvsUE_TE, TEvsUE_UE, TEvsUE_log2FoldChange, TEvsUE_pvalue, TEvsUE_padj
# for every comparison Novogene ran (including unrelated A/B/C/D groups,
# which we simply ignore by not referencing their columns).

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                     header = TRUE, sep = "\t", check.names = FALSE,
                     row.names = "gene_id",
                     quote = "")

cat("Total genes in master file (before subsetting):", nrow(master), "\n")

##############################################################
# STEP 1b: Drop all A/B/C/D group columns (unrelated samples)
##############################################################
# This removes their count/fpkm columns AND any comparison columns
# involving A/B/C/D (DvsC, DvsB, CvsB, DvsA, CvsA, BvsA), so none of
# that data can influence gene totals or filtering below.

cols_to_drop <- grep("^(A|B|C|D)_[0-9]+_(count|fpkm)$|^(DvsC|DvsB|CvsB|DvsA|CvsA|BvsA)_", 
                     colnames(master), value = TRUE)

cat("Dropping", length(cols_to_drop), "columns related to A/B/C/D groups\n")

master <- master[, !(colnames(master) %in% cols_to_drop)]

cat("Remaining columns:", ncol(master), "\n")

##############################################################
# STEP 1c: Drop genes with zero counts across ALL TC/TE/UC/UE samples
##############################################################
# These would be genes only detected in the A/B/C/D samples - they
# shouldn't count toward your gene universe at all. This is a belt-
# and-suspenders check on top of dropping the A-D columns themselves,
# since padj is already NA for these in each comparison, but this
# keeps the master gene set itself clean too.

relevant_count_cols <- c("TC1_count","TC2_count","TC3_count","TC4_count",
                         "TE1_count","TE2_count","TE3_count","TE4_count",
                         "UC1_count","UC2_count","UC3_count","UC4_count",
                         "UE1_count","UE2_count","UE3_count","UE4_count")

has_signal <- rowSums(master[, relevant_count_cols]) > 0

cat("Dropping", sum(!has_signal), "genes with zero counts across all TC/TE/UC/UE samples\n")

master <- master[has_signal, ]

cat("Genes remaining after filtering:", nrow(master), "\n")

##############################################################
# STEP 2: Map each of your 6 comparisons to its column prefix
##############################################################
# Convention confirmed via volcano plot: "AvsB" = A relative to B,
# positive log2FoldChange = higher in A (first-listed group).
# All 5 Novogene comparisons below already have the more-exposed
# group listed first, so positive = higher in the more-exposed
# condition throughout - no sign flips needed.

comparisons <- list(
  "UE vs UC"    = "UEvsUC",   # positive = higher in UE
  "TC vs UC"    = "TCvsUC",   # positive = higher in TC
  "TE vs TC"    = "TEvsTC",   # positive = higher in TE
  "TE vs UC"    = "TEvsUC",   # positive = higher in TE
  "TE vs UE"    = "TEvsUE"    # positive = higher in TE
)

##############################################################
# STEP 3: Loop through the 5 Novogene comparisons
##############################################################

summary_list <- list()

for (comp_name in names(comparisons)) {
  
  prefix   <- comparisons[[comp_name]]
  lfc_col  <- paste0(prefix, "_log2FoldChange")
  padj_col <- paste0(prefix, "_padj")
  
  # Total genes tested = rows with a non-NA padj for this comparison
  total_tested <- sum(!is.na(master[[padj_col]]))
  
  is_sig <- !is.na(master[[padj_col]]) &
    master[[padj_col]] < PADJ_CUTOFF &
    abs(master[[lfc_col]]) > LFC_CUTOFF
  
  n_sig  <- sum(is_sig)
  n_up   <- sum(is_sig & master[[lfc_col]] > 0)
  n_down <- sum(is_sig & master[[lfc_col]] < 0)
  
  pct_transcriptome <- round(100 * n_sig / total_tested, 3)
  
  summary_list[[comp_name]] <- data.frame(
    Comparison         = comp_name,
    Total_DEGs         = n_sig,
    Upregulated        = n_up,
    Downregulated      = n_down,
    Pct_Transcriptome  = pct_transcriptome,
    Genes_Tested       = total_tested
  )
}

##############################################################
# STEP 4: Add the interaction term (separate file, your own DESeq2 run)
##############################################################

interaction_res <- read.csv("results/02_DEGs/interaction_term_all_results_novogene.csv")

# Adjust these column names if your interaction file's headers differ
int_lfc_col  <- "log2FoldChange"
int_padj_col <- "padj"

int_total_tested <- sum(!is.na(interaction_res[[int_padj_col]]))

int_is_sig <- !is.na(interaction_res[[int_padj_col]]) &
  interaction_res[[int_padj_col]] < PADJ_CUTOFF &
  abs(interaction_res[[int_lfc_col]]) > LFC_CUTOFF

int_n_sig  <- sum(int_is_sig)
int_n_up   <- sum(int_is_sig & interaction_res[[int_lfc_col]] > 0)
int_n_down <- sum(int_is_sig & interaction_res[[int_lfc_col]] < 0)

int_pct <- round(100 * int_n_sig / int_total_tested, 3)

summary_list[["Interaction"]] <- data.frame(
  Comparison         = "Interaction",
  Total_DEGs         = int_n_sig,
  Upregulated        = int_n_up,
  Downregulated      = int_n_down,
  Pct_Transcriptome  = int_pct,
  Genes_Tested       = int_total_tested
)

##############################################################
# STEP 5: Combine into one summary table
##############################################################

summary_table <- bind_rows(summary_list)

print(summary_table)


##trouble shooting
dim(master)
class(master$UEvsUC_padj)
class(master$UEvsUC_log2FoldChange)

##############################################################
# STEP 6: Save
##############################################################

write.csv(summary_table, "results/02_DEGs/summary_DEG_table_all_comparisons_v2.csv", row.names = FALSE)

cat("Summary table saved to summary_DEG_table_all_comparisons_v2.csv\n")

