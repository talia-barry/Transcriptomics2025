# ===================================================================
# 02d_venn_diagram.R
# Purpose: Venn diagram of DEG overlap among comparisons + direction table
# Input:   Novogene master file
# Output:  results/02_DEGs/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/02_DEGs", recursive = TRUE, showWarnings = FALSE)

##############################################################
# DEG Overlap Venn Diagram + Direction Companion Table
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

library(dplyr)
library(ggVennDiagram)
library(ggplot2)

# ---- SIGNIFICANCE THRESHOLDS (must match your summary table) ----
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.5

##############################################################
# STEP 1: Load and clean master file (same steps as summary table script)
##############################################################

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                     header = TRUE, sep = "\t", check.names = FALSE,
                     row.names = "gene_id",
                     quote="")

cols_to_drop <- grep("^(A|B|C|D)_[0-9]+_(count|fpkm)$|^(DvsC|DvsB|CvsB|DvsA|CvsA|BvsA)_",
                     colnames(master), value = TRUE)
master <- master[, !(colnames(master) %in% cols_to_drop)]

dim(master)

relevant_count_cols <- c("TC1_count","TC2_count","TC3_count","TC4_count",
                         "TE1_count","TE2_count","TE3_count","TE4_count",
                         "UC1_count","UC2_count","UC3_count","UC4_count",
                         "UE1_count","UE2_count","UE3_count","UE4_count")
has_signal <- rowSums(master[, relevant_count_cols]) > 0
master <- master[has_signal, ]

master$gene_id <- rownames(master)

##############################################################
# STEP 2: Define the 4 core comparisons and their column prefixes
##############################################################
# Positive log2FC = higher in the more-exposed group (confirmed convention)

comparisons <- list(
  "UC vs UE" = "UEvsUC",   # positive = higher in UE
  "UC vs TC" = "TCvsUC",   # positive = higher in TC
  "TC vs TE" = "TEvsTC",   # positive = higher in TE
  "UE vs TE" = "TEvsUE",   # positive = higher in TE
  "UC vs TE" = "TEvsUC"    # positivie = higher in TE
)

##############################################################
# STEP 3: Extract significant gene ID sets for each comparison
##############################################################

sig_sets <- list()
sig_data <- list()   # keep full rows (gene_id, log2FC, padj) for the companion table

for (comp_name in names(comparisons)) {
  
  prefix   <- comparisons[[comp_name]]
  lfc_col  <- paste0(prefix, "_log2FoldChange")
  padj_col <- paste0(prefix, "_padj")
  
  is_sig <- !is.na(master[[padj_col]]) &
    master[[padj_col]] < PADJ_CUTOFF &
    abs(master[[lfc_col]]) > LFC_CUTOFF
  
  sig_genes <- master$gene_id[is_sig]
  
  sig_sets[[comp_name]] <- sig_genes
  
  sig_data[[comp_name]] <- data.frame(
    gene_id = sig_genes,
    log2FC  = master[[lfc_col]][is_sig],
    padj    = master[[padj_col]][is_sig]
  )
  
  cat(comp_name, "-", length(sig_genes), "significant genes\n")
}

##############################################################
# STEP 4: Build the 4-way Venn diagram
##############################################################

venn_plot <- ggVennDiagram(sig_sets, label = "count") +
  scale_fill_gradient(low = "white", high = "#4575B4") +
  labs(title = "DEG Overlap Across Core Comparisons") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))

print(venn_plot)

ggsave("results/02_DEGs/DEG_overlap_venn.png", venn_plot, width = 8, height = 8, dpi = 300)
ggsave("results/02_DEGs/DEG_overlap_venn.pdf", venn_plot, width = 8, height = 8)

cat("\nSaved DEG_overlap_venn.png / .pdf\n")

##############################################################
# STEP 5: Build companion table - genes appearing in 2+ comparisons,
# with their log2FC/direction in EACH comparison they're significant in
##############################################################

# Get every unique significant gene across all 4 comparisons
all_sig_genes <- unique(unlist(sig_sets))

# Count how many of the 4 comparisons each gene is significant in
gene_counts <- sapply(all_sig_genes, function(g) {
  sum(sapply(sig_sets, function(s) g %in% s))
})

shared_genes <- names(gene_counts[gene_counts >= 2])

cat("\n", length(shared_genes), "genes are significant in 2 or more comparisons\n")

if (length(shared_genes) > 0) {
  
  # Build a wide table: one row per gene, one log2FC column per comparison
  companion_table <- data.frame(gene_id = shared_genes)
  
  for (comp_name in names(comparisons)) {
    lookup <- setNames(sig_data[[comp_name]]$log2FC, sig_data[[comp_name]]$gene_id)
    companion_table[[paste0(comp_name, "_log2FC")]] <- lookup[shared_genes]
  }
  
  # Pull gene name/description from master for readability
  companion_table$gene_name <- master$gene_name[match(shared_genes, master$gene_id)]
  companion_table$gene_description <- master$gene_description[match(shared_genes, master$gene_id)]
  
  print(companion_table)
  
  write.csv(companion_table, "results/02_DEGs/DEG_overlap_companion_table.csv", row.names = FALSE)
  cat("\nSaved DEG_overlap_companion_table.csv\n")
  
} else {
  cat("\nNo genes overlap between any of the 4 comparisons - nothing to save.\n")
}


