# ===================================================================
# 02c_volcano_plots.R
# Purpose: Volcano plots for all comparisons including the interaction term
# Input:   Novogene master file + 02a output
# Output:  results/02_DEGs/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/02_DEGs", recursive = TRUE, showWarnings = FALSE)

##############################################################
# Volcano Plots - all comparisons including interaction term
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)   # for combining panels into one figure

# ---- SIGNIFICANCE THRESHOLDS (must match your summary table) ----
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.5

##############################################################
# STEP 1: Load data (reuses same loading/cleaning logic as your
# summary table script - master file + interaction file)
##############################################################

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                     header = TRUE, sep = "\t", check.names = FALSE,
                     row.names = "gene_id",
                      quote = "")
dim(master) 

cols_to_drop <- grep("^(A|B|C|D)_[0-9]+_(count|fpkm)$|^(DvsC|DvsB|CvsB|DvsA|CvsA|BvsA)_",
                     colnames(master), value = TRUE)
master <- master[, !(colnames(master) %in% cols_to_drop)]

relevant_count_cols <- c("TC1_count","TC2_count","TC3_count","TC4_count",
                         "TE1_count","TE2_count","TE3_count","TE4_count",
                         "UC1_count","UC2_count","UC3_count","UC4_count",
                         "UE1_count","UE2_count","UE3_count","UE4_count")
has_signal <- rowSums(master[, relevant_count_cols]) > 0
master <- master[has_signal, ]

interaction_res <- read.csv("results/02_DEGs/interaction_term_all_results_novogene.csv", quote = "")

dim(interaction_res)

##############################################################
# STEP 2: Helper function - classify each gene and build a plot
##############################################################

make_volcano <- function(df, lfc_col, padj_col, title) {
  
  # Build a clean plotting dataframe
  plot_df <- data.frame(
    log2FC = df[[lfc_col]],
    padj   = df[[padj_col]]
  )
  
  # Drop rows with NA (genes excluded by independent filtering)
  plot_df <- plot_df[!is.na(plot_df$padj) & !is.na(plot_df$log2FC), ]
  
  # Classify each gene
  plot_df$category <- "NS"
  plot_df$category[plot_df$padj < PADJ_CUTOFF & plot_df$log2FC >  LFC_CUTOFF] <- "Up"
  plot_df$category[plot_df$padj < PADJ_CUTOFF & plot_df$log2FC < -LFC_CUTOFF] <- "Down"
  plot_df$category <- factor(plot_df$category, levels = c("Up", "Down", "NS"))
  
  # -log10(padj) for y-axis, capped to avoid Inf for padj = 0
  plot_df$neglog10padj <- -log10(plot_df$padj)
  plot_df$neglog10padj[is.infinite(plot_df$neglog10padj)] <-
    max(plot_df$neglog10padj[is.finite(plot_df$neglog10padj)], na.rm = TRUE) + 1
  
  n_up   <- sum(plot_df$category == "Up")
  n_down <- sum(plot_df$category == "Down")
  
  # Split into NS vs significant so they can have different sizes/styling
  ns_df  <- plot_df[plot_df$category == "NS", ]
  sig_df <- plot_df[plot_df$category != "NS", ]
  
  ggplot(plot_df, aes(x = log2FC, y = neglog10padj, fill = category)) +
    geom_point(shape = 21, color = "navy", stroke = 0.3, size = 2, alpha = 0.95) +
    scale_fill_manual(values = c("Up" = "#00E676",     # replace with your own color
                                 "Down" = "#FF1744",   # replace with your own color
                                 "NS" = "#2979FF"),    # replace with your own color
                      labels = c(paste0("Up (", n_up, ")"),
                                 paste0("Down (", n_down, ")"),
                                 "NS")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
    labs(title = title, x = "log2(Fold Change)", y = expression(-log[10]("adj. p-value")),
         fill = "") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}
##############################################################
# STEP 3: Generate each of the 6 volcano plots
##############################################################

p1 <- make_volcano(master, "UEvsUC_log2FoldChange", "UEvsUC_padj", "UE vs UC")
p2 <- make_volcano(master, "TCvsUC_log2FoldChange", "TCvsUC_padj", "TC vs UC")
p3 <- make_volcano(master, "TEvsTC_log2FoldChange", "TEvsTC_padj", "TE vs TC")
p4 <- make_volcano(master, "TEvsUC_log2FoldChange", "TEvsUC_padj", "TE vs UC")
p5 <- make_volcano(master, "TEvsUE_log2FoldChange", "TEvsUE_padj", "TE vs UE")
p6 <- make_volcano(interaction_res, "log2FoldChange", "padj", "Interaction")

#View the plots in R studio First
p1
p2
p3
p4
p5
p6

##############################################################
# STEP 4: Combine into one multi-panel figure (2 rows x 3 cols)
##############################################################

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title = "Differential Expression Across Comparisons",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
  )

print(combined_plot)

##############################################################
# STEP 5: Save
##############################################################

ggsave("results/02_DEGs/volcano_plots_all_comparisons.png", combined_plot,
       width = 14, height = 9, dpi = 300)

ggsave("results/02_DEGs/volcano_plots_all_comparisons.pdf", combined_plot,
       width = 14, height = 9)

cat("Saved volcano_plots_all_comparisons.png and .pdf\n")
