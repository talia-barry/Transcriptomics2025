# ===================================================================
# 05a_compare_Guo2022.R
# Purpose: Compares DEGs with Zeng & Guo 2022 mantle DEGs (same 2017 genome, direct ID match)
# Input:   data/Guo2022/DEGs_Zeng_Guo_2022.txt + Novogene master file
# Output:  results/05_Guo2022_comparison/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/05_Guo2022_comparison", recursive = TRUE, showWarnings = FALSE)


##############################################################
# Comparison to Guo et al. 2022 Mantle Transcriptomics DEGs
# Trans25 - Eastern oyster predator exposure RNA-seq
#
# Both this analysis and Guo et al. 2022 used the 2017 C. virginica
# genome assembly, so gene IDs can be matched directly (no homology/
# name-based matching needed).
#
# IMPORTANT INTERPRETATION NOTE: Guo et al. 2022's comparison is
# L (larger) vs S (smaller) oyster growth phenotype in mantle tissue -
#  this comparison tests whether GROWTH-associated genes
# overlap with PREDATOR-RESPONSE genes, and in what direction:
#   - SAME direction (e.g. up in larger oysters AND up under acute
#     exposure) may indicate a gene tied to general growth/metabolic
#     status rather than defense specifically.
#   - OPPOSITE direction (e.g. up in larger oysters but down under
#     acute exposure) may support a growth-defense tradeoff, where
#     genes promoting growth are suppressed when resources shift
#     toward defense.
##############################################################
rm(list = ls())

library(dplyr)

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.5

# All outputs (tables AND figures) save here
fig_dir <- "results/05_Guo2022_comparison"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

##############################################################
# STEP 1: Load the Guo 2022 DEG file
##############################################################

guo <- read.delim("data/Guo2022/DEGs_Zeng_Guo_2022.txt", header = TRUE, sep = "\t",
                  check.names = FALSE, skip = 1, quote = "")
cat("guo dimensions:", dim(guo), "\n")  # sanity check

# Standardize column names for clarity
guo <- guo %>%
  rename(gene_id = Gene_ID,
         guo_log2FC = log2FoldChange,
         guo_pvalue = pvalue,
         guo_padj = padj,
         guo_significant = `significant(LvsS)`)

guo$gene_id <- as.character(guo$gene_id)

cat(nrow(guo), "genes in Guo et al. 2022 DEG list\n")
cat(sum(guo$guo_significant %in% c("UP", "DOWN")), "significant in Guo et al. 2022\n")

##############################################################
# STEP 2: Load your master comparisons file (same as previous scripts)
##############################################################

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                     header = TRUE, sep = "\t", check.names = FALSE, quote = "")
cat("master dimensions:", dim(master), "\n")  # sanity check - should be ~39,529 rows

dim(master)
master$gene_id <- as.character(master$gene_id)

##############################################################
# STEP 3: Define your six comparisons and their column prefixes
##############################################################

comparisons <- list(
  "UE_vs_UC" = "UEvsUC",
  "TC_vs_UC" = "TCvsUC",
  "TE_vs_TC" = "TEvsTC",
  "TE_vs_UC" = "TEvsUC",
  "TE_vs_UE" = "TEvsUE"
)

##############################################################
# STEP 4: Join Guo gene list with your master file, pull each
# comparison's log2FC/padj, and flag significance + direction match
##############################################################

result <- guo %>%
  dplyr::select(gene_id, guo_log2FC, guo_pvalue, guo_padj, guo_significant,
                `Potentially shell related`,
                `NCBI annotation`, `Blast swiss prot`) %>%
  left_join(master %>% dplyr::select(gene_id), by = "gene_id")

result$in_your_dataset <- result$gene_id %in% master$gene_id

for (comp_name in names(comparisons)) {
  
  prefix   <- comparisons[[comp_name]]
  lfc_col  <- paste0(prefix, "_log2FoldChange")
  padj_col <- paste0(prefix, "_padj")
  
  lookup_lfc  <- setNames(master[[lfc_col]], master$gene_id)
  lookup_padj <- setNames(master[[padj_col]], master$gene_id)
  
  result[[paste0(comp_name, "_log2FC")]] <- lookup_lfc[result$gene_id]
  result[[paste0(comp_name, "_padj")]]   <- lookup_padj[result$gene_id]
  
  result[[paste0(comp_name, "_significant")]] <-
    !is.na(result[[paste0(comp_name, "_padj")]]) &
    result[[paste0(comp_name, "_padj")]] < PADJ_CUTOFF &
    abs(result[[paste0(comp_name, "_log2FC")]]) > LFC_CUTOFF
  
  guo_direction <- ifelse(result$guo_log2FC > 0, "UP", "DOWN")
  your_direction <- ifelse(result[[paste0(comp_name, "_log2FC")]] > 0, "UP",
                           ifelse(result[[paste0(comp_name, "_log2FC")]] < 0, "DOWN", NA))
  
  result[[paste0(comp_name, "_concordant")]] <-
    ifelse(is.na(your_direction), NA, guo_direction == your_direction)
}

##############################################################
# STEP 5: Summary columns
##############################################################

sig_cols <- paste0(names(comparisons), "_significant")
concordant_cols <- paste0(names(comparisons), "_concordant")

result$n_significant_in_your_data <- rowSums(result[sig_cols], na.rm = TRUE)
result$n_concordant_direction <- rowSums(result[concordant_cols], na.rm = TRUE)
result$n_comparisons_tested <- rowSums(!is.na(result[concordant_cols]))

##############################################################
# STEP 6: Save full results
##############################################################

result <- result %>% arrange(desc(n_significant_in_your_data), desc(guo_significant == "UP" | guo_significant == "DOWN"))

write.csv(result, file.path(fig_dir, "Guo2022_comparison_full.csv"), row.names = FALSE)

cat("\nFull comparison table saved to Guo2022_comparison_full.csv\n")
cat(nrow(result), "total genes from Guo et al. 2022 processed\n")
cat(sum(result$in_your_dataset), "of those genes are present in your dataset\n")
cat(sum(result$n_significant_in_your_data > 0),
    "genes are significant in at least one of your comparisons\n")

##############################################################
# STEP 7: Highlight table
##############################################################

highlight_table <- result %>%
  filter(guo_significant %in% c("UP", "DOWN") & n_significant_in_your_data > 0) %>%
  dplyr::select(gene_id, `NCBI annotation`, `Potentially shell related`,
                guo_log2FC, guo_padj, guo_significant,
                all_of(paste0(names(comparisons), "_log2FC")),
                all_of(paste0(names(comparisons), "_padj")),
                all_of(paste0(names(comparisons), "_significant")),
                all_of(paste0(names(comparisons), "_concordant")),
                n_significant_in_your_data, n_concordant_direction, n_comparisons_tested)

write.csv(highlight_table, file.path(fig_dir, "Guo2022_comparison_significant_overlap.csv"), row.names = FALSE)

cat("\nHighlight table (significant in both) saved -", nrow(highlight_table), "genes\n")

##############################################################
# STEP 8: Shell-related genes only
##############################################################

shell_table <- result %>%
  filter(!is.na(`Potentially shell related`) & `Potentially shell related` != "") %>%
  dplyr::select(gene_id, `Potentially shell related`, `NCBI annotation`, `Blast swiss prot`,
                guo_log2FC, guo_padj, guo_significant, in_your_dataset,
                all_of(paste0(names(comparisons), "_log2FC")),
                all_of(paste0(names(comparisons), "_padj")),
                all_of(paste0(names(comparisons), "_significant")),
                all_of(paste0(names(comparisons), "_concordant")),
                n_significant_in_your_data, n_concordant_direction, n_comparisons_tested)

write.csv(shell_table, file.path(fig_dir, "Guo2022_shell_related_genes.csv"), row.names = FALSE)

cat("\nShell-related genes saved -", nrow(shell_table), "total\n")
cat(sum(shell_table$in_your_dataset), "of those are present in your dataset\n")

shell_in_data <- shell_table %>% filter(in_your_dataset)

cat("\nConcordance summary - shell-related genes only (n =", nrow(shell_in_data), "):\n")
for (comp_name in names(comparisons)) {
  col <- paste0(comp_name, "_concordant")
  tested <- shell_in_data[!is.na(shell_in_data[[col]]), ]
  n_concordant <- sum(tested[[col]] == TRUE)
  n_discordant <- sum(tested[[col]] == FALSE)
  pct <- if (nrow(tested) > 0) round(100 * n_concordant / nrow(tested), 1) else NA
  cat(sprintf("  %-12s tested=%d  concordant=%d  discordant=%d  %%concordant=%s%%\n",
              comp_name, nrow(tested), n_concordant, n_discordant, pct))
}

##############################################################
# STEP 9: Binomial tests
##############################################################

general_in_data <- result %>% filter(in_your_dataset)

run_binom_tests <- function(data, comparisons, label) {
  
  cat("\n=== Binomial tests -", label, "(n =", nrow(data), ") ===\n")
  
  test_results <- list()
  
  for (comp_name in names(comparisons)) {
    col <- paste0(comp_name, "_concordant")
    tested <- data[!is.na(data[[col]]), ]
    n_concordant <- sum(tested[[col]] == TRUE)
    n_total <- nrow(tested)
    
    if (n_total == 0) {
      cat(sprintf("  %-12s no genes tested - skipped\n", comp_name))
      next
    }
    
    test <- binom.test(n_concordant, n_total, p = 0.5, alternative = "two.sided")
    
    cat(sprintf("  %-12s concordant=%d/%d (%.1f%%)  p=%.4f\n",
                comp_name, n_concordant, n_total,
                100 * n_concordant / n_total, test$p.value))
    
    test_results[[comp_name]] <- data.frame(
      comparison = comp_name,
      gene_set = label,
      n_tested = n_total,
      n_concordant = n_concordant,
      pct_concordant = round(100 * n_concordant / n_total, 1),
      p_value = round(test$p.value, 4)
    )
  }
  
  bind_rows(test_results)
}

general_binom <- run_binom_tests(general_in_data, comparisons, "General (all overlapping genes)")
shell_binom   <- run_binom_tests(shell_in_data, comparisons, "Shell-related genes")

binom_summary <- bind_rows(general_binom, shell_binom)

write.csv(binom_summary, file.path(fig_dir, "Guo2022_concordance_binomial_tests.csv"), row.names = FALSE)

cat("\n\nBinomial test summary saved to Guo2022_concordance_binomial_tests.csv\n")
print(binom_summary)

##############################################################
# STEP 10: Figures
##############################################################

library(ggplot2)
library(tidyr)

comp_labels <- names(comparisons)

##############################################################
# STEP 10a: Bar chart
##############################################################
# FIX: gene-set labels now compute n dynamically from the actual data
# instead of using hardcoded numbers that can silently go stale.

bar_data <- binom_summary %>%
  mutate(gene_set = factor(gene_set,
                           levels = c("General (all overlapping genes)", "Shell-related genes"),
                           labels = c(paste0("General (n=", nrow(general_in_data), ")"),
                                      paste0("Shell-related (n=", nrow(shell_in_data), ")"))),
         comparison = factor(comparison, levels = comp_labels),
         sig_label = ifelse(p_value < 0.05, "*", ""))

bar_plot <- ggplot(bar_data, aes(x = comparison, y = pct_concordant, fill = gene_set)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sig_label, y = pct_concordant + 4),
            position = position_dodge(width = 0.8), size = 6, fontface = "bold") +
  labs(title = "Directional Concordance with Guo et al. (2022)",
       x = "Comparison", y = "% Concordant with Guo et al. (2022) Direction",
       fill = "Gene Set") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))
bar_plot

ggsave(file.path(fig_dir, "concordance_barchart.png"), bar_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "concordance_barchart.pdf"), bar_plot, width = 8, height = 6)

cat("\nSaved concordance_barchart.png/.pdf\n")

##############################################################
# STEP 10b: Heatmap
##############################################################

heatmap_data <- shell_in_data %>%
  dplyr::select(gene_id, `NCBI annotation`, guo_log2FC,
                all_of(paste0(comp_labels, "_log2FC"))) %>%
  rename_with(~ gsub("_log2FC", "", .x), all_of(paste0(comp_labels, "_log2FC"))) %>%
  rename(`Guo 2022 (L vs S)` = guo_log2FC) %>%
  pivot_longer(cols = c(`Guo 2022 (L vs S)`, all_of(comp_labels)),
               names_to = "comparison", values_to = "log2FC") %>%
  mutate(comparison = factor(comparison, levels = c("Guo 2022 (L vs S)", comp_labels)),
         gene_label = paste0(gene_id, ": ", substr(`NCBI annotation`, 1, 40)))

heatmap_plot <- ggplot(heatmap_data, aes(x = comparison, y = gene_label, fill = log2FC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       na.value = "grey85") +
  labs(title = paste0("Shell-Related Gene Expression (n=", nrow(shell_in_data),
                      "): Guo et al. (2022) vs. Your Comparisons"),
       x = "Comparison", y = "Gene", fill = "log2FC") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8))

heatmap_plot

ggsave(file.path(fig_dir, "shell_genes_heatmap.png"), heatmap_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "shell_genes_heatmap.pdf"), heatmap_plot, width = 10, height = 8)

cat("Saved shell_genes_heatmap.png/.pdf\n")

##############################################################
# STEP 10c: Scatter plots
##############################################################
# Generated for TC vs UC, UE vs UC, and TE vs TC (selection rationale
# established from prior regression/binomial results - re-check after
# this fix in case the strongest comparisons have shifted)

scatter_data <- general_in_data %>%
  mutate(shell_related = ifelse(!is.na(`Potentially shell related`) &
                                  `Potentially shell related` != "",
                                "Shell-related", "Other"))

make_scatter_plot <- function(data, comp_col, comp_label, filename_prefix) {
  
  p <- ggplot(data, aes(x = guo_log2FC, y = .data[[comp_col]], color = shell_related)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Other" = "grey60", "Shell-related" = "#B2182B")) +
    labs(title = paste0(comp_label, " log2FC vs. Guo et al. (2022) log2FC"),
         x = "Guo et al. 2022 log2FC (L vs S)", y = paste0("log2FC (", comp_label, ")"),
         color = "") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(fig_dir, paste0(filename_prefix, "_scatter_vs_Guo2022.png")), p, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_prefix, "_scatter_vs_Guo2022.pdf")), p, width = 8, height = 6)
  
  cat("Saved", filename_prefix, "scatter plot\n")
  
  p
}

scatter_TC_vs_UC <- make_scatter_plot(scatter_data, "TC_vs_UC_log2FC", "TC vs UC", "TCvsUC")
scatter_UE_vs_UC <- make_scatter_plot(scatter_data, "UE_vs_UC_log2FC", "UE vs UC", "UEvsUC")
scatter_TE_vs_TC <- make_scatter_plot(scatter_data, "TE_vs_TC_log2FC", "TE vs TC", "TEvsTC")

scatter_TE_vs_TC
scatter_TC_vs_UC
scatter_UE_vs_UC

cat("\nAll figures saved to:", fig_dir, "\n")

##############################################################
# STEP 11: Regression plots
##############################################################

build_long_data <- function(data, comp_labels) {
  data %>%
    dplyr::select(gene_id, guo_log2FC, all_of(paste0(comp_labels, "_log2FC"))) %>%
    rename_with(~ gsub("_log2FC", "", .x), all_of(paste0(comp_labels, "_log2FC"))) %>%
    pivot_longer(cols = all_of(comp_labels), names_to = "comparison", values_to = "your_log2FC") %>%
    mutate(comparison = factor(comparison, levels = comp_labels)) %>%
    filter(!is.na(your_log2FC) & !is.na(guo_log2FC))
}

get_regression_stats <- function(long_data) {
  long_data %>%
    group_by(comparison) %>%
    summarise(
      n = n(),
      r2 = summary(lm(your_log2FC ~ guo_log2FC))$r.squared,
      slope = coef(lm(your_log2FC ~ guo_log2FC))[2],
      p_value = summary(lm(your_log2FC ~ guo_log2FC))$coefficients[2, 4],
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("R\u00b2 = %.3f\np = %.3f\nn = %d", r2, p_value, n))
}

make_regression_plot <- function(data, comp_labels, title_label, filename_prefix) {
  
  long_data <- build_long_data(data, comp_labels)
  stats_data <- get_regression_stats(long_data)
  
  stats_data <- stats_data %>%
    mutate(x_pos = min(long_data$guo_log2FC, na.rm = TRUE),
           y_pos = max(long_data$your_log2FC, na.rm = TRUE))
  
  p <- ggplot(long_data, aes(x = guo_log2FC, y = your_log2FC)) +
    geom_point(alpha = 0.6, color = "#2166AC", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", fill = "#B2182B", alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_text(data = stats_data, aes(x = x_pos, y = y_pos, label = label),
              hjust = 0, vjust = 1, size = 3, fontface = "italic") +
    facet_wrap(~ comparison, scales = "free") +
    labs(title = paste0(title_label, ": Guo et al. (2022) log2FC vs. log2FC"),
         x = "Guo et al. 2022 log2FC (L vs S)", y = "log2FC") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(fig_dir, paste0(filename_prefix, "_regression.png")), p, width = 12, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_prefix, "_regression.pdf")), p, width = 12, height = 8)
  
  cat("\nSaved", filename_prefix, "regression plot\n")
  print(stats_data %>% dplyr::select(comparison, n, r2, slope, p_value))
  
  p
}

# FIX: title labels now compute n dynamically instead of using hardcoded
# numbers that can silently go stale.

general_regression_plot <- make_regression_plot(general_in_data, comp_labels,
                                                paste0("General Gene Set (n=", nrow(general_in_data), ")"),
                                                "general")
general_regression_plot

shell_regression_plot <- make_regression_plot(shell_in_data, comp_labels,
                                              paste0("Shell-Related Genes (n=", nrow(shell_in_data), ")"),
                                              "shell_related")
shell_regression_plot

cat("\nAll regression plots saved to:", fig_dir, "\n")

##############################################################
# STEP 12: Standalone regression results table
##############################################################

general_stats <- get_regression_stats(build_long_data(general_in_data, comp_labels)) %>%
  mutate(gene_set = paste0("General (n=", nrow(general_in_data), ")"))

shell_stats <- get_regression_stats(build_long_data(shell_in_data, comp_labels)) %>%
  mutate(gene_set = paste0("Shell-related (n=", nrow(shell_in_data), ")"))

regression_summary <- bind_rows(general_stats, shell_stats) %>%
  dplyr::select(comparison, gene_set, n, slope, r2, p_value) %>%
  rename(R2 = r2) %>%
  mutate(slope = round(slope, 3),
         R2 = round(R2, 3),
         p_value = round(p_value, 4)) %>%
  arrange(gene_set, comparison)

write.csv(regression_summary, file.path(fig_dir, "Guo2022_regression_results.csv"), row.names = FALSE)

cat("\nRegression results table saved to Guo2022_regression_results.csv\n")
print(regression_summary)
