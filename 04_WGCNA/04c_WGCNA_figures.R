# ===================================================================
# 04c_WGCNA_figures.R
# Purpose: WGCNA figures and module summary table (run locally from 04b outputs)
# Input:   results/04_WGCNA/WGCNA3/
# Output:  results/04_WGCNA/WGCNA3/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/04_WGCNA/WGCNA3", recursive = TRUE, showWarnings = FALSE)

##############################################################
# WGCNA - Local plotting/summary script
# Reads the .RData + CSVs pulled down from Amarel - no
# recomputation, so this is fast and light on memory.
##############################################################

library(WGCNA)
library(dplyr)
library(ggplot2)

# ---- Set this to wherever you saved the files locally ----
results_dir <- "results/04_WGCNA/WGCNA3"

load(file.path(results_dir, "WGCNA_results3.RData"))
# brings in: net, MEs, moduleColors, gene_module_map, kME_all,
#            coldata, traits, moduleTraitCor, moduleTraitPvalue

combined <- read.csv(file.path(results_dir, "WGCNA_significant_modules_combined3.csv"),
                     colClasses = c(gene_id = "character"))
all_genes <- read.csv(file.path(results_dir, "WGCNA_all_gene_module_assignments3.csv"),
                      colClasses = c(gene_id = "character"))

sig_modules <- unique(paste0("ME", combined$module_number))

##############################################################
# FIGURE 1: Legible module-trait heatmap (significant modules only)
##############################################################

plotCor <- moduleTraitCor[sig_modules, , drop = FALSE]
plotPvalue <- moduleTraitPvalue[sig_modules, , drop = FALSE]

textMatrix <- paste0(signif(plotCor, 2), "\n(", signif(plotPvalue, 1), ")")
dim(textMatrix) <- dim(plotCor)

n_rows <- nrow(plotCor)
label_cex <- if (n_rows > 25) 0.6 else if (n_rows > 12) 0.75 else 0.95

pdf(file.path(results_dir, "FIG_module_trait_heatmap.pdf"), width = 8, height = max(5, 1 + n_rows * 0.4))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = plotCor,
               xLabels = colnames(plotCor),
               yLabels = rownames(plotCor),
               ySymbols = rownames(plotCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = label_cex,
               cex.lab.y = label_cex,
               main = "Module-trait relationships (significant modules only)")
dev.off()


##############################################################
# FIGURE 2: Module size barplot (all modules, significant ones
# highlighted) - shows how big/small each module is at a glance
##############################################################

module_sizes <- as.data.frame(table(net$colors))
colnames(module_sizes) <- c("module_number", "size")
module_sizes$module_number <- as.numeric(as.character(module_sizes$module_number))
module_sizes$is_significant <- paste0("ME", module_sizes$module_number) %in% sig_modules
module_sizes <- module_sizes %>% arrange(desc(size))
module_sizes$module_number <- factor(module_sizes$module_number, levels = module_sizes$module_number)

p_sizes <- ggplot(module_sizes, aes(x = module_number, y = size, fill = is_significant)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
                    name = NULL) +
  labs(title = "Module sizes", x = "Module", y = "Number of genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
p_sizes

ggsave(file.path(results_dir, "FIG_module_sizes.pdf"), p_sizes, width = 12, height = 5)

##############################################################
# FIGURE 3: Eigengene relatedness (module-module relationships)
# Shows which modules cluster together / are redundant with each
# other - useful for describing overall network structure.
##############################################################

pdf(file.path(results_dir, "FIG_eigengene_network.pdf"), width = 8, height = 8)
plotEigengeneNetworks(MEs, "Eigengene network", marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(5, 4, 1, 2), plotDendrograms = TRUE,
                      xLabelsAngle = 90)
dev.off()

##############################################################
# FIGURE 4: Gene significance vs module membership scatterplots
# Classic WGCNA validation figure - for each significant module,
# plots each gene's kME (module membership) against its individual
# correlation with the relevant trait. A strong positive relationship
# supports that the module is biologically coherent, not noise.
##############################################################

gs_mm_dir <- file.path(results_dir, "FIG_gene_significance_vs_MM")
if (!dir.exists(gs_mm_dir)) dir.create(gs_mm_dir)

expr_present <- exists("expr_matrix")  # only possible if expr_matrix was also saved

if (expr_present) {
  for (me_name in sig_modules) {
    mod_num <- as.numeric(gsub("ME", "", me_name))
    mod_color <- unique(gene_module_map$module_color[gene_module_map$module_number == mod_num])
    mod_genes <- gene_module_map$gene_id[gene_module_map$module_number == mod_num]
    
    # Trait this module was significant for (pick the strongest one)
    best_trait <- colnames(moduleTraitCor)[which.max(abs(moduleTraitCor[me_name, ]))]
    
    gene_sig <- as.numeric(cor(expr_matrix[, mod_genes], traits[[best_trait]], use = "p"))
    module_membership <- kME_all[mod_genes, paste0("kME", mod_num)]
    
    df <- data.frame(MM = module_membership, GS = gene_sig)
    
    p <- ggplot(df, aes(x = MM, y = GS)) +
      geom_point(alpha = 0.5, color = mod_color) +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
      labs(title = paste0("Module ", mod_num, " (", mod_color, ") vs ", best_trait),
           x = "Module Membership (kME)", y = paste0("Gene Significance (", best_trait, ")")) +
      theme_minimal()
    
    ggsave(file.path(gs_mm_dir, paste0("MM_vs_GS_module_", mod_num, ".pdf")), p, width = 5, height = 5)
  }
  cat("Saved gene significance vs module membership plots to", gs_mm_dir, "\n")
} else {
  cat("expr_matrix not found in .RData - skipping Figure 4.",
      "If you want this plot, add expr_matrix to the save() call on Amarel and re-download.\n")
}

p
##############################################################
# TABLE 1: Significant modules summary - one row per module,
# ready to drop into a thesis results table
##############################################################

summary_rows <- lapply(sig_modules, function(me_name) {
  mod_num <- as.numeric(gsub("ME", "", me_name))
  mod_color <- unique(gene_module_map$module_color[gene_module_map$module_number == mod_num])
  size <- sum(gene_module_map$module_number == mod_num)
  
  best_trait <- colnames(moduleTraitCor)[which.max(abs(moduleTraitCor[me_name, ]))]
  best_r <- moduleTraitCor[me_name, best_trait]
  best_p <- moduleTraitPvalue[me_name, best_trait]
  
  top_hub <- combined %>%
    filter(module_number == mod_num) %>%
    arrange(desc(abs(kME))) %>%
    dplyr::slice(1)
  
  data.frame(module_number = mod_num,
             module_color = mod_color,
             n_genes = size,
             strongest_trait = best_trait,
             r = round(best_r, 3),
             p = round(best_p, 4),
             top_hub_gene = top_hub$gene_id,
             top_hub_kME = round(top_hub$kME, 3))
})

summary_table <- bind_rows(summary_rows) %>% arrange(desc(abs(r)))

write.csv(summary_table, file.path(results_dir, "TABLE_significant_modules_summary.csv"), row.names = FALSE)
cat("\nSaved TABLE_significant_modules_summary.csv\n")

##############################################################
# NOTE: figures already generated on Amarel and downloaded
# (reuse as-is, no need to remake locally):
#   - sample_clustering2.png  (outlier check dendrogram)
#   - module_dendrogram2.png  (gene clustering + module colors)
#   - module_eigengene_sample_heatmap2.png (ME across samples)
##############################################################

cat("\nDone. New outputs in:", results_dir, "\n")
cat("FIG_module_trait_heatmap.pdf, FIG_module_sizes.pdf,\n")
cat("FIG_eigengene_network.pdf, FIG_gene_significance_vs_MM/ (if expr_matrix available),\n")
cat("TABLE_significant_modules_summary.csv\n")

##########################################################
#Check number of modules
################################################

length(unique(net$colors))                 # total, including grey
length(setdiff(unique(net$colors), 0))     # total, excluding grey
length(sig_modules)

