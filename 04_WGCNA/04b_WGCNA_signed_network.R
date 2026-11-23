# ===================================================================
# 04b_WGCNA_signed_network.R
# Purpose: Signed WGCNA network on all filtered genes (26,172); run on Amarel with 04b_run_WGCNA_amarel.sh
# Input:   Novogene master file
# Output:  results/04_WGCNA/WGCNA3/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/04_WGCNA/WGCNA3", recursive = TRUE, showWarnings = FALSE)

##############################################################
# WGCNA - Weighted Gene Co-expression Network Analysis
# Trans25 - Eastern oyster predator exposure RNA-seq
# GROUP-LEVEL VERSION: traits = TC/TE/UC/UE (each vs. other three)
#THIS ONE HAS CORRECT SIGNED NETWORK TYPE
#THIS ONE USES ALL GENS NO VARIABILITY RESTRICTION
##############################################################

# .libPaths("/scratch/tab307/Rlibs")   # Amarel-specific R library location

library(WGCNA)
library(DESeq2)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

wgcna_dir <- "results/04_WGCNA/WGCNA3"
if (!dir.exists(wgcna_dir)) dir.create(wgcna_dir)

##############################################################
# STEP 1: Load raw counts
##############################################################

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                      header = TRUE, sep = "\t", check.names = FALSE, quote = "")
cat("master dimensions:", dim(master), "\n")

count_cols <- c("TC1_count","TC2_count","TC3_count","TC4_count",
                 "TE1_count","TE2_count","TE3_count","TE4_count",
                 "UC1_count","UC2_count","UC3_count","UC4_count",
                 "UE1_count","UE2_count","UE3_count","UE4_count")

counts <- master[, count_cols]
rownames(counts) <- master$gene_id
colnames(counts) <- gsub("_count", "", colnames(counts))

cat("Counts matrix dimensions:", dim(counts), "\n")

##############################################################
# STEP 2: Sample metadata + DESeq2 object
##############################################################

sample_names <- colnames(counts)
coldata <- data.frame(
  sample       = sample_names,
  longterm     = substr(sample_names, 1, 1),   # "T"/"U" - kept for DESeq2 design
  shortterm    = substr(sample_names, 2, 2),   # "C"/"E" - kept for DESeq2 design
  sample_group = substr(sample_names, 1, 2),   # "TC"/"TE"/"UC"/"UE" - used for WGCNA traits
  row.names    = sample_names
)
coldata$longterm  <- factor(coldata$longterm,  levels = c("U", "T"))
coldata$shortterm <- factor(coldata$shortterm, levels = c("C", "E"))

stopifnot(all(colnames(counts) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ longterm + shortterm)

##############################################################
# STEP 3: Filter low-count genes, variance-stabilize
##############################################################

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
cat("Genes remaining after low-count filter:", nrow(dds), "\n")

vsd <- vst(dds, blind = TRUE)
expr_matrix <- t(assay(vsd))

cat("Expression matrix for WGCNA:", dim(expr_matrix), "\n")

##############################################################
# STEP 4: Outlier check
##############################################################

gsg <- goodSamplesGenes(expr_matrix, verbose = 3)
cat("All genes/samples good:", gsg$allOK, "\n")

if (!gsg$allOK) {
  expr_matrix <- expr_matrix[gsg$goodSamples, gsg$goodGenes]
  cat("Filtered to good samples/genes:", dim(expr_matrix), "\n")
}

sampleTree <- hclust(dist(expr_matrix), method = "average")

png(file.path(wgcna_dir, "sample_clustering3.png"), width = 1200, height = 800, res = 150)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex = 0.9)
dev.off()

##############################################################
# STEP 5: Soft-thresholding power
##############################################################

powers <- c(1:30)
sft <- pickSoftThreshold(expr_matrix, powerVector = powers, networkType = "signed", verbose = 5)

cat("\nSoft-thresholding fit indices by power:\n")

png(file.path(wgcna_dir, "soft_threshold_selection3.png"), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

soft_power <- 11   # confirmed correct for this dataset

##############################################################
# STEP 6: Construct network, detect modules
# (same network as before - only the traits in Step 7 change)
##############################################################

cor <- WGCNA::cor

net <- blockwiseModules(expr_matrix,
                         power = soft_power,
			 networkType = "signed",
                         TOMType = "signed",
                         maxBlockSize = 40000,
                         minModuleSize = 30,
                         reassignThreshold = 0,
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)

cat("\nModules detected:\n")
print(table(net$colors))

moduleColors <- labels2colors(net$colors)

png(file.path(wgcna_dir, "module_dendrogram3.png"), width = 1600, height = 900, res = 150)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                     "Module colors", dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
dev.off()

##############################################################
# STEP 7: Correlate module eigengenes with GROUP traits
# Each column = that group (4 samples) vs. the other three (12
# samples) pooled - e.g. TC = 1 for TC samples, 0 for TE/UC/UE.
# Note: with only 4 samples per group these are noisier/less
# powered than the pooled longterm/shortterm (8 vs 8) traits -
# treat alongside, not instead of, those results.
##############################################################

MEs <- net$MEs
MEs <- orderMEs(MEs)

traits <- data.frame(
  TC = as.numeric(coldata$sample_group == "TC"),
  TE = as.numeric(coldata$sample_group == "TE"),
  UC = as.numeric(coldata$sample_group == "UC"),
  UE = as.numeric(coldata$sample_group == "UE")
)
rownames(traits) <- rownames(coldata)

moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expr_matrix))

cat("\nModule-trait correlations:\n")
print(round(moduleTraitCor, 3))
cat("\nModule-trait p-values:\n")
print(round(moduleTraitPvalue, 4))

##############################################################
# STEP 8: Identify significant modules (excludes grey/unassigned)
##############################################################

sig_threshold_r <- 0.5
sig_threshold_p <- 0.05

sig_modules <- character()
for (trait in colnames(traits)) {
  sig <- rownames(moduleTraitCor)[abs(moduleTraitCor[, trait]) > sig_threshold_r &
                                     moduleTraitPvalue[, trait] < sig_threshold_p]
  sig_modules <- union(sig_modules, sig)
}

sig_modules <- setdiff(sig_modules, "ME0")

cat("\nSignificant modules (|r|>0.5, p<0.05, excluding grey):\n")
print(sig_modules)

##############################################################
# STEP 9: Module-trait heatmap - significant modules only
##############################################################

if (length(sig_modules) == 0) {
  cat("\nNo significant modules - skipping heatmap.\n")
} else {
  plotCor <- moduleTraitCor[sig_modules, , drop = FALSE]
  plotPvalue <- moduleTraitPvalue[sig_modules, , drop = FALSE]

  textMatrix <- paste0(signif(plotCor, 2), "\n(", signif(plotPvalue, 1), ")")
  dim(textMatrix) <- dim(plotCor)

  n_rows <- nrow(plotCor)
  plot_height <- max(600, 150 + n_rows * 45)
  label_cex <- if (n_rows > 25) 0.55 else if (n_rows > 12) 0.7 else 0.9

  png(file.path(wgcna_dir, "module_trait_heatmap3.png"),
      width = 1200, height = plot_height, res = 150)
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
}

##############################################################
# STEP 9b: Module eigengene heatmap across samples
##############################################################

library(pheatmap)

annotation_col <- data.frame(
  sample_group = coldata$sample_group,
  row.names = rownames(coldata)
)

png(file.path(wgcna_dir, "module_eigengene_sample_heatmap3.png"),
    width = 2200, height = 1800, res = 200)
pheatmap(
  t(MEs),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize = 8,
  main = "Module Eigengenes Across Samples",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

##############################################################
# STEP 10: Combined gene list with kME + rank
##############################################################

gene_module_map <- data.frame(
  gene_id = colnames(expr_matrix),
  module_number = net$colors,
  module_color = moduleColors
)

kME_all <- signedKME(expr_matrix, MEs)

combined_list <- list()

for (me_name in sig_modules) {
  mod_num <- as.numeric(gsub("ME", "", me_name))
  mod_color <- unique(gene_module_map$module_color[gene_module_map$module_number == mod_num])

  mod_genes <- gene_module_map %>% filter(module_number == mod_num)
  kME_col <- paste0("kME", mod_num)

  stopifnot(kME_col %in% colnames(kME_all))

  mod_genes$kME <- kME_all[mod_genes$gene_id, kME_col]

  mod_genes <- mod_genes %>%
    arrange(desc(abs(kME))) %>%
    mutate(rank_in_module = row_number())

  combined_list[[me_name]] <- mod_genes
  cat("Module", mod_num, "(", mod_color, "):", nrow(mod_genes), "genes\n")
}

combined_df <- bind_rows(combined_list)

write.csv(combined_df,
          file.path(wgcna_dir, "WGCNA_significant_modules_combined3.csv"),
          row.names = FALSE)

cat("\nSaved combined CSV with", nrow(combined_df), "genes across",
    length(sig_modules), "significant modules.\n")

write.csv(gene_module_map, file.path(wgcna_dir, "WGCNA_all_gene_module_assignments3.csv"), row.names = FALSE)

cat("\nWGCNA analysis complete. All outputs saved.\n")
###################
#SAVE ALL OUTPUT
######################
save(net, MEs, moduleColors, gene_module_map, kME_all, coldata, traits,
     moduleTraitCor, moduleTraitPvalue, expr_matrix,
     file = file.path(wgcna_dir, "WGCNA_results3.RData"))
