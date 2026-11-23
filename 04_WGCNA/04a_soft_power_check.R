# ===================================================================
# 04a_soft_power_check.R
# Purpose: Soft-threshold power selection for the signed WGCNA network (run locally)
# Input:   Novogene master file
# Output:  results/04_WGCNA/WGCNA_signed_check/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/04_WGCNA/WGCNA_signed_check", recursive = TRUE, showWarnings = FALSE)

###############################################################
# WGCNA - Soft threshold check ONLY (SIGNED, full gene set)
# Run LOCALLY in RStudio - displays plot on screen.
#   *** networkType = "signed" added to pickSoftThreshold() -
#      matches the signed network you'll build with
#      blockwiseModules(). Without this, pickSoftThreshold()
#      fits the curve as if unsigned by default, and the power
#      it recommends won't be correct for a signed network.
##############################################################

library(WGCNA)
library(DESeq2)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

wgcna_dir <- "results/04_WGCNA/WGCNA_signed_check"
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
  sample_group = substr(sample_names, 1, 2),   # "TC"/"TE"/"UC"/"UE" - each is its own group
  row.names    = sample_names
)
coldata$sample_group <- factor(coldata$sample_group, levels = c("TC", "TE", "UC", "UE"))

stopifnot(all(colnames(counts) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ sample_group)

##############################################################
# STEP 3: Filter low-count genes, variance-stabilize
# (no variance filter here - all genes that pass this
# low-count filter go into the soft power check)
##############################################################

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
cat("Genes remaining after low-count filter:", nrow(dds), "\n")

vsd <- vst(dds, blind = TRUE)
expr_matrix <- t(assay(vsd))

cat("Expression matrix for soft power check (all genes):", dim(expr_matrix), "\n")

##############################################################
# STEP 4: Outlier check
##############################################################

gsg <- goodSamplesGenes(expr_matrix, verbose = 3)
cat("All genes/samples good:", gsg$allOK, "\n")

if (!gsg$allOK) {
  expr_matrix <- expr_matrix[gsg$goodSamples, gsg$goodGenes]
  cat("Filtered to good samples/genes:", dim(expr_matrix), "\n")
}

##############################################################
# STEP 5: Soft-thresholding power - SIGNED network
# networkType = "signed" matches the signed blockwiseModules()
# network you'll build next - signed networks typically need a
# HIGHER power than unsigned to reach R^2 ~ 0.85 (often in the
# 12-20 range rather than 6-12), so don't assume 6 still applies.
##############################################################

powers <- c(1:30)   # widened range since signed networks often need higher power
sft <- pickSoftThreshold(expr_matrix, powerVector = powers,
                         networkType = "signed", ##ADDED HEREEEEE
                         verbose = 5)

cat("\nSoft-thresholding fit indices by power (signed network):\n")
print(sft$fitIndices[, c("Power", "SFT.R.sq", "slope", "mean.k.", "median.k.", "max.k.")])

dev.new(width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence (signed, all genes)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity (signed, all genes)", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

cat("\nInspect the plot window and the table above to pick soft_power.\n")
cat("Remember: this power is for networkType = \"signed\" - use the same\n")
cat("networkType = \"signed\" in your blockwiseModules() call along with it.\n")

# Save expr_matrix + sft so you don't have to redo Steps 1-4 later
save(expr_matrix, sft, coldata,
     file = file.path(wgcna_dir, "WGCNA_signed_pre_network.RData"))
cat("Saved expr_matrix + sft to WGCNA_signed_pre_network.RData\n")