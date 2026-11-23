# ===================================================================
# 02a_interaction_term_DESeq2.R
# Purpose: DESeq2 interaction model (~ longterm + shortterm + longterm:shortterm) on Novogene raw counts
# Input:   data/novogene/gene_count_novogene_2025.txt
# Output:  results/02_DEGs/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/02_DEGs", recursive = TRUE, showWarnings = FALSE)

##############################################################
# Interaction term (longterm:shortterm) analysis
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

library(DESeq2)

# ---- 1. SET YOUR WORKING DIRECTORY ----
# (no setwd(): run from the RStudio Project at the repository root)

# ---- 2. LOAD RAW COUNT MATRIX ----
# This is the combined file that includes your 16 samples PLUS unrelated
# A/B/C/D group samples and a bunch of annotation columns (gene_name, gene_chr,
# gene_start, etc). We only want gene_id as rownames and our 16 sample columns.
# Adjust filename/path as needed.
raw <- read.table("data/novogene/gene_count_novogene_2025.txt",
                  header = TRUE, sep = "\t", check.names = FALSE,
                  row.names = "gene_id")

# ---- 2b. KEEP ONLY YOUR 16 SAMPLES (drop A/B/C/D groups + annotation columns) ----
my_samples <- c("TC1","TC2","TC3","TC4",
                "TE1","TE2","TE3","TE4",
                "UC1","UC2","UC3","UC4",
                "UE1","UE2","UE3","UE4")

counts <- raw[, my_samples]

# Quick sanity check - should be 16 columns, all integer-looking counts
stopifnot(ncol(counts) == 16)
str(counts)

# ---- 3. BUILD SAMPLE METADATA AUTOMATICALLY FROM COLUMN NAMES ----
# Parses "TC1" -> longterm = T, shortterm = C, etc, directly from sample names
# so there's no separate file to keep in sync with the counts matrix.
coldata <- data.frame(
  sample    = my_samples,
  longterm  = substr(my_samples, 1, 1),   # "T" or "U"
  shortterm = substr(my_samples, 2, 2),   # "C" or "E"
  row.names = my_samples
)

# Make sure factor levels are set with the CONTROL group as reference
coldata$longterm  <- factor(coldata$longterm,  levels = c("U", "T"))
coldata$shortterm <- factor(coldata$shortterm, levels = c("C", "E"))

print(coldata)

# ---- 4. CHECK SAMPLE ORDER MATCHES ----
# CRITICAL: counts columns and coldata rows must be in identical order
stopifnot(all(colnames(counts) == rownames(coldata)))

# ---- 5. BUILD DESeq DATASET WITH FULL FACTORIAL MODEL ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ longterm + shortterm + longterm:shortterm
)

# ---- 6. PRE-FILTER LOW COUNT GENES (recommended, speeds up + cleans results) ----
#keep <- rowSums(counts(dds) >= 10) >= 4   # at least 10 counts in at least 4 samples
#dds <- dds[keep, ]

# ---- 7. RUN DESeq2 ----
dds <- DESeq(dds)

# ---- 8. CHECK AVAILABLE COEFFICIENT NAMES ----
# This will show you the exact name of the interaction term to pull
resultsNames(dds)
# You're looking for something like: "longtermT.shorttermE"

# ---- 9. EXTRACT THE INTERACTION TERM RESULTS ----
# Replace "longtermT.shorttermE" below with whatever resultsNames(dds) printed
res_interaction <- results(dds, name = "longtermT.shorttermE")

# Order by adjusted p-value
res_interaction <- res_interaction[order(res_interaction$padj), ]

# ---- 10. SUMMARY ----
summary(res_interaction)

# ---- 11. APPLY YOUR SIGNIFICANCE THRESHOLD ----
sig_interaction <- subset(res_interaction, padj < 0.05 & abs(log2FoldChange) > 1.5)
cat("Number of genes with significant interaction effect:", nrow(sig_interaction), "\n")

# ---- 12. SAVE RESULTS ----
write.csv(as.data.frame(res_interaction), "results/02_DEGs/interaction_term_all_results_novogene.csv")
write.csv(as.data.frame(sig_interaction), "results/02_DEGs/interaction_term_significant_DEGs_novogene.csv")

cat("Done. Full results saved to interaction_term_all_results.csv\n")
cat("Significant interaction DEGs saved to interaction_term_significant_DEGs.csv\n")
