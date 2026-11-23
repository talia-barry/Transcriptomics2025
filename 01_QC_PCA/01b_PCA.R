# ===================================================================
# 01b_PCA.R
# Purpose: PCA of samples (VST counts) from the Novogene master file
# Input:   data/novogene/MasterFile_novogene_DEGs_2025.xls
# Output:  results/01_QC_PCA/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/01_QC_PCA", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PCA of DEGs from Novogene Masterfile
# Crassostrea virginica - Predator Exposure RNA-seq (Chapter 4)
# ============================================================
# Run inside your RStudio Project (no setwd needed)

library(DESeq2)
library(dplyr)      # loaded after DESeq2 so dplyr::select() takes priority
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# 1. Load the masterfile
# ------------------------------------------------------------
masterfile <- read.delim(
  "data/novogene/MasterFile_novogene_DEGs_2025.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE, 
  quote = "")

dim(masterfile)

# ------------------------------------------------------------
# 2. Flag DEGs using only the TC/TE/UC/UE-relevant comparisons
# ------------------------------------------------------------
comparisons <- c("TvsU", "EvsN", "TEvsUC", "TEvsUE",
                 "TCvsUC", "UEvsUC", "TEvsTC")

padj_cols   <- paste0(comparisons, "_padj")
log2fc_cols <- paste0(comparisons, "_log2FoldChange")

stopifnot(all(padj_cols %in% colnames(masterfile)))
stopifnot(all(log2fc_cols %in% colnames(masterfile)))

deg_flags <- sapply(seq_along(comparisons), function(i) {
  padj   <- masterfile[[padj_cols[i]]]
  log2fc <- masterfile[[log2fc_cols[i]]]
  (!is.na(padj) & padj < 0.05) & (!is.na(log2fc) & abs(log2fc) > 1.5)
})

masterfile$is_DEG <- apply(deg_flags, 1, any)

deg_table <- masterfile %>% dplyr::filter(is_DEG)
cat("Number of DEGs (union across TC/TE/UC/UE comparisons):",
    nrow(deg_table), "\n")

# ------------------------------------------------------------
# 3. Build raw count matrix (TC/TE/UC/UE samples only, DEGs only)
# ------------------------------------------------------------
count_cols <- grep("^(TC|TE|UC|UE)[0-9]_count$", colnames(deg_table), value = TRUE)

count_mat <- deg_table %>%
  dplyr::select(all_of(count_cols)) %>%
  as.matrix()

rownames(count_mat) <- deg_table$gene_id
mode(count_mat) <- "integer"

# ------------------------------------------------------------
# 4. Build sample metadata
# ------------------------------------------------------------
sample_names <- sub("_count$", "", count_cols)   # e.g. "TC1"

col_data <- data.frame(
  row.names      = sample_names,
  group          = substr(sample_names, 1, nchar(sample_names) - 1),  # TC, TE, UC, UE
  exposure_6wk   = ifelse(grepl("^T", sample_names), "Trained", "Untrained"),
  reexposure_24h = ifelse(grepl("E[0-9]$", sample_names), "Exposed", "Control")
)

colnames(count_mat) <- sample_names
count_mat <- count_mat[, rownames(col_data)]  # ensure order matches

# ------------------------------------------------------------
# 5. Variance-stabilizing transform + PCA
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = col_data,
  design    = ~ group
)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

pca <- prcomp(t(assay(vsd)), scale. = TRUE, center = TRUE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- as.data.frame(pca$x) %>%
  mutate(sample = rownames(.)) %>%
  dplyr::left_join(
    col_data %>% mutate(sample = rownames(col_data)),
    by = "sample"
  )

# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.68, color = NA) +
  geom_point(aes(shape = exposure_6wk), size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  labs(
    title = "PCA of DEGs — C. virginica predator exposure",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    color = "Group", fill = "Group", shape = "6-wk Exposure"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(p)
ggsave("results/01_QC_PCA/PCA_DEGs_oldGenome.png", p, width = 7, height = 5.5, dpi = 300)
