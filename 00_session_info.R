# ===================================================================
# 00_session_info.R
# Purpose: Record the R version and package versions used in this
#          analysis, for reproducibility. Run once from the repository
#          root (open Transcriptomics2025.Rproj); it writes session_info.txt.
# Note: 04b (WGCNA network) was run on the Rutgers Amarel cluster inside
#       the R 4.5.1 container listed in 04_WGCNA/04b_run_WGCNA_amarel.sh;
#       its WGCNA/DESeq2 versions are recorded in the 04b job log.
# ===================================================================

pkgs <- c(
  # Bioconductor
  "DESeq2", "WGCNA", "fgsea", "rrvgo", "clusterProfiler",
  "AnnotationHub", "AnnotationDbi", "GO.db", "org.Hs.eg.db",
  # CRAN
  "dplyr", "tidyr", "stringr", "ggplot2", "ggpubr", "ggrepel",
  "ggVennDiagram", "patchwork", "pheatmap", "treemap"
)

# Load each package quietly so it appears in sessionInfo()
missing <- character(0)
for (p in pkgs) {
  ok <- suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE))
  if (ok) suppressPackageStartupMessages(library(p, character.only = TRUE)) else missing <- c(missing, p)
}
if (length(missing) > 0) cat("Not installed (skipped):", paste(missing, collapse = ", "), "\n")

# Simple version table: one line per package
versions <- data.frame(
  package = pkgs,
  version = vapply(pkgs, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_,
    character(1)),
  row.names = NULL
)

sink("session_info.txt")
cat("Session info for Transcriptomics2025 analyses\n")
cat("Recorded:", format(Sys.time(), "%Y-%m-%d"), "\n")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  cat("RStudio:", as.character(rstudioapi::versionInfo()$version), "\n")
cat("\n## Package versions\n")
print(versions, row.names = FALSE)
cat("\n## Full sessionInfo()\n")
print(sessionInfo())
sink()

print(versions, row.names = FALSE)
cat("\nSaved session_info.txt\n")
