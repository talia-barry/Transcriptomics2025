# ===================================================================
# 04d_WGCNA_GO_enrichment.R
# Purpose: GO enrichment (clusterProfiler) of significant WGCNA modules
# Input:   results/04_WGCNA/WGCNA3/
# Output:  results/04_WGCNA/WGCNA3/Go_enrichment_signed/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/04_WGCNA/WGCNA3/Go_enrichment_signed", recursive = TRUE, showWarnings = FALSE)

##############################################################
# GO enrichment for WGCNA significant modules for SIGNED network
# Crassostrea virginica - Trans25
##############################################################

# BiocManager::install(c("AnnotationHub", "clusterProfiler"))

library(AnnotationHub)
library(clusterProfiler)
library(dplyr)

##############################################################
# STEP 1: Find and load the C. virginica OrgDb from AnnotationHub
##############################################################

ah <- AnnotationHub()
hits <- query(ah, c("Crassostrea virginica", "OrgDb"))
hits   # inspect - note the AH#### ID for the OrgDb entry

orgdb <- ah[["AH117535"]]   # <-- use with the actual ID from `hits` above

keytypes(orgdb)   # check what ID types this OrgDb supports (GID, ENTREZID, SYMBOL, etc.)

##############################################################
# STEP 2: Check whether your gene IDs actually match this OrgDb
# before running full enrichment - if overlap is near 0%, this
# OrgDb won't work for your annotation and you'll need a custom
# GO mapping (e.g. from eggNOG-mapper/InterProScan) with topGO instead.
##############################################################

all_genes <- read.csv("results/04_WGCNA/WGCNA3/WGCNA_all_gene_module_assignments3.csv", colClasses = c(gene_id = "character"))
my_ids <- unique(all_genes$gene_id)

orgdb_ids <- keys(orgdb, keytype = "GID")   # try "ENTREZID" instead if GID errors

overlap_pct <- 100 * length(intersect(my_ids, orgdb_ids)) / length(my_ids)
cat("Percent of your gene IDs found in this OrgDb:", round(overlap_pct, 1), "%\n")

# STOP HERE and check overlap_pct before proceeding.
# If it's high (>70-80%), continue below.

##############################################################
# STEP 3: Run GO enrichment per significant module
##############################################################

combined <- read.csv("results/04_WGCNA/WGCNA3/WGCNA_significant_modules_combined3.csv", colClasses = c(gene_id = "character"))
modules <- sort(unique(combined$module_number))

background_genes <- my_ids   # universe = ALL genes tested in WGCNA, not just significant ones

go_dir <- "results/04_WGCNA/WGCNA3/Go_enrichment_signed"
if (!dir.exists(go_dir)) dir.create(go_dir)

go_results <- list()

for (m in modules) {
  mod_genes <- combined$gene_id[combined$module_number == m]
  mod_color <- unique(combined$module_color[combined$module_number == m])
  
  cat("\nRunning GO enrichment for module", m, "(", mod_color, ") -", length(mod_genes), "genes\n")
  
  ego <- tryCatch({
    enrichGO(gene = mod_genes,
             universe = background_genes,
             OrgDb = orgdb,
             keyType = "GID",       # match whatever worked in Step 2
             ont = "BP",             # Biological Process - change to "MF"/"CC"/"ALL" if wanted
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) {
    cat("  No enrichment result / error for module", m, ":", conditionMessage(e), "\n")
    NULL
  })
  
  go_results[[paste0("module_", m)]] <- ego
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    out_df <- as.data.frame(ego)
    write.csv(out_df, file.path(go_dir, paste0("GO_module_", m, "_", mod_color, ".csv")), row.names = FALSE)
    cat("  Saved", nrow(out_df), "significant GO terms\n")
  } else {
    cat("  No significant GO terms for this module\n")
  }
}

##############################################################
# STEP 4: Combined summary across all modules (for a quick scan)
##############################################################

summary_list <- lapply(names(go_results), function(nm) {
  ego <- go_results[[nm]]
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  df <- as.data.frame(ego)
  df$module <- nm
  df
})

summary_df <- bind_rows(summary_list)
write.csv(summary_df, file.path(go_dir, "GO_enrichment_all_modules_summary_signed.csv"), row.names = FALSE)

cat("\nGO enrichment complete. Results in", go_dir, "\n")






