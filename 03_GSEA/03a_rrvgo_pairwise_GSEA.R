# ===================================================================
# 03a_rrvgo_pairwise_GSEA.R
# Purpose: Clusters significant GO terms (GSEA, FDR < 0.25) with rrvgo for the pairwise comparisons; makes Table 2.2
# Input:   data/novogene/GSEA/<comparison>/gsea_report_*.tsv (GSEA reports from Novogene)
# Output:  results/03_GSEA/rrvgo_outputs/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/03_GSEA/rrvgo_outputs", recursive = TRUE, showWarnings = FALSE)

##############################################################
# GO Term Clustering (rrvgo) - all GSEA comparisons, both directions
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

# ---- INSTALL (run once if not already installed) ----
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("rrvgo", "GOSemSim", "org.Hs.eg.db"))
# install.packages(c("dplyr", "stringr", "treemap"))

library(rrvgo)
library(dplyr)
library(stringr)

##############################################################
# STEP 1: Settings
##############################################################

FDR_CUTOFF <- 0.25   # matches your GSEA significance convention
SIM_THRESHOLD <- 0.4 # matches your validated REVIGO clustering threshold

base_dir   <- "data/novogene/GSEA"
output_dir <- "results/03_GSEA/rrvgo_outputs"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Each comparison and its two phenotype labels (matches your file naming:
# gsea_report_for_{comparison}_{phenotype}.tsv)
comparisons <- list(
  "UEvsUC" = c("UE", "UC"),
  "TCvsUC" = c("TC", "UC"),
  "TEvsTC" = c("TE", "TC"),
  "TEvsUC" = c("TE", "UC"),
  "TEvsUE" = c("TE", "UE")
)

ontologies <- c("BP", "MF", "CC")

##############################################################
# STEP 2: Helper function - extract GO IDs + FDR from one GSEA report
##############################################################

load_gsea_report <- function(filepath) {
  
  df <- read.delim(filepath, header = TRUE, check.names = FALSE)
  
  # Extract GO ID from the NAME column, e.g. "RIBOSOME_BIOGENESIS(GO:0042254)"
  df$GOID <- str_extract(df$NAME, "GO:[0-9]+")
  
  # Drop rows where no GO ID could be parsed (shouldn't normally happen)
  df <- df[!is.na(df$GOID), ]
  
  # Keep just what we need
  df <- df %>% select(GOID, `FDR q-val`, NES) %>%
    rename(fdr = `FDR q-val`)
  
  df
}

##############################################################
# STEP 3: Helper function - run rrvgo for one GO set + one ontology
##############################################################
# UPDATED: now accepts nes_values so NES is available inside the function
# (previously NES was only merged in later, in Step 5).

run_rrvgo <- function(go_ids, fdr_values, nes_values, ont, label) {
  
  if (length(go_ids) < 2) {
    cat("  Skipping", ont, "- fewer than 2 terms after filtering\n")
    return(NULL)
  }
  
  simMatrix <- tryCatch({
    calculateSimMatrix(go_ids,
                       orgdb = "org.Hs.eg.db",
                       ont = ont,
                       method = "Wang")   # graph-topology based - no species
    # annotation frequency needed, so
    # choice of orgdb here barely matters
  }, error = function(e) {
    cat("  No terms found for ontology", ont, "- skipping\n")
    return(NULL)
  })
  
  # Safely determine row count, regardless of what calculateSimMatrix returned
  n_rows <- tryCatch({
    r <- nrow(simMatrix)
    if (is.null(r) || is.na(r)) 0 else r
  }, error = function(e) 0)
  
  if (is.null(simMatrix) || n_rows < 2) {
    cat("  Skipping", ont, "- not enough valid terms in this ontology\n")
    return(NULL)
  }
  
  # UPDATED: Score boxes by strength of enrichment (magnitude of NES), not
  # FDR. NES can be negative (enriched toward the opposite phenotype), so we
  # use abs() for box sizing - treemapPlot()/reduceSimMatrix() need positive
  # values for area and representative-term selection. The real, signed NES
  # is still shown in the label and saved in the CSV below.
  scores_named <- setNames(abs(nes_values), go_ids)
  scores_named <- scores_named[rownames(simMatrix)]
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores_named,
                                  threshold = SIM_THRESHOLD,
                                  orgdb = "org.Hs.eg.db")
  
  # Save reduced term table - keep the real FDR and NES alongside the
  # abs(NES) score used for box sizing
  reducedTerms$fdr <- fdr_values[match(reducedTerms$go, go_ids)]
  reducedTerms$NES <- nes_values[match(reducedTerms$go, go_ids)]
  
  # Build a plot-only copy with the real (signed) NES appended to the label
  # (keeps the CSV's "term" column clean for write.csv below)
  treemapTerms <- reducedTerms
  treemapTerms$term <- paste0(treemapTerms$term,
                              "\n(NES=", round(treemapTerms$NES, 2), ")")
  
  write.csv(reducedTerms,
            file.path(output_dir, paste0(label, "_", ont, "_reduced_terms.csv")),
            row.names = FALSE)
  
  # Save treemap - plots treemapTerms (with NES in label) instead of reducedTerms
  png(file.path(output_dir, paste0(label, "_", ont, "_treemap.png")),
      width = 1600, height = 1200, res = 150)
  treemapPlot(treemapTerms, title = paste0(label, " - ", ont))
  dev.off()
  
  cat("  Saved", ont, "treemap + table for", label, "\n")
  
  reducedTerms
}

##############################################################
# STEP 4: Loop through every comparison x phenotype x ontology
##############################################################

all_reduced <- list()

for (comp_name in names(comparisons)) {
  
  phenotypes <- comparisons[[comp_name]]
  
  for (pheno in phenotypes) {
    
    filepath <- file.path(base_dir, comp_name,
                          paste0("gsea_report_for_", comp_name, "_", pheno, ".tsv"))
    
    if (!file.exists(filepath)) {
      cat("File not found, skipping:", filepath, "\n")
      next
    }
    
    cat("\n=== Processing", comp_name, "-", pheno, "===\n")
    
    report <- load_gsea_report(filepath)
    
    # Filter to your significance threshold
    sig <- report[!is.na(report$fdr) & report$fdr < FDR_CUTOFF, ]
    
    cat("  ", nrow(sig), "significant terms at FDR <", FDR_CUTOFF, "\n")
    
    if (nrow(sig) < 2) {
      cat("  Not enough significant terms to cluster - skipping this phenotype\n")
      next
    }
    
    label <- paste0(comp_name, "_", pheno)
    
    for (ont in ontologies) {
      cat(" Ontology:", ont, "\n")
      # UPDATED: now passes sig$NES so run_rrvgo can size boxes by abs(NES)
      # and show real NES on labels
      result <- run_rrvgo(sig$GOID, sig$fdr, sig$NES, ont, label)
      if (!is.null(result)) {
        # Tag with metadata so the combined table (Step 5) can be sorted/filtered
        result$comparison <- comp_name
        result$phenotype  <- pheno
        result$ontology   <- ont
        all_reduced[[paste0(label, "_", ont)]] <- result
      }
    }
  }
}

cat("\n\nDone. All treemaps and reduced term tables saved to:\n", output_dir, "\n")

##############################################################
# STEP 5: Combine every comparison/phenotype/ontology's reduced terms
# into one master table, with NES merged back in from the original
# GSEA reports, sorted by significance (lowest FDR first).
#
# NOTE: This step is self-contained - it re-reads the individual
# "*_reduced_terms.csv" files already saved to output_dir, plus the
# original GSEA report TSVs, so it can be re-run on its own without
# repeating Steps 1-4, AS LONG AS base_dir, output_dir, comparisons,
# and load_gsea_report() are still defined in your R session.
##############################################################

reduced_files <- list.files(output_dir, pattern = "_reduced_terms\\.csv$", full.names = TRUE)

if (length(reduced_files) == 0) {
  cat("\nNo reduced term files found in", output_dir, "- nothing to combine.\n")
} else {
  
  combined_list <- list()
  
  for (f in reduced_files) {
    
    fname <- basename(f)
    
    # Filename format: {comparison}_{phenotype}_{ontology}_reduced_terms.csv
    # e.g. "UEvsUC_UE_BP_reduced_terms.csv"
    parts <- str_match(fname, "^(.+)_([A-Z]{2})_reduced_terms\\.csv$")
    label_part <- parts[2]   # e.g. "UEvsUC_UE"
    ont         <- parts[3]  # e.g. "BP"
    
    label_parts <- str_match(label_part, "^(.+)_([A-Za-z]+)$")
    comp_name <- label_parts[2]  # e.g. "UEvsUC"
    pheno     <- label_parts[3]  # e.g. "UE"
    
    reduced_df <- read.csv(f)
    
    # Load the matching original GSEA report to pull NES
    # (reduced_df already has NES from Step 3/4, but this re-pull from the
    # source report is left as-is and simply overwrites with the same values)
    report_path <- file.path(base_dir, comp_name,
                             paste0("gsea_report_for_", comp_name, "_", pheno, ".tsv"))
    
    if (!file.exists(report_path)) {
      cat("  Could not find original report for", fname, "- NES will be NA\n")
      reduced_df$NES <- NA
    } else {
      report <- load_gsea_report(report_path)
      nes_lookup <- setNames(report$NES, report$GOID)
      reduced_df$NES <- nes_lookup[reduced_df$go]
    }
    
    reduced_df$comparison <- comp_name
    reduced_df$phenotype  <- pheno
    reduced_df$ontology   <- ont
    
    combined_list[[fname]] <- reduced_df
  }
  
  combined <- bind_rows(combined_list)
  
  # Reorder columns so the useful ones are up front
  combined <- combined %>%
    select(comparison, phenotype, ontology, go, term, parentTerm,
           fdr, NES, size, cluster, everything()) %>%
    arrange(fdr)
  
  write.csv(combined,
            file.path(output_dir, "ALL_reduced_terms_combined.csv"),
            row.names = FALSE)
  
  cat("\nCombined table saved to:",
      file.path(output_dir, "ALL_reduced_terms_combined.csv"), "\n")
  cat("Sorted by FDR (most significant first) -", nrow(combined), "total rows.\n")
  cat("Open this file and sort/filter by 'fdr' or 'NES' to identify your top",
      "terms across every comparison at once.\n")
  
  ##############################################################
  # STEP 6: Supplementary table - full clean list, NA-term rows dropped
  ##############################################################
  # Some GO IDs don't resolve to a readable term name (often obsolete/
  # updated IDs) - these are uninformative on their own, so they're
  # excluded here. Check them manually via QuickGO/AmiGO if you want
  # them back in.
  
  supp_table <- combined %>%
    filter(!is.na(term) & term != "") %>%
    select(comparison, phenotype, ontology, go, term, parentTerm,
           NES, fdr, size) %>%
    arrange(fdr)
  
  write.csv(supp_table,
            file.path(output_dir, "Supplementary_Table_GSEA_GO_terms.csv"),
            row.names = FALSE)
  
  cat("\nSupplementary table saved -", nrow(supp_table), "rows (NA-term rows excluded)\n")
  
  ##############################################################
  # STEP 7: Curated table for main text - specific GO IDs grouped by theme
  ##############################################################
  # EDIT THIS LIST as your themes/results change - add or remove GO IDs
  # and reassign theme labels to match whatever your results text says.
  
  curated_ids_themes <- c(
    "GO:0042254" = "Ribosome biogenesis / RNA processing",
    "GO:0022613" = "Ribosome biogenesis / RNA processing",
    "GO:0006364" = "Ribosome biogenesis / RNA processing",
    "GO:0019725" = "Redox homeostasis / oxidoreductase activity",
    "GO:0045454" = "Redox homeostasis / oxidoreductase activity",
    "GO:0016614" = "Redox homeostasis / oxidoreductase activity",
    "GO:0016616" = "Redox homeostasis / oxidoreductase activity",
    "GO:0016853" = "Redox homeostasis / oxidoreductase activity",
    "GO:0016866" = "Redox homeostasis / oxidoreductase activity",
    "GO:0005839" = "Secondary: protein turnover / mitochondrial",
    "GO:0004298" = "Secondary: protein turnover / mitochondrial",
    "GO:0070003" = "Secondary: protein turnover / mitochondrial",
    "GO:0005739" = "Secondary: protein turnover / mitochondrial",
    "GO:0098798" = "Secondary: protein turnover / mitochondrial"
  )
  
  theme_order <- c("Ribosome biogenesis / RNA processing",
                   "Redox homeostasis / oxidoreductase activity",
                   "Secondary: protein turnover / mitochondrial")
  
  curated_table <- supp_table %>%
    filter(go %in% names(curated_ids_themes)) %>%
    mutate(theme = curated_ids_themes[go],
           theme = factor(theme, levels = theme_order)) %>%
    select(theme, go, term, ontology, NES, fdr, size) %>%
    arrange(theme, fdr)
  
  write.csv(curated_table,
            file.path(output_dir, "Table_2_2_GSEA_GO_terms_curated.csv"),
            row.names = FALSE)
  
  cat("\nCurated table saved -", nrow(curated_table), "rows across",
      length(unique(curated_table$theme)), "themes\n")
  print(curated_table)

}  # closes the if/else from Step 5 (was missing; added for GitHub version)

  ##############################################################
  # STEP 8: Curated table for TE (TEvsUC comparison) only
  ##############################################################
  # GO ID curation like Table 2.2, but restricted to the
  # TEvsUC comparison, TE phenotype side only.

te_table <- supp_table %>%
  filter(comparison == "TEvsUC", phenotype == "TE") %>%
  select(go, term, ontology, NES, fdr, size) %>%
  arrange(fdr)

te_table

write.csv(te_table,
          file.path(output_dir, "Table_GSEA_GO_terms_TEvsUC_TE.csv"),
          row.names = FALSE)
