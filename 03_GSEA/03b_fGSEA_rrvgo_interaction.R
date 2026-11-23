# ===================================================================
# 03b_fGSEA_rrvgo_interaction.R
# Purpose: fgsea on the interaction term (GO gene sets built from AnnotationHub C. virginica OrgDb) + rrvgo clustering; makes Table 3.2
# Input:   Novogene master file + 02a output
# Output:  results/03_GSEA/
# Genome: NCBI C. virginica 2017 assembly (C_virginica-3.0); read mapping + counts by Novogene
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/03_GSEA", recursive = TRUE, showWarnings = FALSE)

##############################################################
# Build a custom GO .gmt file using AnnotationHub's pre-built
# C. virginica OrgDb (no manual NCBI download needed)
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

# install.packages("BiocManager")
# BiocManager::install(c("AnnotationHub", "AnnotationDbi"))

library(AnnotationHub)
library(AnnotationDbi)
library(dplyr)

##############################################################
# STEP 1: Connect to AnnotationHub and pull the C. virginica OrgDb
##############################################################
# This queries Bioconductor's AnnotationHub for the pre-built
# organism database - first run will download and cache it locally
# (~few hundred MB, one-time download to ~/.AnnotationHub)

ah <- AnnotationHub()

query_result <- query(ah, c("OrgDb", "Crassostrea virginica"))
print(query_result)

# Should show one record - typically "AH96120" as of this writing,
# but confirm from the printed list above in case IDs have shifted
orgdb <- ah[["AH117535"]]

print(orgdb)

##############################################################
# STEP 2: Load your gene IDs (to restrict the gmt to genes in your data)
##############################################################

master <- read.table("data/novogene/MasterFile_novogene_DEGs_2025.xls",
                     header = TRUE, sep = "\t", check.names = FALSE,
                     quote = "", comment.char = "")   # see 02a: apostrophes otherwise drop genes

your_gene_ids <- as.character(master$gene_id)

cat(length(your_gene_ids), "genes in your dataset\n")

##############################################################
# STEP 3: Pull GO annotations for your genes from the OrgDb
##############################################################
# keytype "GID" is standard for AnnotationHub-built NCBI OrgDbs -
# if this errors, run keytypes(orgdb) to see the correct key name

go_annot <- AnnotationDbi::select(orgdb,
                                  keys = your_gene_ids,
                                  keytype = "GID",
                                  columns = c("GO", "ONTOLOGY"))

cat("GO annotation rows retrieved:", nrow(go_annot), "\n")

# Drop rows with no GO annotation
go_annot <- go_annot[!is.na(go_annot$GO), ]

##############################################################
# STEP 4: Get readable GO term names (GO IDs alone aren't descriptive)
##############################################################
# install.packages("BiocManager"); BiocManager::install("GO.db")

library(GO.db)

go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(go_annot$GO),
                                  keytype = "GOID",
                                  columns = "TERM")

go_annot <- go_annot %>%
  left_join(go_terms, by = c("GO" = "GOID"))

##############################################################
# STEP 5: Build the .gmt structure
##############################################################

gmt_data <- go_annot %>%
  distinct(GID, GO, TERM) %>%
  group_by(GO, TERM) %>%
  summarise(genes = list(unique(GID)), .groups = "drop") %>%
  filter(lengths(genes) >= 15, lengths(genes) <= 5000)   # matches your set size convention

cat("GO terms after size filtering (min=15, max=5000):", nrow(gmt_data), "\n")

##############################################################
# STEP 6: Write out as a .gmt file
##############################################################

gmt_lines <- apply(gmt_data, 1, function(row) {
  paste(c(row$GO, row$TERM, unlist(row$genes)), collapse = "\t")
})

writeLines(gmt_lines, "results/03_GSEA/Cvirginica_GO.gmt")

cat("\nSaved Cvirginica_GO.gmt with", length(gmt_lines), "gene sets\n")
cat("Ready to use with fgsea's gmtPathways() function.\n")


##############################################################
# STEP 7: Run pre-ranked GSEA (fgsea) on the interaction term
# using the .gmt file just built above
##############################################################

# BiocManager::install("fgsea")
library(fgsea)

FDR_CUTOFF <- 0.25   # matches your GSEA significance convention elsewhere
N_PERM     <- 1000   # matches Novogene's GSEA settings (nperm: 1000)
SEED       <- 42     # for reproducibility of permutation-based p-values

set.seed(SEED)

# ---- Load interaction term results (full list, not just significant) ----
interaction_res <- read.csv("results/02_DEGs/interaction_term_all_results_novogene.csv")

# ---- Load the GO gene sets just built ----
pathways <- gmtPathways("results/03_GSEA/Cvirginica_GO.gmt")

cat("\nLoaded", length(pathways), "GO gene sets for fgsea\n")
cat(nrow(interaction_res), "genes in interaction results\n")

# ---- Build the ranked gene list ----
# DESeq2's Wald statistic ("stat" column) is the standard ranking metric
# for preranked GSEA on DESeq2 output. Adjust the column name below if
# your file's header differs.

stopifnot("stat" %in% colnames(interaction_res))

if ("gene_id" %in% colnames(interaction_res)) {
  gene_ids <- as.character(interaction_res$gene_id)
} else {
  gene_ids <- rownames(interaction_res)
}

ranks <- interaction_res$stat
names(ranks) <- gene_ids
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

cat(length(ranks), "genes with a valid test statistic, ready for GSEA\n")

# ---- Run fgsea ----
fgsea_res <- fgsea(pathways = pathways,
                   stats = ranks,
                   minSize = 15,
                   maxSize = 5000,
                   nPermSimple = N_PERM)

fgsea_res <- fgsea_res %>% arrange(padj)

cat("\nGSEA complete.\n")
cat(sum(fgsea_res$padj < FDR_CUTOFF), "gene sets significant at FDR <", FDR_CUTOFF, "\n")

##check the results
# How many gene sets and ranked genes do we actually have?
length(pathways)
length(ranks)

# Peek at the ID formats on both sides
head(names(ranks))
head(pathways[[1]])

# Check actual overlap
all_pathway_genes <- unique(unlist(pathways))
length(all_pathway_genes)
sum(names(ranks) %in% all_pathway_genes)

# ---- Save results ----
fgsea_out <- fgsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ",")))

write.csv(fgsea_out, "results/03_GSEA/interaction_term_fgsea_results.csv", row.names = FALSE)
cat("Full results saved to interaction_term_fgsea_results.csv\n")

fgsea_sig <- fgsea_out %>% filter(padj < FDR_CUTOFF)
write.csv(fgsea_sig, "results/03_GSEA/interaction_term_fgsea_significant.csv", row.names = FALSE)
cat(nrow(fgsea_sig), "significant gene sets saved to interaction_term_fgsea_significant.csv\n")

# ---- Quick look at top results ----
sum(fgsea_res$padj < FDR_CUTOFF)
fgsea_res %>% dplyr::select(pathway, NES, pval, padj, size) %>% arrange(pval) %>% head(15)

##############################################################
# GO Term Clustering (rrvgo) - Interaction Term fgsea results
# Trans25 - Eastern oyster predator exposure RNA-seq
##############################################################

##############################################################
# RELOAD BLOCK - run this if starting a fresh R session/environment
# and you only want to run Step 3 onward below. Re-loads everything
# Step 3+ depends on: libraries, output_dir, and the fgsea results
# (re-read from the CSV saved earlier, not recomputed).
##############################################################

library(rrvgo)
library(dplyr)

SIM_THRESHOLD <- 0.4   # matches your validated clustering threshold
output_dir <- "results/03_GSEA/rrvgo_interaction_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

fgsea_sig <- read.csv("results/03_GSEA/interaction_term_fgsea_significant.csv")

cat(nrow(fgsea_sig), "significant gene sets loaded\n")

go_ids <- fgsea_sig$pathway
padj_values <- fgsea_sig$padj
nes_values <- setNames(fgsea_sig$NES, fgsea_sig$pathway)

ontologies <- c("BP", "MF", "CC")

##############################################################
# STEP 3: Run rrvgo separately per ontology (BP/MF/CC)
##############################################################
# UPDATED: box size now driven by abs(NES) instead of -log10(padj),
# and treemap labels now show the real (signed) NES value.

all_reduced <- list()

for (ont in ontologies) {
  
  cat("\n=== Ontology:", ont, "===\n")
  
  simMatrix <- tryCatch({
    calculateSimMatrix(go_ids,
                       orgdb = "org.Hs.eg.db",
                       ont = ont,
                       method = "Wang")
  }, error = function(e) {
    cat("  No terms found for this ontology - skipping\n")
    return(NULL)
  })
  
  n_rows <- tryCatch({
    r <- nrow(simMatrix)
    if (is.null(r) || is.na(r)) 0 else r
  }, error = function(e) 0)
  
  if (is.null(simMatrix) || n_rows < 2) {
    cat("  Not enough valid terms in this ontology to cluster\n")
    next
  }
  
  # UPDATED: score boxes by strength of enrichment (magnitude of NES),
  # not padj. NES can be negative, so abs() is needed for box area and
  # for reduceSimMatrix()'s representative-term selection. The real,
  # signed NES is kept separately below for the CSV and the label.
  scores_named <- setNames(abs(nes_values[go_ids]), go_ids)
  scores_named <- scores_named[rownames(simMatrix)]
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores_named,
                                  threshold = SIM_THRESHOLD,
                                  orgdb = "org.Hs.eg.db")
  
  # Attach real padj and NES values back in (not the abs(NES) score used for sizing)
  reducedTerms$padj <- padj_values[match(reducedTerms$go, go_ids)]
  reducedTerms$NES  <- nes_values[reducedTerms$go]
  reducedTerms$ontology <- ont
  
  # Build a plot-only copy with the real (signed) NES appended to the label
  # (keeps the CSV's "term" column clean for write.csv below)
  treemapTerms <- reducedTerms
  treemapTerms$term <- paste0(treemapTerms$term,
                              "\n(NES=", round(treemapTerms$NES, 2), ")")
  
  write.csv(reducedTerms,
            file.path(output_dir, paste0("interaction_", ont, "_reduced_terms.csv")),
            row.names = FALSE)
  
  png(file.path(output_dir, paste0("interaction_", ont, "_treemap.png")),
      width = 1600, height = 1200, res = 150)
  treemapPlot(treemapTerms, title = paste0("Interaction Term - ", ont))
  dev.off()
  
  cat("  Saved", ont, "treemap + table -", nrow(reducedTerms), "terms\n")
  
  all_reduced[[ont]] <- reducedTerms
}

##############################################################
# STEP 4: Combine into one table, sorted by significance
##############################################################

if (length(all_reduced) > 0) {
  
  combined <- bind_rows(all_reduced) %>%
    dplyr::select(ontology, go, term, parentTerm, NES, padj, size, cluster, everything()) %>%
    arrange(padj)
  
  write.csv(combined,
            file.path(output_dir, "interaction_ALL_reduced_terms.csv"),
            row.names = FALSE)
  
  cat("\n\nCombined table saved -", nrow(combined), "total rows\n")
  print(combined %>% dplyr::select(ontology, go, term, NES, padj))
  
} else {
  cat("\nNo ontologies had enough terms to cluster.\n")
}

##############################################################
# STEP 5: Supplementary table - all significant terms
# (Run this independently any time, as long as 'combined' exists
# in your session from Step 4 above)
##############################################################

supp_table <- combined %>%
  dplyr::select(ontology, go, term, NES, padj, size) %>%
  arrange(padj)

write.csv(supp_table,
          file.path(output_dir, "Supplementary_Table_Interaction_fgsea_terms.csv"),
          row.names = FALSE)

cat("\nSupplementary table saved -", nrow(supp_table), "rows\n")

##############################################################
# STEP 6: Curated table for main text - grouped by theme
# (Run this independently any time, as long as 'supp_table' exists
# in your session from Step 5 above)
##############################################################
# EDIT THIS LIST if your interpretation of the themes changes

curated_ids_themes <- c(
  "GO:0006457" = "Protein folding / ER-mitochondrial",
  "GO:0005789" = "Protein folding / ER-mitochondrial",
  "GO:0005783" = "Protein folding / ER-mitochondrial",
  "GO:0005743" = "Protein folding / ER-mitochondrial",
  "GO:0005634" = "Protein folding / ER-mitochondrial",
  "GO:0005829" = "Protein folding / ER-mitochondrial",
  "GO:0016491" = "Redox / proteolytic remodeling",
  "GO:0006508" = "Redox / proteolytic remodeling",
  "GO:0006869" = "Redox / proteolytic remodeling",
  "GO:0031012" = "Redox / proteolytic remodeling",
  "GO:0005764" = "Redox / proteolytic remodeling",
  "GO:0009986" = "Redox / proteolytic remodeling",
  "GO:0005524" = "Secondary / general cellular machinery",
  "GO:0000398" = "Secondary / general cellular machinery",
  "GO:0003924" = "Secondary / general cellular machinery",
  "GO:0003682" = "Secondary / general cellular machinery",
  "GO:0003729" = "Secondary / general cellular machinery",
  "GO:0005525" = "Secondary / general cellular machinery",
  "GO:0005730" = "Secondary / general cellular machinery",
  "GO:0003723" = "Secondary / general cellular machinery",
  "GO:0008017" = "Secondary / general cellular machinery",
  "GO:0000139" = "Secondary / general cellular machinery"
)

theme_order <- c("Protein folding / ER-mitochondrial",
                 "Redox / proteolytic remodeling",
                 "Secondary / general cellular machinery")

curated_table <- supp_table %>%
  filter(go %in% names(curated_ids_themes)) %>%
  mutate(theme = curated_ids_themes[go],
         theme = factor(theme, levels = theme_order)) %>%
  dplyr::select(theme, go, term, ontology, NES, padj, size) %>%
  arrange(theme, padj)

write.csv(curated_table,
          file.path(output_dir, "Table_3_2_Interaction_fgsea_curated.csv"),
          row.names = FALSE)

cat("\nCurated table saved -", nrow(curated_table), "rows across",
    length(unique(curated_table$theme)), "themes\n")
print(curated_table)