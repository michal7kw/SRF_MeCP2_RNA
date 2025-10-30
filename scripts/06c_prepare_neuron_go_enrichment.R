#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 6c: Prepare Neuron GO/KEGG Enrichment Files
# ====================================================================
# This is a ONE-TIME SETUP script that creates GO and KEGG enrichment
# files for neuron data from the previous analysis.
#
# PURPOSE:
#   - Uses the SAME methodology as NPC enrichment (06_functional_enrichment.R)
#   - Ensures comparable results between NPCs and Neurons
#   - Creates both GO and KEGG enrichment for up/downregulated genes
#
# PREREQUISITES:
#   - Neuron differential expression file must exist:
#     /beegfs/.../SRF_MeCP2_rna_old/results/05_deseq2/Neuron_differential_expression.csv
#
# OUTPUT:
#   - GO_enrichment_upregulated.csv
#   - GO_enrichment_downregulated.csv
#   - KEGG_enrichment_upregulated.csv
#   - KEGG_enrichment_downregulated.csv
#   - (All saved in neuron GO directory)
#
# NOTE: This script mirrors the NPC enrichment analysis (step 06) to ensure
#       results are directly comparable across cell types.

suppressPackageStartupMessages({
    library(tidyverse)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

cat("=========================================================\n")
cat("Neuron GO/KEGG Enrichment Analysis\n")
cat("(Mirroring NPC Analysis Parameters)\n")
cat("=========================================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ============================================================
# Configuration - MATCHING NPC ANALYSIS
# ============================================================

# Input paths (neuron data from old analysis)
NEURON_BASE_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results"
neuron_deseq_file <- file.path(NEURON_BASE_DIR, "05_deseq2", "Neuron_differential_expression.csv")
neuron_go_dir <- file.path(NEURON_BASE_DIR, "go_analysis_neu")

# Parameters - SAME AS NPC ANALYSIS (from 06_functional_enrichment.R)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1              # |log2FC| > 1
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.2

cat("Analysis Parameters (matching NPC analysis):\n")
cat("  - Adjusted p-value cutoff:", PADJ_CUTOFF, "\n")
cat("  - Log2 fold change cutoff:", LFC_CUTOFF, "\n")
cat("  - GO/KEGG p-value cutoff:", PVALUE_CUTOFF, "\n")
cat("  - GO/KEGG q-value cutoff:", QVALUE_CUTOFF, "\n\n")

# ============================================================
# Validation
# ============================================================

cat("Validating inputs...\n")

if (!file.exists(neuron_deseq_file)) {
    stop("ERROR: Neuron DESeq2 file not found: ", neuron_deseq_file)
}

if (!dir.exists(neuron_go_dir)) {
    cat("  Creating output directory:", neuron_go_dir, "\n")
    dir.create(neuron_go_dir, recursive = TRUE)
}

cat("  ✓ DESeq2 results found:", basename(neuron_deseq_file), "\n")
cat("  ✓ Output directory ready\n\n")

# ============================================================
# 1. Load DESeq2 Results
# ============================================================

cat("============================================================\n")
cat("STEP 1: Loading Neuron DESeq2 Results\n")
cat("============================================================\n")

deseq_results <- read.csv(neuron_deseq_file)
cat("  Loaded", nrow(deseq_results), "genes\n")

# Check if gene_id column exists
if (!"gene_id" %in% colnames(deseq_results)) {
    stop("ERROR: 'gene_id' column not found in DESeq2 results")
}

# Set gene_id as rownames
rownames(deseq_results) <- deseq_results$gene_id

# Remove version numbers from Ensembl IDs (e.g., ENSG00000000003.15 -> ENSG00000000003)
rownames(deseq_results) <- gsub("\\..*", "", rownames(deseq_results))

cat("  Column names:", paste(colnames(deseq_results), collapse = ", "), "\n\n")

# ============================================================
# 2. Prepare Gene Lists
# ============================================================

cat("============================================================\n")
cat("STEP 2: Filtering Significant Genes\n")
cat("============================================================\n")

# Get significant genes
sig_genes <- deseq_results %>%
    filter(padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_CUTOFF) %>%
    rownames()

# Upregulated genes
up_genes <- deseq_results %>%
    filter(padj < PADJ_CUTOFF, log2FoldChange > LFC_CUTOFF) %>%
    rownames()

# Downregulated genes
down_genes <- deseq_results %>%
    filter(padj < PADJ_CUTOFF, log2FoldChange < -LFC_CUTOFF) %>%
    rownames()

cat("  Significant genes (padj <", PADJ_CUTOFF, ", |log2FC| >", LFC_CUTOFF, "):", length(sig_genes), "\n")
cat("    - Upregulated:", length(up_genes), "\n")
cat("    - Downregulated:", length(down_genes), "\n\n")

if (length(sig_genes) == 0) {
    stop("ERROR: No significant genes found with current thresholds!")
}

# Convert Ensembl IDs to Entrez IDs
cat("  Converting Ensembl IDs to Entrez IDs...\n")
sig_genes_entrez <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
up_genes_entrez <- bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_genes_entrez <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("  Converted", nrow(sig_genes_entrez), "genes to Entrez IDs\n")
cat("    - Upregulated:", nrow(up_genes_entrez), "\n")
cat("    - Downregulated:", nrow(down_genes_entrez), "\n\n")

# ============================================================
# 3. Gene Ontology Enrichment Analysis
# ============================================================

cat("============================================================\n")
cat("STEP 3: Running GO Enrichment Analysis\n")
cat("============================================================\n")

# GO enrichment for upregulated genes
cat("  Analyzing upregulated genes...\n")
go_up <- enrichGO(
    gene = up_genes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF,
    readable = TRUE
)

# GO enrichment for downregulated genes
cat("  Analyzing downregulated genes...\n")
go_down <- enrichGO(
    gene = down_genes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF,
    readable = TRUE
)

# Save GO results
cat("\n  Saving GO enrichment results...\n")

if (!is.null(go_up) && nrow(go_up) > 0) {
    go_up_file <- file.path(neuron_go_dir, "GO_enrichment_upregulated.csv")
    write.csv(as.data.frame(go_up), go_up_file, row.names = FALSE)
    cat("    ✓ GO upregulated:", nrow(go_up), "terms →", basename(go_up_file), "\n")
} else {
    cat("    ⚠ No enriched GO terms for upregulated genes\n")
}

if (!is.null(go_down) && nrow(go_down) > 0) {
    go_down_file <- file.path(neuron_go_dir, "GO_enrichment_downregulated.csv")
    write.csv(as.data.frame(go_down), go_down_file, row.names = FALSE)
    cat("    ✓ GO downregulated:", nrow(go_down), "terms →", basename(go_down_file), "\n")
} else {
    cat("    ⚠ No enriched GO terms for downregulated genes\n")
}

# ============================================================
# 4. KEGG Pathway Enrichment Analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 4: Running KEGG Enrichment Analysis\n")
cat("============================================================\n")

# KEGG enrichment for upregulated genes
cat("  Analyzing upregulated genes...\n")
kegg_up <- enrichKEGG(
    gene = up_genes_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF
)

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# KEGG enrichment for downregulated genes
cat("  Analyzing downregulated genes...\n")
kegg_down <- enrichKEGG(
    gene = down_genes_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF
)

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# Save KEGG results
cat("\n  Saving KEGG enrichment results...\n")

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up_file <- file.path(neuron_go_dir, "KEGG_enrichment_upregulated.csv")
    write.csv(as.data.frame(kegg_up), kegg_up_file, row.names = FALSE)
    cat("    ✓ KEGG upregulated:", nrow(kegg_up), "pathways →", basename(kegg_up_file), "\n")
} else {
    cat("    ⚠ No enriched KEGG pathways for upregulated genes\n")
}

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down_file <- file.path(neuron_go_dir, "KEGG_enrichment_downregulated.csv")
    write.csv(as.data.frame(kegg_down), kegg_down_file, row.names = FALSE)
    cat("    ✓ KEGG downregulated:", nrow(kegg_down), "pathways →", basename(kegg_down_file), "\n")
} else {
    cat("    ⚠ No enriched KEGG pathways for downregulated genes\n")
}

# ============================================================
# Final Summary
# ============================================================

cat("\n=========================================================\n")
cat("NEURON ENRICHMENT ANALYSIS COMPLETE\n")
cat("=========================================================\n\n")

cat("Output directory:", neuron_go_dir, "\n\n")

cat("Generated files:\n")

# List all output files
output_files <- c(
    "GO_enrichment_upregulated.csv",
    "GO_enrichment_downregulated.csv",
    "KEGG_enrichment_upregulated.csv",
    "KEGG_enrichment_downregulated.csv"
)

for (f in output_files) {
    full_path <- file.path(neuron_go_dir, f)
    if (file.exists(full_path)) {
        file_size <- file.size(full_path)
        cat("  ✓", f, sprintf("(%.1f KB)\n", file_size / 1024))
    } else {
        cat("  ✗", f, "(not created)\n")
    }
}

cat("\nThese files are now ready for cross-cell-type comparison:\n")
cat("  → scripts/07_compare_celltype_enrichment.R\n\n")

cat("Parameters used (matching NPC analysis):\n")
cat("  - padj <", PADJ_CUTOFF, ", |log2FC| >", LFC_CUTOFF, "\n")
cat("  - Organism: Homo sapiens (org.Hs.eg.db)\n")
cat("  - GO ontologies: All (BP, MF, CC)\n")
cat("  - KEGG organism: hsa (Homo sapiens)\n\n")

cat("End time:", as.character(Sys.time()), "\n")
cat("=========================================================\n")
