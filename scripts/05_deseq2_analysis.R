#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 5: Differential Expression Analysis with DESeq2
# This script performs differential expression analysis between Control and Mutant conditions

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
})

cat("==========================================\n")
cat("Starting DESeq2 Differential Expression Analysis\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths
count_file <- "results/counts/gene_counts.txt"
metadata_file <- "metadata/sample_info.csv"
output_dir <- "results/DESeq2"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load and prepare data
# ============================================================
cat("Loading count data...\n")
count_data <- read.table(count_file, header = TRUE, row.names = 1, skip = 1)

# Remove annotation columns (Chr, Start, End, Strand, Length)
count_data <- count_data[, -c(1:5)]

# Clean sample names
colnames(count_data) <- gsub("results.alignment.", "", colnames(count_data))
colnames(count_data) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(count_data))

cat("Loaded counts for", nrow(count_data), "genes across", ncol(count_data), "samples\n")

# Load metadata
cat("Loading sample metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1)
metadata$condition <- factor(metadata$condition, levels = c("Control", "Mutant"))
metadata$replicate <- factor(metadata$replicate)

cat("Sample information:\n")
print(metadata)
cat("\n")

# Verify sample order matches
if (!all(rownames(metadata) == colnames(count_data))) {
    stop("Sample names in count data and metadata do not match!")
}

# ============================================================
# 2. Create DESeq2 object
# ============================================================
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = metadata,
    design = ~ condition
)

# Pre-filtering: keep genes with at least 10 reads total
cat("Pre-filtering low count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Retained", sum(keep), "genes after filtering\n\n")

# ============================================================
# 3. Run DESeq2 analysis
# ============================================================
cat("Running DESeq2 differential expression analysis...\n")
dds <- DESeq(dds)

# Get results
cat("Extracting results...\n")
res <- results(dds, contrast = c("condition", "Mutant", "Control"))

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Summary of results
cat("\nDifferential Expression Summary:\n")
cat("================================\n")
summary(res)
cat("\n")

# Count significant genes
sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat("Significant genes (padj < 0.05):", sig_genes, "\n")
cat("  - Upregulated in Mutant:", sig_up, "\n")
cat("  - Downregulated in Mutant:", sig_down, "\n\n")

# ============================================================
# 4. Save results
# ============================================================
cat("Saving results...\n")

# Save all results
write.csv(
    as.data.frame(res_ordered),
    file = file.path(output_dir, "DESeq2_results_all.csv"),
    row.names = TRUE
)

# Save significant genes only (padj < 0.05)
res_sig <- subset(res_ordered, padj < 0.05)
write.csv(
    as.data.frame(res_sig),
    file = file.path(output_dir, "DESeq2_results_significant.csv"),
    row.names = TRUE
)

# Save significant genes with stricter criteria (padj < 0.05 and |log2FC| > 1)
res_sig_strict <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(
    as.data.frame(res_sig_strict),
    file = file.path(output_dir, "DESeq2_results_significant_FC2.csv"),
    row.names = TRUE
)

cat("Significant genes (padj < 0.05, |log2FC| > 1):", nrow(res_sig_strict), "\n")

# ============================================================
# 5. Normalized counts
# ============================================================
cat("\nExtracting normalized counts...\n")
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(
    normalized_counts,
    file = file.path(output_dir, "normalized_counts.csv"),
    row.names = TRUE
)

# ============================================================
# 6. Variance stabilized transformation for visualization
# ============================================================
cat("Performing variance stabilizing transformation...\n")
vsd <- vst(dds, blind = FALSE)
write.csv(
    assay(vsd),
    file = file.path(output_dir, "vst_counts.csv"),
    row.names = TRUE
)

# Save DESeq2 object
cat("Saving DESeq2 object...\n")
saveRDS(dds, file = file.path(output_dir, "dds.rds"))
saveRDS(vsd, file = file.path(output_dir, "vsd.rds"))
saveRDS(res, file = file.path(output_dir, "results.rds"))

cat("\nDESeq2 analysis complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nKey output files:\n")
cat("  - DESeq2_results_all.csv: All genes with statistics\n")
cat("  - DESeq2_results_significant.csv: Significant genes (padj < 0.05)\n")
cat("  - DESeq2_results_significant_FC2.csv: Significant genes (padj < 0.05, |log2FC| > 1)\n")
cat("  - normalized_counts.csv: DESeq2 normalized counts\n")
cat("  - vst_counts.csv: Variance stabilized counts\n")
cat("  - *.rds: R objects for downstream analysis\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
