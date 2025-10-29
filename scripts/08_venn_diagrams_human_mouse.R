#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 8: Venn Diagrams - Human vs Mouse DEGs
# This script creates Venn diagrams comparing upregulated genes between
# human and mouse neurons/NPCs with ortholog conversion

suppressPackageStartupMessages({
    library(tidyverse)
    library(VennDiagram)
    library(biomaRt)
    library(RColorBrewer)
})

cat("==========================================\n")
cat("Creating Venn Diagrams: Human vs Mouse\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ============================================================
# Configuration
# ============================================================

# File paths
HUMAN_NEURON_FILE <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results/05_deseq2/Neuron_differential_expression_upregulated.csv"
HUMAN_NPC_FILE <- "results/DESeq2/DESeq2_results_significant.csv"
MOUSE_NEURON_FILE <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_CUTandTAG/DATA/DEA_NEU.csv"
MOUSE_NPC_FILE <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_CUTandTAG/DATA/DEA_NSC.csv"

# Output directory
OUTPUT_DIR <- "results/venn_diagrams"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Thresholds (consistent across all datasets)
LOG2FC_THRESHOLD <- 0.5
PADJ_THRESHOLD <- 0.05

# ============================================================
# Function: Convert mouse genes to human orthologs
# ============================================================
convert_mouse_to_human <- function(mouse_genes) {
    cat("  Converting", length(mouse_genes), "mouse genes to human orthologs...\n")

    # Method 0: Check for local ortholog mapping file first
    local_ortholog_file <- file.path(OUTPUT_DIR, "mouse_human_orthologs.csv")
    if (file.exists(local_ortholog_file)) {
        cat("    Found local ortholog mapping file:", local_ortholog_file, "\n")
        tryCatch({
            all_orthologs <- read.csv(local_ortholog_file, stringsAsFactors = FALSE)

            # Filter to requested genes
            orthologs <- all_orthologs %>%
                filter(mouse_gene %in% mouse_genes) %>%
                select(mouse_gene, human_gene) %>%
                distinct()

            cat("    Found", nrow(orthologs), "ortholog mappings from local file\n")
            cat("    Unique human orthologs:", length(unique(orthologs$human_gene)), "\n")

            return(orthologs)
        }, error = function(e) {
            cat("    ERROR reading local file:", conditionMessage(e), "\n")
            cat("    Trying other methods...\n")
        })
    }

    # Method 1: Try biomaRt with multiple mirror attempts
    cat("    Attempting biomaRt connection...\n")

    mirrors <- c(
        "https://www.ensembl.org",
        "https://useast.ensembl.org",
        "https://uswest.ensembl.org",
        "https://asia.ensembl.org"
    )

    for (mirror in mirrors) {
        cat("      Trying mirror:", mirror, "\n")
        orthologs <- tryCatch({
            mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = mirror)
            human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = mirror)

            # Query for orthologs
            result <- getLDS(
                attributes = c("mgi_symbol"),
                filters = "mgi_symbol",
                values = mouse_genes,
                mart = mouse_mart,
                attributesL = c("hgnc_symbol"),
                martL = human_mart,
                uniqueRows = TRUE
            )

            colnames(result) <- c("mouse_gene", "human_gene")

            # Remove empty human gene symbols
            result <- result %>%
                filter(human_gene != "", !is.na(human_gene))

            cat("      SUCCESS! Found", nrow(result), "ortholog mappings\n")
            cat("      Unique human orthologs:", length(unique(result$human_gene)), "\n")

            return(result)

        }, error = function(e) {
            cat("      FAILED:", conditionMessage(e), "\n")
            return(NULL)
        })

        if (!is.null(orthologs) && nrow(orthologs) > 0) {
            return(orthologs)
        }
    }

    # Method 2: Try downloading ortholog data directly from Ensembl
    cat("    biomaRt failed. Attempting to download ortholog data from Ensembl FTP...\n")
    orthologs <- tryCatch({
        download_ensembl_orthologs(mouse_genes)
    }, error = function(e) {
        cat("      FAILED:", conditionMessage(e), "\n")
        return(NULL)
    })

    if (!is.null(orthologs) && nrow(orthologs) > 0) {
        return(orthologs)
    }

    # Method 3: Simple capitalization heuristic for common genes
    cat("    All methods failed. Using capitalization heuristic (mouse lowercase -> human uppercase)...\n")
    cat("    WARNING: This is approximate and may not be accurate for all genes!\n")
    cat("    For better results, manually download ortholog data from Ensembl BioMart\n")

    # Create ortholog mapping using simple case conversion
    # Most mouse genes: Mecp2, Gapdh, etc.
    # Corresponding human genes: MECP2, GAPDH, etc.
    orthologs <- data.frame(
        mouse_gene = mouse_genes,
        human_gene = toupper(mouse_genes),
        stringsAsFactors = FALSE
    ) %>%
        # Remove empty entries
        filter(mouse_gene != "", human_gene != "")

    cat("    Converted", nrow(orthologs), "genes using capitalization\n")
    cat("    Mouse genes (first 10):", paste(head(orthologs$mouse_gene, 10), collapse = ", "), "\n")
    cat("    Human genes (first 10):", paste(head(orthologs$human_gene, 10), collapse = ", "), "\n")

    if (nrow(orthologs) == 0) {
        warning("Capitalization method returned 0 genes!")
    }

    return(orthologs)
}

# ============================================================
# Function: Download ortholog data from Ensembl FTP
# ============================================================
download_ensembl_orthologs <- function(mouse_genes) {
    cache_file <- file.path(OUTPUT_DIR, "ensembl_mouse_human_orthologs_cache.rds")

    # Check if we have a cached version
    if (file.exists(cache_file)) {
        cat("      Found cached ortholog data\n")
        all_orthologs <- readRDS(cache_file)
    } else {
        cat("      Downloading ortholog data...\n")

        # Try to download from Ensembl BioMart web service directly
        url <- "https://www.ensembl.org/biomart/martservice"

        # XML query for mouse-human orthologs
        query <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="mmusculus_gene_ensembl" interface="default">
<Attribute name="mgi_symbol"/>
<Attribute name="hsapiens_homolog_associated_gene_name"/>
</Dataset>
</Query>'

        # Download
        temp_file <- tempfile()
        download.file(
            url = paste0(url, "?query=", URLencode(query)),
            destfile = temp_file,
            quiet = FALSE,
            method = "auto"
        )

        # Read the downloaded file
        all_orthologs <- read.delim(temp_file, stringsAsFactors = FALSE)
        colnames(all_orthologs) <- c("mouse_gene", "human_gene")

        # Remove empty entries
        all_orthologs <- all_orthologs %>%
            filter(mouse_gene != "", human_gene != "", !is.na(human_gene))

        # Cache for future use
        saveRDS(all_orthologs, cache_file)
        cat("      Cached ortholog data for future use\n")

        unlink(temp_file)
    }

    # Filter to requested genes
    orthologs <- all_orthologs %>%
        filter(mouse_gene %in% mouse_genes)

    cat("      Found", nrow(orthologs), "ortholog mappings\n")

    return(orthologs)
}

# ============================================================
# Function: Load and filter human genes
# ============================================================
load_human_genes <- function(file_path, source_name, prefiltered = FALSE) {
    cat("\nLoading", source_name, "data...\n")
    cat("  File:", file_path, "\n")

    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    data <- read.csv(file_path)
    cat("  Loaded", nrow(data), "genes\n")

    if (prefiltered) {
        # File is already filtered (upregulated only)
        genes <- data$gene_symbol
        cat("  Pre-filtered file - using all", length(genes), "genes\n")
    } else {
        # Need to filter for upregulated genes
        genes <- data %>%
            filter(log2FoldChange > LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD) %>%
            pull(gene_symbol)
        cat("  Filtered to", length(genes), "upregulated genes (log2FC >", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, ")\n")
    }

    # Remove empty or NA gene symbols
    genes <- genes[!is.na(genes) & genes != ""]
    cat("  Final gene count:", length(genes), "\n")

    return(genes)
}

# ============================================================
# Function: Load and filter mouse genes, convert to human
# ============================================================
load_mouse_genes <- function(file_path, source_name) {
    cat("\nLoading", source_name, "data...\n")
    cat("  File:", file_path, "\n")

    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(data), "total genes\n")
    cat("  Column names:", paste(colnames(data), collapse = ", "), "\n")

    # Debug: Check data structure
    cat("  First few log2FoldChange values:", paste(head(data$log2FoldChange, 5), collapse = ", "), "\n")
    cat("  First few padj values:", paste(head(data$padj, 5), collapse = ", "), "\n")

    # Filter for upregulated genes
    filtered_data <- data %>%
        filter(!is.na(log2FoldChange), !is.na(padj)) %>%
        filter(log2FoldChange > LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD)

    cat("  After filtering: ", nrow(filtered_data), "genes meet criteria\n")

    if (nrow(filtered_data) == 0) {
        warning("No mouse genes pass the filtering criteria!")
        return(character())
    }

    mouse_genes <- filtered_data %>% pull(gene)
    cat("  Mouse genes passing filter:", length(mouse_genes), "\n")
    cat("  First 10 mouse genes:", paste(head(mouse_genes, 10), collapse = ", "), "\n")

    # Convert to human orthologs
    orthologs <- convert_mouse_to_human(mouse_genes)

    if (is.null(orthologs) || nrow(orthologs) == 0) {
        warning("No orthologs found for ", source_name)
        return(character())
    }

    cat("  Ortholog conversion returned", nrow(orthologs), "mappings\n")

    # Get unique human gene symbols
    human_genes <- unique(orthologs$human_gene)
    # Remove empty or NA genes
    human_genes <- human_genes[!is.na(human_genes) & human_genes != ""]

    cat("  Final count:", length(human_genes), "unique human gene symbols\n")
    cat("  First 10 human orthologs:", paste(head(human_genes, 10), collapse = ", "), "\n")

    # Save ortholog mapping
    mapping_file <- file.path(OUTPUT_DIR, paste0(gsub(" ", "_", source_name), "_ortholog_mapping.csv"))
    write.csv(orthologs, mapping_file, row.names = FALSE)
    cat("  Saved ortholog mapping to:", mapping_file, "\n")

    return(human_genes)
}

# ============================================================
# Function: Create Venn diagram
# ============================================================
create_venn_diagram <- function(gene_list1, gene_list2,
                                label1, label2,
                                output_file, title) {

    cat("\nCreating Venn diagram:", title, "\n")
    cat("  ", label1, "genes:", length(gene_list1), "\n")
    cat("  ", label2, "genes:", length(gene_list2), "\n")

    # Calculate overlaps
    overlap <- intersect(gene_list1, gene_list2)
    unique1 <- setdiff(gene_list1, gene_list2)
    unique2 <- setdiff(gene_list2, gene_list1)

    cat("  Overlap:", length(overlap), "genes\n")
    cat("  Unique to", label1, ":", length(unique1), "genes\n")
    cat("  Unique to", label2, ":", length(unique2), "genes\n")

    # Save overlap genes to file
    overlap_file <- sub("\\.png$", "_overlap.txt", output_file)
    writeLines(sort(overlap), overlap_file)
    cat("  Saved overlapping genes to:", overlap_file, "\n")

    # Save unique genes to files
    unique1_file <- sub("\\.png$", paste0("_unique_", gsub(" ", "_", label1), ".txt"), output_file)
    writeLines(sort(unique1), unique1_file)
    cat("  Saved unique", label1, "genes to:", unique1_file, "\n")

    unique2_file <- sub("\\.png$", paste0("_unique_", gsub(" ", "_", label2), ".txt"), output_file)
    writeLines(sort(unique2), unique2_file)
    cat("  Saved unique", label2, "genes to:", unique2_file, "\n")

    # Create Venn diagram
    venn_plot <- venn.diagram(
        x = list(gene_list1, gene_list2),
        category.names = c(label1, label2),
        filename = NULL,
        output = TRUE,

        # Appearance
        fill = c("#e74c3c", "#3498db"),
        alpha = 0.5,

        # Numbers
        cex = 1.5,
        fontfamily = "sans",

        # Category labels
        cat.cex = 1.3,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.fontfamily = "sans",

        # Main title
        main = title,
        main.cex = 1.5,
        main.fontface = "bold",
        main.fontfamily = "sans",

        # Other
        margin = 0.1
    )

    # Save as PNG
    png(output_file, width = 2400, height = 2400, res = 300)
    grid.draw(venn_plot)
    dev.off()
    cat("  Saved Venn diagram to:", output_file, "\n")

    # Save as PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    pdf(pdf_file, width = 8, height = 8)
    grid.draw(venn_plot)
    dev.off()
    cat("  Saved Venn diagram to:", pdf_file, "\n")

    # Return summary statistics
    return(list(
        total1 = length(gene_list1),
        total2 = length(gene_list2),
        overlap = length(overlap),
        unique1 = length(unique1),
        unique2 = length(unique2),
        overlap_genes = overlap
    ))
}

# ============================================================
# Main Analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 1: Loading Human Gene Data\n")
cat("============================================================\n")

# Load human neuron genes (pre-filtered file)
human_neuron_genes <- load_human_genes(
    HUMAN_NEURON_FILE,
    "Human Neurons (old data)",
    prefiltered = TRUE
)

# Load human NPC genes (need filtering)
human_npc_genes <- load_human_genes(
    HUMAN_NPC_FILE,
    "Human NPCs (current data)",
    prefiltered = FALSE
)

cat("\n============================================================\n")
cat("STEP 2: Loading Mouse Gene Data and Converting to Human Orthologs\n")
cat("============================================================\n")

# Load and convert mouse neuron genes
mouse_neuron_genes_human <- load_mouse_genes(
    MOUSE_NEURON_FILE,
    "Mouse Neurons"
)

# Load and convert mouse NPC genes
mouse_npc_genes_human <- load_mouse_genes(
    MOUSE_NPC_FILE,
    "Mouse NPCs"
)

cat("\n============================================================\n")
cat("STEP 3: Creating Venn Diagrams\n")
cat("============================================================\n")

# Venn diagram 1: Human Neurons vs Mouse Neurons
neuron_stats <- create_venn_diagram(
    gene_list1 = human_neuron_genes,
    gene_list2 = mouse_neuron_genes_human,
    label1 = "Human Neurons",
    label2 = "Mouse Neurons",
    output_file = file.path(OUTPUT_DIR, "venn_neurons_human_vs_mouse.png"),
    title = "Upregulated Genes: Human vs Mouse Neurons\n(log2FC > 0.5, padj < 0.05)"
)

# Venn diagram 2: Human NPCs vs Mouse NPCs
npc_stats <- create_venn_diagram(
    gene_list1 = human_npc_genes,
    gene_list2 = mouse_npc_genes_human,
    label1 = "Human NPCs",
    label2 = "Mouse NPCs",
    output_file = file.path(OUTPUT_DIR, "venn_npcs_human_vs_mouse.png"),
    title = "Upregulated Genes: Human vs Mouse NPCs\n(log2FC > 0.5, padj < 0.05)"
)

# ============================================================
# Summary Report
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY REPORT\n")
cat("============================================================\n\n")

cat("NEURONS (Human vs Mouse):\n")
cat("  Human Neurons:", neuron_stats$total1, "upregulated genes\n")
cat("  Mouse Neurons:", neuron_stats$total2, "upregulated genes (human orthologs)\n")
cat("  Overlap:", neuron_stats$overlap, "genes (",
    round(neuron_stats$overlap / min(neuron_stats$total1, neuron_stats$total2) * 100, 1),
    "% of smaller set)\n")
cat("  Unique to Human:", neuron_stats$unique1, "genes\n")
cat("  Unique to Mouse:", neuron_stats$unique2, "genes\n\n")

cat("NPCs (Human vs Mouse):\n")
cat("  Human NPCs:", npc_stats$total1, "upregulated genes\n")
cat("  Mouse NPCs:", npc_stats$total2, "upregulated genes (human orthologs)\n")
cat("  Overlap:", npc_stats$overlap, "genes (",
    round(npc_stats$overlap / min(npc_stats$total1, npc_stats$total2) * 100, 1),
    "% of smaller set)\n")
cat("  Unique to Human:", npc_stats$unique1, "genes\n")
cat("  Unique to Mouse:", npc_stats$unique2, "genes\n\n")

# Save summary to file
summary_file <- file.path(OUTPUT_DIR, "venn_summary_report.txt")
sink(summary_file)
cat("VENN DIAGRAM ANALYSIS SUMMARY\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("==========================================\n\n")

cat("FILTERING CRITERIA:\n")
cat("  log2FoldChange >", LOG2FC_THRESHOLD, "\n")
cat("  padj <", PADJ_THRESHOLD, "\n\n")

cat("NEURONS (Human vs Mouse):\n")
cat("  Human Neurons:", neuron_stats$total1, "upregulated genes\n")
cat("  Mouse Neurons:", neuron_stats$total2, "upregulated genes (human orthologs)\n")
cat("  Overlap:", neuron_stats$overlap, "genes (",
    round(neuron_stats$overlap / min(neuron_stats$total1, neuron_stats$total2) * 100, 1),
    "% of smaller set)\n")
cat("  Unique to Human:", neuron_stats$unique1, "genes\n")
cat("  Unique to Mouse:", neuron_stats$unique2, "genes\n\n")

cat("NPCs (Human vs Mouse):\n")
cat("  Human NPCs:", npc_stats$total1, "upregulated genes\n")
cat("  Mouse NPCs:", npc_stats$total2, "upregulated genes (human orthologs)\n")
cat("  Overlap:", npc_stats$overlap, "genes (",
    round(npc_stats$overlap / min(npc_stats$total1, npc_stats$total2) * 100, 1),
    "% of smaller set)\n")
cat("  Unique to Human:", npc_stats$unique1, "genes\n")
cat("  Unique to Mouse:", npc_stats$unique2, "genes\n\n")

if (length(neuron_stats$overlap_genes) > 0) {
    cat("OVERLAPPING NEURON GENES (", length(neuron_stats$overlap_genes), "):\n", sep = "")
    cat(paste(sort(neuron_stats$overlap_genes), collapse = ", "), "\n\n")
}

if (length(npc_stats$overlap_genes) > 0 && length(npc_stats$overlap_genes) <= 100) {
    cat("OVERLAPPING NPC GENES (first 100):\n")
    cat(paste(sort(npc_stats$overlap_genes)[1:min(100, length(npc_stats$overlap_genes))], collapse = ", "), "\n\n")
}

sink()
cat("Saved summary report to:", summary_file, "\n")

cat("\n==========================================\n")
cat("Venn diagram analysis complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Generated files:\n")
cat("  - venn_neurons_human_vs_mouse.png/pdf\n")
cat("  - venn_npcs_human_vs_mouse.png/pdf\n")
cat("  - *_overlap.txt (overlapping genes)\n")
cat("  - *_unique_*.txt (unique genes)\n")
cat("  - *_ortholog_mapping.csv (mouse-to-human mappings)\n")
cat("  - venn_summary_report.txt (full summary)\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
