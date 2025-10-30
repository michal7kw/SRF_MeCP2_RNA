#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 9: Venn Diagrams with Proper MGI Orthologs
# This script uses the official MGI (Mouse Genome Informatics) ortholog database
# for accurate mouse-to-human gene conversion

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggVennDiagram)
    library(RColorBrewer)
    library(gridExtra)
    library(grid)
})

cat("==========================================\n")
cat("Creating Venn Diagrams with MGI Orthologs\n")
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

# MGI ortholog database
MGI_ORTHOLOG_FILE <- "results/venn_diagrams_MGI/ortholog_db/HOM_MouseHumanSequence.rpt"

# Output directory
OUTPUT_DIR <- "results/venn_diagrams_MGI"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Thresholds
LOG2FC_THRESHOLD <- 0.5
PADJ_THRESHOLD <- 0.05

# ============================================================
# Function: Load and parse MGI ortholog database
# ============================================================
load_mgi_orthologs <- function() {
    cat("Loading MGI ortholog database...\n")
    cat("  File:", MGI_ORTHOLOG_FILE, "\n")

    if (!file.exists(MGI_ORTHOLOG_FILE)) {
        stop("MGI ortholog file not found: ", MGI_ORTHOLOG_FILE,
             "\nPlease download from: http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
    }

    # Read the MGI file
    mgi_data <- read.delim(MGI_ORTHOLOG_FILE, stringsAsFactors = FALSE, quote = "")

    cat("  Loaded", nrow(mgi_data), "records\n")

    # Extract mouse genes
    mouse_genes <- mgi_data %>%
        filter(Common.Organism.Name == "mouse, laboratory") %>%
        select(DB.Class.Key, mouse_gene = Symbol)

    # Extract human genes
    human_genes <- mgi_data %>%
        filter(Common.Organism.Name == "human") %>%
        select(DB.Class.Key, human_gene = Symbol)

    # Join mouse and human genes by orthology group
    orthologs <- mouse_genes %>%
        inner_join(human_genes, by = "DB.Class.Key") %>%
        select(mouse_gene, human_gene) %>%
        distinct() %>%
        filter(mouse_gene != "", human_gene != "")

    cat("  Created", nrow(orthologs), "mouse-human ortholog pairs\n")
    cat("  Unique mouse genes:", length(unique(orthologs$mouse_gene)), "\n")
    cat("  Unique human genes:", length(unique(orthologs$human_gene)), "\n")

    # Save the mapping for reference
    mapping_file <- file.path(OUTPUT_DIR, "MGI_mouse_human_orthologs.csv")
    write.csv(orthologs, mapping_file, row.names = FALSE)
    cat("  Saved ortholog mapping to:", mapping_file, "\n\n")

    return(orthologs)
}

# ============================================================
# Function: Convert mouse genes to human using MGI database
# ============================================================
convert_mouse_to_human_mgi <- function(mouse_genes, mgi_orthologs) {
    cat("  Converting", length(mouse_genes), "mouse genes to human orthologs using MGI database...\n")

    # Match with MGI orthologs
    orthologs <- mgi_orthologs %>%
        filter(mouse_gene %in% mouse_genes)

    # Get genes without orthologs
    no_ortholog <- setdiff(mouse_genes, orthologs$mouse_gene)

    cat("    Found orthologs for:", nrow(orthologs), "genes\n")
    cat("    No ortholog found for:", length(no_ortholog), "genes\n")

    if (length(no_ortholog) > 0) {
        cat("    Genes without orthologs (first 20):\n      ",
            paste(head(no_ortholog, 20), collapse = ", "), "\n")
    }

    return(orthologs)
}

# ============================================================
# Function: Load and filter human genes
# ============================================================
load_human_genes <- function(file_path, source_name, prefiltered = FALSE, direction = "up") {
    cat("\nLoading", source_name, "data...\n")
    cat("  File:", file_path, "\n")
    cat("  Direction:", direction, "\n")

    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    data <- read.csv(file_path)
    cat("  Loaded", nrow(data), "genes\n")

    if (prefiltered) {
        genes <- data$gene_symbol
        cat("  Pre-filtered file - using all", length(genes), "genes\n")
    } else {
        if (direction == "up") {
            genes <- data %>%
                filter(log2FoldChange > LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD) %>%
                pull(gene_symbol)
            cat("  Filtered to", length(genes), "upregulated genes (log2FC >", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, ")\n")
        } else if (direction == "down") {
            genes <- data %>%
                filter(log2FoldChange < -LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD) %>%
                pull(gene_symbol)
            cat("  Filtered to", length(genes), "downregulated genes (log2FC < -", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, ")\n")
        } else {
            stop("Invalid direction: ", direction, ". Must be 'up' or 'down'")
        }
    }

    genes <- genes[!is.na(genes) & genes != ""]
    cat("  Final gene count:", length(genes), "\n")

    return(genes)
}

# ============================================================
# Function: Load and filter mouse genes, convert to human
# ============================================================
load_mouse_genes_mgi <- function(file_path, source_name, mgi_orthologs, direction = "up") {
    cat("\nLoading", source_name, "data...\n")
    cat("  File:", file_path, "\n")
    cat("  Direction:", direction, "\n")

    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(data), "total genes\n")

    # Filter for up/downregulated genes
    if (direction == "up") {
        filtered_data <- data %>%
            filter(!is.na(log2FoldChange), !is.na(padj)) %>%
            filter(log2FoldChange > LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD)
        cat("  Filtered to", nrow(filtered_data), "upregulated genes (log2FC >", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, ")\n")
    } else if (direction == "down") {
        filtered_data <- data %>%
            filter(!is.na(log2FoldChange), !is.na(padj)) %>%
            filter(log2FoldChange < -LOG2FC_THRESHOLD, padj < PADJ_THRESHOLD)
        cat("  Filtered to", nrow(filtered_data), "downregulated genes (log2FC < -", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, ")\n")
    } else {
        stop("Invalid direction: ", direction, ". Must be 'up' or 'down'")
    }

    if (nrow(filtered_data) == 0) {
        warning("No mouse genes pass the filtering criteria!")
        return(list(human_genes = character(), ortholog_mapping = data.frame()))
    }

    mouse_genes <- filtered_data %>% pull(gene)
    cat("  Mouse genes:", length(mouse_genes), "\n")

    # Convert to human orthologs using MGI
    orthologs <- convert_mouse_to_human_mgi(mouse_genes, mgi_orthologs)

    if (nrow(orthologs) == 0) {
        warning("No orthologs found for ", source_name)
        return(list(human_genes = character(), ortholog_mapping = data.frame()))
    }

    # Get unique human gene symbols
    human_genes <- unique(orthologs$human_gene)
    human_genes <- human_genes[!is.na(human_genes) & human_genes != ""]

    cat("  Converted to", length(human_genes), "unique human orthologs\n")

    # Calculate statistics
    percent_mapped <- (nrow(orthologs) / length(mouse_genes)) * 100
    cat("  Mapping success rate:", sprintf("%.1f%%\n", percent_mapped))

    # Save ortholog mapping
    mapping_file <- file.path(OUTPUT_DIR, paste0(gsub(" ", "_", source_name), "_MGI_ortholog_mapping.csv"))
    write.csv(orthologs, mapping_file, row.names = FALSE)
    cat("  Saved ortholog mapping to:", mapping_file, "\n")

    return(list(human_genes = human_genes, ortholog_mapping = orthologs))
}

# ============================================================
# Function: Create Venn diagram with ggVennDiagram
# ============================================================
create_venn_diagram <- function(gene_list1, gene_list2,
                                label1, label2,
                                output_file, title) {

    cat("\nCreating Venn diagram:", title, "\n")
    cat("  ", label1, "genes:", length(gene_list1), "\n")
    cat("  ", label2, "genes:", length(gene_list2), "\n")

    overlap <- intersect(gene_list1, gene_list2)
    unique1 <- setdiff(gene_list1, gene_list2)
    unique2 <- setdiff(gene_list2, gene_list1)

    cat("  Overlap:", length(overlap), "genes\n")
    cat("  Unique to", label1, ":", length(unique1), "genes\n")
    cat("  Unique to", label2, ":", length(unique2), "genes\n")

    # Calculate overlap percentage
    if (length(gene_list1) > 0 && length(gene_list2) > 0) {
        overlap_pct1 <- (length(overlap) / length(gene_list1)) * 100
        overlap_pct2 <- (length(overlap) / length(gene_list2)) * 100
        cat("  Overlap as % of", label1, ":", sprintf("%.1f%%\n", overlap_pct1))
        cat("  Overlap as % of", label2, ":", sprintf("%.1f%%\n", overlap_pct2))
    }

    # Save gene lists
    overlap_file <- sub("\\.png$", "_overlap.txt", output_file)
    writeLines(sort(overlap), overlap_file)
    cat("  Saved overlapping genes to:", overlap_file, "\n")

    unique1_file <- sub("\\.png$", paste0("_unique_", gsub(" ", "_", label1), ".txt"), output_file)
    writeLines(sort(unique1), unique1_file)

    unique2_file <- sub("\\.png$", paste0("_unique_", gsub(" ", "_", label2), ".txt"), output_file)
    writeLines(sort(unique2), unique2_file)

    # Prepare data for ggVennDiagram
    venn_list <- list(gene_list1, gene_list2)
    names(venn_list) <- c(label1, label2)

    # Create Venn diagram with ggVennDiagram (without set labels)
    venn_plot <- ggVennDiagram(
        venn_list,
        label = "count",
        label_alpha = 0,
        label_size = 5,
        edge_size = 1.5,
        set_size = 0  # Hide default set labels
    ) +
        scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
        scale_color_manual(values = c("#E74C3C", "#3498DB")) +
        scale_x_continuous(expand = expansion(mult = 0.2)) +
        scale_y_continuous(expand = expansion(mult = 0.35)) +
        # Add set labels above and below - outside the circles
        annotate("text", x = 0, y = 10, label = label1,
                 size = 5, fontface = "bold", color = "#E74C3C") +
        annotate("text", x = 0, y = -8, label = label2,
                 size = 5, fontface = "bold", color = "#3498DB") +
        labs(title = title) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 2)),
            plot.margin = margin(t = 2, r = 20, b = 2, l = 20),
            legend.position = "none"
        ) +
        coord_fixed(clip = "off")

    # Create statistics table
    stats_df <- data.frame(
        Metric = c(
            paste(label1, "total"),
            paste(label2, "total"),
            "Overlap",
            paste("Unique to", label1),
            paste("Unique to", label2),
            paste("Overlap % of", label1),
            paste("Overlap % of", label2)
        ),
        Value = c(
            length(gene_list1),
            length(gene_list2),
            length(overlap),
            length(unique1),
            length(unique2),
            sprintf("%.1f%%", overlap_pct1),
            sprintf("%.1f%%", overlap_pct2)
        )
    )

    # Create table grob
    stats_table <- tableGrob(
        stats_df,
        rows = NULL,
        theme = ttheme_minimal(
            core = list(
                fg_params = list(hjust = 0, x = 0.1, fontsize = 10),
                bg_params = list(fill = c("#F0F0F0", "#FFFFFF"))
            ),
            colhead = list(
                fg_params = list(fontface = "bold", fontsize = 11),
                bg_params = list(fill = "#4981BF", col = "white")
            )
        )
    )

    # Combine plot and table
    combined_plot <- grid.arrange(
        venn_plot,
        stats_table,
        ncol = 1,
        heights = c(0.65, 0.35)
    )

    # Save PNG
    ggsave(
        output_file,
        plot = combined_plot,
        width = 10,
        height = 10,
        dpi = 300,
        bg = "white"
    )
    cat("  Saved PNG:", output_file, "\n")

    # Save PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(
        pdf_file,
        plot = combined_plot,
        width = 10,
        height = 10,
        bg = "white"
    )
    cat("  Saved PDF:", pdf_file, "\n")

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
cat("STEP 1: Loading MGI Ortholog Database\n")
cat("============================================================\n")

mgi_orthologs <- load_mgi_orthologs()

cat("\n============================================================\n")
cat("STEP 2: Loading Human Gene Data - UPREGULATED\n")
cat("============================================================\n")

human_neuron_genes_up <- load_human_genes(
    HUMAN_NEURON_FILE,
    "Human Neurons",
    prefiltered = TRUE,
    direction = "up"
)

human_npc_genes_up <- load_human_genes(
    HUMAN_NPC_FILE,
    "Human NPCs",
    prefiltered = FALSE,
    direction = "up"
)

cat("\n============================================================\n")
cat("STEP 2b: Loading Human Gene Data - DOWNREGULATED\n")
cat("============================================================\n")

# Note: The pre-filtered neuron file only contains upregulated genes
# We need to load the raw data for downregulated genes
# Check if a downregulated file exists, otherwise we'll need the raw data
human_neuron_genes_down <- tryCatch({
    load_human_genes(
        "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results/05_deseq2/Neuron_differential_expression_downregulated.csv",
        "Human Neurons",
        prefiltered = TRUE,
        direction = "down"
    )
}, error = function(e) {
    cat("  No pre-filtered downregulated file found for neurons. Trying raw data...\n")
    tryCatch({
        load_human_genes(
            "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results/05_deseq2/Neuron_differential_expression.csv",
            "Human Neurons",
            prefiltered = FALSE,
            direction = "down"
        )
    }, error = function(e2) {
        cat("  WARNING: Could not load downregulated neuron genes:", conditionMessage(e2), "\n")
        return(character())
    })
})

human_npc_genes_down <- load_human_genes(
    HUMAN_NPC_FILE,
    "Human NPCs",
    prefiltered = FALSE,
    direction = "down"
)

cat("\n============================================================\n")
cat("STEP 3: Loading Mouse Data and Converting with MGI Orthologs - UPREGULATED\n")
cat("============================================================\n")

mouse_neuron_result_up <- load_mouse_genes_mgi(
    MOUSE_NEURON_FILE,
    "Mouse Neurons",
    mgi_orthologs,
    direction = "up"
)

mouse_npc_result_up <- load_mouse_genes_mgi(
    MOUSE_NPC_FILE,
    "Mouse NPCs",
    mgi_orthologs,
    direction = "up"
)

cat("\n============================================================\n")
cat("STEP 3b: Loading Mouse Data and Converting with MGI Orthologs - DOWNREGULATED\n")
cat("============================================================\n")

mouse_neuron_result_down <- load_mouse_genes_mgi(
    MOUSE_NEURON_FILE,
    "Mouse Neurons",
    mgi_orthologs,
    direction = "down"
)

mouse_npc_result_down <- load_mouse_genes_mgi(
    MOUSE_NPC_FILE,
    "Mouse NPCs",
    mgi_orthologs,
    direction = "down"
)

cat("\n============================================================\n")
cat("STEP 4: Creating Venn Diagrams - UPREGULATED GENES\n")
cat("============================================================\n")

# Venn diagram 1: Neurons - Upregulated
neuron_stats_up <- create_venn_diagram(
    gene_list1 = human_neuron_genes_up,
    gene_list2 = mouse_neuron_result_up$human_genes,
    label1 = "Human Neurons",
    label2 = "Mouse Neurons",
    output_file = file.path(OUTPUT_DIR, "venn_neurons_human_vs_mouse_MGI_upregulated.png"),
    title = "Upregulated Genes: Human vs Mouse Neurons\n(MGI Orthologs, log2FC > 0.5, padj < 0.05)"
)

# Venn diagram 2: NPCs - Upregulated
npc_stats_up <- create_venn_diagram(
    gene_list1 = human_npc_genes_up,
    gene_list2 = mouse_npc_result_up$human_genes,
    label1 = "Human NPCs",
    label2 = "Mouse NPCs",
    output_file = file.path(OUTPUT_DIR, "venn_npcs_human_vs_mouse_MGI_upregulated.png"),
    title = "Upregulated Genes: Human vs Mouse NPCs\n(MGI Orthologs, log2FC > 0.5, padj < 0.05)"
)

cat("\n============================================================\n")
cat("STEP 4b: Creating Venn Diagrams - DOWNREGULATED GENES\n")
cat("============================================================\n")

# Venn diagram 3: Neurons - Downregulated
neuron_stats_down <- create_venn_diagram(
    gene_list1 = human_neuron_genes_down,
    gene_list2 = mouse_neuron_result_down$human_genes,
    label1 = "Human Neurons",
    label2 = "Mouse Neurons",
    output_file = file.path(OUTPUT_DIR, "venn_neurons_human_vs_mouse_MGI_downregulated.png"),
    title = "Downregulated Genes: Human vs Mouse Neurons\n(MGI Orthologs, log2FC < -0.5, padj < 0.05)"
)

# Venn diagram 4: NPCs - Downregulated
npc_stats_down <- create_venn_diagram(
    gene_list1 = human_npc_genes_down,
    gene_list2 = mouse_npc_result_down$human_genes,
    label1 = "Human NPCs",
    label2 = "Mouse NPCs",
    output_file = file.path(OUTPUT_DIR, "venn_npcs_human_vs_mouse_MGI_downregulated.png"),
    title = "Downregulated Genes: Human vs Mouse NPCs\n(MGI Orthologs, log2FC < -0.5, padj < 0.05)"
)

# ============================================================
# Summary Report
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY REPORT (Using MGI Orthologs)\n")
cat("============================================================\n\n")

cat("UPREGULATED GENES\n")
cat("=================\n\n")

cat("NEURONS:\n")
cat("  Human:", neuron_stats_up$total1, "upregulated genes\n")
cat("  Mouse:", neuron_stats_up$total2, "upregulated genes (MGI orthologs)\n")
cat("  Overlap:", neuron_stats_up$overlap, "genes (",
    sprintf("%.1f%%", neuron_stats_up$overlap / min(neuron_stats_up$total1, neuron_stats_up$total2) * 100),
    "of smaller set)\n")
cat("  Unique to Human:", neuron_stats_up$unique1, "genes\n")
cat("  Unique to Mouse:", neuron_stats_up$unique2, "genes\n\n")

cat("NPCs:\n")
cat("  Human:", npc_stats_up$total1, "upregulated genes\n")
cat("  Mouse:", npc_stats_up$total2, "upregulated genes (MGI orthologs)\n")
cat("  Overlap:", npc_stats_up$overlap, "genes (",
    sprintf("%.1f%%", npc_stats_up$overlap / min(npc_stats_up$total1, npc_stats_up$total2) * 100),
    "of smaller set)\n")
cat("  Unique to Human:", npc_stats_up$unique1, "genes\n")
cat("  Unique to Mouse:", npc_stats_up$unique2, "genes\n\n")

cat("DOWNREGULATED GENES\n")
cat("===================\n\n")

cat("NEURONS:\n")
cat("  Human:", neuron_stats_down$total1, "downregulated genes\n")
cat("  Mouse:", neuron_stats_down$total2, "downregulated genes (MGI orthologs)\n")
cat("  Overlap:", neuron_stats_down$overlap, "genes (",
    sprintf("%.1f%%", neuron_stats_down$overlap / min(neuron_stats_down$total1, neuron_stats_down$total2) * 100),
    "of smaller set)\n")
cat("  Unique to Human:", neuron_stats_down$unique1, "genes\n")
cat("  Unique to Mouse:", neuron_stats_down$unique2, "genes\n\n")

cat("NPCs:\n")
cat("  Human:", npc_stats_down$total1, "downregulated genes\n")
cat("  Mouse:", npc_stats_down$total2, "downregulated genes (MGI orthologs)\n")
cat("  Overlap:", npc_stats_down$overlap, "genes (",
    sprintf("%.1f%%", npc_stats_down$overlap / min(npc_stats_down$total1, npc_stats_down$total2) * 100),
    "of smaller set)\n")
cat("  Unique to Human:", npc_stats_down$unique1, "genes\n")
cat("  Unique to Mouse:", npc_stats_down$unique2, "genes\n\n")

# Save summary
summary_file <- file.path(OUTPUT_DIR, "venn_summary_MGI_orthologs.txt")
sink(summary_file)
cat("VENN DIAGRAM ANALYSIS - MGI ORTHOLOGS\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("==========================================\n\n")
cat("FILTERING CRITERIA:\n")
cat("  Upregulated: log2FoldChange >", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, "\n")
cat("  Downregulated: log2FoldChange < -", LOG2FC_THRESHOLD, ", padj <", PADJ_THRESHOLD, "\n")
cat("  Ortholog database: MGI (Mouse Genome Informatics)\n\n")

cat("===============================\n")
cat("UPREGULATED GENES\n")
cat("===============================\n\n")

cat("NEURONS:\n")
cat("  Human:", neuron_stats_up$total1, "\n")
cat("  Mouse:", neuron_stats_up$total2, "(MGI orthologs)\n")
cat("  Overlap:", neuron_stats_up$overlap, sprintf("(%.1f%%)\n", neuron_stats_up$overlap / min(neuron_stats_up$total1, neuron_stats_up$total2) * 100))
if (length(neuron_stats_up$overlap_genes) > 0) {
    cat("  Overlapping genes:\n    ", paste(neuron_stats_up$overlap_genes, collapse = ", "), "\n\n")
}

cat("NPCs:\n")
cat("  Human:", npc_stats_up$total1, "\n")
cat("  Mouse:", npc_stats_up$total2, "(MGI orthologs)\n")
cat("  Overlap:", npc_stats_up$overlap, sprintf("(%.1f%%)\n", npc_stats_up$overlap / min(npc_stats_up$total1, npc_stats_up$total2) * 100))
if (length(npc_stats_up$overlap_genes) <= 100) {
    cat("  Overlapping genes:\n    ", paste(npc_stats_up$overlap_genes, collapse = ", "), "\n\n")
} else {
    cat("  Overlapping genes (first 100):\n    ", paste(head(npc_stats_up$overlap_genes, 100), collapse = ", "), "\n\n")
}

cat("\n===============================\n")
cat("DOWNREGULATED GENES\n")
cat("===============================\n\n")

cat("NEURONS:\n")
cat("  Human:", neuron_stats_down$total1, "\n")
cat("  Mouse:", neuron_stats_down$total2, "(MGI orthologs)\n")
cat("  Overlap:", neuron_stats_down$overlap, sprintf("(%.1f%%)\n", neuron_stats_down$overlap / min(neuron_stats_down$total1, neuron_stats_down$total2) * 100))
if (length(neuron_stats_down$overlap_genes) > 0) {
    cat("  Overlapping genes:\n    ", paste(neuron_stats_down$overlap_genes, collapse = ", "), "\n\n")
}

cat("NPCs:\n")
cat("  Human:", npc_stats_down$total1, "\n")
cat("  Mouse:", npc_stats_down$total2, "(MGI orthologs)\n")
cat("  Overlap:", npc_stats_down$overlap, sprintf("(%.1f%%)\n", npc_stats_down$overlap / min(npc_stats_down$total1, npc_stats_down$total2) * 100))
if (length(npc_stats_down$overlap_genes) <= 100) {
    cat("  Overlapping genes:\n    ", paste(npc_stats_down$overlap_genes, collapse = ", "), "\n\n")
} else {
    cat("  Overlapping genes (first 100):\n    ", paste(head(npc_stats_down$overlap_genes, 100), collapse = ", "), "\n\n")
}
sink()

cat("Saved summary report to:", summary_file, "\n")

cat("\n==========================================\n")
cat("Analysis complete! (Using MGI Orthologs)\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Generated files:\n")
cat("  UPREGULATED GENES:\n")
cat("    - venn_neurons_human_vs_mouse_MGI_upregulated.png/pdf\n")
cat("    - venn_npcs_human_vs_mouse_MGI_upregulated.png/pdf\n")
cat("  DOWNREGULATED GENES:\n")
cat("    - venn_neurons_human_vs_mouse_MGI_downregulated.png/pdf\n")
cat("    - venn_npcs_human_vs_mouse_MGI_downregulated.png/pdf\n")
cat("  OTHER FILES:\n")
cat("    - MGI_mouse_human_orthologs.csv (full MGI database)\n")
cat("    - *_MGI_ortholog_mapping.csv (genes used in this analysis)\n")
cat("    - *_overlap.txt, *_unique_*.txt (gene lists)\n")
cat("    - venn_summary_MGI_orthologs.txt\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
