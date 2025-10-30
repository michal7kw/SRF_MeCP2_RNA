#!/usr/bin/env Rscript

# Script to validate the quality of mouse-to-human ortholog mappings
# generated using the capitalization method

suppressPackageStartupMessages({
    library(tidyverse)
})

cat("==========================================\n")
cat("Validating Ortholog Mapping Quality\n")
cat("==========================================\n\n")

INPUT_DIR <- "results/venn_diagrams"

# ============================================================
# Load ortholog mappings
# ============================================================

cat("Loading ortholog mapping files...\n")
neuron_orthologs <- read.csv(file.path(INPUT_DIR, "Mouse_Neurons_ortholog_mapping.csv"))
npc_orthologs <- read.csv(file.path(INPUT_DIR, "Mouse_NPCs_ortholog_mapping.csv"))

cat("  Neurons:", nrow(neuron_orthologs), "mappings\n")
cat("  NPCs:", nrow(npc_orthologs), "mappings\n\n")

# ============================================================
# Function to assess mapping quality
# ============================================================
assess_mapping_quality <- function(orthologs, dataset_name) {
    cat("==========================================\n")
    cat("Assessing:", dataset_name, "\n")
    cat("==========================================\n")

    total <- nrow(orthologs)

    # Mouse-specific gene patterns
    gm_genes <- orthologs %>% filter(grepl("^Gm[0-9]", mouse_gene))
    rik_genes <- orthologs %>% filter(grepl("Rik$", mouse_gene))
    predicted_genes <- orthologs %>% filter(grepl("^[0-9].*Rik$", mouse_gene))
    pseudogenes <- orthologs %>% filter(grepl("-ps[0-9]?$", mouse_gene))

    mouse_specific <- nrow(gm_genes) + nrow(rik_genes) + nrow(predicted_genes) + nrow(pseudogenes)

    # Likely real orthologs
    likely_real <- total - mouse_specific

    cat("\nTotal mappings:", total, "\n\n")

    cat("Mouse-specific genes (likely NO human ortholog):\n")
    cat("  Gm* genes:", nrow(gm_genes), sprintf("(%.1f%%)\n", nrow(gm_genes)/total*100))
    cat("  *Rik genes:", nrow(rik_genes), sprintf("(%.1f%%)\n", nrow(rik_genes)/total*100))
    cat("  Predicted genes:", nrow(predicted_genes), sprintf("(%.1f%%)\n", nrow(predicted_genes)/total*100))
    cat("  Pseudogenes:", nrow(pseudogenes), sprintf("(%.1f%%)\n", nrow(pseudogenes)/total*100))
    cat("  Total mouse-specific:", mouse_specific, sprintf("(%.1f%%)\n\n", mouse_specific/total*100))

    cat("Likely real orthologs:", likely_real, sprintf("(%.1f%%)\n\n", likely_real/total*100))

    # Known gene validation
    cat("Known gene validation (checking if these exist):\n")
    known_genes <- c("Mecp2", "Gapdh", "Actb", "Tubb3", "Pax6", "Sox2", "Nestin")

    for (gene in known_genes) {
        exists <- gene %in% orthologs$mouse_gene
        if (exists) {
            human <- orthologs$human_gene[orthologs$mouse_gene == gene]
            cat("  ", gene, "->", human, "(FOUND)\n")
        } else {
            cat("  ", gene, "(NOT FOUND in this dataset)\n")
        }
    }

    cat("\n")

    # Return statistics
    return(list(
        total = total,
        mouse_specific = mouse_specific,
        likely_real = likely_real,
        percent_real = likely_real/total*100
    ))
}

# ============================================================
# Assess both datasets
# ============================================================

neuron_stats <- assess_mapping_quality(neuron_orthologs, "Neuron Genes")
npc_stats <- assess_mapping_quality(npc_orthologs, "NPC Genes")

# ============================================================
# Check overlap genes for validation
# ============================================================

cat("==========================================\n")
cat("Overlap Gene Validation\n")
cat("==========================================\n\n")

# Load overlap genes
neuron_overlap <- readLines(file.path(INPUT_DIR, "venn_neurons_human_vs_mouse_overlap.txt"))
npc_overlap <- readLines(file.path(INPUT_DIR, "venn_npcs_human_vs_mouse_overlap.txt"))

cat("Neuron overlapping genes (", length(neuron_overlap), "):\n", sep="")
cat("  ", paste(neuron_overlap, collapse=", "), "\n\n")

cat("NPC overlapping genes (first 50 of ", length(npc_overlap), "):\n", sep="")
cat("  ", paste(head(npc_overlap, 50), collapse=", "), "\n\n")

# Check how many overlap genes are from mouse-specific genes
neuron_overlap_mouse_specific <- sum(grepl("^GM[0-9]|RIK$", neuron_overlap))
npc_overlap_mouse_specific <- sum(grepl("^GM[0-9]|RIK$", npc_overlap))

cat("Mouse-specific genes in overlap:\n")
cat("  Neurons:", neuron_overlap_mouse_specific, "out of", length(neuron_overlap),
    sprintf("(%.1f%%)\n", neuron_overlap_mouse_specific/length(neuron_overlap)*100))
cat("  NPCs:", npc_overlap_mouse_specific, "out of", length(npc_overlap),
    sprintf("(%.1f%%)\n", npc_overlap_mouse_specific/length(npc_overlap)*100))

cat("\n")

# ============================================================
# Summary Report
# ============================================================

cat("==========================================\n")
cat("SUMMARY\n")
cat("==========================================\n\n")

cat("Overall Mapping Quality:\n\n")

cat("NEURONS:\n")
cat("  Total mappings:", neuron_stats$total, "\n")
cat("  Likely real orthologs:", neuron_stats$likely_real,
    sprintf("(%.1f%%)\n", neuron_stats$percent_real))
cat("  Mouse-specific genes:", neuron_stats$mouse_specific,
    sprintf("(%.1f%%)\n", 100-neuron_stats$percent_real))
cat("\n")

cat("NPCs:\n")
cat("  Total mappings:", npc_stats$total, "\n")
cat("  Likely real orthologs:", npc_stats$likely_real,
    sprintf("(%.1f%%)\n", npc_stats$percent_real))
cat("  Mouse-specific genes:", npc_stats$mouse_specific,
    sprintf("(%.1f%%)\n", 100-npc_stats$percent_real))
cat("\n")

cat("CONCLUSIONS:\n")
cat("  - The capitalization method worked for ~", round(mean(c(neuron_stats$percent_real, npc_stats$percent_real))), "%% of genes\n", sep="")
cat("  - Mouse-specific genes (Gm*, *Rik) do NOT have true human orthologs\n")
cat("  - These should be excluded for accurate comparison\n")
cat("  - The overlap genes appear to be predominantly real orthologs\n")
cat("  - For publication-quality results, use proper ortholog database (e.g., Ensembl BioMart)\n\n")

cat("RECOMMENDATION:\n")
cat("  Consider filtering out mouse-specific genes before creating Venn diagrams:\n")
cat("  - Remove genes matching: ^Gm[0-9], *Rik$, ^[0-9].*Rik$, -ps[0-9]?\n")
cat("  - This will give more accurate overlap statistics\n\n")

cat("==========================================\n")
