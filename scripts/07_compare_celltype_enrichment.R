#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 7: Cross-Cell-Type GO/KEGG Enrichment Comparison
# =======================================================================
# This script creates comparative visualization of GO/KEGG enrichment
# between NPCs and Neurons, showing upregulated and downregulated genes.
#
# PREREQUISITES:
#   1. NPC GO/KEGG enrichment files must exist:
#      - results/GO_KEGG/GO_enrichment_upregulated.csv
#      - results/GO_KEGG/GO_enrichment_downregulated.csv
#
#   2. Neuron GO enrichment files must exist (run 06c first if needed):
#      - [neuron_dir]/go_analysis_neu/GO_enrichment_upregulated.csv
#      - [neuron_dir]/go_analysis_neu/GO_enrichment_downregulated.csv
#
# OUTPUT:
#   - GO_comparative_NPCs_Neurons.png/pdf
#   - KEGG_comparative_NPCs_Neurons.png/pdf (if available)
#
# VISUALIZATION:
#   Creates dot plots with 4 columns: NPCs UP | NPCs Down | Ns UP | Ns Down
#   Similar to the format in DATA/90-1244168066/go_plots.png

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

cat("=========================================================\n")
cat("Cross-Cell-Type Enrichment Comparison\n")
cat("=========================================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ============================================================
# Configuration
# ============================================================

# NPC paths (current analysis)
NPC_GO_DIR <- "results/GO_KEGG"

# Neuron paths (from previous analysis)
NEURON_BASE_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results"
NEURON_GO_DIR <- file.path(NEURON_BASE_DIR, "go_analysis_neu")

# Output directory
OUTPUT_DIR <- "results/GO_KEGG/comparative_plots"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Plotting parameters
TOP_N_TERMS <- 20         # Number of top terms to show
MIN_COUNT <- 5            # Minimum gene count to include
MAX_PVAL <- 0.05          # Maximum adjusted p-value to include

# ============================================================
# Validation
# ============================================================

cat("Validating prerequisites...\n\n")

validate_file <- function(file_path, description) {
    if (file.exists(file_path)) {
        cat("  ✓", description, "\n")
        cat("    →", file_path, "\n")
        return(TRUE)
    } else {
        cat("  ✗", description, "NOT FOUND\n")
        cat("    →", file_path, "\n")
        return(FALSE)
    }
}

cat("NPC GO Enrichment Files:\n")
npc_go_up_exists <- validate_file(
    file.path(NPC_GO_DIR, "GO_enrichment_upregulated.csv"),
    "NPC upregulated"
)
npc_go_down_exists <- validate_file(
    file.path(NPC_GO_DIR, "GO_enrichment_downregulated.csv"),
    "NPC downregulated"
)

cat("\nNeuron GO Enrichment Files:\n")
neuron_go_up_exists <- validate_file(
    file.path(NEURON_GO_DIR, "GO_enrichment_upregulated.csv"),
    "Neuron upregulated"
)
neuron_go_down_exists <- validate_file(
    file.path(NEURON_GO_DIR, "GO_enrichment_downregulated.csv"),
    "Neuron downregulated"
)

cat("\n")

# Check if we have minimum required files
if (!(npc_go_up_exists || npc_go_down_exists)) {
    stop("ERROR: No NPC GO enrichment files found!\n",
         "       Please run: scripts/06_functional_enrichment.R")
}

if (!(neuron_go_up_exists || neuron_go_down_exists)) {
    stop("ERROR: No Neuron GO enrichment files found!\n",
         "       Please run: scripts/06c_prepare_neuron_go_enrichment.R")
}

cat("✓ All required files present\n\n")

# ============================================================
# Function: Load and Prepare Combined Data
# ============================================================

prepare_combined_celltype_data <- function(npc_up_file, npc_down_file,
                                           neuron_up_file, neuron_down_file,
                                           type = "GO") {

    cat("Loading", type, "enrichment data for NPCs and Neurons...\n")

    # Initialize empty data frame
    combined_data <- data.frame()

    # Helper function to load and process a single file
    load_condition <- function(file_path, cell_type, direction) {
        if (!file.exists(file_path)) {
            cat("  ⊘", cell_type, direction, "- file not found, skipping\n")
            return(NULL)
        }

        data <- read.csv(file_path) %>%
            filter(Count >= MIN_COUNT, p.adjust <= MAX_PVAL) %>%
            arrange(p.adjust) %>%
            head(TOP_N_TERMS) %>%
            mutate(
                CellType = cell_type,
                Direction = direction,
                Term = Description
            )

        # Select relevant columns
        if ("ONTOLOGY" %in% colnames(data)) {
            data <- data %>% select(Term, Count, p.adjust, CellType, Direction, ONTOLOGY)
        } else {
            data <- data %>% select(Term, Count, p.adjust, CellType, Direction)
        }

        cat("  ✓", cell_type, direction, "-", nrow(data), "terms\n")
        return(data)
    }

    # Load all four conditions
    npc_up <- load_condition(npc_up_file, "NPCs", "UP")
    npc_down <- load_condition(npc_down_file, "NPCs", "Down")
    neuron_up <- load_condition(neuron_up_file, "Ns", "UP")
    neuron_down <- load_condition(neuron_down_file, "Ns", "Down")

    # Combine all data
    combined_data <- bind_rows(npc_up, npc_down, neuron_up, neuron_down)

    if (nrow(combined_data) == 0) {
        cat("  ⚠ No data loaded!\n")
        return(NULL)
    }

    # Get unique terms
    unique_terms <- unique(combined_data$Term)
    cat("  → Total unique terms:", length(unique_terms), "\n\n")

    # Create complete data frame with all combinations
    all_combinations <- expand.grid(
        Term = unique_terms,
        CellType = c("NPCs", "Ns"),
        Direction = c("UP", "Down"),
        stringsAsFactors = FALSE
    )

    # Merge with actual data
    plot_data <- all_combinations %>%
        left_join(combined_data, by = c("Term", "CellType", "Direction")) %>%
        mutate(Condition = paste(CellType, Direction, sep = " "))

    return(plot_data)
}

# ============================================================
# Function: Create Comparative Dot Plot
# ============================================================

create_combined_dotplot <- function(data, title, output_file) {

    cat("Creating comparative plot:", title, "\n")

    # Filter out rows with no data (NA p-values)
    data_present <- data %>% filter(!is.na(p.adjust))

    if (nrow(data_present) == 0) {
        cat("  ⚠ No data to plot\n\n")
        return(NULL)
    }

    cat("  Plotting", length(unique(data_present$Term)), "unique terms\n")

    # Order terms by their best (lowest) p-value across all conditions
    term_order <- data %>%
        group_by(Term) %>%
        summarize(best_pval = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
        arrange(best_pval) %>%
        pull(Term)

    data$Term <- factor(data$Term, levels = rev(term_order))

    # Set factor levels for Condition (column order in plot)
    data$Condition <- factor(data$Condition,
                            levels = c("NPCs UP", "NPCs Down", "Ns UP", "Ns Down"))
    data_present$Condition <- factor(data_present$Condition,
                                    levels = c("NPCs UP", "NPCs Down", "Ns UP", "Ns Down"))

    # Simplify term names by removing GO: prefix
    data$Term_display <- gsub("^GO:\\d+\\s*", "", data$Term)
    data_present$Term_display <- gsub("^GO:\\d+\\s*", "", data_present$Term)

    # Update term factor levels
    term_order_display <- gsub("^GO:\\d+\\s*", "", term_order)
    data$Term_display <- factor(data$Term_display, levels = rev(term_order_display))
    data_present$Term_display <- factor(data_present$Term_display,
                                       levels = rev(term_order_display))

    # Create the plot
    p <- ggplot(data_present, aes(x = Condition, y = Term_display)) +
        geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
        scale_color_gradient(
            low = "#4575b4",      # Blue for more significant
            high = "#d73027",     # Red for less significant
            name = "Adjusted p-value",
            trans = "log10",
            breaks = c(1e-7, 1e-6, 1e-5, 1e-4),
            labels = c("1e-7", "1e-6", "1e-5", "1e-4")
        ) +
        scale_size_continuous(
            name = "Count",
            breaks = c(10, 25, 50, 100),
            range = c(3, 12)
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 11, face = "bold", color = "black"),
            axis.title = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "bottom",
            legend.box = "horizontal",
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        ) +
        labs(title = title)

    # Calculate plot dimensions
    n_terms <- length(unique(data_present$Term_display))
    plot_height <- max(8, n_terms * 0.35)
    plot_width <- 10

    # Save plot
    ggsave(output_file, p, width = plot_width, height = plot_height, dpi = 300)
    cat("  ✓ Saved PNG:", output_file, "\n")

    # Also save PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(pdf_file, p, width = plot_width, height = plot_height)
    cat("  ✓ Saved PDF:", pdf_file, "\n\n")

    return(p)
}

# ============================================================
# Load Exclusion Terms (if available)
# ============================================================

EXCLUSION_FILE <- "to_hide.txt"

excluded_terms <- character()
if (file.exists(EXCLUSION_FILE)) {
    cat("Loading term exclusion list...\n")
    excluded_terms <- readLines(EXCLUSION_FILE, warn = FALSE)
    excluded_terms <- excluded_terms[nchar(excluded_terms) > 0]  # Remove empty lines
    cat("  ✓ Loaded", length(excluded_terms), "terms to exclude\n")
    cat("  Terms to hide:\n")
    for (term in excluded_terms) {
        cat("    -", term, "\n")
    }
    cat("\n")
} else {
    cat("ℹ No exclusion file found (", EXCLUSION_FILE, ") - showing all terms\n\n")
}

# ============================================================
# Main Analysis: GO Enrichment Comparison
# ============================================================

cat("=========================================================\n")
cat("GO Enrichment Comparison\n")
cat("=========================================================\n\n")

npc_go_up_file <- file.path(NPC_GO_DIR, "GO_enrichment_upregulated.csv")
npc_go_down_file <- file.path(NPC_GO_DIR, "GO_enrichment_downregulated.csv")
neuron_go_up_file <- file.path(NEURON_GO_DIR, "GO_enrichment_upregulated.csv")
neuron_go_down_file <- file.path(NEURON_GO_DIR, "GO_enrichment_downregulated.csv")

go_data <- prepare_combined_celltype_data(
    npc_up_file = npc_go_up_file,
    npc_down_file = npc_go_down_file,
    neuron_up_file = neuron_go_up_file,
    neuron_down_file = neuron_go_down_file,
    type = "GO"
)

if (!is.null(go_data)) {
    # Create full plot (all terms)
    cat("Creating full GO plot (all terms)...\n")
    create_combined_dotplot(
        data = go_data,
        title = "GO Enrichment: NPCs and Neurons (All Terms)",
        output_file = file.path(OUTPUT_DIR, "GO_comparative_NPCs_Neurons_full.png")
    )

    # Create filtered plot (excluding specified terms)
    if (length(excluded_terms) > 0) {
        cat("Creating filtered GO plot (excluding specified terms)...\n")

        # Filter out excluded terms
        go_data_filtered <- go_data %>%
            filter(!Term %in% excluded_terms)

        # Count how many terms were removed
        n_removed <- length(unique(go_data$Term)) - length(unique(go_data_filtered$Term))
        cat("  Excluded", n_removed, "terms\n")

        if (nrow(go_data_filtered %>% filter(!is.na(p.adjust))) > 0) {
            create_combined_dotplot(
                data = go_data_filtered,
                title = "GO Enrichment: NPCs and Neurons (Filtered)",
                output_file = file.path(OUTPUT_DIR, "GO_comparative_NPCs_Neurons_filtered.png")
            )
        } else {
            cat("  ⚠ No terms remaining after filtering\n\n")
        }
    }
} else {
    cat("⚠ Skipping GO plots - no data available\n\n")
}

# ============================================================
# Main Analysis: KEGG Enrichment Comparison (Optional)
# ============================================================

cat("=========================================================\n")
cat("KEGG Pathway Enrichment Comparison\n")
cat("=========================================================\n\n")

npc_kegg_up_file <- file.path(NPC_GO_DIR, "KEGG_enrichment_upregulated.csv")
npc_kegg_down_file <- file.path(NPC_GO_DIR, "KEGG_enrichment_downregulated.csv")
neuron_kegg_up_file <- file.path(NEURON_GO_DIR, "KEGG_enrichment_upregulated.csv")
neuron_kegg_down_file <- file.path(NEURON_GO_DIR, "KEGG_enrichment_downregulated.csv")

# Check if KEGG files exist
kegg_available <- any(file.exists(c(npc_kegg_up_file, npc_kegg_down_file,
                                    neuron_kegg_up_file, neuron_kegg_down_file)))

if (kegg_available) {
    kegg_data <- prepare_combined_celltype_data(
        npc_up_file = npc_kegg_up_file,
        npc_down_file = npc_kegg_down_file,
        neuron_up_file = neuron_kegg_up_file,
        neuron_down_file = neuron_kegg_down_file,
        type = "KEGG"
    )

    if (!is.null(kegg_data)) {
        create_combined_dotplot(
            data = kegg_data,
            title = "KEGG Pathway Enrichment: NPCs and Neurons",
            output_file = file.path(OUTPUT_DIR, "KEGG_comparative_NPCs_Neurons.png")
        )
    } else {
        cat("⚠ Skipping KEGG plot - no data available\n\n")
    }
} else {
    cat("⊘ KEGG enrichment files not found - skipping KEGG comparison\n\n")
}

# ============================================================
# Final Summary
# ============================================================

cat("=========================================================\n")
cat("CROSS-CELL-TYPE COMPARISON COMPLETE\n")
cat("=========================================================\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("Generated files:\n")
output_files <- list.files(OUTPUT_DIR, full.names = FALSE, pattern = "\\.(png|pdf)$")
if (length(output_files) > 0) {
    # Separate by type
    go_full <- grep("GO.*full", output_files, value = TRUE)
    go_filtered <- grep("GO.*filtered", output_files, value = TRUE)
    kegg <- grep("KEGG", output_files, value = TRUE)

    if (length(go_full) > 0) {
        cat("\n  GO Plots (All Terms):\n")
        for (f in go_full) {
            cat("    ✓", f, "\n")
        }
    }

    if (length(go_filtered) > 0) {
        cat("\n  GO Plots (Filtered - excluding", length(excluded_terms), "terms):\n")
        for (f in go_filtered) {
            cat("    ✓", f, "\n")
        }
    }

    if (length(kegg) > 0) {
        cat("\n  KEGG Plots:\n")
        for (f in kegg) {
            cat("    ✓", f, "\n")
        }
    }
} else {
    cat("  (No files generated)\n")
}

if (length(excluded_terms) > 0) {
    cat("\n📝 Note: Filtered plots exclude these", length(excluded_terms), "terms:\n")
    for (term in excluded_terms) {
        cat("    -", term, "\n")
    }
}

cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("=========================================================\n")
