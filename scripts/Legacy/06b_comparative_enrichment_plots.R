#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 6b: Comparative Enrichment Plots - NPCs vs Neurons
# This script creates comparative dot plots showing GO/KEGG enrichment
# for both NPCs and Neurons (Ns) with upregulated and downregulated genes

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

cat("==========================================\n")
cat("Creating Comparative Enrichment Plots\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths
# NPCs (current analysis)
npc_input_dir <- "results/GO_KEGG"

# Neurons (from old analysis)
neuron_base_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_rna_old/results"
neuron_deg_dir <- file.path(neuron_base_dir, "05_deseq2")
neuron_go_dir <- file.path(neuron_base_dir, "go_analysis_neu")

output_dir <- "results/GO_KEGG/comparative_plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set parameters
top_n_terms <- 20  # Number of top terms to show
min_count <- 5     # Minimum gene count to include a term
max_pval <- 0.05   # Maximum adjusted p-value to include

# ============================================================
# Function to prepare combined data for NPCs and Neurons plotting
# ============================================================
prepare_combined_celltype_data <- function(npc_up_file, npc_down_file,
                                           neuron_up_file, neuron_down_file,
                                           type = "GO") {

    cat("Processing", type, "enrichment results for NPCs and Neurons...\n")

    # Initialize empty data frame
    combined_data <- data.frame()

    # Load NPC upregulated results
    if (file.exists(npc_up_file)) {
        npc_up_data <- read.csv(npc_up_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                CellType = "NPCs",
                Direction = "UP",
                Term = Description
            )

        # Select columns - check if ONTOLOGY exists
        if ("ONTOLOGY" %in% colnames(npc_up_data)) {
            npc_up_data <- npc_up_data %>% select(Term, Count, p.adjust, CellType, Direction, ONTOLOGY)
        } else {
            npc_up_data <- npc_up_data %>% select(Term, Count, p.adjust, CellType, Direction)
        }

        cat("  - Loaded", nrow(npc_up_data), "NPC upregulated terms\n")
        combined_data <- bind_rows(combined_data, npc_up_data)
    }

    # Load NPC downregulated results
    if (file.exists(npc_down_file)) {
        npc_down_data <- read.csv(npc_down_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                CellType = "NPCs",
                Direction = "Down",
                Term = Description
            )

        # Select columns - check if ONTOLOGY exists
        if ("ONTOLOGY" %in% colnames(npc_down_data)) {
            npc_down_data <- npc_down_data %>% select(Term, Count, p.adjust, CellType, Direction, ONTOLOGY)
        } else {
            npc_down_data <- npc_down_data %>% select(Term, Count, p.adjust, CellType, Direction)
        }

        cat("  - Loaded", nrow(npc_down_data), "NPC downregulated terms\n")
        combined_data <- bind_rows(combined_data, npc_down_data)
    }

    # Load Neuron upregulated results
    if (file.exists(neuron_up_file)) {
        neuron_up_data <- read.csv(neuron_up_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                CellType = "Ns",
                Direction = "UP",
                Term = Description
            )

        # Select columns - check if ONTOLOGY exists
        if ("ONTOLOGY" %in% colnames(neuron_up_data)) {
            neuron_up_data <- neuron_up_data %>% select(Term, Count, p.adjust, CellType, Direction, ONTOLOGY)
        } else {
            neuron_up_data <- neuron_up_data %>% select(Term, Count, p.adjust, CellType, Direction)
        }

        cat("  - Loaded", nrow(neuron_up_data), "Neuron upregulated terms\n")
        combined_data <- bind_rows(combined_data, neuron_up_data)
    }

    # Load Neuron downregulated results
    if (file.exists(neuron_down_file)) {
        neuron_down_data <- read.csv(neuron_down_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                CellType = "Ns",
                Direction = "Down",
                Term = Description
            )

        # Select columns - check if ONTOLOGY exists
        if ("ONTOLOGY" %in% colnames(neuron_down_data)) {
            neuron_down_data <- neuron_down_data %>% select(Term, Count, p.adjust, CellType, Direction, ONTOLOGY)
        } else {
            neuron_down_data <- neuron_down_data %>% select(Term, Count, p.adjust, CellType, Direction)
        }

        cat("  - Loaded", nrow(neuron_down_data), "Neuron downregulated terms\n")
        combined_data <- bind_rows(combined_data, neuron_down_data)
    }

    # Get unique terms from all conditions
    unique_terms <- unique(combined_data$Term)
    cat("  - Total unique terms:", length(unique_terms), "\n")

    # Create a complete data frame with all combinations
    # This ensures we have entries for terms that appear in only some conditions
    all_combinations <- expand.grid(
        Term = unique_terms,
        CellType = c("NPCs", "Ns"),
        Direction = c("UP", "Down"),
        stringsAsFactors = FALSE
    )

    # Merge with actual data
    plot_data <- all_combinations %>%
        left_join(combined_data, by = c("Term", "CellType", "Direction"))

    # Create combined column for plotting (NPCs UP, NPCs Down, Ns UP, Ns Down)
    plot_data <- plot_data %>%
        mutate(Condition = paste(CellType, Direction, sep = " "))

    return(plot_data)
}

# ============================================================
# Function to create combined NPCs and Neurons comparative dot plot
# ============================================================
create_combined_dotplot <- function(data, title, output_file) {

    # Filter out rows with no data (NA p-values)
    data_present <- data %>% filter(!is.na(p.adjust))

    if (nrow(data_present) == 0) {
        cat("  No data to plot for", title, "\n")
        return(NULL)
    }

    # Order terms by their best (lowest) p-value across all conditions
    term_order <- data %>%
        group_by(Term) %>%
        summarize(best_pval = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
        arrange(best_pval) %>%
        pull(Term)

    data$Term <- factor(data$Term, levels = rev(term_order))

    # Set factor levels for Condition to control column order: NPCs UP, NPCs Down, Ns UP, Ns Down
    data$Condition <- factor(data$Condition, levels = c("NPCs UP", "NPCs Down", "Ns UP", "Ns Down"))
    data_present$Condition <- factor(data_present$Condition, levels = c("NPCs UP", "NPCs Down", "Ns UP", "Ns Down"))

    # Simplify term names by removing prefixes like "GO:"
    data$Term_display <- gsub("^GO:\\d+\\s*", "", data$Term)
    data_present$Term_display <- gsub("^GO:\\d+\\s*", "", data_present$Term)

    # Update term factor levels
    term_order_display <- gsub("^GO:\\d+\\s*", "", term_order)
    data$Term_display <- factor(data$Term_display, levels = rev(term_order_display))
    data_present$Term_display <- factor(data_present$Term_display, levels = rev(term_order_display))

    # Create the plot with facets for cell types
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
            axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 0),
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

    # Also save PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(pdf_file, p, width = plot_width, height = plot_height)

    cat("  Saved:", output_file, "\n")
    cat("  Saved:", pdf_file, "\n")

    return(p)
}

# ============================================================
# Process GO enrichment for NPCs and Neurons
# ============================================================
cat("\n============================================================\n")
cat("Processing GO Enrichment Results\n")
cat("============================================================\n")

# NPC GO enrichment files (current analysis)
npc_go_up_file <- file.path(npc_input_dir, "GO_enrichment_upregulated.csv")
npc_go_down_file <- file.path(npc_input_dir, "GO_enrichment_downregulated.csv")

# Check if we need to create Neuron GO enrichment files
# First, check if they already exist in the neuron GO directory
neuron_go_up_file <- file.path(neuron_go_dir, "GO_enrichment_upregulated.csv")
neuron_go_down_file <- file.path(neuron_go_dir, "GO_enrichment_downregulated.csv")

# If neuron GO files don't exist, create them from the neuron DEG files
if (!file.exists(neuron_go_up_file) || !file.exists(neuron_go_down_file)) {
    cat("\nNeuron GO enrichment files not found. Creating them from DEG files...\n")

    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(org.Hs.eg.db)
    })

    # Load neuron up/downregulated gene lists
    neuron_up_deg_file <- file.path(neuron_deg_dir, "Neuron_differential_expression_upregulated.csv")
    neuron_down_deg_file <- file.path(neuron_deg_dir, "Neuron_differential_expression_downregulated.csv")

    # Process upregulated genes
    if (file.exists(neuron_up_deg_file)) {
        neuron_up_genes <- read.csv(neuron_up_deg_file) %>%
            pull(gene_symbol) %>%
            unique() %>%
            na.omit()

        cat("  Neuron upregulated genes:", length(neuron_up_genes), "\n")

        if (length(neuron_up_genes) > 0) {
            # Convert to Entrez IDs
            neuron_up_entrez <- bitr(neuron_up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db, drop = TRUE)

            # Run GO enrichment
            go_neuron_up <- enrichGO(
                gene = neuron_up_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE
            )

            if (!is.null(go_neuron_up) && nrow(go_neuron_up) > 0) {
                write.csv(as.data.frame(go_neuron_up), neuron_go_up_file, row.names = FALSE)
                cat("  Created:", neuron_go_up_file, "\n")
            } else {
                neuron_go_up_file <- NULL
            }
        } else {
            neuron_go_up_file <- NULL
        }
    } else {
        cat("  Warning: Neuron upregulated DEG file not found\n")
        neuron_go_up_file <- NULL
    }

    # Process downregulated genes
    if (file.exists(neuron_down_deg_file)) {
        neuron_down_genes <- read.csv(neuron_down_deg_file) %>%
            pull(gene_symbol) %>%
            unique() %>%
            na.omit()

        cat("  Neuron downregulated genes:", length(neuron_down_genes), "\n")

        if (length(neuron_down_genes) > 0) {
            # Convert to Entrez IDs
            neuron_down_entrez <- bitr(neuron_down_genes, fromType = "SYMBOL", toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db, drop = TRUE)

            # Run GO enrichment
            go_neuron_down <- enrichGO(
                gene = neuron_down_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE
            )

            if (!is.null(go_neuron_down) && nrow(go_neuron_down) > 0) {
                write.csv(as.data.frame(go_neuron_down), neuron_go_down_file, row.names = FALSE)
                cat("  Created:", neuron_go_down_file, "\n")
            } else {
                neuron_go_down_file <- NULL
            }
        } else {
            neuron_go_down_file <- NULL
        }
    } else {
        cat("  Warning: Neuron downregulated DEG file not found\n")
        neuron_go_down_file <- NULL
    }
}

# Check if we have the required files
if ((file.exists(npc_go_up_file) || file.exists(npc_go_down_file)) &&
    (!is.null(neuron_go_up_file) && !is.null(neuron_go_down_file))) {

    # Prepare combined data
    go_data <- prepare_combined_celltype_data(
        npc_up_file = npc_go_up_file,
        npc_down_file = npc_go_down_file,
        neuron_up_file = neuron_go_up_file,
        neuron_down_file = neuron_go_down_file,
        type = "GO"
    )

    # Create combined comparative plot
    cat("\nCreating combined GO comparative plot (NPCs and Neurons)...\n")
    create_combined_dotplot(
        go_data,
        title = "GO Enrichment: NPCs and Neurons",
        output_file = file.path(output_dir, "GO_comparative_NPCs_Neurons.png")
    )

} else {
    cat("Required GO enrichment files not found. Skipping GO comparative plots.\n")
    if (!file.exists(npc_go_up_file)) cat("  Missing:", npc_go_up_file, "\n")
    if (!file.exists(npc_go_down_file)) cat("  Missing:", npc_go_down_file, "\n")
    if (is.null(neuron_go_up_file)) {
        cat("  Missing: Neuron upregulated GO file\n")
    } else if (!file.exists(neuron_go_up_file)) {
        cat("  Missing:", neuron_go_up_file, "\n")
    }
    if (is.null(neuron_go_down_file)) {
        cat("  Missing: Neuron downregulated GO file\n")
    } else if (!file.exists(neuron_go_down_file)) {
        cat("  Missing:", neuron_go_down_file, "\n")
    }
}

# ============================================================
# Process KEGG enrichment for NPCs and Neurons
# ============================================================
cat("\n============================================================\n")
cat("Processing KEGG Enrichment Results\n")
cat("============================================================\n")

# NPC KEGG enrichment files (current analysis)
npc_kegg_up_file <- file.path(npc_input_dir, "KEGG_enrichment_upregulated.csv")
npc_kegg_down_file <- file.path(npc_input_dir, "KEGG_enrichment_downregulated.csv")

# Neuron KEGG enrichment files
neuron_kegg_up_file <- file.path(neuron_go_dir, "KEGG_enrichment_upregulated.csv")
neuron_kegg_down_file <- file.path(neuron_go_dir, "KEGG_enrichment_downregulated.csv")

# If neuron KEGG files don't exist, try alternative locations
if (!file.exists(neuron_kegg_up_file) || !file.exists(neuron_kegg_down_file)) {
    cat("\nNeuron KEGG enrichment files not found in expected location.\n")
    cat("Skipping KEGG comparative plots (neuron data not available).\n")
    neuron_kegg_up_file <- NULL
    neuron_kegg_down_file <- NULL
}

# Check if we have the required files
if ((file.exists(npc_kegg_up_file) || file.exists(npc_kegg_down_file)) &&
    (!is.null(neuron_kegg_up_file) && !is.null(neuron_kegg_down_file))) {

    # Prepare combined data
    kegg_data <- prepare_combined_celltype_data(
        npc_up_file = npc_kegg_up_file,
        npc_down_file = npc_kegg_down_file,
        neuron_up_file = neuron_kegg_up_file,
        neuron_down_file = neuron_kegg_down_file,
        type = "KEGG"
    )

    # Create combined comparative plot
    cat("\nCreating combined KEGG comparative plot (NPCs and Neurons)...\n")
    create_combined_dotplot(
        kegg_data,
        title = "KEGG Pathway Enrichment: NPCs and Neurons",
        output_file = file.path(output_dir, "KEGG_comparative_NPCs_Neurons.png")
    )

} else {
    cat("Required KEGG enrichment files not found. Skipping KEGG comparative plots.\n")
}

cat("\n==========================================\n")
cat("Comparative enrichment plots complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nGenerated plots:\n")
cat("  - GO_comparative_NPCs_Neurons.png: Combined GO comparison for NPCs and Neurons\n")
cat("  - KEGG_comparative_NPCs_Neurons.png: Combined KEGG comparison for NPCs and Neurons\n")
cat("\nAll plots saved in both PNG and PDF formats\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
