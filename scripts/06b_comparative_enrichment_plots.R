#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 6b: Comparative Enrichment Plots
# This script creates comparative dot plots showing GO/KEGG enrichment
# for upregulated and downregulated genes side by side

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

cat("==========================================\n")
cat("Creating Comparative Enrichment Plots\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths
input_dir <- "results/GO_KEGG"
output_dir <- "results/GO_KEGG/comparative_plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set parameters
top_n_terms <- 20  # Number of top terms to show
min_count <- 5     # Minimum gene count to include a term
max_pval <- 0.05   # Maximum adjusted p-value to include

# ============================================================
# Function to prepare data for comparative plotting
# ============================================================
prepare_comparative_data <- function(up_file, down_file, type = "GO") {

    cat("Processing", type, "enrichment results...\n")

    # Initialize empty data frame
    combined_data <- data.frame()

    # Load upregulated results
    if (file.exists(up_file)) {
        up_data <- read.csv(up_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                Direction = "UP",
                Term = Description
            ) %>%
            select(Term, Count, p.adjust, Direction, ONTOLOGY)

        cat("  - Loaded", nrow(up_data), "upregulated terms\n")
        combined_data <- bind_rows(combined_data, up_data)
    }

    # Load downregulated results
    if (file.exists(down_file)) {
        down_data <- read.csv(down_file) %>%
            filter(Count >= min_count, p.adjust <= max_pval) %>%
            arrange(p.adjust) %>%
            head(top_n_terms) %>%
            mutate(
                Direction = "Down",
                Term = Description
            ) %>%
            select(Term, Count, p.adjust, Direction, ONTOLOGY)

        cat("  - Loaded", nrow(down_data), "downregulated terms\n")
        combined_data <- bind_rows(combined_data, down_data)
    }

    # Get unique terms from both directions
    unique_terms <- unique(combined_data$Term)
    cat("  - Total unique terms:", length(unique_terms), "\n")

    # Create a complete data frame with all combinations
    # This ensures we have entries for terms that appear in only one direction
    all_combinations <- expand.grid(
        Term = unique_terms,
        Direction = c("UP", "Down"),
        stringsAsFactors = FALSE
    )

    # Merge with actual data
    plot_data <- all_combinations %>%
        left_join(combined_data, by = c("Term", "Direction"))

    return(plot_data)
}

# ============================================================
# Function to create comparative dot plot
# ============================================================
create_comparative_dotplot <- function(data, title, output_file, show_ontology = TRUE) {

    # Filter out rows with no data (NA p-values)
    data_present <- data %>% filter(!is.na(p.adjust))

    if (nrow(data_present) == 0) {
        cat("  No data to plot for", title, "\n")
        return(NULL)
    }

    # Order terms by their best (lowest) p-value across all directions
    term_order <- data %>%
        group_by(Term) %>%
        summarize(best_pval = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
        arrange(best_pval) %>%
        pull(Term)

    data$Term <- factor(data$Term, levels = rev(term_order))
    data$Direction <- factor(data$Direction, levels = c("UP", "Down"))

    # Add ontology labels if available and requested
    if (show_ontology && "ONTOLOGY" %in% colnames(data)) {
        # Get ontology for each term (take first non-NA value)
        term_ontology <- data %>%
            filter(!is.na(ONTOLOGY)) %>%
            group_by(Term) %>%
            summarize(ont = first(ONTOLOGY), .groups = "drop")

        # Add abbreviated ontology to term labels
        data <- data %>%
            left_join(term_ontology, by = "Term") %>%
            mutate(Term_label = if_else(!is.na(ont), paste0(Term, " (", ont, ")"), as.character(Term)))

        data$Term_label <- factor(data$Term_label,
                                   levels = rev(paste0(term_order, " (",
                                                      term_ontology$ont[match(term_order, term_ontology$Term)], ")")))
        y_var <- "Term_label"
    } else {
        y_var <- "Term"
    }

    # Create the plot
    p <- ggplot(data_present, aes_string(x = "Direction", y = y_var)) +
        geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
        scale_color_gradient(
            low = "blue",
            high = "red",
            name = "Adjusted p-value",
            trans = "log10",
            breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
            labels = c("1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2"),
            limits = c(1e-7, max(data_present$p.adjust, na.rm = TRUE))
        ) +
        scale_size_continuous(
            name = "Count",
            breaks = c(10, 25, 50, 100),
            range = c(2, 10)
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 11, face = "bold"),
            axis.title = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right"
        ) +
        labs(title = title)

    # Save plot
    ggsave(output_file, p, width = 12, height = max(8, nrow(data_present) * 0.3), dpi = 300)

    # Also save PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(pdf_file, p, width = 12, height = max(8, nrow(data_present) * 0.3))

    cat("  Saved:", output_file, "\n")
    cat("  Saved:", pdf_file, "\n")

    return(p)
}

# ============================================================
# Function to create simplified comparative plot (like the example image)
# ============================================================
create_simplified_comparative_plot <- function(data, title, output_file,
                                               selected_terms = NULL,
                                               term_labels = NULL) {

    # If specific terms are provided, filter to those
    if (!is.null(selected_terms)) {
        data <- data %>% filter(Term %in% selected_terms)
    }

    # Filter out rows with no data
    data_present <- data %>% filter(!is.na(p.adjust))

    if (nrow(data_present) == 0) {
        cat("  No data to plot for", title, "\n")
        return(NULL)
    }

    # Apply custom term labels if provided
    if (!is.null(term_labels)) {
        data <- data %>%
            mutate(Term_display = term_labels[match(Term, names(term_labels))])
        data$Term_display[is.na(data$Term_display)] <- data$Term[is.na(data$Term_display)]
    } else {
        data$Term_display <- data$Term
    }

    # Order terms by their best (lowest) p-value
    term_order <- data %>%
        group_by(Term_display) %>%
        summarize(best_pval = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(best_pval)) %>%  # Reverse order for plotting
        pull(Term_display)

    data$Term_display <- factor(data$Term_display, levels = term_order)
    data$Direction <- factor(data$Direction, levels = c("UP", "Down"))

    # Create the plot with styling similar to the example image
    p <- ggplot(data, aes(x = Direction, y = Term_display)) +
        geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
        scale_color_gradient(
            low = "#4575b4",      # Blue for significant
            high = "#d73027",     # Red for less significant
            name = "Adjusted p-value",
            trans = "log10"
        ) +
        scale_size_continuous(
            name = "Count",
            range = c(3, 12)
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.y = element_text(size = 11, color = "black"),
            axis.text.x = element_text(size = 12, face = "bold", color = "black"),
            axis.title = element_blank(),
            panel.grid.major.x = element_line(colour = "grey85"),
            panel.grid.major.y = element_line(colour = "grey85"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "bottom",
            legend.box = "horizontal",
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        ) +
        labs(title = title)

    # Save plot
    ggsave(output_file, p, width = 8, height = max(6, length(unique(data$Term_display)) * 0.4), dpi = 300)

    # Also save PDF
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(pdf_file, p, width = 8, height = max(6, length(unique(data$Term_display)) * 0.4))

    cat("  Saved:", output_file, "\n")
    cat("  Saved:", pdf_file, "\n")

    return(p)
}

# ============================================================
# Process GO enrichment
# ============================================================
cat("\n============================================================\n")
cat("Processing GO Enrichment Results\n")
cat("============================================================\n")

go_up_file <- file.path(input_dir, "GO_enrichment_upregulated.csv")
go_down_file <- file.path(input_dir, "GO_enrichment_downregulated.csv")

if (file.exists(go_up_file) || file.exists(go_down_file)) {
    go_data <- prepare_comparative_data(go_up_file, go_down_file, type = "GO")

    # Create full comparative plot with ontology labels
    cat("\nCreating full GO comparative plot...\n")
    create_comparative_dotplot(
        go_data,
        title = "GO Enrichment: Upregulated vs Downregulated",
        output_file = file.path(output_dir, "GO_comparative_full.png"),
        show_ontology = TRUE
    )

    # Create simplified plot (top terms, simplified labels)
    cat("\nCreating simplified GO comparative plot...\n")

    # Select top terms by p-value from each direction
    top_go_up <- go_data %>%
        filter(Direction == "UP", !is.na(p.adjust)) %>%
        arrange(p.adjust) %>%
        head(15) %>%
        pull(Term)

    top_go_down <- go_data %>%
        filter(Direction == "Down", !is.na(p.adjust)) %>%
        arrange(p.adjust) %>%
        head(15) %>%
        pull(Term)

    selected_go_terms <- unique(c(top_go_up, top_go_down))

    create_simplified_comparative_plot(
        go_data,
        title = "GO Enrichment: Mutant vs Control",
        output_file = file.path(output_dir, "GO_comparative_simplified.png"),
        selected_terms = selected_go_terms
    )

} else {
    cat("GO enrichment files not found. Skipping GO comparative plots.\n")
}

# ============================================================
# Process KEGG enrichment
# ============================================================
cat("\n============================================================\n")
cat("Processing KEGG Enrichment Results\n")
cat("============================================================\n")

kegg_up_file <- file.path(input_dir, "KEGG_enrichment_upregulated.csv")
kegg_down_file <- file.path(input_dir, "KEGG_enrichment_downregulated.csv")

if (file.exists(kegg_up_file) || file.exists(kegg_down_file)) {
    kegg_data <- prepare_comparative_data(kegg_up_file, kegg_down_file, type = "KEGG")

    # For KEGG, use Description as Term
    if ("Description" %in% colnames(read.csv(kegg_up_file))) {
        # Re-read and process KEGG data
        combined_kegg <- data.frame()

        if (file.exists(kegg_up_file)) {
            kegg_up <- read.csv(kegg_up_file) %>%
                filter(Count >= min_count, p.adjust <= max_pval) %>%
                arrange(p.adjust) %>%
                head(top_n_terms) %>%
                mutate(
                    Direction = "UP",
                    Term = Description
                ) %>%
                select(Term, Count, p.adjust, Direction)
            combined_kegg <- bind_rows(combined_kegg, kegg_up)
        }

        if (file.exists(kegg_down_file)) {
            kegg_down <- read.csv(kegg_down_file) %>%
                filter(Count >= min_count, p.adjust <= max_pval) %>%
                arrange(p.adjust) %>%
                head(top_n_terms) %>%
                mutate(
                    Direction = "Down",
                    Term = Description
                ) %>%
                select(Term, Count, p.adjust, Direction)
            combined_kegg <- bind_rows(combined_kegg, kegg_down)
        }

        # Get unique terms
        unique_kegg_terms <- unique(combined_kegg$Term)

        # Create complete data frame
        all_kegg_combinations <- expand.grid(
            Term = unique_kegg_terms,
            Direction = c("UP", "Down"),
            stringsAsFactors = FALSE
        )

        kegg_data <- all_kegg_combinations %>%
            left_join(combined_kegg, by = c("Term", "Direction"))
    }

    # Create KEGG comparative plots
    cat("\nCreating KEGG comparative plot...\n")
    create_comparative_dotplot(
        kegg_data,
        title = "KEGG Pathway Enrichment: Upregulated vs Downregulated",
        output_file = file.path(output_dir, "KEGG_comparative_full.png"),
        show_ontology = FALSE
    )

    # Create simplified plot
    cat("\nCreating simplified KEGG comparative plot...\n")

    top_kegg_up <- kegg_data %>%
        filter(Direction == "UP", !is.na(p.adjust)) %>%
        arrange(p.adjust) %>%
        head(15) %>%
        pull(Term)

    top_kegg_down <- kegg_data %>%
        filter(Direction == "Down", !is.na(p.adjust)) %>%
        arrange(p.adjust) %>%
        head(15) %>%
        pull(Term)

    selected_kegg_terms <- unique(c(top_kegg_up, top_kegg_down))

    create_simplified_comparative_plot(
        kegg_data,
        title = "KEGG Pathway Enrichment: Mutant vs Control",
        output_file = file.path(output_dir, "KEGG_comparative_simplified.png"),
        selected_terms = selected_kegg_terms
    )

} else {
    cat("KEGG enrichment files not found. Skipping KEGG comparative plots.\n")
}

cat("\n==========================================\n")
cat("Comparative enrichment plots complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nGenerated plots:\n")
cat("  - GO_comparative_full.png: Full GO comparison with ontology labels\n")
cat("  - GO_comparative_simplified.png: Simplified GO comparison (top terms)\n")
cat("  - KEGG_comparative_full.png: Full KEGG pathway comparison\n")
cat("  - KEGG_comparative_simplified.png: Simplified KEGG comparison (top terms)\n")
cat("\nAll plots saved in both PNG and PDF formats\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
