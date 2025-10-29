#!/usr/bin/env Rscript

# Helper script to download mouse-human ortholog mappings
# This can be run separately if the main script has network issues

suppressPackageStartupMessages({
    library(tidyverse)
})

cat("==========================================\n")
cat("Downloading Mouse-Human Ortholog Mappings\n")
cat("==========================================\n\n")

OUTPUT_FILE <- "results/venn_diagrams/mouse_human_orthologs.csv"
dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

# Method 1: Try downloading from Ensembl BioMart
cat("Attempting to download from Ensembl BioMart...\n")

tryCatch({
    url <- "https://www.ensembl.org/biomart/martservice"

    query <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="mmusculus_gene_ensembl" interface="default">
<Attribute name="mgi_symbol"/>
<Attribute name="hsapiens_homolog_associated_gene_name"/>
<Attribute name="hsapiens_homolog_orthology_type"/>
<Attribute name="hsapiens_homolog_orthology_confidence"/>
</Dataset>
</Query>'

    temp_file <- tempfile()
    cat("  Downloading...\n")
    download.file(
        url = paste0(url, "?query=", URLencode(query)),
        destfile = temp_file,
        quiet = FALSE
    )

    orthologs <- read.delim(temp_file, stringsAsFactors = FALSE)
    colnames(orthologs) <- c("mouse_gene", "human_gene", "orthology_type", "confidence")

    # Filter for high-confidence orthologs
    orthologs <- orthologs %>%
        filter(
            mouse_gene != "",
            human_gene != "",
            !is.na(human_gene),
            confidence == 1  # High confidence orthologs
        ) %>%
        select(mouse_gene, human_gene, orthology_type) %>%
        distinct()

    cat("  Downloaded", nrow(orthologs), "ortholog mappings\n")
    cat("  Unique mouse genes:", length(unique(orthologs$mouse_gene)), "\n")
    cat("  Unique human genes:", length(unique(orthologs$human_gene)), "\n")

    # Save
    write.csv(orthologs, OUTPUT_FILE, row.names = FALSE)
    cat("\nSUCCESS! Saved ortholog mappings to:", OUTPUT_FILE, "\n")

    # Show some examples
    cat("\nExample mappings:\n")
    print(head(orthologs, 10))

    unlink(temp_file)

}, error = function(e) {
    cat("\nERROR downloading from Ensembl:", conditionMessage(e), "\n")
    cat("\nAlternative: You can manually download ortholog data from:\n")
    cat("  https://www.ensembl.org/biomart/martview\n")
    cat("  1. Select 'Mouse genes (GRCm39)' database\n")
    cat("  2. Add attributes: MGI symbol, Human ortholog gene name\n")
    cat("  3. Download as CSV and save to:", OUTPUT_FILE, "\n")
})

cat("\n==========================================\n")
