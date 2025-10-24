#!/bin/bash

# Pipeline Status Checker
# This script provides a quick overview of pipeline progress

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=========================================="
echo "RNA-seq Pipeline Status Check"
echo -e "==========================================${NC}\n"

# Function to check if files exist
check_files() {
    local pattern=$1
    local count=$(find $pattern -type f 2>/dev/null | wc -l)
    echo $count
}

# Function to print status
print_status() {
    local step=$1
    local expected=$2
    local found=$3
    local desc=$4

    if [ $found -eq 0 ]; then
        echo -e "${RED}✗${NC} Step $step: $desc (0/$expected files)"
    elif [ $found -lt $expected ]; then
        echo -e "${YELLOW}⚠${NC} Step $step: $desc ($found/$expected files)"
    else
        echo -e "${GREEN}✓${NC} Step $step: $desc ($found/$expected files)"
    fi
}

echo -e "${BLUE}Pipeline Steps:${NC}"
echo "---------------"

# Step 1: FastQC
fastqc_count=$(check_files "results/fastqc/*.html")
print_status "1" "12" "$fastqc_count" "FastQC - Quality Control"

# Step 2: STAR Alignment
bam_count=$(check_files "results/alignment/*_Aligned.sortedByCoord.out.bam")
print_status "2" "6" "$bam_count" "STAR - Read Alignment"

# Step 3: featureCounts
if [ -f "results/counts/gene_counts.txt" ]; then
    echo -e "${GREEN}✓${NC} Step 3: featureCounts - Gene Quantification (1/1 files)"
else
    echo -e "${RED}✗${NC} Step 3: featureCounts - Gene Quantification (0/1 files)"
fi

# Step 4: MultiQC
if [ -f "results/multiqc/multiqc_report.html" ]; then
    echo -e "${GREEN}✓${NC} Step 4: MultiQC - QC Report (1/1 files)"
else
    echo -e "${RED}✗${NC} Step 4: MultiQC - QC Report (0/1 files)"
fi

# Step 5: DESeq2
deseq2_count=$(check_files "results/DESeq2/*.csv")
print_status "5" "3" "$deseq2_count" "DESeq2 - Differential Expression"

# Step 6: GO/KEGG
enrichment_count=$(check_files "results/GO_KEGG/*.csv")
if [ $enrichment_count -gt 0 ]; then
    echo -e "${GREEN}✓${NC} Step 6: GO/KEGG - Functional Enrichment ($enrichment_count files)"
else
    echo -e "${RED}✗${NC} Step 6: GO/KEGG - Functional Enrichment (0 files)"
fi

# Step 7: Visualizations
plot_count=$(check_files "results/plots/*.png")
print_status "7" "8" "$plot_count" "Visualizations - Plots"

echo ""
echo -e "${BLUE}Current Jobs:${NC}"
echo "-------------"
squeue -u $USER 2>/dev/null || echo "No jobs running"

echo ""
echo -e "${BLUE}Recent Job History:${NC}"
echo "-------------------"
sacct -u $USER --format=JobID,JobName%20,State,Elapsed,MaxRSS -S $(date -d '1 day ago' +%Y-%m-%d) 2>/dev/null | tail -15 || echo "sacct not available"

echo ""
echo -e "${BLUE}Disk Usage:${NC}"
echo "-----------"
echo "Results directory: $(du -sh results/ 2>/dev/null | cut -f1)"
echo "Logs directory:    $(du -sh logs/ 2>/dev/null | cut -f1)"
echo "Total:             $(du -sh . 2>/dev/null | cut -f1)"

echo ""
echo -e "${BLUE}Key Results Summary:${NC}"
echo "--------------------"

if [ -f "results/DESeq2/DESeq2_results_significant.csv" ]; then
    sig_genes=$(tail -n +2 results/DESeq2/DESeq2_results_significant.csv | wc -l)
    echo -e "Significant genes (padj < 0.05): ${GREEN}$sig_genes${NC}"
fi

if [ -f "results/DESeq2/DESeq2_results_significant_FC2.csv" ]; then
    sig_genes_fc=$(tail -n +2 results/DESeq2/DESeq2_results_significant_FC2.csv | wc -l)
    echo -e "Significant genes (padj < 0.05, |FC| > 2): ${GREEN}$sig_genes_fc${NC}"
fi

if [ -f "results/GO_KEGG/GO_enrichment_all_genes.csv" ]; then
    go_terms=$(tail -n +2 results/GO_KEGG/GO_enrichment_all_genes.csv | wc -l)
    echo -e "Enriched GO terms: ${GREEN}$go_terms${NC}"
fi

if [ -f "results/GO_KEGG/KEGG_enrichment_all_genes.csv" ]; then
    kegg_pathways=$(tail -n +2 results/GO_KEGG/KEGG_enrichment_all_genes.csv | wc -l)
    echo -e "Enriched KEGG pathways: ${GREEN}$kegg_pathways${NC}"
fi

echo ""
echo -e "${YELLOW}Tip: Run 'bash run_pipeline.sh' to start the pipeline${NC}"
echo -e "${YELLOW}Tip: Run 'watch -n 5 bash check_pipeline_status.sh' for live monitoring${NC}"
echo -e "${BLUE}==========================================${NC}"
