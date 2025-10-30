#!/bin/bash
#SBATCH --job-name=07_compare_celltype_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="./logs/07_compare_celltype_enrichment.err"
#SBATCH --output="./logs/07_compare_celltype_enrichment.out"

# RNA-seq Pipeline Step 6: Run Functional Enrichment Analysis
# This script runs the R script for GO and KEGG enrichment analysis

echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-gsea

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting Functional Enrichment Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Run the R script
Rscript scripts/07_compare_celltype_enrichment.R

echo ""
echo "End time: $(date)"
echo "=========================================="
