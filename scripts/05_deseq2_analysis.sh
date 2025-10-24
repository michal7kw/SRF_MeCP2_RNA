#!/bin/bash
#SBATCH --job-name=05_DESeq2
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="./logs/05_deseq2.err"
#SBATCH --output="./logs/05_deseq2.out"

# RNA-seq Pipeline Step 5: Run DESeq2 Analysis
# This script runs the R script for differential expression analysis

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting DESeq2 Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Run the R script
Rscript scripts/05_deseq2_analysis.R

echo ""
echo "End time: $(date)"
echo "=========================================="
