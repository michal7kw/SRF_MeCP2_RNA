#!/bin/bash
#SBATCH --job-name=08_venn_diagrams_human_mouse
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="./logs/08_venn_diagrams_human_mouse.err"
#SBATCH --output="./logs/08_venn_diagrams_human_mouse.out"

# This script runs the R script to create all plots

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-data-vis

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting Visualization Generation"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Run the R script
Rscript scripts/08_venn_diagrams_human_mouse.R

echo ""
echo "End time: $(date)"
echo "=========================================="
