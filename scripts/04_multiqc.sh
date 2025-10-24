#!/bin/bash
#SBATCH --job-name=04_multiQC
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error="./logs/04_multiqc.err"
#SBATCH --output="./logs/04_multiqc.out"

# RNA-seq Pipeline Step 4: Aggregate QC Reports with MultiQC
# This script aggregates all QC reports into a single interactive HTML report

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/multiqc

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting MultiQC Report Generation"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Define paths
OUTPUT_DIR="results/multiqc"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "Aggregating QC reports from:"
echo "  - FastQC results"
echo "  - STAR alignment logs"
echo "  - featureCounts summary"
echo ""

# Run MultiQC
echo "Running MultiQC..."
multiqc \
    . \
    --outdir ${OUTPUT_DIR} \
    --filename multiqc_report \
    --title "RNA-seq Analysis: Control vs Mutant NPCs" \
    --comment "Human Neural Progenitor Cells - Differential Expression Analysis" \
    --force \
    --verbose

echo ""
echo "MultiQC report generated successfully!"
echo "Report location: ${OUTPUT_DIR}/multiqc_report.html"
echo ""
echo "To view the report, download it to your local machine and open in a web browser."
echo ""
echo "End time: $(date)"
echo "=========================================="
