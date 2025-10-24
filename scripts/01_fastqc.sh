#!/bin/bash
#SBATCH --job-name=01_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="./logs/01_fastqc.err"
#SBATCH --output="./logs/01_fastqc.out"

# RNA-seq Pipeline Step 1: Quality Control with FastQC
# This script runs FastQC on all raw FASTQ files

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/quality

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting FastQC Quality Control"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Define paths
FASTQ_DIR="DATA/90-1244168066/00_fastq"
OUTPUT_DIR="results/fastqc"
THREADS=32

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Get list of all FASTQ files
FASTQ_FILES=(${FASTQ_DIR}/*.fastq.gz)

echo "Found ${#FASTQ_FILES[@]} FASTQ files to process"
echo ""

# Run FastQC on all files
echo "Running FastQC..."
fastqc \
    --outdir ${OUTPUT_DIR} \
    --threads ${THREADS} \
    --quiet \
    ${FASTQ_FILES[@]}

echo ""
echo "FastQC analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "End time: $(date)"
echo "=========================================="
