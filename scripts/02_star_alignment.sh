#!/bin/bash
#SBATCH --job-name=02_STAR_align
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="./logs/02_star_alignment_%a.err"
#SBATCH --output="./logs/02_star_alignment_%a.out"
#SBATCH --array=1-6

# RNA-seq Pipeline Step 2: STAR Alignment
# This script aligns reads to the human reference genome using STAR
# Uses SLURM job arrays to process samples in parallel

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Change to project directory
PROJECT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_top/SRF_MeCP2_RNA"
cd ${PROJECT_DIR}

echo "=========================================="
echo "Starting STAR Alignment"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Processing array task: ${SLURM_ARRAY_TASK_ID}"
echo ""

# Define paths
FASTQ_DIR="DATA/90-1244168066/00_fastq"
STAR_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCh38"
GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/downloads_grch38_gencode_v44/gencode.v44.annotation.gtf"
OUTPUT_DIR="results/alignment"
THREADS=16

# # Create output directory
# mkdir -p ${OUTPUT_DIR}

# Sample names array
SAMPLES=(G1 G2 G3 M1 M2 M3)

# Get sample name for this array task
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Processing sample: ${SAMPLE}"
echo ""

# Define input files
R1="${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

# Define output prefix
OUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}_"

# Check if input files exist
if [[ ! -f ${R1} ]] || [[ ! -f ${R2} ]]; then
    echo "ERROR: Input files not found for sample ${SAMPLE}"
    exit 1
fi

echo "Input R1: ${R1}"
echo "Input R2: ${R2}"
echo "Output prefix: ${OUT_PREFIX}"
echo ""

# Run STAR alignment
echo "Running STAR alignment..."
STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${STAR_INDEX} \
    --sjdbGTFfile ${GTF} \
    --readFilesIn ${R1} ${R2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUT_PREFIX} \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --limitBAMsortRAM 60000000000

echo ""
echo "Creating BAM index..."
samtools index -@ ${THREADS} ${OUT_PREFIX}Aligned.sortedByCoord.out.bam

echo ""
echo "Alignment complete for sample: ${SAMPLE}"
echo "Output files:"
ls -lh ${OUT_PREFIX}*
echo ""
echo "End time: $(date)"
echo "=========================================="
