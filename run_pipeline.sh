#!/bin/bash

# Master RNA-seq Analysis Pipeline
# This script orchestrates the complete RNA-seq analysis workflow
# Usage: bash run_pipeline.sh [step]
#   - If no step is provided, runs the entire pipeline
#   - If a step number (1-7) is provided, runs only that step

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    local color=$1
    shift
    echo -e "${color}$@${NC}"
}

print_message $BLUE "=========================================="
print_message $BLUE "RNA-seq Analysis Pipeline"
print_message $BLUE "Control vs Mutant NPCs"
print_message $BLUE "=========================================="
echo ""

# Check if a specific step is requested
STEP=${1:-all}

# Function to submit a job and return job ID
submit_job() {
    local script=$1
    local job_id=$(sbatch --parsable ${script})
    echo ${job_id}
}

# Function to submit a job with dependency
submit_job_with_dependency() {
    local script=$1
    local dependency=$2
    local job_id=$(sbatch --parsable --dependency=afterok:${dependency} ${script})
    echo ${job_id}
}

# Initialize job IDs
JOB_FASTQC=""
JOB_STAR=""
JOB_COUNTS=""
JOB_MULTIQC=""
JOB_DESEQ2=""
JOB_ENRICHMENT=""
JOB_VIS=""

# ============================================================
# Step 1: Quality Control with FastQC
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "1" ]]; then
    print_message $GREEN "\n[Step 1/7] Submitting FastQC job..."
    JOB_FASTQC=$(submit_job scripts/01_fastqc.sh)
    print_message $YELLOW "  Job ID: ${JOB_FASTQC}"

    if [[ ${STEP} == "1" ]]; then
        print_message $BLUE "\nSubmitted Step 1 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 2: STAR Alignment
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "2" ]]; then
    print_message $GREEN "\n[Step 2/7] Submitting STAR alignment job array..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_FASTQC} ]]; then
        # Run after FastQC completes (optional, can run in parallel)
        JOB_STAR=$(submit_job scripts/02_star_alignment.sh)
    else
        JOB_STAR=$(submit_job scripts/02_star_alignment.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_STAR}"

    if [[ ${STEP} == "2" ]]; then
        print_message $BLUE "\nSubmitted Step 2 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 3: Gene Quantification with featureCounts
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "3" ]]; then
    print_message $GREEN "\n[Step 3/7] Submitting featureCounts job..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_STAR} ]]; then
        # Must run after STAR completes
        JOB_COUNTS=$(submit_job_with_dependency scripts/03_featurecounts.sh ${JOB_STAR})
    else
        JOB_COUNTS=$(submit_job scripts/03_featurecounts.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_COUNTS}"
    print_message $YELLOW "  Depends on: ${JOB_STAR}"

    if [[ ${STEP} == "3" ]]; then
        print_message $BLUE "\nSubmitted Step 3 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 4: MultiQC Report
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "4" ]]; then
    print_message $GREEN "\n[Step 4/7] Submitting MultiQC job..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_COUNTS} ]]; then
        # Run after all QC steps complete
        JOB_MULTIQC=$(submit_job_with_dependency scripts/04_multiqc.sh ${JOB_COUNTS})
    else
        JOB_MULTIQC=$(submit_job scripts/04_multiqc.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_MULTIQC}"
    print_message $YELLOW "  Depends on: ${JOB_COUNTS}"

    if [[ ${STEP} == "4" ]]; then
        print_message $BLUE "\nSubmitted Step 4 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 5: DESeq2 Differential Expression Analysis
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "5" ]]; then
    print_message $GREEN "\n[Step 5/7] Submitting DESeq2 analysis job..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_COUNTS} ]]; then
        # Must run after featureCounts completes
        JOB_DESEQ2=$(submit_job_with_dependency scripts/05_deseq2_analysis.sh ${JOB_COUNTS})
    else
        JOB_DESEQ2=$(submit_job scripts/05_deseq2_analysis.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_DESEQ2}"
    print_message $YELLOW "  Depends on: ${JOB_COUNTS}"

    if [[ ${STEP} == "5" ]]; then
        print_message $BLUE "\nSubmitted Step 5 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 6: Functional Enrichment Analysis (GO and KEGG)
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "6" ]]; then
    print_message $GREEN "\n[Step 6/7] Submitting functional enrichment job..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_DESEQ2} ]]; then
        # Must run after DESeq2 completes
        JOB_ENRICHMENT=$(submit_job_with_dependency scripts/06_functional_enrichment.sh ${JOB_DESEQ2})
    else
        JOB_ENRICHMENT=$(submit_job scripts/06_functional_enrichment.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_ENRICHMENT}"
    print_message $YELLOW "  Depends on: ${JOB_DESEQ2}"

    if [[ ${STEP} == "6" ]]; then
        print_message $BLUE "\nSubmitted Step 6 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Step 7: Generate Visualizations
# ============================================================
if [[ ${STEP} == "all" ]] || [[ ${STEP} == "7" ]]; then
    print_message $GREEN "\n[Step 7/7] Submitting visualization job..."

    if [[ ${STEP} == "all" ]] && [[ -n ${JOB_DESEQ2} ]]; then
        # Must run after DESeq2 completes
        JOB_VIS=$(submit_job_with_dependency scripts/07_visualizations.sh ${JOB_DESEQ2})
    else
        JOB_VIS=$(submit_job scripts/07_visualizations.sh)
    fi

    print_message $YELLOW "  Job ID: ${JOB_VIS}"
    print_message $YELLOW "  Depends on: ${JOB_DESEQ2}"

    if [[ ${STEP} == "7" ]]; then
        print_message $BLUE "\nSubmitted Step 7 only. Exiting."
        exit 0
    fi
fi

# ============================================================
# Summary
# ============================================================
print_message $BLUE "\n=========================================="
print_message $GREEN "Pipeline submitted successfully!"
print_message $BLUE "=========================================="
echo ""
print_message $YELLOW "Submitted jobs:"
[[ -n ${JOB_FASTQC} ]] && echo "  Step 1 (FastQC):       ${JOB_FASTQC}"
[[ -n ${JOB_STAR} ]] && echo "  Step 2 (STAR):         ${JOB_STAR}"
[[ -n ${JOB_COUNTS} ]] && echo "  Step 3 (featureCounts): ${JOB_COUNTS}"
[[ -n ${JOB_MULTIQC} ]] && echo "  Step 4 (MultiQC):      ${JOB_MULTIQC}"
[[ -n ${JOB_DESEQ2} ]] && echo "  Step 5 (DESeq2):       ${JOB_DESEQ2}"
[[ -n ${JOB_ENRICHMENT} ]] && echo "  Step 6 (Enrichment):   ${JOB_ENRICHMENT}"
[[ -n ${JOB_VIS} ]] && echo "  Step 7 (Visualization): ${JOB_VIS}"
echo ""
print_message $YELLOW "Monitor job status with: squeue -u \$USER"
print_message $YELLOW "Check logs in: ./logs/"
echo ""
print_message $GREEN "Expected runtime: 8-12 hours for complete pipeline"
print_message $BLUE "=========================================="
