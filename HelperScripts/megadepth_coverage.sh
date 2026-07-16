#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64gb
#SBATCH --tmp=100gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

# This script calculates coverage in fixed-size windows for a list of BAM files
# using megadepth. Run as a Slurm job array where each task
# processes one BAM file.
#
# Submit: sbatch --array=0-$(( $(wc -l < bam_list.txt) - 1 )) coverage_windows.sh

# Dependencies
# megadepth version 1.2.x
export PATH=${PATH}:/users/6/pmorrell/Apps/megadepth

#------------------
# User provided input arguments

# File containing one BAM (or CRAM) path per line, no header
BAM_LIST="/scratch.global/pmorrell/Cowpea_Diversity/bam_files.txt"

# Output directory
OUT_DIR="/scratch.global/pmorrell/Cowpea_Diversity/coverage_output"

# Window size in bp for coverage estimation
WIN_SIZE="100"

#------------------
mkdir -p ${OUT_DIR}

if [[ ! -f "${BAM_LIST}" ]]; then
    echo "ERROR: BAM list file not found: ${BAM_LIST}" >&2
    exit 1
fi

if ! command -v megadepth >/dev/null 2>&1; then
    echo "ERROR: megadepth not found on PATH." >&2
    echo "Current PATH: ${PATH}" >&2
    exit 1
fi

# Prepare array for Slurm job array
mapfile -t BAM_ARR < <(grep -v '^[[:space:]]*$' "${BAM_LIST}")

if [[ ${#BAM_ARR[@]} -eq 0 ]]; then
    echo "ERROR: BAM_LIST is empty: ${BAM_LIST}" >&2
    exit 1
fi

# Determine maximum array limit
MAX_ARRAY_LIMIT=$((${#BAM_ARR[@]} - 1))
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Not running as array job - get index from command-line argument or default to 0
    if [[ $# -gt 0 ]]; then
        SLURM_ARRAY_TASK_ID="$1"
        echo "Using provided BAM index: ${SLURM_ARRAY_TASK_ID}"
    else
        echo "WARNING: SLURM_ARRAY_TASK_ID is not set and no index provided."
        echo "Option 1 - Submit as array job:"
        echo "  sbatch --array=0-${MAX_ARRAY_LIMIT} megadepth_coverage.sh"
        echo "Option 2 - Process a single BAM file:"
        echo "  sbatch megadepth_coverage.sh 0  # (or any index 0-${MAX_ARRAY_LIMIT})"
        exit 1
    fi
fi

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID > MAX_ARRAY_LIMIT )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is out of range 0-${MAX_ARRAY_LIMIT}." >&2
    exit 1
fi

# Get the current BAM file we are processing
CURR_BAM=${BAM_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing BAM file: ${CURR_BAM}"

if [[ ! -f "${CURR_BAM}" ]]; then
    echo "ERROR: BAM file not found: ${CURR_BAM}" >&2
    exit 1
fi

#------------------
function run_megadepth() {
    local bam_file="$1"
    local win_size="$2"
    local threads="${SLURM_NTASKS_PER_NODE:-1}"
    local sample_prefix
    local output_prefix
    local output_file
    sample_prefix=$(basename "${bam_file}" .bam)
    sample_prefix=$(basename "${sample_prefix}" .cram)
    sample_prefix=$(basename "${sample_prefix}" .fastq.gz)
    output_prefix="${sample_prefix}.megadepth.${win_size}bp"
    output_file="${output_prefix}.annotation.bed"
    set -x # For debugging
    # --annotation <bp>: coverage sums over fixed-size windows.
    # --no-annotation-stdout writes to an annotation file with the given prefix.
    megadepth \
        "${bam_file}" \
        --annotation "${win_size}" \
        --threads "${threads}" \
        --prefix "${output_prefix}" \
        --no-annotation-stdout
    set +x
    if [[ ! -s "${output_file}" ]]; then
        echo "ERROR: Output file is missing or empty: ${output_file}" >&2
        exit 1
    fi
    echo "Wrote output: ${OUT_DIR}/${output_file}"
}
export -f run_megadepth

#------------------
# Go into output directory
cd "${OUT_DIR}"

# Run coverage calculation per sample
run_megadepth "${CURR_BAM}" "${WIN_SIZE}"