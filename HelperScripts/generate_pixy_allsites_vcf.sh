#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=2
#SBATCH --mem=32g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Generate all-sites VCFs suitable for pixy, one output per chromosome.
# Chromosome names are discovered from BAM index stats.

set -euo pipefail

# Load tools if available via modules
if command -v module >/dev/null 2>&1; then
    if ! command -v bcftools >/dev/null 2>&1; then
        module load bcftools/1.21 2>/dev/null || module load bcftools
    fi
    if ! command -v samtools >/dev/null 2>&1; then
        module load samtools/1.21 2>/dev/null || module load samtools
    fi
fi

function usage() {
    cat <<'EOF'
Usage:
    generate_pixy_allsites_vcf.sh --method bcftools \
    --bam-list BAM_LIST \
    --reference REF_FASTA \
    --output-dir OUT_DIR \
    [--output-prefix PREFIX] \
    [--threads 2]
    
Notes:
    - Creates one VCF per chromosome found in BAM index stats.
    - --reference expects a FASTA file (for example, ref.fa), not a .fai path.
    - Sample names are normalized by default (for example, SAMPLE_1.fastq.gz -> SAMPLE).
    - Requires BAM index files (.bai) to be present.
EOF
}

function require_cmd() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    }
}

function ensure_readable_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || {
        echo "[ERROR] ${label} not found or empty: ${path}" >&2
        exit 1
    }
}

function ensure_readable_lines() {
    local list_file="$1"
    local label="$2"
    ensure_readable_file "${list_file}" "${label} list"
    while IFS= read -r line || [[ -n "${line}" ]]; do
        [[ -z "${line}" ]] && continue
        [[ -s "${line}" ]] || {
            echo "[ERROR] ${label} path from list is not readable: ${line}" >&2
            exit 1
        }
    done < "${list_file}"
}

function index_vcf_if_needed() {
    local vcf_path="$1"
    local threads_arg="$2"
    if [[ -e "${vcf_path}.tbi" || -e "${vcf_path}.csi" ]]; then
        return 0
    fi
    bcftools index --threads "${threads_arg}" -t "${vcf_path}"
}

function get_chromosomes_from_bam_list() {
    local list_file="$1"
    while IFS= read -r bam_path || [[ -n "${bam_path}" ]]; do
        [[ -z "${bam_path}" ]] && continue
        samtools idxstats "${bam_path}" | awk '$1 != "*" {print $1}'
    done < "${list_file}" | awk '!seen[$1]++'
}

function normalize_sample_names_in_vcf() {
    local vcf_path="$1"
    local tmp_dir="$2"
    local threads_arg="$3"
    local old_samples_file="${tmp_dir}/samples.old.txt"
    local map_file="${tmp_dir}/samples.map.txt"
    local dup_file="${tmp_dir}/samples.duplicate.txt"
    local out_vcf="${tmp_dir}/reheadered.$(basename "${vcf_path}")"

    bcftools query -l "${vcf_path}" > "${old_samples_file}"
    awk 'BEGIN{OFS="\t"} {
        old=$0
        new=$0
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", new)
        sub(/\.fastq\.gz$/, "", new)
        sub(/\.fq\.gz$/, "", new)
        sub(/_R?[12]$/, "", new)
        print old, new
    }' "${old_samples_file}" > "${map_file}"

    cut -f2 "${map_file}" | sort | uniq -d > "${dup_file}"
    if [[ -s "${dup_file}" ]]; then
        echo "[ERROR] Normalized sample names are not unique for ${vcf_path}" >&2
        echo "[ERROR] Duplicate normalized names (first 10):" >&2
        head -n 10 "${dup_file}" >&2
        exit 1
    fi

    if awk -F '\t' '$1 != $2 {changed=1} END {exit changed ? 0 : 1}' "${map_file}"; then
        echo "[pixy-allsites] Normalizing sample names in ${vcf_path}"
        bcftools reheader -s "${map_file}" -o "${out_vcf}" "${vcf_path}"
        mv -f "${out_vcf}" "${vcf_path}"
        bcftools index --threads "${threads_arg}" -f -t "${vcf_path}"
    fi
}

# --- Default Variables ---
method="bcftools"
threads=2
min_mapq=20
min_baseq=20
min_qual=20
bam_list=""
reference=""
output_dir=""
output_prefix="pixy_allsites"

# --- Parse Arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-list) bam_list="$2"; shift 2 ;;
        --reference) reference="$2"; shift 2 ;;
        --output-dir) output_dir="$2"; shift 2 ;;
        --output-prefix) output_prefix="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --method) method="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1"; usage; exit 1 ;;
    esac
done

# --- Validation ---
[[ -n "${reference}" ]] || { echo "[ERROR] --reference is required"; exit 1; }
[[ -n "${output_dir}" ]] || { echo "[ERROR] --output-dir is required"; exit 1; }
ensure_readable_file "${reference}" "Reference FASTA"
[[ -s "${reference}.fai" ]] || {
    echo "[ERROR] Missing FASTA index: ${reference}.fai" >&2
    echo "[ERROR] Create it with: samtools faidx ${reference}" >&2
    exit 1
}
mkdir -p "${output_dir}"

# --- Execute bcftools Pipeline ---
if [[ "${method}" == "bcftools" ]]; then
    [[ -n "${bam_list}" ]] || { echo "[ERROR] --bam-list is required"; exit 1; }
    ensure_readable_lines "${bam_list}" "BAM"
    require_cmd bcftools
    require_cmd samtools

    mapfile -t chromosomes < <(get_chromosomes_from_bam_list "${bam_list}")
    [[ ${#chromosomes[@]} -gt 0 ]] || {
        echo "[ERROR] No chromosome names found from BAM list: ${bam_list}" >&2
        exit 1
    }

    echo "[pixy-allsites] Chromosome source: BAM index stats from ${bam_list}"
    echo "[pixy-allsites] Chromosomes to process: ${#chromosomes[@]}"
    echo "[pixy-allsites] Reference: ${reference}"
    echo "[pixy-allsites] Output directory: ${output_dir}"

    tmp_root="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/pixy_allsites_${SLURM_JOB_ID:-$$}"
    mkdir -p "${tmp_root}"
    trap 'rm -rf "${tmp_root}"' EXIT

    for chrom in "${chromosomes[@]}"; do
        output="${output_dir}/${output_prefix}.${chrom}.vcf.gz"
        echo "[pixy-allsites] Processing ${chrom} -> ${output}"

        bcftools mpileup \
            --threads "${threads}" \
            -r "${chrom}" \
            -f "${reference}" \
            -b "${bam_list}" \
            -a AD,DP \
            -q "${min_mapq}" \
            -Q "${min_baseq}" \
            -Ou | \
        bcftools call \
            --threads "${threads}" \
            -m \
            -A \
            -Ou | \
        bcftools filter \
            --threads "${threads}" \
            -i "QUAL>=${min_qual}" \
            -Ou | \
        bcftools norm \
            --threads "${threads}" \
            -f "${reference}" \
            -Oz \
            -o "${output}"

        index_vcf_if_needed "${output}" "${threads}"
        normalize_sample_names_in_vcf "${output}" "${tmp_root}" "${threads}"
    done
else
    echo "[ERROR] This simplified script only supports --method bcftools"
    exit 1
fi

echo "[pixy-allsites] Completed processing all chromosomes from BAM list: ${bam_list}"