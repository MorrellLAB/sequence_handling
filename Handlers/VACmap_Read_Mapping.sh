#!/bin/env bash

#   This script maps long reads (ONT/PacBio CLR/PacBio HiFi) using VACmap,
#   a structural-variation-aware aligner, and outputs sorted/indexed BAM files.
#   VACmap: https://github.com/micahvista/VACmap

set -euo pipefail

#   shellcheck disable=SC2034
declare -a VACmap_Read_Mapping_Dependencies=(vacmap samtools)

#   Infer VACmap alignment mode from sequencing platform when VACMAP_MODE is not set.
#   Modes: H (high-error ONT/CLR), L (low-error HiFi), S (small-variant sensitivity),
#          R (translocation detection), asm (assembly alignment).
function resolve_vacmap_mode() {
    local seq_platform="$1"
    local explicit_mode="${VACMAP_MODE:-}"

    if [[ -n "${explicit_mode}" ]]; then
        printf '%s\n' "${explicit_mode}"
        return 0
    fi

    case "${seq_platform}" in
        ONT)         printf 'H\n' ;;
        PACBIO)      printf 'H\n' ;;
        PACBIO_HIFI) printf 'L\n' ;;
        *)
            echo "[ERROR] Cannot infer VACmap mode from SEQ_PLATFORM='${seq_platform}'." \
                 "Set VACMAP_MODE explicitly (H, L, S, R, or asm)." >&2
            return 1
            ;;
    esac
}

export -f resolve_vacmap_mode

#   A function to map long reads using VACmap
function VACmap_Read_Mapping() {
    local sample_list="$1"   # Full path to a text file listing read files (one per line)
    local mode="$2"          # VACmap alignment mode: H, L, S, R, or asm
    local project="$3"       # Project name used for output naming
    local seq_platform="$4"  # Sequencing platform (ONT, PACBIO, PACBIO_HIFI)
    local out_dir="$5"       # Parent output directory
    local reference="$6"     # Reference genome FASTA path
    local threads="$7"       # Number of threads
    local kmer="${8:-15}"    # k-mer size (default: 15; HiFi recommended: 19)
    local window="${9:-10}"  # Minimizer window size (default: 10)
    local workdir="${10:-}"  # Temp dir required only for mode=asm

    local out="${out_dir}/VACmap_Read_Mapping"
    mkdir -p "${out}"

    echo "[VACmap] Starting VACmap_Read_Mapping" >&2
    echo "[VACmap] sample_list  = ${sample_list}" >&2
    echo "[VACmap] mode         = ${mode}" >&2
    echo "[VACmap] project      = ${project}" >&2
    echo "[VACmap] seq_platform = ${seq_platform}" >&2
    echo "[VACmap] reference    = ${reference}" >&2
    echo "[VACmap] threads      = ${threads}" >&2
    echo "[VACmap] kmer         = ${kmer}" >&2
    echo "[VACmap] window       = ${window}" >&2
    echo "[VACmap] out          = ${out}" >&2
    echo "[VACmap] VACmap version: $(vacmap --version 2>&1 || true)" >&2
    echo "[VACmap] samtools version: $(samtools --version | head -1)" >&2

    # Validate inputs
    if [[ ! -f "${reference}" ]]; then
        echo "[ERROR] Reference genome not found: ${reference}" >&2
        exit 1
    fi

    if [[ ! -f "${sample_list}" ]]; then
        echo "[ERROR] Sample list not found: ${sample_list}" >&2
        exit 1
    fi

    # Validate alignment mode
    case "${mode}" in
        H|L|S|R|asm) ;;
        *)
            echo "[ERROR] Invalid VACmap mode '${mode}'. Must be one of: H L S R asm" >&2
            exit 1
            ;;
    esac

    # Assembly mode requires a working directory
    if [[ "${mode}" == "asm" && -z "${workdir}" ]]; then
        echo "[ERROR] VACmap assembly mode (asm) requires a working directory (VACMAP_WORKDIR)." >&2
        exit 1
    fi

    # Load sample array (FASTA/FASTQ, plain or gzipped, and BAM)
    local -a sample_array
    mapfile -t sample_array < <(grep -E '\.(fastq|fastq\.gz|fasta|fasta\.gz|fa|fq|fa\.gz|fq\.gz|bam)$' "${sample_list}")

    if [[ ${#sample_array[@]} -eq 0 ]]; then
        echo "[ERROR] No valid read files found in sample list: ${sample_list}" >&2
        exit 1
    fi

    # Determine which sample to process (Slurm array vs. sequential)
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local sample="${sample_array[${SLURM_ARRAY_TASK_ID}]}"
        _vacmap_process_sample \
            "${sample}" "${out}" "${mode}" "${project}" "${seq_platform}" \
            "${reference}" "${threads}" "${kmer}" "${window}" "${workdir}"
    elif [[ -n "${PBS_ARRAYID:-}" ]]; then
        local sample="${sample_array[${PBS_ARRAYID}]}"
        _vacmap_process_sample \
            "${sample}" "${out}" "${mode}" "${project}" "${seq_platform}" \
            "${reference}" "${threads}" "${kmer}" "${window}" "${workdir}"
    else
        echo "[VACmap] No array scheduler detected; processing all ${#sample_array[@]} samples sequentially." >&2
        for sample in "${sample_array[@]}"; do
            _vacmap_process_sample \
                "${sample}" "${out}" "${mode}" "${project}" "${seq_platform}" \
                "${reference}" "${threads}" "${kmer}" "${window}" "${workdir}"
        done
    fi

    # Emit a manifest of output BAM files for downstream handlers
    find "${out}" -maxdepth 1 -name "*.sorted.bam" ! -name "*.bai" | sort \
        > "${out}/${project}_vacmap_mapped.txt"
    echo "[VACmap] BAM manifest written: ${out}/${project}_vacmap_mapped.txt" >&2
}

#   Process a single read file through VACmap
function _vacmap_process_sample() {
    local sample="$1"
    local out="$2"
    local mode="$3"
    local project="$4"
    local seq_platform="$5"
    local reference="$6"
    local threads="$7"
    local kmer="$8"
    local window="$9"
    local workdir="${10:-}"

    local temp_dir
    temp_dir=$(mktemp -d)
    # shellcheck disable=SC2064
    trap "rm -rf '${temp_dir}'" EXIT

    if [[ ! -f "${sample}" ]]; then
        echo "[ERROR] Read file not found: ${sample}" >&2
        exit 1
    fi

    # Derive a clean sample name (strip path, .gz, and one extension layer)
    local sample_name
    sample_name=$(basename "${sample}" .gz)
    local base_name="${sample_name%.*}"
    # Strip common trimming suffixes added by Fastplong/fastp
    base_name="${base_name%_trimmed}"

    local out_bam="${out}/${base_name}_${project}.sorted.bam"
    local out_flagstat="${out}/${base_name}_${project}.flagstat.txt"
    local sample_log="${out}/${base_name}_${project}.vacmap.log"

    echo "[VACmap] Sample: ${base_name} | Input: ${sample} | Output: ${out_bam}" \
        | tee -a "${sample_log}" >&2

    # Build the vacmap invocation
    local -a vacmap_cmd=(
        vacmap
        -ref   "${reference}"
        -read  "${sample}"
        -mode  "${mode}"
        -t     "${threads}"
        -k     "${kmer}"
        -w     "${window}"
        -o     "${out_bam}"
        --rg-id "${base_name}"
        --rg-sm "${base_name}"
        --rg-pl "${seq_platform}"
    )

    # Assembly mode requires an explicit working directory
    if [[ "${mode}" == "asm" ]]; then
        mkdir -p "${workdir}"
        vacmap_cmd+=(-workdir "${workdir}")
    fi

    echo "[VACmap] Command: ${vacmap_cmd[*]}" | tee -a "${sample_log}" >&2

    # Run VACmap; exit with informative message on failure
    "${vacmap_cmd[@]}" >> "${sample_log}" 2>&1 || {
        echo "[ERROR] VACmap failed for sample ${base_name}. See log: ${sample_log}" >&2
        exit 1
    }

    # Index the BAM file (VACmap outputs sorted BAM when -o ends in .sorted.bam)
    samtools index -@ "${threads}" "${out_bam}" || {
        echo "[ERROR] Failed to index ${out_bam}" >&2
        exit 1
    }

    # Alignment statistics
    samtools flagstat -@ "${threads}" "${out_bam}" > "${out_flagstat}" || {
        echo "[ERROR] samtools flagstat failed for ${out_bam}" >&2
        exit 1
    }

    echo "[VACmap] Completed: ${base_name}" | tee -a "${sample_log}" >&2
    echo "[VACmap] Flagstat summary:" >&2
    cat "${out_flagstat}" >&2
}

export -f VACmap_Read_Mapping
export -f _vacmap_process_sample
