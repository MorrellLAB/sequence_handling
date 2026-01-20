#!/bin/env bash

#   This script maps long reads (ONT/PacBio) using minimap2
#   and outputs sorted/indexed BAM files

set -o pipefail

#   shellcheck disable=SC2034
declare -a Long_Read_Mapping_Dependencies=(minimap2 samtools)

#   A function to run read mapping using minimap2 for long reads
function Long_Read_Mapping() {
    local sample_list="$1"      # The list of paths to samples to be mapped
    local preset="$2"           # minimap2 preset: map-ont, map-hifi, or map-pb
    local project="$3"          # Project name
    local seq_platform="$4"     # Sequencing platform (ONT or PACBIO)
    local out_dir="$5"          # Output directory
    local reference="$6"        # The reference genome to map to
    local threads="$7"          # Number of threads to use
    
    local out="${out_dir}/Long_Read_Mapping"
    mkdir -p "${out}"
    
    echo "Starting Long_Read_Mapping..." >&2
    echo "sample_list = ${sample_list}" >&2
    echo "preset = ${preset}" >&2
    echo "project = ${project}" >&2
    echo "seq_platform = ${seq_platform}" >&2
    echo "out = ${out}" >&2
    echo "reference = ${reference}" >&2
    echo "threads = ${threads}" >&2
    
    # Verify reference genome exists
    if [[ ! -f "${reference}" ]]; then
        echo "ERROR: Reference genome not found: ${reference}" >&2
        exit 1
    fi
    
    # Verify sample list exists
    if [[ ! -f "${sample_list}" ]]; then
        echo "ERROR: Sample list not found: ${sample_list}" >&2
        exit 1
    fi
    
    # Get the sample array (handle various extensions)
    local -a sample_array
    mapfile -t sample_array < <(grep -E ".fastq|.fastq.gz|.fasta|.fasta.gz|.fa|.fq|.fa.gz|.fq.gz" "${sample_list}")
    
    # Get which sample in the list we are working on
    # Use PBS_ARRAYID or SLURM_ARRAY_TASK_ID if available, otherwise process all
    if [[ -n "${PBS_ARRAYID}" ]]; then
        local sample="${sample_array[${PBS_ARRAYID}]}"
    elif [[ -n "${SLURM_ARRAY_TASK_ID}" ]]; then
        local sample="${sample_array[${SLURM_ARRAY_TASK_ID}]}"
    else
        # Process all samples if not running as array job
        echo "Processing all samples sequentially..." >&2
        for sample in "${sample_array[@]}"; do
            process_sample "${sample}" "${out}" "${preset}" "${project}" "${seq_platform}" "${reference}" "${threads}"
        done
        # Create a list of mapped BAM files
        find "${out}" -name "*.bam" ! -name "*.bai" | sort > "${out}/${project}_long_reads_mapped.txt"
        return 0
    fi
    
    # Process single sample (array job)
    process_sample "${sample}" "${out}" "${preset}" "${project}" "${seq_platform}" "${reference}" "${threads}"
}

#   Helper function to process a single sample
function process_sample() {
    local sample="$1"
    local out="$2"
    local preset="$3"
    local project="$4"
    local seq_platform="$5"
    local reference="$6"
    local threads="$7"
    
    # Get the name of the sample without the path
    # If it is gzipped, this will strip off the .gz
    local sample_name
    sample_name=$(basename "${sample}" .gz)
    # Get the name without the file extension
    local base_name
    base_name="${sample_name%.*}"
    # Remove _trimmed suffix if present
    base_name="${base_name%_trimmed}"
    
    local out_bam="${out}/${base_name}_${project}.sorted.bam"
    local out_flagstat="${out}/${base_name}_${project}.flagstat.txt"
    
    echo "Processing sample: ${base_name}" >&2
    echo "Input: ${sample}" >&2
    echo "Output: ${out_bam}" >&2
    
    # Build the read group string
    local rg_string="@RG\tID:${base_name}\tSM:${base_name}\tPL:${seq_platform}"
    
    # Run minimap2 with appropriate preset and pipe to samtools for sorting
    (set -x; minimap2 \
        -ax "${preset}" \
        -t "${threads}" \
        -R "${rg_string}" \
        -I 6g \
        "${reference}" \
        "${sample}" \
        | samtools sort -@ 4 -o "${out_bam}" -)
    
    # Index the BAM file
    samtools index -@ "${threads}" "${out_bam}"
    
    # Generate alignment statistics
    samtools flagstat -@ "${threads}" "${out_bam}" > "${out_flagstat}"
    
    echo "Completed: ${out_bam}" >&2
}

export -f Long_Read_Mapping
export -f process_sample
