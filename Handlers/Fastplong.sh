#!/bin/env bash

#   Handler for long-read quality filtering using fastplong
#   For ONT and PacBio reads

set -o pipefail

#   shellcheck disable=SC2034
declare -a Fastplong_Handler_Dependencies=(fastplong parallel)

function runFastplong() {
    local sampleName="$1"
    local input="$2"
    local outdir="$3"
    local adapters="$4"
    local out="${outdir}/${sampleName}"
    
    mkdir -p "${out}"

    echo "Running fastplong on: $input"
    
    # fastplong parameters optimized for long reads
    fastplong \
        -i "${input}" \
        -o "${out}/${sampleName}_trimmed.fastq.gz" \
        --thread 4 \
        --adapter_fasta "${adapters}" \
        --html "${out}/${sampleName}_fastp.html" \
        --json "${out}/${sampleName}_fastp.json" \
        --report_title "${sampleName} fastp report"

}

export -f runFastplong

#   The main Fastplong_Handler function
function Fastplong() {
    local sampleList="$1"    # List of samples (full paths to long-read fastq files)
    local outPrefix="$2"/Fastplong # Output directory
    local project="$3"       # Project name
    local adapters="$4"      # Path to adapter FASTA file

    echo "Starting Fastplong function..." >&2
    echo "sampleList = $sampleList" >&2
    echo "outPrefix = $outPrefix" >&2
    echo "project = $project" >&2
    echo "adapters = $adapters" >&2

    if [[ ! -f "$sampleList" ]]; then
        echo "ERROR: sampleList file not found: $sampleList" >&2
        exit 1
    fi

    mkdir -p "$outPrefix"

    # Process all samples (all long reads are single-end)
    while read -r sample; do
        sampleName=$(basename "${sample}" .fastq.gz)
        sampleName=$(basename "${sampleName}" .fq.gz)
        runFastplong "${sampleName}" "${sample}" "${outPrefix}" "${adapters}"
    done < "${sampleList}"

    find "${outPrefix}" -name "*_trimmed.fastq.gz" | sort > "${outPrefix}/${project}_fastplong_trimmed.txt"
}

