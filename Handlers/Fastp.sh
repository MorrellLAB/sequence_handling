#!/bin/env bash

#   This script performs quality trimming and adapter removal
#   using fastp.
#   Please install fastp before use.

set -o pipefail

#   What are the dependencies for Fastp_Handler?
declare -a Fastp_Handler_Dependencies=(fastp parallel)

#   Function to run fastp on paired-end samples
function runFastpPaired() {
    local sampleName="$1"
    local forward="$2"
    local reverse="$3"
    local outdir="$4"
    local out="${outdir}/${sampleName}"

    mkdir -p "${out}" # Create the output directory

    echo "Running fastp on: $forward and $reverse"

    fastp \
    -i "${forward}" \
    -I "${reverse}" \
    -o "${out}/${sampleName}_R1_trimmed.fastq.gz" \
    -O "${out}/${sampleName}_R2_trimmed.fastq.gz" \
    --thread 4 \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html "${out}/${sampleName}_fastp.html" \
    --json "${out}/${sampleName}_fastp.json" \
    --report_title "${sampleName} fastp report"
}

export -f runFastpPaired

#   Function to run fastp on single-end samples
function runFastpSingle() {
    local sampleName="$1"
    local single="$2"
    local out="$3"/"${sampleName}"
    mkdir -p "${out}"

    fastp \
        -i "${single}" \
        -o "${out}/${sampleName}_single_trimmed.fastq.gz" \
        --thread 4 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --html "${out}/${sampleName}_fastp.html" \
        --json "${out}/${sampleName}_fastp.json" \
        --report_title "${sampleName} fastp report"
}

export -f runFastpSingle

#   The main Fastp_Handler function
function Fastp() {
    local sampleList="$1"    # List of samples
    local forwardNaming="$2" # Forward naming
    local reverseNaming="$3" # Reverse naming
    local singleNaming="$4"  # Singles naming
    local outPrefix="$5"/Fastp # Outdirectory
    local project="$6"       # Project name

    echo "Starting Fastp function..." >&2
    echo "sampleList = $sampleList" >&2
    echo "forwardNaming = $forwardNaming" >&2
    echo "reverseNaming = $reverseNaming" >&2
    echo "singleNaming = $singleNaming" >&2
    echo "outPrefix = $outPrefix" >&2
    echo "project = $project" >&2

    if [[ ! -f "$sampleList" ]]; then
        echo "ERROR: sampleList file not found: $sampleList" >&2
        exit 1
    fi

    mkdir -p "$outPrefix"

    local -a forwardSamples=($(grep -E "${forwardNaming}" "${sampleList}"))
    local -a reverseSamples=($(grep -E "${reverseNaming}" "${sampleList}"))
    local -a singleSamples=()

    if [[ "$singleNaming" != "NONE" ]]; then
        singleSamples=($(grep -E "${singleNaming}" "${sampleList}"))
    fi

    # Paired-end trimming
    if [[ "${#forwardSamples[@]}" -ne "${#reverseSamples[@]}" ]]; then
        echo "Unequal numbers of forward and reverse reads, exiting..." >&2
        exit 1	
    fi

    if [[ "${#forwardSamples[@]}" -gt 0 ]]; then
        for ((i=0; i<${#forwardSamples[@]}; i++)); do
            sampleName=$(basename "${forwardSamples[$i]}" "$forwardNaming")
            runFastpPaired "${sampleName}" "${forwardSamples[$i]}" "${reverseSamples[$i]}" "${outPrefix}"
        done
    fi

    # Single-end trimming
    if [[ "${#singleSamples[@]}" -gt 0 ]]; then
        for singleSample in "${singleSamples[@]}"; do
            sampleName=$(basename "${singleSample}" "${singleNaming}")
            runFastpSingle "${sampleName}" "${singleSample}" "${outPrefix}"
        done
    fi

    find "${outPrefix}" -name "*.fastq.gz" | sort > "${outPrefix}/${project}_fastp_trimmed.txt"
}

