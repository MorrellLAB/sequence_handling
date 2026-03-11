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
    local skipAdapter="$5"
    local out="${outdir}/${sampleName}"

    mkdir -p "${out}"

    echo "Running fastplong on sample: ${sampleName} (skipAdapter=${skipAdapter})"

    local -a adapterArgs=()
    if [[ "${skipAdapter}" == "true" ]]; then
        adapterArgs=(--disable_adapter_trimming)
    else
        adapterArgs=(--adapter_fasta "${adapters}")
    fi

    # fastplong parameters optimized for long reads
    fastplong \
        -i "${input}" \
        -o "${out}/${sampleName}_trimmed.fastq.gz" \
        --thread 4 \
        "${adapterArgs[@]}" \
        --html "${out}/${sampleName}_fastp.html" \
        --json "${out}/${sampleName}_fastp.json" \
        --report_title "${sampleName} fastp report"

}

export -f runFastplong

#   The main Fastplong_Handler function
function Fastplong() {
    local sampleList="$1"    # List of samples (full paths to long-read fastq files)
    local outPrefix="$2"         # Output directory (caller provides full path)
    local project="$3"       # Project name
    local adapters="$4"      # Path to adapter FASTA file
    local configPath="${5:-${CONFIG_FASTP:-}}" # Optional path to config for auto-detect

    echo "Starting Fastplong function..." >&2
    echo "sampleList = $sampleList" >&2
    echo "outPrefix = $outPrefix" >&2
    echo "project = $project" >&2
    echo "adapters = $adapters" >&2
    echo "configPath = $configPath" >&2

    if [[ ! -f "$sampleList" ]]; then
        echo "ERROR: sampleList file not found: $sampleList" >&2
        exit 1
    fi

    mkdir -p "$outPrefix"

    # Determine whether to skip adapter trimming
    local skipAdapter="false"

    if [[ "${FORCE_ADAPTER_TRIM:-}" == "true" ]]; then
        skipAdapter="false"
    elif [[ "${SKIP_ADAPTER_TRIM:-}" == "true" ]]; then
        skipAdapter="true"
    elif [[ -n "$configPath" ]] && [[ -f "$configPath" ]]; then
        local seqPlatform
        local minimapType
        seqPlatform=$(grep -E "^SEQ_PLATFORM=" "$configPath" | tail -n1 | cut -d= -f2 | tr '[:lower:]' '[:upper:]')
        minimapType=$(grep -E "^MINIMAP2_READ_TYPE=" "$configPath" | tail -n1 | cut -d= -f2 | tr '[:lower:]' '[:upper:]')

        if [[ "$seqPlatform" == "PACBIO" ]] && [[ "$minimapType" =~ HIFI ]]; then
            skipAdapter="true"
        fi
    fi

    echo "skipAdapter = $skipAdapter" >&2

    # Stream-concatenate per-sample fastq files via FIFO to avoid large temp files
    function process_sample_stream() {
        local sampleName="$1"
        shift
        local files=("$@")

        [[ ${#files[@]} -eq 0 ]] && return

        local sampleDir="${outPrefix}/${sampleName}"
        local fifoPath="${sampleDir}/${sampleName}_concat.fastq.gz"

        mkdir -p "${sampleDir}"

        echo "Streaming ${#files[@]} files for ${sampleName} into fastplong..."

        mkfifo "${fifoPath}"
        runFastplong "${sampleName}" "${fifoPath}" "${outPrefix}" "${adapters}" "${skipAdapter}" &
        local fastplong_pid=$!

        cat "${files[@]}" > "${fifoPath}"

        wait "${fastplong_pid}"
        rm -f "${fifoPath}"
    }

    # Process all samples (all long reads are single-end)
    # Supports two formats:
    #   Flat list:    one fastq.gz path per line (sample name derived from filename)
    #   Named groups: sample_name line, then one or more fastq.gz paths
    local currentSample=""
    local -a sampleFiles=()

    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" ]] && continue
        if [[ "$line" =~ \.(fastq|fq)\.gz$ ]]; then
            if [[ -z "$currentSample" ]]; then
                # Flat list: each fastq.gz is its own sample
                local sName
                sName=$(basename "$line")
                sName="${sName%.fastq.gz}"
                sName="${sName%.fq.gz}"
                process_sample_stream "$sName" "$line"
            else
                sampleFiles+=("$line")
            fi
        else
            # Named format: process previous sample, start new one
            if [[ -n "$currentSample" ]] && [[ ${#sampleFiles[@]} -gt 0 ]]; then
                process_sample_stream "$currentSample" "${sampleFiles[@]}"
            fi
            currentSample="$line"
            sampleFiles=()
        fi
    done < "${sampleList}"

    # Process the last named sample (named format only)
    if [[ -n "$currentSample" ]] && [[ ${#sampleFiles[@]} -gt 0 ]]; then
        process_sample_stream "$currentSample" "${sampleFiles[@]}"
    fi

    find "${outPrefix}" -name "*_trimmed.fastq.gz" | sort > "${outPrefix}/${project}_fastplong_trimmed.txt"
}

