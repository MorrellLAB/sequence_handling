#!/bin/bash

#   This script performs long-read processing using fastplong
#   which removes adapter sequences and filters by quality
#   Please install fastplong before use.

set -o pipefail

#   What are the dependencies for Fastplong?
declare -a Fastplong_Dependencies=(fastplong parallel)

#   A function to detect if we should skip adapter trimming for PacBio HiFi reads
function detect_pacbio_hifi() {
    local sample_file="$1"
    #   Check if this is a PacBio HiFi file (typically indicated by filename patterns or read structure)
    #   PacBio HiFi reads are typically labeled with 'hifi' or 'ccs' in the filename
    if [[ "${sample_file}" =~ hifi|ccs|HiFi|CCS ]]; then
        echo "true"
    else
        echo "false"
    fi
}

#   A function to process a sample using direct file concatenation
function process_sample_direct() {
    local sampleName="$1" # Name of the sample
    local sampleFiles="$2" # Comma-separated list of sample files
    local outDirectory="$3" # Output directory
    local skipAdapter="$4" # Whether to skip adapter trimming (for PacBio HiFi)
    local seqHand="$5" # The sequence_handling directory
    
    #   Make the output directories
    mkdir -p "${outDirectory}/${sampleName}"
    local out="${outDirectory}/${sampleName}"
    local stats="${out}/stats"
    mkdir -p "${stats}"
    
    #   Convert comma-separated list to array
    IFS=',' read -ra fileArray <<< "${sampleFiles}"
    
    #   Create a temporary concatenated file
    local tmpConcatFile="${out}/${sampleName}_concat_temp.fastq"
    
    #   Concatenate all input files
    echo "Processing sample: ${sampleName}"
    echo "Concatenating ${#fileArray[@]} file(s)..."
    
    #   Clear any existing temp file
    true > "${tmpConcatFile}"
    
    #   Concatenate files, handling different compression formats
    for file in "${fileArray[@]}"; do
        if [[ ! -f "${file}" ]]; then
            echo "Error: Input file not found: ${file}" >&2
            exit 1
        fi
        
        if [[ "${file}" =~ \.gz$ ]]; then
            gzip -cd "${file}" >> "${tmpConcatFile}"
        elif [[ "${file}" =~ \.bz2$ ]]; then
            bzip2 -cd "${file}" >> "${tmpConcatFile}"
        else
            cat "${file}" >> "${tmpConcatFile}"
        fi
    done
    
    #   Verify concatenated file was created
    if [[ ! -s "${tmpConcatFile}" ]]; then
        echo "Error: Failed to create concatenated file for ${sampleName}" >&2
        rm -f "${tmpConcatFile}"
        exit 1
    fi
    
    #   Prepare output file names
    local outputFile="${out}/${sampleName}_filtered.fastq.gz"
    
    #   Run fastplong with appropriate options
    echo "Running fastplong on ${sampleName}..."
    
    if [[ "${skipAdapter}" == "true" ]]; then
        #   Skip adapter trimming for PacBio HiFi reads
        if ! fastplong --skip-adapters --input "${tmpConcatFile}" --output "${outputFile}"; then
            echo "Error: fastplong failed for sample ${sampleName}" >&2
            rm -f "${tmpConcatFile}"
            exit 1
        fi
    else
        #   Standard processing with adapter trimming
        if ! fastplong --input "${tmpConcatFile}" --output "${outputFile}"; then
            echo "Error: fastplong failed for sample ${sampleName}" >&2
            rm -f "${tmpConcatFile}"
            exit 1
        fi
    fi
    
    #   Verify output file was created
    if [[ ! -f "${outputFile}" ]]; then
        echo "Error: Output file not created for ${sampleName}: ${outputFile}" >&2
        rm -f "${tmpConcatFile}"
        exit 1
    fi
    
    #   Clean up temporary concatenated file
    rm -f "${tmpConcatFile}"
    
    echo "Successfully processed ${sampleName}"
}

#   Export the function
export -f process_sample_direct
export -f detect_pacbio_hifi

#   Main handler function for Fastplong processing
function Fastplong() {
    local sampleList="$1" # List of samples (one sample per line, with files separated by commas)
    local outPrefix="$2" # Output directory prefix
    local project="$3" # Project name
    local seqHand="$4" # The sequence_handling directory
    
    #   Create the output directory
    local outDirectory="${outPrefix}/Fastplong"
    mkdir -p "${outDirectory}"
    
    #   Check if helper scripts directory exists
    if [[ ! -d "${seqHand}"/HelperScripts ]]; then
        echo "Cannot find directory with helper scripts, exiting..." >&2
        exit 1
    fi
    
    #   Read the sample list
    if [[ ! -f "${sampleList}" ]]; then
        echo "Error: Sample list file not found: ${sampleList}" >&2
        exit 1
    fi
    
    #   Process each sample
    declare -a sampleNames=()
    declare -a sampleFileLists=()
    declare -a skipAdapterFlags=()
    
    #   Parse the sample list
    while IFS= read -r line; do
        #   Skip empty lines and comments
        [[ -z "${line}" || "${line}" =~ ^# ]] && continue
        
        #   Extract sample name (first file basename without extension)
        local firstFile
        local sampleName
        local skipAdapter
        firstFile=$(echo "${line}" | cut -d',' -f1)
        #   Remove common sequencing file extensions (.fastq.gz, .fq.gz, .fastq, .fq, etc.)
        sampleName=$(basename "${firstFile}" | sed -e 's/\.fastq\.gz$//' -e 's/\.fq\.gz$//' -e 's/\.fastq\.bz2$//' -e 's/\.fq\.bz2$//' -e 's/\.fastq$//' -e 's/\.fq$//')
        
        #   Detect if we should skip adapter trimming
        skipAdapter=$(detect_pacbio_hifi "${firstFile}")
        
        sampleNames+=("${sampleName}")
        sampleFileLists+=("${line}")
        skipAdapterFlags+=("${skipAdapter}")
    done < "${sampleList}"
    
    #   Check if we have any samples to process
    if [[ ${#sampleNames[@]} -eq 0 ]]; then
        echo "Error: No samples found in sample list" >&2
        exit 1
    fi
    
    #   Process samples in parallel
    echo "Processing ${#sampleNames[@]} sample(s)..."
    if ! parallel --verbose --xapply process_sample_direct {1} {2} "${outDirectory}" {3} "${seqHand}" \
        ::: "${sampleNames[@]}" \
        ::: "${sampleFileLists[@]}" \
        ::: "${skipAdapterFlags[@]}"; then
        echo "Error: Parallel processing failed" >&2
        exit 1
    fi
    
    #   Create a list of output files
    find "${outDirectory}" -name "*_filtered.fastq.gz" -type f | sort > "${outDirectory}/${project}_fastplong_filtered.txt"
    
    #   Verify the output list was created
    if [[ ! -f "${outDirectory}/${project}_fastplong_filtered.txt" ]]; then
        echo "Error: Failed to create output file list" >&2
        exit 1
    fi
    
    echo "Fastplong processing complete. Output list: ${outDirectory}/${project}_fastplong_filtered.txt"
}

#   Export the function
export -f Fastplong
