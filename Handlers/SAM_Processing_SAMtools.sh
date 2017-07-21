#!/bin/bash

#   This script proceses SAM files
#   including sorting and deduplicating

set -o pipefail

#   What are the dependencies for SAM_Processing
declare -a SAM_Processing_Dependencies=(samtools parallel)

#   A function to see if our reference FASTA is indexed
function checkFaidx() {
    local reference="$1" # What is our reference FASTA file?
    if ! [[ -f "${reference}" ]]; then echo "Cannot find the reference genome, exiting..." >&2; exit 1; fi # Make sure it exists
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference directory
    if ! [[ $(find "${referenceDirectory}" -maxdepth 1 -name "${referenceName}.fai") ]]; then return 1; fi # Check to make sure we have the index file for our reference FASTA file
    if ! [[ -r "${referenceDirectory}"/"${referenceName}.fai" ]]; then echo "Reference index does not have read permissions, exiting..." >&2; exit 1; fi # Make sure we can read the index
}

#   Export the function
export -f checkFaidx

#   A function to index the FASTA and exit
function faidxReference() {
    local reference="$1" # What is our reference FASTA file?
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    if ! [[ -w "${referenceDirectory}" ]]; then echo "Cannot create reference index because you do not have write permissions for ${referenceDirectory}" >&2; exit 1; fi # Make sure we can create the index files
    echo "Indexing reference for SAM Processing, will quit upon completion..." >&2
    samtools faidx "${reference}" # Index our reference FASTA file
    echo "Please re-run sequence_handling to process SAM files" >&2
    exit 10 # Exit the script with a unique exit status
}

#   Export the function
export -f faidxReference

#   A function to make our outdirectories
function makeOutDirectories() {
    local outBase="$1"/SAMtools
    #mkdir -p "${outBase}"/Fixed_Header_SAM "${outBase}"/Raw_BAM/stats "${outBase}"/Sorted_BAM/stats "${outBase}"/Finished/stats
    mkdir -p "${outBase}"/Statistics/Raw_SAM_Stats "${outBase}"/Statistics/Sorted_BAM_Stats "${outBase}"/Statistics/Finished_BAM_Stats
}

#   Export the function
export -f makeOutDirectories

#   A function to process the SAM files using SAMTools
function SAMToolsProcessing() {
    local SAMFile="$1"
    local reference="$2"
    local out="$3/SAMtools"
    #   Sample name, taken from full name of SAM file
    sampleName=$(basename "${SAMFile}" .sam)
    #   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFile}" > "${out}"/"${sampleName}"_fixed_header.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${reference}" "${out}"/"${sampleName}"_fixed_header.sam > "${out}/${sampleName}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${out}/${sampleName}_raw.bam" > "${out}/Statistics/Raw_SAM_Stats/${sampleName}_raw.txt"
    #   Sort the raw BAM file
    samtools sort "${out}/${sampleName}_raw.bam" > "${out}/${sampleName}_sorted.bam"
    #   Create alignment statistics for the sorted BAM file
    samtools stats "${out}/${sampleName}_sorted.bam" > "${out}/Statistics/Sorted_BAM_Stats/${sampleName}_sorted.txt"
    #   Deduplicate the sorted BAM file
    samtools rmdup "${out}/${sampleName}_sorted.bam" "${out}/${sampleName}.bam"
    #   Create alignment statistics using SAMTools
    samtools flagstat "${out}/${sampleName}.bam" > "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt"
    #   Create an index for our BAM file
    samtools index "${out}/${sampleName}.bam"
    #   Rename the index file
    mv "${out}/${sampleName}.bam.bai" "${out}/${sampleName}.bai"
    #   Remove intermediate files - comment out these lines if you need to troubleshoot
    rm "${out}"/"${sampleName}"_fixed_header.sam
    rm "${out}/${sampleName}_raw.bam"
    rm "${out}/${sampleName}_sorted.bam"
}

#   Export the function
export -f SAMToolsProcessing

#   A function to run the SAM processing
function SAM_Processing() {
    local SAMList="$1" # What is our list of samples?
    local outDirectory="$2"/SAM_Processing # Where are we storing our results?
    local referenceSequence="$3" # What is our reference sequence?
    local project="$4" # What do we call our results?
    makeOutDirectories "${outDirectory}" # Make our outdirectories
    parallel SAMToolsProcessing {} "${referenceSequence}" "${outDirectory}" :::: "${SAMList}" # Process our SAM files using SAMTools
    find "${outDirectory}/SAMtools" -name "*.bam" | sort > "${outDirectory}"/SAMtools/"${project}"_BAM_list.txt # Create a list of finished files
}

#   Export the function
export -f SAM_Processing
