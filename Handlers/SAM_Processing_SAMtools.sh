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
}

#   Export the function
export -f checkFaidx

#   A function to index the FASTA and exit
function faidxReference() {
    local reference="$1" # What is our reference FASTA file?
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
    mkdir -p "${outBase}"/Fixed_Header_SAM "${outBase}"/Raw_BAM/stats "${outBase}"/Sorted_BAM/stats "${outBase}"/Finished/stats
}

#   Export the function
export -f makeOutDirectories

#   A function to process the SAM files using SAMTools
function SAMToolsProcessing() {
    local SAMFile="$1"
    local reference="$2"
    local out="$3/SAMtools"
    local indexOpts="$4"
    #   Sample name, taken from full name of SAM file
    sampleName=$(basename "${SAMFile}" .sam)
    #   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFile}" > "${out}"/Fixed_Header_SAM/"${sampleName}"_fixed_header.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${reference}" "${out}"/Fixed_Header_SAM/"${sampleName}"_fixed_header.sam > "${out}/Raw_BAM/${sampleName}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${out}/Raw_BAM/${sampleName}_raw.bam" > "${out}/Raw_BAM/stats/${sampleName}_raw.stats"
    #   Sort the raw BAM file
    samtools sort "${out}/Raw_BAM/${sampleName}_raw.bam" > "${out}/Sorted_BAM/${sampleName}_sorted.bam"
    #   Create alignment statistics for the sorted BAM file
    samtools stats "${out}/Sorted_BAM/${sampleName}_sorted.bam" > "${out}/Sorted_BAM/stats/${sampleName}_sorted.stats"
    #   Deduplicate the sorted BAM file
    samtools rmdup "${out}/Sorted_BAM/${sampleName}_sorted.bam" "${out}/Finished/${sampleName}_finished.bam"
    #   Create alignment statistics using SAMTools
    samtools flagstat "${out}/Finished/${sampleName}_finished.bam" > "${out}/Finished/stats/${sampleName}_finished.stats"
    #   Create a CSI index for our BAM file
    samtools index "${indexOpts}" "${out}/Finished/${sampleName}_finished.bam"
}

#   Export the function
export -f SAMToolsProcessing

#   A function to run the SAM processing
function SAM_Processing() {
    local SAMList="$1" # What is our list of samples?
    local outDirectory="$2"/SAM_Processing # Where are we storing our results?
    local referenceSequence="$3" # What is our reference sequence?
    local project="$4" # What do we call our results?
    local indexType="$5"
    makeOutDirectories "${outDirectory}" # Make our outdirectories
    #   Get index options; ask if we're making BAI or CSI indecies
    [[ "${indexType}" == 'BAI' ]] && local indexOpts='-b' || local indexOpts='-c'
    parallel SAMToolsProcessing {} "${referenceSequence}" "${outDirectory}" "${indexOpts}" :::: "${SAMList}" # Process our SAM files using SAMTools
    find "${outDirectory}/SAMtools/Finished" -name "*_finished.bam" | sort > "${outDirectory}"/SAMtools/"${project}"_Finished_BAM_list.txt # Create a list of finished files
}

#   Export the function
export -f SAM_Processing
