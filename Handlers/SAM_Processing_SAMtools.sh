#!/bin/bash

#   This script proceses SAM files
#   including sorting and deduplicating

set -e
set -o pipefail

#   What are the dependencies for SAM_Processing
declare -a SAM_Processing_Dependencies=(samtools parallel)

#   Find our mapped samples in the outdirectory
function findMappedSAM() {
    local readMappedDirectory="$1" # Where are the mapped reads stored?
    local project="$2" # What is our project called?
    local SAMList="${readMappedDirectory}"/"${project}"_Mapped.txt # Create a name for our sample list
    find "${readMappedDirectory}" -name "*.sam" | sort > "${SAMList}" # Create our list
    echo "${SAMList}" # Return the name of our sample list
}

#   Export the function
export -f findMappedSAM

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
    local outBase="$1"
    mkdir -p "${outBase}"/fixedHeader "${outBase}"/rawBAM/stats "${outBase}"/sorted/stats "${outBase}"/finished/stats
}

#   Export the function
export -f makeOutDirectories

#   A function to process the SAM files using SAMTools
function SAMToolsProcessing() {
    #   Today's date
    local YMD=`date +%Y-%m-%d`
    local SAMFile="$1"
    local reference="$2"
    local out="$3"
    #   Sample name, taken from full name of SAM file
    sampleName=`basename "${SAMFile}" .sam`
    #   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFile}" > "${out}"/fixedHeader/"${sampleName}"_FixedHeader.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${reference}" "${out}"/fixedHeader/"${sampleName}"_FixedHeader.sam > "${out}/rawBAM/${sampleName}_${YMD}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${out}/rawBAM/${sampleName}_${YMD}_raw.bam" > "${out}/rawBAM/stats/${sampleName}_${YMD}_raw_stats.out"
    #   Sort the raw BAM file
    samtools sort "${out}/rawBAM/${sampleName}_${YMD}_raw.bam" "${out}/sorted/${sampleName}_${YMD}_sorted"
    #   Create alignment statistics for the sorted BAM file
    samtools stats "${out}/sorted/{sampleName}_${YMD}_sorted.bam" > "${out}/out/sorted/stats/${sampleName}_${YMD}_sorted_stats.out"
    #   Deduplicate the sorted BAM file
    samtools rmdup "${out}/sorted/${sampleName}_${YMD}_sorted.bam" "${out}/finished/${sampleName}_${YMD}_deduped.bam"
    #   Create alignment statistics using SAMTools
    samtools flagstat "${out}/finished/${sampleName}_${YMD}_deduped.bam" > "${out}/finished/stats/${sampleName}_${YMD}_finished_stats.out"
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
    find "${outDirectory}/finished" -name "*_deduped.bam" | sort > "${outDirectory}"/"${project}"_Processed_BAM.txt # Create a list of finished files
    # samtools merge -r "${outDirectory}"/"${project}"_Merged.bam $(find "${outDirectory}/finished" -name "*_deduped.bam") # Merge the finished BAM files into one BAM file
}

#   Export the function
export -f SAM_Processing
