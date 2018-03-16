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
    mkdir -p "${outBase}"/Statistics/Raw_SAM_Stats "${outBase}"/Statistics/Sorted_BAM_Stats "${outBase}"/Statistics/Finished_BAM_Stats "${outBase}"/Intermediates/Sorted "${outBase}"/Intermediates/Fixed_Header "${outBase}"/Intermediates/Raw_BAM
}

#   Export the function
export -f makeOutDirectories

#   A function to process the SAM files using SAMTools
function SAMToolsProcessing() {
    local SAMFile="$1"
    local reference="$2"
    local out="$3/SAMtools"
    local project="$4"
    #   Sample name, taken from full name of SAM file
    sampleName=$(basename "${SAMFile}" .sam)
    #   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFile}" > "${out}"/Intermediates/Fixed_Header/"${sampleName}"_fixed_header.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${reference}" "${out}"/Intermediates/Fixed_Header/"${sampleName}"_fixed_header.sam > "${out}/Intermediates/Raw_BAM/${sampleName}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${out}/Intermediates/Raw_BAM/${sampleName}_raw.bam" > "${out}/Statistics/Raw_SAM_Stats/${sampleName}_raw.txt"
    #   Sort the raw BAM file
    samtools sort "${out}/Intermediates/Raw_BAM/${sampleName}_raw.bam" > "${out}/Intermediates/Sorted/${sampleName}_sorted.bam"
    #   Create alignment statistics for the sorted BAM file
    samtools stats "${out}/Intermediates/Sorted/${sampleName}_sorted.bam" > "${out}/Statistics/Sorted_BAM_Stats/${sampleName}_sorted.txt"
    #   Deduplicate the sorted BAM file
    samtools rmdup "${out}/Intermediates/Sorted/${sampleName}_sorted.bam" "${out}/${sampleName}.bam"
    #   Create alignment statistics using SAMTools
    samtools flagstat "${out}/${sampleName}.bam" > "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt"
    #   Add the data from flagstat to the summary file
    local num_reads=$(head -n 1 "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt" | cut -f 1 -d " ")
    local percent_mapped=$(grep "%" "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt" | head -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local percent_paired=$(grep "%" "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt" | head -n 2 | tail -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local percent_singleton=$(grep "%" "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt" | tail -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local num_split_chr=$(tail -n 2 "${out}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt" | head -n 1 | cut -f 1 -d " ")
    local percent_split_chr=$(echo "${num_split_chr}/${num_reads}" | bc -l)
    echo -e "${sampleName}\t${num_reads}\t${percent_mapped}\t${percent_paired}\t${percent_singleton}\t${percent_split_chr}" >> "${out}/Statistics/${project}_mapping_summary_unfinished.txt"
    #   Create an index for our BAM file
    samtools index "${out}/${sampleName}.bam"
    #   Rename the index file
    mv "${out}/${sampleName}.bam.bai" "${out}/${sampleName}.bai"
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
    #   Create the header for the mapping stats summary file
    echo -e "Sample name\tTotal reads\tPercent mapped\tPercent paired\tPercent singletons\tFraction with mate mapped to different chr" > "${outDirectory}/SAMtools/Statistics/${project}_mapping_summary_unfinished.txt"
    #   Process our SAM files using SAMTools
    parallel SAMToolsProcessing {} "${referenceSequence}" "${outDirectory}" "${project}" :::: "${SAMList}"
    #   Sort the mapping stats summary file
    echo -e "Sample name\tTotal reads\tPercent mapped\tPercent paired\tPercent singletons\tFraction with mate mapped to different chr" > "${outDirectory}/SAMtools/Statistics/${project}_mapping_summary.txt"
    tail -n +2 "${outDirectory}/SAMtools/Statistics/${project}_mapping_summary_unfinished.txt" | sort >> "${outDirectory}/SAMtools/Statistics/${project}_mapping_summary.txt"
    rm "${outDirectory}/SAMtools/Statistics/${project}_mapping_summary_unfinished.txt"
    #   Create a list of finished files
    find "${outDirectory}/SAMtools" -name "*.bam" | sort > "${outDirectory}"/SAMtools/"${project}"_BAM_list.txt 
    #   Remove intermediate files
    rm -rf "${outDirectory}/SAMtools/Intermediates" 
}

#   Export the function
export -f SAM_Processing
