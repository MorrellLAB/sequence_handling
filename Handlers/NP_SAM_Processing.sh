#!/bin/sh

set -o pipefail

#   What are the dependencies for NP_SAM_Processing?
declare -a NP_SAM_Processing_Dependencies=(samtools)

#   A function to run SAM_Processing
function NP_SAM_Processing() {
    local sample_list="$1" # The list of paths to samples
    local outDirectory="$2/SAM_Processing" # The location to output BAM files
    local reference="$3" # The reference sequence
    local project="$4" # The name of the project
    #   Make sure the output directory exists
    mkdir -p "${out}/Statistics/Finished_BAM_Stats" 
    #   Make an array of samples
    declare -a sample_array=($(grep -E ".sam" "${sample_list}"))
    #   Get which sample in the list we are working on
    local sample="${sample_array[${PBS_ARRAYID}]}"
    #   Get the name of the sample without the path in front
    local sampleName=$(basename "${sample}" .sam)
    #   Convert to BAM format
    (set -x; samtools view -b -T "${reference}" "${sample}" > "${outDirectory}/${sampleName}.bam")
    #   Run flagstat on the finished file to assess quality
    (set -x; samtools flagstat "${outDirectory}/${sampleName}.bam" > "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt")
    #   Add the data from flagstat to the summary file
    local num_reads=$(head -n 1 "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt" | cut -f 1 -d " ")
    local percent_mapped=$(grep "%" "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt" | head -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local percent_paired=$(grep "%" "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt" | head -n 2 | tail -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local percent_singleton=$(grep "%" "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt" | tail -n 1 | cut -f 2 -d "(" | cut -f 1 -d " ")
    local num_split_chr=$(tail -n 2 "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_stats.txt" | head -n 1 | cut -f 1 -d " ")
    local percent_split_chr=$(echo "${num_split_chr}/${num_reads}" | bc -l)
    echo -e "${sampleName}\t${num_reads}\t${percent_mapped}\t${percent_paired}\t${percent_singleton}\t${percent_split_chr}" >> "${outDirectory}/Statistics/${project}_mapping_summary.txt"
    #   Index the finished BAM file
    samtools index "${outDirectory}/${sampleName}.bam"
    #   Rename the index file
    mv "${outDirectory}/${sampleName}.bam.bai" "${outDirectory}/${sampleName}.bai"
}

export -f NP_SAM_Processing