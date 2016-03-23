#!/bin/env bash

#   This script performs quality trimming
#   on a series of FASTQ samples using sickle and seqqs.
#   Please install these before use.

set -e
set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Quality_Trimming_Dependencies=(sickle seqqs Rscript)

#   A function to perform the paired-end trimming and plotting
#       Adapted from Tom Kono and Peter Morrell
function trimAutoplotPaired() {
    #   Set the arguments for trimming
    local sampleName="$1" #   Name of the sample
    local forward="$2" #  Forward file
    local reverse="$3" #  Reverse file
    local out="$4"/"${sampleName}" #  Outdirectory
    local threshold="$5" #    Threshold Value
    local encoding="$6" # Platform for sequencing
    local seqHand="$7" #  The sequence_handling directory
    if [[ -d "${seqHand}"/HelperScripts ]] #   Check to see if helper scripts directory exists
    then
        local helper="${seqHand}"/HelperScripts #   The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts!"
        exit 1
    fi
    #   Make the out directory
    local stats="${out}"/stats
    local plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Run seqqs on the raw samples
    seqqs -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R1 "${forward}"
    seqqs -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R2 "${reverse}"
    #   Trim the sequences based on quality
    sickle pe -t "${encoding}" -q "${threshold}" \
        -f "${forward}" \
        -r "${reverse}" \
        -o "${out}"/"${sampleName}"_R1_trimmed.fastq \
        -p "${out}"/"${sampleName}"_R2_trimmed.fastq \
        -s "${out}"/"${sampleName}"_singles_trimmed.fastq
    #   Run seqqs on the trimmed samples
    seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R1 "${out}"/"${sampleName}"_R1_trimmed.fastq
    seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R2 "${out}"/"${sampleName}"_R2_trimmed.fastq
    #   Gzip the trimmed files
    gzip "${out}"/"${sampleName}"_R1_trimmed.fastq
    gzip "${out}"/"${sampleName}"_R2_trimmed.fastq
    gzip "${out}"/"${sampleName}"_singles_trimmed.fastq
    #   Fix the quality scores
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R1_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R2_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R1_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R2_qual.txt
    #   Make the forward plots
    Rscript "${helper}"/plot_seqqs.R \
    "${stats}/raw_${sampleName}_R1_nucl.txt" "${stats}/raw_${sampleName}_R1_len.txt" "${stats}/raw_${sampleName}_R1_qual.txt_adj" \
    "${stats}/trimmed_${sampleName}_R1_nucl.txt" "${stats}/trimmed_${sampleName}_R1_len.txt" "${stats}/trimmed_${sampleName}_R1_qual.txt_adj" \
    "${sampleName}" "forward"
    #   Make the reverse plots
    Rscript "${helper}"/plot_seqqs.R \
    "${stats}/raw_${sampleName}_R2_nucl.txt" "${stats}/raw_${sampleName}_R2_len.txt" "${stats}/raw_${sampleName}_R2_qual.txt_adj" \
    "${stats}/trimmed_${sampleName}_R2_nucl.txt" "${stats}/trimmed_${sampleName}_R2_len.txt" "${stats}/trimmed_${sampleName}_R2_qual.txt_adj" \
    "${sampleName}" "reverse"
}

#   Export the function
export -f trimAutoplotPaired

#   A function to perform the single-end trimming and plotting
#       Adapted from Tom Kono and Peter Morrell
function trimAutoplotSingle() {
    #   Set the arguments for trimming
    local sampleName="$1" #   Name of the sample
    local single="$2" #  Single file
    local out="$4"/"${sampleName}" #  Outdirectory
    local threshold="$5" #    Threshold Value
    local encoding="$6" # Platform for sequencing
    local seqHand="$7" #  The sequence_handling directory
    if [[ -d "${seqHand}"/HelperScripts ]] #   Check to see if helper scripts directory exists
    then
        local helper="${seqHand}"/HelperScripts #   The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts!"
        exit 1
    fi
    #   Make the out directories
    local stats="${out}"/stats
    local plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Run seqqs on the raw samples
    seqqs -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_single "${single}"
    #   Trim the sequences based on quality
    sickle se -t "${encoding}" -q "${threshold}" \
        -f "${single}" \
        -o "${out}"/"${sampleName}"_single_trimmed.fastq \
    #   Run seqqs on the trimmed samples
    seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_single "${out}"/"${sampleName}"_single_trimmed.fastq
    #   Gzip the trimmed files
    gzip "${out}"/"${sampleName}"_single_trimmed.fastq
    #   Fix the quality scores
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_single_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_single_qual.txt
    #   Make the single plots
    Rscript "${helper}"/plot_seqqs.R \
    "${stats}/raw_${sampleName}_single_nucl.txt" "${stats}/raw_${sampleName}_single_len.txt" "${stats}/raw_${sampleName}_single_qual.txt_adj" \
    "${stats}/trimmed_${sampleName}_single_nucl.txt" "${stats}/trimmed_${sampleName}_single_len.txt" "${stats}/trimmed_${sampleName}_single_qual.txt_adj" \
    "${sampleName}" "single"
}

#   A function to run the quality trimming
function Quality_Trimming() {
    local sampleList="$1" #   List of samples
    local forwardNaming="$2" #  Forward naming
    local reverseNaming="$3" #  Reverse naming
    local singleNaming="$4" #  Singles naming
    local outPrefix="$5"/"Quality_Trimming" #  Outdirectory
    local threshold="$6" #  Threshold Value
    local encoding="$7" #  Platform for sequencing
    local seqHand="$8" #  The sequence_handling directory
    #   Create arrays of forward and reverse samples
    local -a forwardSamples=(`grep -E "${forwardNaming}" "${sampleList}"`)
    local -a reverseSamples=(`grep -E "${reverseNaming}" "${sampleList}"`)
    local -s singleSamples=(`grep -E "${singleNaming}" "${sampleList}"`)
    #   Check to see whether we have paired-end or single samples
    if [[ ! -z "${forwardSamples[@]}" && ! -z "${reverseSamples[@]}" ]] # If we have paired-end samples
    then
        # Make sure we have equal numbers of forward and reverse samples
        if [[ "${#forwardSamples[@]}" -ne "${#reverseSamples[@]}" ]] 
            then echo "Unequal numbers of forward and reverse reads." >&2 
            exit 1
        fi 
        #   Create an array of sample names
        declare -a pairedNames=($(parallel basename {} "${forwardNaming}" ::: "${forwardSamples[@]}"))
        #   Run the paired trimmer in parallel
        parallel --xapply trimAutoplotPaired {1} {2} {3} ${outPrefix} ${threshold} ${encoding} ${seqHand} ::: ${pairedNames[@]} ::: ${forwardSamples[@]} ::: ${reverseSamples[@]}
    fi
    if ! [[ -z "${singleSamples[@]}" ]] # If we have single-end samples
    then
        #   Create an array of sample names
        declare -a singleNames=($(parallel basename {} "${singleNaming}" ::: "${singleSamples[@]}"))
        #   Run the single trimmer in parallel
        parallel --xapply trimAutoplotSingle {1} {2} ${outPrefix} ${threshold} ${encoding} ${seqHand} ::: ${singleNames[@]} ::: ${singleSamples[@]}
    fi
}

#   Export the function
export -f Quality_Trimming