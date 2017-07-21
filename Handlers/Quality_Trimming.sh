#!/bin/env bash

#   This script performs quality trimming
#   on a series of FASTQ samples using sickle and seqqs.
#   Please install these before use.

set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Quality_Trimming_Dependencies=(sickle seqqs Rscript parallel)

#   A function to perform the paired-end trimming and plotting
#       Adapted from Tom Kono and Peter Morrell
function trimAutoplotPaired() {
    local sampleName="$1" # Name of the sample
    local forward="$2" # Forward file
    local reverse="$3" # Reverse file
    local out="$4"/"${sampleName}" # Outdirectory
    local threshold="$5" # Threshold Value
    local encoding="$6" # Platform for sequencing
    local seqHand="$7" # The sequence_handling directory
    if [[ -d "${seqHand}"/HelperScripts ]] # Check to see if helper scripts directory exists
    then
        local helper="${seqHand}"/HelperScripts # The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts, exiting..."
        exit 1
    fi
    #   Make the out directory
    local stats="${out}"/stats
    local plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Check compression of forward sample
    if [[ $( echo "${forward}" | rev | cut -f 1 -d '.' | rev) == 'gz' ]] # If gzipped...
    then
        echo "gzip"
        local forwardTrim="${out}/${sampleName}_Forward_PIPE" # Create a name for the pipe
        rm -f "${forwardTrim}" # Remove any existing pipes
        mkfifo "${forwardTrim}" # Make the pipe
        gzip -cd "${forward}" | seqqs -e -p "${stats}/raw_${sampleName}_R1" - > "${forwardTrim}" & # Uncompress the file, run seqqs, and write to pipe
    elif [[ $( echo "${forward}" | rev | cut -f 1 -d '.' | rev) == 'bz2' ]] # If bzipped...
    then
        echo 'bzip'
        local forwardTrim="${out}/${sampleName}_Forward_PIPE" # Create a name for the pipe
        rm -f "${forwardTrim}" # Remove any existing pipes
        mkfifo "${forwardTrim}" # Make the pipe
        bzip2 -cd "${forward}" | seqqs -e -p "${stats}/raw_${sampleName}_R1" - > "${forwardTrim}" & # Uncompress the file, run seqqs, and write to pipe
    else # Otherwise
        echo 'nope'
        local forwardTrim="${forward}" # Use the name of the sample as 'toTrim'
        seqqs -p "${stats}/raw_${sampleName}_R1" "${forwardTrim}" # Run seqqs
    fi
    #   Check compression on the reverse samples
    if [[ $( echo "${reverse}" | rev | cut -f 1 -d '.' | rev) == 'gz' ]] # If gzipped...
    then
        echo "gzip"
        local reverseTrim="${out}/${sampleName}_Reverse_PIPE" # Create a name for the pipe
        rm -f "${reverseTrim}" # Remove any existing pipes
        mkfifo "${reverseTrim}" # Make the pipe
        gzip -cd "${reverse}" | seqqs -e -p "${stats}/raw_${sampleName}_R2" - > "${reverseTrim}" & # Uncompress the file, run seqqs, and write to pipe
    elif [[ $( echo "${reverse}" | rev | cut -f 1 -d '.' | rev) == 'bz2' ]] # If bzipped...
    then
        echo 'bzip'
        local reverseTrim="${out}/${sampleName}_Reverse_PIPE" # Create a name for the pipe
        rm -f "${reverseTrim}" # Remove any existing pipes
        mkfifo "${reverseTrim}" # Make the pipe
        bzip2 -cd "${reverse}" | seqqs -e -p "${stats}/raw_${sampleName}_R2" - > "${reverseTrim}" & # Uncompress the file, run seqqs, and write to pipe
    else # Otherwise
        echo 'nope'
        local reverseTrim="${reverse}" # Use the name of the sample as 'toTrim'
        seqqs -p "${stats}/raw_${sampleName}_R2" "${toTrim}" # Run seqqs
    fi
    #   Trim the sequences based on quality
    sickle pe -t "${encoding}" -q "${threshold}" --gzip-output \
        -f "${forwardTrim}" \
        -r "${reverseTrim}" \
        -o "${out}"/"${sampleName}"_R1_trimmed.fastq.gz \
        -p "${out}"/"${sampleName}"_R2_trimmed.fastq.gz \
        -s "${out}"/"${sampleName}"_singles_trimmed.fastq.gz
    #   Run seqqs on the trimmed samples
    gzip -cd "${out}"/"${sampleName}"_R1_trimmed.fastq.gz | seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R1 -
    gzip -cd "${out}"/"${sampleName}"_R2_trimmed.fastq.gz | seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R2 -
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
    local sampleName="$1" # Name of the sample
    local single="$2" # Single file
    local out="$3"/"${sampleName}" # Outdirectory
    local threshold="$4" # Threshold Value
    local encoding="$5" # Platform for sequencing
    local seqHand="$6" # The sequence_handling directory
    if [[ -d "${seqHand}"/HelperScripts ]] # Check to see if helper scripts directory exists
    then
        local helper="${seqHand}"/HelperScripts # The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts, exiting..."
        exit 1
    fi
    #   Make the out directories
    local stats="${out}"/stats
    local plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Check compression type and run seqqs on raw samples
    if [[ $( echo "${single}" | rev | cut -f 1 -d '.' | rev) == 'gz' ]] # If gzipped...
    then
        echo "gzip"
        local toTrim="${out}/${sampleName}_single_PIPE" # Create a name for the pipe
        rm -f "${toTrim}" # Remove any existing pipes
        mkfifo "${toTrim}" # Make the pipe
        gzip -cd "${single}" | seqqs -e -p "${stats}/raw_${sampleName}_single" - > "${toTrim}" & # Uncompress the file, run seqqs, and write to pipe
    elif [[ $( echo "${single}" | rev | cut -f 1 -d '.' | rev) == 'bz2' ]] # If bzipped...
    then
        echo 'bzip'
        local toTrim="${out}/${sampleName}_single_PIPE" # Create a name for the pipe
        rm -f "${toTrim}" # Remove any existing pipes
        mkfifo "${toTrim}" # Make the pipe
        bzip2 -cd "${single}" | seqqs -e -p "${stats}/raw_${sampleName}_single" - > "${toTrim}" & # Uncompress the file, run seqqs, and write to pipe
    else # Otherwise
        echo 'nope'
        local toTrim="${single}" # Use the name of the sample as 'toTrim'
        seqqs -p "${stats}/raw_${sampleName}_single" "${toTrim}" # Run seqqs
    fi
    #   Trim the sequences based on quality
    sickle se -t "${encoding}" -q "${threshold}" --gzip-output \
        -f "${toTrim}" \
        -o "${out}"/"${sampleName}"_single_trimmed.fastq.gz \
    #   Run seqqs on the trimmed samples
    gzip -cd "${out}"/"${sampleName}"_single_trimmed.fastq.gz | seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_single -
    #   Fix the quality scores
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_single_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_single_qual.txt
    #   Make the single plots
    Rscript "${helper}"/plot_seqqs.R \
    "${stats}/raw_${sampleName}_single_nucl.txt" "${stats}/raw_${sampleName}_single_len.txt" "${stats}/raw_${sampleName}_single_qual.txt_adj" \
    "${stats}/trimmed_${sampleName}_single_nucl.txt" "${stats}/trimmed_${sampleName}_single_len.txt" "${stats}/trimmed_${sampleName}_single_qual.txt_adj" \
    "${sampleName}" "single"
}

#   Export the function
export -f trimAutoplotSingle

#   A function to run the quality trimming
function Quality_Trimming() {
    local sampleList="$1" # List of samples
    local forwardNaming="$2" # Forward naming
    local reverseNaming="$3" # Reverse naming
    local singleNaming="$4" # Singles naming
    local outPrefix="$5"/"Quality_Trimming" # Outdirectory
    local threshold="$6" # Threshold Value
    local encoding="$7" # Platform for sequencing
    local seqHand="$8" # The sequence_handling directory
    local project="$9" # The name of our project
    #   Create arrays of forward, reverse, and single samples
    local -a forwardSamples=(`grep -E "${forwardNaming}" "${sampleList}"`)
    local -a reverseSamples=(`grep -E "${reverseNaming}" "${sampleList}"`)
    local -a singleSamples=(`grep -E "${singleNaming}" "${sampleList}"`)
    #   Check to see whether we have paired-end or single samples
    if [[ ! -z "${forwardSamples[@]}" && ! -z "${reverseSamples[@]}" ]] # If we have paired-end samples
    then
        # Make sure we have equal numbers of forward and reverse samples
        if [[ "${#forwardSamples[@]}" -ne "${#reverseSamples[@]}" ]] 
            then echo "Unequal numbers of forward and reverse reads, exiting..." >&2 
            exit 1
        fi 
        #   Create an array of sample names
        declare -a pairedNames=($(parallel --verbose basename {} "${forwardNaming}" ::: "${forwardSamples[@]}"))
        #   Run the paired trimmer in parallel
        parallel --verbose --xapply trimAutoplotPaired {1} {2} {3} ${outPrefix} ${threshold} ${encoding} ${seqHand} ::: ${pairedNames[@]} ::: ${forwardSamples[@]} ::: ${reverseSamples[@]}
    fi
    if ! [[ -z "${singleSamples[@]}" ]] # If we have single-end samples
    then
        #   Create an array of sample names
        declare -a singleNames=($(parallel basename {} "${singleNaming}" ::: "${singleSamples[@]}"))
        #   Run the single trimmer in parallel
        parallel --verbose --xapply trimAutoplotSingle {1} {2} ${outPrefix} ${threshold} ${encoding} ${seqHand} ::: ${singleNames[@]} ::: ${singleSamples[@]}
    fi
    find "${outPrefix}" -type p -exec rm {} \; # Clean up all pipes
    find "${outPrefix}" -name "*singles*" -prune -o -name "*.fastq.gz" -print | sort > "${outPrefix}"/"${project}"_trimmed_quality.txt
}

#   Export the function
export -f Quality_Trimming
