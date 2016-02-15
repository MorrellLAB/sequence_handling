#!/bin/env bash

#   This script performs quality and adapter trimming
#   on a series of FASTQ samples using Sickle, Seqqs,
#   and Scythe. Please install these before use.

set -e
set -o pipefail

#   What are the dependencies for Quality_Trimming?
declare -a Quality_Trimming_Dependencies=(sickle seqqs scythe Rscript)

#   Test to see if there are equal numbers of forward and reverse reads
function checkCounts() {
    local rawSamples="$1" # Sample list
    local forwardNaming="$2" # Forward naming scheme
    local reverseNaming="$3" # Reverse naminc scheme
    local fwdCount="`grep -ce ${forwardNaming} ${rawSamples}`" # Count the number of forward samples
    local revCount="`grep -ce ${reverseNaming} ${rawSamples}`" # Count the number of reverse samples
    if ! [[ "${fwdCount}" -eq "${revCount}" ]]; then return 1; fi # If our numbers are not equal, return an exit status of 1
}

#   Export the function
export -f checkCounts

#   A function to perform the trimming and plotting
# #       Adapted from Tom Kono and Peter Morrell
# function trimAutoplot() {
#     #   Set the arguments for trimming
#     local name="$1" #   Name of the sample
#     local forward="$2" #  Forward file
#     local reverse="$3" #  Reverse file
#     local out="$4"/"${name}" #  Outdirectory
#     local adapters="$5"
#     local prior="$6"
#     local threshold="$7"
#     local platform="$8"
#     local helperScripts="$9"
#     #   Make the out directories
#     local stats="${out}"/stats # Directory for stats
#     local plots="${stats}"/plots # Directory for plots
#     mkdir -p "${plots}" # Make the directories
#     #   Trim the sequences based on quality and adapters
#     sickle pe -t "${platform}" -q "${threshold}" \
#         -f <(seqqs -e -q "${platform}" -p "${stats}"/raw_"${name}"_R1 "${forward}" | scythe -a "${adapters}" -p "${prior}" - 2> "${stats}"/"${name}"_R1_scythe.stderr) \
#         -r <(seqqs -e -q "${platform}" -p "${stats}"/raw_"${name}"_R2 "${reverse}" | scythe -a "${adapters}" -p "${prior}" - 2> "${stats}"/"${name}"_R2_scythe.stderr) \
#         -o >(seqqs -e -q "${platform}" -p "${stats}"/trimmed_"${name}"_R1 - | gzip > "${out}"/"${name}"_R1_trimmed.fq.gz) \
#         -p >(seqqs -e -q "${platform}" -p "${stats}"/trimmed_"${name}"_R2 - | gzip > "${out}"/"${name}"_R2_trimmed.fq.gz) \
#         -s >(seqqs -e -q "${platform}" -p "${stats}"/trimmed_"${name}"_singles - | gzip > "${out}"/"${name}"_singles_trimmed.fq.gz) > "${stats}"/"${name}"_sickle.stderr
#     #   Fix the quality scores
#     "${helperScripts}"/fix_quality.sh "${stats}"/raw_"${name}"_R1_qual.txt 33
#     "${helperScripts}"/fix_quality.sh "${stats}"/raw_"${name}"_R2_qual.txt 33
#     "${helperScripts}"/fix_quality.sh "${stats}"/trimmed_"${name}"_R1_qual.txt 33
#     "${helperScripts}"/fix_quality.sh "${stats}"/trimmed_"${name}"_R2_qual.txt 33
#     #   Make the plots
#     Rscript "${helperScripts}"/plot_seqqs.R "${stats}" "${name}"
# }

#   A function to help make an array for
function makeOutdirs() {
    local name="$1" # Name of the sample
    local out="$2"/"${name}" # Outdirectory
    echo "${out}" "${out}" # Echo the name of the outdirectory twice for forward and reverse samples
}

#   Export the function
export -f makeOutdirs

function checkCompressed() {
    local names="$1" # Name of the sample
    local forward="$2" # Forward file
    local reverse="$3" # Reverse file
    local out="$4" # Directory with pipes
    if [[ "$(echo ${forward} | rev | cut -f 1 -d '.' | rev)" == "gz" ]] # If we have gzipped files
    then
        local F_PIPE="${out}/${name}_F" # Make a name for a pipe
        mkfio "${F_PIPE}" # Make the named pipe
        gzip -cd "${forward}" > "${F_PIPE}" & # Unzip the forward file to the pipe
        echo "${F_PIPE}" # Echo the name of the pipe
    else # Otherwise
        echo "${forward}" # Echo the name of the forward file
    fi
    if [[ "$(echo ${reverse} | rev | cut -f 1 -d '.' | rev)" == "gz" ]] # Echo the name of the pipe
    then
        local R_PIPE="${out}/${name}_R" # Make a name for the pipe
        mkfio "${R_PIPE}" # Make the named pipe
        gzip -cd "${reverse}" > "${R_PIPE}" & # Unzip the reverse file to the pipe
        echo "${R_PIPE}" # Echo the name of the pipe
    else # Otherwise
        echo "${reverse}" # Echo the name of the reverse file
    fi
}

#   A function to fix the name for forward and reverse samples
function fixNames() {
    local name="$1" # Name of the sample
    local forwardNaming="$(echo $2 | cut -f 1 -d '.')" # Forward naming scheme, without the extension
    local reverseNaming="$(echo $3 | cut -f 1 -d '.')" # Reverse naming scheme, without the extension
    echo "${name}${forwardNaming}" "${name}${reverseNaming}" # Echo the forward and reverse naming scheme
}

#   Export the function
export -f fixNames

#       Adapted from Tom Kono
function trimAutoplot() {
    #   Set the arguments for the trimming
    local name="$1" # Name for the sample
    local sample="$2" # What sample are we working with
    local out="$3" #  Outdirectory
    local adapters="$4" # Adapter file
    local prior="$5" # Prior
    local platform="$6" # Quality encoding platform
    local helperScripts="$7" # Directory with helper scripts
    local matches="${out}/${name}.match"
    local trimmed="${out}/${name}_ScytheTrimmed.fastq.gz"
    scythe -a "${adapters}" -p "${prior}" -m "${matches}" --quiet | gzip -c > "${trimmed}"
}

#   Export the function
export -f trimAutoplot

#   A function to run the quality trimming
function Quality_Trimming() {
    local rawSamples="$1" # What is our list of samples?
    local outDirectory="$2"/Quality_Trimming # Where are we storing our results?
    local project="$3" # What do we call our results?
    local forwardNaming="$4" # What is the extension indicating a forward read?
    local reverseNaming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local threshold="$8" # What is the threshold for trimming?
    local platform="$9" # What platform did we sequence on?
    local helperScripts="${10}"/Helper_Scripts # Where are the helper scripts stored?
    checkCounts "${rawSamples}" "${forwardNaming}" "${reverseNaming}" # Check the counts of forward and reverse reads
    if [[ "$?" -ne 0 ]]; then echo "Unbalanced forward and reverse reads" >&2; exit 1; fi # If not an equal amount, exit out with error
    mkdir -p outDirectory # Make our out directory
    local -a forwardSamples=($(grep -E "${forwardNaming}" "${rawSamples}")) # Get the forward samples
    local -a reverseSamples=($(grep -E "${reverseNaming}" "${rawSamples}")) # Get the reverse samples
    local -a sampleNames=($(parallel basename {} "${forwardNaming}" ::: "${forwardSamples}")) # Get our sample names
    local -a fullSamples=($(parallel --xapply checkCompressed {1} {2} {3} "${outDirectory}" ::: "${sampleNames[@]}" ::: "${forwardSamples[@]}" ::: "${reverseSamples[@]}")) # Create a list of all samples or pipes if gzipped samples
    local -a outdirectories=($(parallel makeOutdirs {} "${outDirectory}" ::: "${sampleNames[@]}")) # Create a list of outdirectories
    local -a fixedNames=($(parallel fixNames {} "${forwardNaming}" "${reverseNaming}" ::: "${sampleNames[@]}"))
    # parallel --xapply trimAutoplot {1} {2} {3} "${outDirectory}" "${adapters}" "${prior}" "${threshold}" "${platform}" "${helperScripts}" ::: "${sampleNames[@]}" ::: "${forwardSamples[@]}" ::: "${reverseSamples[@]}" # Run trimAutoplot in parallel
    # parallel --xapply trimAutoplot {1} {2} {3} "${outDirectory}" "${adapters}" "${prior}" "${platform}" "${helperScripts}" ::: "${sampleNames[@]}" ::: "${forwardSamples[@]}" ::: "${reverseSamples[@]}"
    # parallel --xapply trimAutoplot {1} {2} {3} "${adapters}" "${prior}" "${platform}" "${helperScripts}" ::: "${fixedNames[@]}" ::: "${fullSamples[@]}" ::: "${outdirectories}"
    echo $adapters $prior $platform $helperScripts
    find "${outDirectory}" -regex ".*_R[1-2]_trimmed.fq.gz" | sort > "${outDirectory}"/"${project}"_trimmed.txt # Create our list of trimmmed files
}
