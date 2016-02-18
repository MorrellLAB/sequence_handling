#!/bin/bash

#   This script generates a series of QSub submissions for read mapping
#   The Burrows-Wheeler Aligner (BWA) and the Portable Batch System (PBS)
#   are required to use this script

set -e
set -o pipefail

#   What are the dependencies for Read_Mapping?
declare -a Read_Mapping_Dependencies=(bwa)

#   A function to parse BWA settings
function ParseBWASettings() {
	local POSITIONALS='' # Create a string of positional arguments for BWA mem
	if [[ "${PAIRED}" == true ]]; then POSITIONALS="${POSITIONALS}"'-P '; fi # Add paired to the positionals
	if [[ "${INTERLEAVED}" == true ]]; then POSITIONALS="${POSITIONALS}"'-p '; fi # Add interleaved to the positionals
	if [[ "${SECONDARY}" == true ]]; then POSITIONALS='-a'; else SECONDARY=''; fi # Add secondary to the positionals
	if [[ "${APPEND}" == true ]]; then POSITIONALS="${POSITIONALS}"'-C '; fi # Add append to the positionals
	if [[ "${HARD}" == true ]]; then POSITIONALS="${POSITIONALS}"'-H '; fi # Add hard to the positionals
	if [[ "${SPLIT}" == true ]]; then POSITIONALS="${POSITIONALS}"'-M '; fi # Add split to the positionals
	if [[ "${VERBOSITY}" == 'disabled' ]]; then VERBOSITY=0; elif [[ "${VERBOSITY}" == 'errors' ]]; then VERBOSITY=1; elif [[ "${VERBOSITY}" == 'warnings' ]]; then VERBOSITY=2; elif [[ "${VERBOSITY}" == 'all' ]]; then VERBOSITY=3; elif [[ "${VERBOSITY}" == 'debug' ]]; then VERBOSITY=4; else echo "Failed to recognize verbosity level, exiting..."; exit 1; fi # Set the verbosity level
	MEM_SETTINGS=$(echo "-t ${THREADS} -k ${SEED} -w ${WIDTH} -d ${DROPOFF} -r ${RE_SEED} -A ${MATCH} -B ${MISMATCH} -O ${GAP} -E ${EXTENSION} -L ${CLIP} -U ${UNPAIRED} -T ${THRESHOLD} -v ${VERBOSITY} ${POSITIONALS}")
    echo "${MEM_SETTINGS}"
}

#   Export the function
export -f ParseBWASettings

#   A function to see if our referenced FASTA is indexed
function checkIndex() {
    local reference="$1"
    if ! [[ -f "${reference}" ]]; then echo "Cannot find reference genome, exiting..." >&2; exit 1; fi
    local referenceDirectory=$(dirname "${reference}")
    local referenceName=$(basename "${reference}")
    if ! [[ $(find "${referenceDirectory}" -maxdepth 1 -name "${referenceName}.bwt") ]]; then return 1; fi
}

#   Export the function
export -f checkIndex

#   A function to index the FASTA and exit
function indexReference() {
    local reference="$1"
    echo "Indexing reference, will quit upon completion..." >&2
    bwa index "${reference}"
    echo "Please re-run sequence_handling to map reads" >&2
    exit 10
}

#   Export the function
export -f indexReference

function createReadGroupID() {
    local sample="$1"
    local project="$2"
    local platform="$3"
    local readGroupID="@RG\tID:${sample}\tLB${project}_${sample}\tPL:${platform}\tPU${sample}\tSM:${sample}"
    echo "${sample}"
}

function Read_Mapping() {
    local sampleName="$1"
    local forwardSample="$2"
    local reverseSample="$3"
    local project="$4"
    local platform="$5"
    local memSettings=$(ParseBWASettings)
    local readGroupID=$(createReadGroupID "${sample}" "${project}" "${platform}")
    # local -a forwardSamples=($(grep -E "${forwardNaming}" "${sampleList}"))
    # local -a reverseSamples=($(grep -E "${reverseNaming}" "${sampleList}"))
    # if ! [[ "${#forwardSamples}" -ne "${reverseSamples}" ]]; then echo "Unequal numbers of forward and reverse reads, exiting..." >&2; exit 1; fi
    exit 3
}
