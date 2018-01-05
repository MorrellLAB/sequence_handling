#!/bin/bash

#   This script realigns mapped regions
#   around insertions and deletions.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_IndelRealigner.job

set -o pipefail

#   What are the dependencies for Indel_Realigner?
declare -a Indel_Realigner_Dependencies=(java)

#   A function to create the target file
function Indel_Realigner() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Indel_Realigner # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local memory="$5" # How much memory can java use?
	local targets="$6" # Where is the targets file?
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
		GATK_IN+=(-I $s)
	done
    #   Make sure the out directory exists
    mkdir -p "${out}"
    #   Change into the output directory
    cd "${out}"
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T IndelRealigner \
	    -R "${reference}" \
        --entropyThreshold 0.10 \
	    --LODThresholdForCleaning 3.0 \
	    --targetIntervals "${targets}" \
	    "${GATK_IN[@]}" \
	    --nWayOut '_realigned.bam')
}

#   Export the function
export -f Indel_Realigner
