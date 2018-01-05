#!/bin/bash

#   This script creates a targets file
#   for use with the GATK Indel Realigner.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_RTC.job

set -o pipefail

#   What are the dependencies for Realigner_Target_Creator?
declare -a Realigner_Target_Creator_Dependencies=(java)

#   A function to create the target file
function Realigner_Target_Creator() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Realigner_Target_Creator # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local memory="$5" # How much memory can java use?
	local project="$6" # What is the name of this project?
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Put the sample list into array format
    #	Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
		GATK_IN+=(-I $s)
	done
    #	Make sure the out directory exists
    mkdir -p "${out}"
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T RealignerTargetCreator \
	    -R "${reference}" \
        -nt 1 \
	    "${GATK_IN[@]}" \
	    -o "${out}/${project}.intervals")
}

#   Export the function
export -f Realigner_Target_Creator
