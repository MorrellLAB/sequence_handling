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
    local memory="$5" # How much memory can Java use?
	local intervals_list="$6" # Where is the list of intervals files?
    local lod="$7" # What is the LOD threshold?
    local entropy="$8" # What is the entropy threshold?
    local qscores="$9" # Do we fix quality scores?
    declare -a intervals_array=($(grep -E ".intervals" "${intervals_list}")) # Put the intervals list into array format
    local current_intervals="${intervals_array[${PBS_ARRAYID}]}" # Get the intervals file for the current sample
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Put the sample list into array format
    local current="${sample_array[${PBS_ARRAYID}]}" # Pull out one sample to work on
    local name=$(basename ${current} .bam) # Get the name of the sample without the extension
    #   Make sure the out directory exists
    mkdir -p "${out}"
    #   Change into the output directory
    cd "${out}"
    if [[ "${qscores}" == true ]]
    then
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T IndelRealigner \
	    -R "${reference}" \
        --entropyThreshold "${entropy}" \
	    --LODThresholdForCleaning "${lod}" \
	    --targetIntervals "${current_intervals}" \
	    -I "${current}" \
        --fix_misencoded_quality_scores \
	    -o "${out}/${name}_realigned.bam")
    else
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T IndelRealigner \
	    -R "${reference}" \
        --entropyThreshold "${entropy}" \
	    --LODThresholdForCleaning "${lod}" \
	    --targetIntervals "${current_intervals}" \
	    -I "${current}" \
	    -o "${out}/${name}_realigned.bam")
    fi
}

#   Export the function
export -f Indel_Realigner
