#!/bin/bash

#   This script realigns mapped regions
#   around insertions and deletions.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_IndelRealigner.job

set -e
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
    local gatkVer="${10}" # Either 3 or 4
    local max_reads_in_mem="${11}"
    declare -a intervals_array=($(grep -E ".intervals" "${intervals_list}")) # Put the intervals list into array format
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Put the sample list into array format
    if [[ "$USE_PBS" == "true" ]]; then
        local current_intervals="${intervals_array[${PBS_ARRAYID}]}" # Get the intervals file for the current sample
        local current="${sample_array[${PBS_ARRAYID}]}" # Pull out one sample to work on
    elif [[ "${USE_SLURM}" == true ]]; then
        local current_intervals="${intervals_array[${SLURM_ARRAY_TASK_ID}]}" # Get the intervals file for the current sample
        local current="${sample_array[${SLURM_ARRAY_TASK_ID}]}" # Pull out one sample to work on
    fi
    local name=$(basename ${current} .bam) # Get the name of the sample without the extension
    #   Make sure the out directory exists
    mkdir -p "${out}"
    #   Change into the output directory
    cd "${out}"
    #   To prevent potential memory issues, the -Xmx value should be a little less than the total amount of physical memory
    #       by at least a few GB.
    #   So, we will subtract a few GB of mem from user provided memory
    mem_num=$(basename ${memory} g)
    # Check that we have requested the minimum required memory (>4g)
    if [ ${mem_num} -le 4 ]; then
        echo "You have 4g or less requested memory, please request more memory for this handler to work. Exiting..."
        exit 31
    fi
    new_mem_num="$((${mem_num}-4))"
    mem=$(printf "${new_mem_num}g")
    #   Start indel realignment
    if [[ "${gatkVer}" == 3 ]]; then
        if [[ "${qscores}" == true ]]
        then
        #   Run GATK using the parameters given
        java -Xmx"${mem}" -jar "${gatk}" \
            -T IndelRealigner \
            -R "${reference}" \
            --entropyThreshold "${entropy}" \
            --LODThresholdForCleaning "${lod}" \
            --targetIntervals "${current_intervals}" \
            --maxReadsInMemory ${max_reads_in_mem} \
            -I "${current}" \
            --fix_misencoded_quality_scores \
            -o "${out}/${name}_realigned.bam"
        else
        #   Run GATK using the parameters given
        java -Xmx"${mem}" -jar "${gatk}" \
            -T IndelRealigner \
            -R "${reference}" \
            --entropyThreshold "${entropy}" \
            --LODThresholdForCleaning "${lod}" \
            --targetIntervals "${current_intervals}" \
            --maxReadsInMemory ${max_reads_in_mem} \
            -I "${current}" \
            -o "${out}/${name}_realigned.bam"
        fi
    fi
}

#   Export the function
export -f Indel_Realigner
