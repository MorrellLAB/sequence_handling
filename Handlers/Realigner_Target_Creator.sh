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
    local memory="$5" # How much memory can Java use?
    local qscores="$6" # Do we fix quality scores?
    local gatkVer="$7" # Either 3 or 4
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Put the sample list into array format
    local current="${sample_array[${PBS_ARRAYID}]}" # Pull out one sample to work on
    local name=$(basename ${current} .bam) # Get the name of the sample without the extension

    #	Make sure the out directory exists
    mkdir -p "${out}"
    if [[ "${gatkVer}" == 3 ]]; then
        if [[ "${qscores}" == true ]]
        then
        #   Run GATK3 using the parameters given
        java -Xmx"${memory}" -jar "${gatk}" \
            -T RealignerTargetCreator \
            -R "${reference}" \
            -nt 1 \
            -I "${current}" \
            --fix_misencoded_quality_scores \
            -o "${out}/${name}.intervals"
        else
        #   Run GATK3 using the parameters given
        java -Xmx"${memory}" -jar "${gatk}" \
            -T RealignerTargetCreator \
            -R "${reference}" \
            -nt 1 \
            -I "${current}" \
            -o "${out}/${name}.intervals"
        fi
    else
        echo "Indel realignment functionality is no longer available in GATK 4. Please use GATK 3."
    fi
}

#   Export the function
export -f Realigner_Target_Creator
