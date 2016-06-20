#!/bin/bash

#   This script realignes BAM files around indels.
#   It uses GATK to do so

set -e
set -o pipefail

#   What are the dependencies for Indel_Realignment?
declare -a Indel_Realignment_Dependencies=(java)

#   A function to check to make sure GATK is where it actually is
function checkGATK() {
    local GATK="$1" # Where is GATK?
    if ! [[ -f "${GATK}" ]]; then echo "Failed to find GATK" >&2; return 1; fi # If we can't find GATK, exit with error
}

#   Export the function
export -f checkGATK

#   Find our mapped samples in the outdirectory
#function findProcessedBAM() {
    #local processedDirectory="$1" # Where are the processed BAM files stored?
    #local project="$2" # What is our project called?
    #local BAMList="${processedDirectory}"/"${project}"_Processed.list # Create a name for our sample list
    #find "${processedDirectory}" -name "*.bam" | sort > "${BAMList}" # Create our list
    #echo "${BAMList}" # Return the name of our sample list
#}

#   Export the function
#export -f findProcessedBAM

#   A function to generate realignment targets
function createRealignmentTargets() {
    local GATK="$1" # Where is GATK?
    local out="$2" # Where are we storing our output files?
    local reference="$3" # What is our reference genome
    local knownIndels="$4" # What are our known indels?
    local maxMem="$5" # What is the maximum amount of memory given to us?
    local project="$6" # What should we call our output file?
    #   Run the realigner target creator
    java -Xmx"${maxMem}" -jar "${GATK}" \
        -T RealignerTargetCreator \
        -R "${reference}" \
        -o "${out}"/"${project}".intervals \
        -known "${knownIndels}"
    echo "${out}/${project}.intervals"
}

#   Export the function
export -f createRealignmentTargets

#   A function to run the indel realignment
function Indel_Realignment() {
    local sampleList="$1" # Where are our samples located?
    local outDirectory="$2"/Indel_Realignment # Where are we storing our results?
    local GATK="$3" # Where is GATK?
    local reference="$4" # What is our reference sequence?
    local knownIndels="$5" # What are our known indels?
    local model="$6" # What model are we working on?
    local LODScore="$7" # What is our log of odds score?
    local maxMem="$8" # How much memory are we allowed to use?
    local project="$9" # What do we call our results?
    local targetIntervals="$10" # What, if any, are our target intervals?
    #   Make the out directory
    mkdir -p "${outDirectory}"
    #   Check to see if we have our target intervals
    if ! [[ -f "${targetIntervals}" ]]; then targetIntervals=$(createRealignmentTargets "${GATK}" "${outDirectory}" "${reference}" "${knownIndels}" "${maxMem}" "${project}"); fi
    cd ${outDirectory}
    #   Run the indel realigner
    java -Xmx"${maxMem}" -jar "${GATK}" \
        -T IndelRealigner \
        -I "${sampleList}" \
        -R "${reference}" \
        -targetIntervals "${targetIntervals}" \
        -known "${knownIndels}" \
        --consensusDeterminationModel "${model}" \
        -LOD "${LODScore}" \
        -nWayOut ".realigned.bam"
}

#   Export the function
export -f Indel_Realignment
