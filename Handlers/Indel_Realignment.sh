#!/bin/bash

set -e
set -o pipefail

declare -a Indel_Realignment_Dependencies(java)

function checkGATK() {
    local GATK="$1"
    if ! [[ -f "${GATK}" ]]; then echo "Failed to find GATK" >&2; return 1; fi
}

function createRealignmentTargets() {
    local GATK="$1"
    local out="$2"
    local reference="$3"
    local knownIndels="$4"
    local maxMem="$5"
    local project="$6"
    java -Xmx"${maxMem}" -jar "${GATK}" \
        -T RealignerTargetCreator \
        -R "${reference}" \
        -o "${out}"/"${project}".intervals \
        -known "${knownIndels}"
    echo "${out}/${project}.intervals"
}

function Indel_Realignment() {
    local sampleList="$1" # What is our list of samples?
    local outDirectory="$2"/Indel_Realignment # Where are we storing our results?
    local GATK="$3" # Where is GATK?
    local reference="$4" # What is our reference sequence?
    local knownIndels="$5" # What are our known indels?
    local model="$6" # What model are we working on?
    local LODScore="$7" # What is our log of odds score?
    local maxMem="$8" # How much memory are we allowed to use?
    local project="$9" # What do we call our results?
    local targetIntervals="$10" # What, if any, are our target intervals?
    #   Check to see if we have our target intervals
    if ! [[ -f "${targetIntervals}" ]]; then targetIntervals=$(createRealignmentTargets "${GATK}" "${outDirectory}" "${reference}" "${knownIndels}" "${maxMem}" "${project}"); fi
    # java -Xmx"${maxMem}" -jar "${GATK}" -T IndelRealigner -nWayOut
}