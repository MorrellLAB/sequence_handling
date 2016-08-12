#!/bin/bash

#   This script proceses SAM files
#   including sorting and deduplicating

set -e
set -o pipefail

#   What are the dependencies for SAM_Processing
declare -a SAM_Processing_Dependencies=(java samtools)

#   A function to check to make sure Picard is where it actually is
function checkPicard() {
    local Picard="$1" # Where is Picard?
    if ! [[ -f "${Picard}" ]]; then echo "Failed to find Picard" >&2; return 1; fi # If we can't find Picard, exit with error
}

#   Export the function
export -f checkPicard

#   A function to make our outdirectories
function makeOutDirectories() {
    local outBase="$1"
    mkdir -p "${outBase}"/Raw_SAM_Stats "${outBase}"/Sorted_BAM "${outBase}"/Deduped_BAM/stats "${outBase}"/Finished/stats
}

#   Export the function
export -f makeOutDirectories

#   A function to check to see if we have a directory to store temp files
function checkTemp() {
    local tmp=${1}
    if [[ -z $tmp ]]
    then
        local arg=''
    else
        local arg="TMP_DIR=$tmp"
    fi
    mkdir -p "${tmp}"
    echo ${arg}
}

#   Export the function
export -f checkTemp

#    A function to process the SAM files using Picard
function SAM_Processing(){
    local SAMFile="$1" # What is our SAM file?
    local outDirectory="$2"/SAM_Processing/Picard # Where do we store our results?
    local picardJar="$3" # Where is our JAR for Picard?
    local platform="$4" # What platform were our samples sequenced on?
    local maxMem="$5" # What is the most amount of memory that we can use?
    local maxFiles="$6" # What is the maximum number of file handles that we can use?
    local tmp="$7" # Where should Picard store temporary files?
    local sampleName=$(basename "${SAMFile}" .sam)
    #    Check if temp directory exists
    local tempArg=$(checkTemp $tmp)
    #   Make the out directories
    makeOutDirectories "${outDirectory}"
    #    Generate metrics on the input SAM file
    samtools flagstat "${SAMFile}" > "${outDirectory}/Raw_SAM_Stats/${sampleName}_raw.stats"
    #   Sort the SAM files and convert to BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} SortSam \
        INPUT="${SAMFile}" \
        OUTPUT="${outDirectory}/Sorted_BAM/${sampleName}_sorted.bam" \
        SO="coordinate" \
        CREATE_INDEX="true" \
        VALIDATION_STRINGENCY="SILENT" \
        "${tempArg}"
    #   Deduplicate the BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} MarkDuplicates \
        INPUT="${outDirectory}/Sorted_BAM/${sampleName}_sorted.bam" \
        OUTPUT="${outDirectory}/Deduped_BAM/${sampleName}_deduped.bam" \
        METRICS_FILE="${outDirectory}/Deduped_BAM/stats/${sampleName}_Duplication_Metrics.txt" \
        REMOVE_DUPLICATES="true" \
        ASSUME_SORTED="true" \
        CREATE_INDEX="true" \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${maxFiles} \
        "${tempArg}"
    #   Add read group information to the BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} AddOrReplaceReadGroups \
        INPUT="${outDirectory}/Deduped_BAM/${sampleName}_deduped.bam" \
        OUTPUT="${outDirectory}/Finished/${sampleName}_finished.bam" \
        RGID="${sampleName}" \
        RGLB="${sampleName}" \
        RGPL="${platform}" \
        RGPU="${sampleName}" \
        RGSM="${sampleName}" \
        CREATE_INDEX="true" \
        "${tempArg}"
    #    Generate metrics on the finished BAM files    
    samtools flagstat "${outDirectory}/Finished/${sampleName}_finished.bam" > "${outDirectory}/Finished/stats/${sampleName}_finished.stats"
}

#    Export the function
export -f SAM_Processing