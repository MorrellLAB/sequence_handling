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
    mkdir -p "${outBase}"/raw_SAM_stats "${outBase}"/sorted_BAM "${outBase}"/deduped_BAM "${outBase}"/deduped_BAM_stats "${outBase}"/finished_BAM "${outBase}"/finished_BAM_stats
}

#   Export the function
export -f makeOutDirectories

function checkTemp() {
    local tmp=${1}
    if [[ -z $temp ]]
    then
        local arg="TMP_DIR=$tmp"
    else
        local arg=""
    fi
    echo ${arg}
}

#   Export the function
export -f checkTemp

#    A function to process the SAM files using Picard
function SAM_Processing(){
    local SAMFile="$1" # What is our SAM file?
    local outDirectory="$2"/SAM_Processing # Where do we store our results?
    local picardJar="$3" # Where is our JAR for Picard?
    local platform="$4" # What platform were our samples sequenced on?
    local maxMem="$5" # What is the most amount of memory that we can use?
    local maxFiles="$6" # What is the maximum number of file handles that we can use?
    local tmp="$7" # Where should Picard store temporary files?
    local sampleName=$(basename "${SAMFile}" .sam)
    #    Test to make sure the input SAM file is valid
    local tempArg=$(checkTemp $tmp)
    makeOutDirectories "${outDirectory}"
    #samtools quickcheck -vvvvv "${SAMFile}"
    #if [ $? -ne 0 ]; then
        #echo "Samtools quickcheck failed. Check input SAM file formatting."
        #exit 1
    #fi
    #    Generate metrics on the input SAM file
    samtools flagstat "${SAMFile}" > "${outDirectory}/raw_SAM_stats/${sampleName}_raw_stats.out"
    #   Sort the SAM files and convert to BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} SortSam \
        INPUT="${SAMFile}" \
        OUTPUT="${outDirectory}/sorted_BAM/${sampleName}_sorted.bam" \
        SO="coordinate" \
        CREATE_INDEX="true" \
        VALIDATION_STRINGENCY="SILENT" \
        ${tempArg}
    #   Deduplicate the BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} MarkDuplicates \
        INPUT="${outDirectory}/sorted_BAM/${sampleName}_sorted.bam" \
        OUTPUT="${outDirectory}/deduped_BAM/${sampleName}_deduped.bam" \
        METRICS_FILE="${outDirectory}/deduped_BAM_stats/${sampleName}_DuplicationMetrics.txt" \
        REMOVE_DUPLICATES="true" \
        ASSUME_SORTED="true" \
        CREATE_INDEX="true" \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${maxFiles} \
        TMP_DIR="${tmp}"
    #   Add read group information to the BAM files
    java -Xmx"${maxMem}" -jar ${picardJar} AddOrReplaceReadGroups \
        INPUT="${outDirectory}/deduped_BAM/${sampleName}_deduped.bam" \
        OUTPUT="${outDirectory}/finished_BAM/${sampleName}_finished.bam" \
        RGID="${sampleName}" \
        RGLB="${sampleName}" \
        RGPL="${platform}" \
        RGPU="${sampleName}" \
        RGSM="${sampleName}" \
        CREATE_INDEX="true" \
        TMP_DIR="${tmp}"
    #    Generate metrics on the finished BAM files    
    samtools flagstat "${outDirectory}/finished_BAM/${sampleName}_finished.bam" > "${outDirectory}/finished_BAM_stats/${sampleName}_finished_stats.out"
    #    Check validity of output SAM files
    #samtools quickcheck -vvvvv "${outDirectory}/finished_BAM/${sampleName}_finished.bam"
    #if [ $? -ne 0 ]; then
        #echo "Samtools quickcheck failed. Check output BAM file formatting."
        #exit 1
    #fi
}

#    Export the function
export -f SAM_Processing

#Picard_Processing /panfs/roc/groups/9/morrellp/hoffmanp/downsample/SAM_Raw/2012BM7F068_ATTCCT_L001_2015-10-12.sam /home/morrellp/wyant008/picard_testing/Downsample /panfs/roc/itascasoft/picard/2.1.1/picard.jar sanger 15g 1000

#   A function to run the SAM processing
#function SAM_Processing() {
    #local SAMList="$1" # What is our list of samples?
    #local outDirectory="$2"/SAM_Processing # Where are we storing our results?
    #local referenceSequence="$3" # What is our reference sequence?
    #local project="$4" # What do we call our results?
    #makeOutDirectories "${outDirectory}" # Make our outdirectories
    #parallel SAMToolsProcessing {} "${referenceSequence}" "${outDirectory}" :::: "${SAMList}" # Process our SAM files using SAMTools
    #find "${outDirectory}/finished" -name "*_deduped.bam" | sort > "${outDirectory}"/"${project}"_Processed_BAM.txt # Create a list of finished files
    # samtools merge -r "${outDirectory}"/"${project}"_Merged.bam $(find "${outDirectory}/finished" -name "*_deduped.bam") # Merge the finished BAM files into one BAM file
#}