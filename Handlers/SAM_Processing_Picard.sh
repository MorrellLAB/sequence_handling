#!/bin/bash

#   This script proceses SAM files
#   including sorting and deduplicating

set -o pipefail

#   What are the dependencies for SAM_Processing
declare -a SAM_Processing_Dependencies=(java samtools)

#   A function to make our outdirectories
function makeOutDirectories() {
    local outBase="$1"
    mkdir -p "${outBase}"/Statistics/Raw_SAM_Stats "${outBase}"/Statistics/Deduplicated_BAM_Stats "${outBase}"/Statistics/Finished_BAM_Stats "${outBase}"/Intermediates/Sorted "${outBase}"/Intermediates/Deduplicated
}

#   Export the function
export -f makeOutDirectories

#    A function to process the SAM files using Picard
function SAM_Processing(){
    local SAMFile="$1" # What is our SAM file?
    local outDirectory="$2"/SAM_Processing/Picard # Where do we store our results?
    local picardJar="$3" # Where is our JAR for Picard?
    local platform="$4" # What platform were our samples sequenced on?
    local maxMem="$5" # What is the most amount of memory that we can use?
    local maxFiles="$6" # What is the maximum number of file handles that we can use?
    local tmp="$7" # Where is the temp directory?
    local sampleName=$(basename "${SAMFile}" .sam)
    #   Make the out directories
    makeOutDirectories "${outDirectory}"
    #   Generate metrics on the input SAM file
    samtools flagstat "${SAMFile}" > "${outDirectory}/Statistics/Raw_SAM_Stats/${sampleName}_raw.txt"
    #   Sort the SAM files and convert to BAM file
    if [[ -z ${tmp} ]] # If tmp is left blank
    then
        #   Sort the SAM file and convert to BAM
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} SortSam \
            INPUT="${SAMFile}" \
            OUTPUT="${outDirectory}/Intermediates/Sorted/${sampleName}_sorted.bam" \
            SO="coordinate" \
            VALIDATION_STRINGENCY="SILENT")
        #   Deduplicate the BAM file
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} MarkDuplicates \
            INPUT="${outDirectory}/Intermediates/Sorted/${sampleName}_sorted.bam" \
            OUTPUT="${outDirectory}/Intermediates/Deduplicated/${sampleName}_deduped.bam" \
            METRICS_FILE="${outDirectory}/Statistics/Deduplicated_BAM_Stats/${sampleName}_deduplicated.txt" \
            REMOVE_DUPLICATES="true" \
            ASSUME_SORTED="true" \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${maxFiles})
        #   Add read group information to the BAM file
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} AddOrReplaceReadGroups \
            INPUT="${outDirectory}/Intermediates/Deduplicated/${sampleName}_deduped.bam" \
            OUTPUT="${outDirectory}/${sampleName}.bam" \
            RGID="${sampleName}" \
            RGLB="${sampleName}" \
            RGPL="${platform}" \
            RGPU="${sampleName}" \
            RGSM="${sampleName}")
    else    # If a tmp is provided
        #   Make sure tmp exists
        mkdir -p ${tmp}
        #   Make sure we have write permissions for tmp
        if ! [[ -w "${tmp}" ]]; then echo "You do not have write permissions for ${tmp}, exiting..." >&2; exit 1; fi
        #   Sort the SAM file and convert to BAM
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} SortSam \
            INPUT="${SAMFile}" \
            OUTPUT="${outDirectory}/Intermediates/Sorted/${sampleName}_sorted.bam" \
            SO="coordinate" \
            VALIDATION_STRINGENCY="SILENT" \
            TMP_DIR="${tmp}")
        #   Deduplicate the BAM file
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} MarkDuplicates \
            INPUT="${outDirectory}/Intermediates/Sorted/${sampleName}_sorted.bam" \
            OUTPUT="${outDirectory}/Intermediates/Deduplicated/${sampleName}_deduped.bam" \
            METRICS_FILE="${outDirectory}/Statistics/Deduplicated_BAM_Stats/${sampleName}_deduplicated.txt" \
            REMOVE_DUPLICATES="true" \
            ASSUME_SORTED="true" \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${maxFiles} \
            TMP_DIR="${tmp}")
        #   Add read group information to the BAM file
        (set -x; java -Xmx"${maxMem}" -jar ${picardJar} AddOrReplaceReadGroups \
            INPUT="${outDirectory}/Intermediates/Deduplicated/${sampleName}_deduped.bam" \
            OUTPUT="${outDirectory}/${sampleName}.bam" \
            RGID="${sampleName}" \
            RGLB="${sampleName}" \
            RGPL="${platform}" \
            RGPU="${sampleName}" \
            RGSM="${sampleName}" \
            TMP_DIR="${tmp}")
    fi
    #   Generate metrics on the finished BAM file
    samtools flagstat "${outDirectory}/${sampleName}.bam" > "${outDirectory}/Statistics/Finished_BAM_Stats/${sampleName}_finished.txt"
    #   Index the finished BAM file
    samtools index "${outDirectory}/${sampleName}.bam"
    #   Rename the index file
    mv "${outDirectory}/${sampleName}.bam.bai" "${outDirectory}/${sampleName}.bai"
    #   Remove intermediate files - comment out these lines if you need to troubleshoot
    rm "${outDirectory}/Intermediates/Sorted/${sampleName}_sorted.bam"
    rm "${outDirectory}/Intermediates/Deduplicated/${sampleName}_deduped.bam"
}

#    Export the function
export -f SAM_Processing
