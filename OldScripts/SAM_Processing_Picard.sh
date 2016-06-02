#!/bin/sh

#   This script proceses SAM files
#   including sorting and deduplicating.
#   Final outputs are BAM files.

set -e
set -o pipefail

#   What are the dependencies for SAM_Processing_Picard
declare -a SAM_Processing_Picard_Dependencies=(samtools, picard, java1.8)

#PBS -l mem=22gb,nodes=1:ppn=8,walltime=36:00:00
#   Define a path to the SAMTools in the 'SAMTOOLS' field on line 65
#       If using MSI, leave the definition as is to use their installation of SAMTools
#       Otherwise, it should look like this:
#           SAMTOOLS=${HOME}/software/samtools
#       Please be sure to comment out (put a '#' symbol in front of) the 'module load samtools' on line 64
#       And to uncomment (remove the '#' symbol) from lines 65 and 66
#   Define paths to the Picard installation, but not to 'picard.jar' itself, on line 69
#       If using MSI's Picard module, leave this blank, otherwise it should look like this:
#           PICARD_DIR=${HOME}/directory/picard
#       Where the 'picard.jar' file is located within '${HOME}/directory/picard'
#   Set the sequencing platform on line 72
#       This should look like:
#           PLATFORM=Illumina
#       This is set to Illumina by default; please see Picard's documentaiton for more details
#SAMPLE_INFO=
#REF_GEN=
#SCRATCH=
#PROJECT=
#	Load the SAMTools module for MSI, else define path to SAMTools installation
#module load samtools
#SAMTOOLS=
#export PATH=$PATH:${SAMTOOLS}
#   Define path to the directory containing 'picard.jar' if not using MSI's Picard module
#PICARD_DIR=
#   Sequencing platform
#PLATFORM=Illumina

#   Find our mapped samples in the outdirectory
function findMappedSAM() {
    local readMappedDirectory="$1" # Where are the mapped reads stored?
    local project="$2" # What is our project called?
    local SAMList="${readMappedDirectory}"/"${project}"_Mapped.txt # Create a name for our sample list
    find "${readMappedDirectory}" -name "*.sam" | sort > "${SAMList}" # Create our list
    echo "${SAMList}" # Return the name of our sample list
}

#   Export the function
export -f findMappedSAM

#   A function to see if our reference FASTA is indexed
function checkFaidx() {
    local reference="$1" # What is our reference FASTA file?
    if ! [[ -f "${reference}" ]]; then echo "Cannot find the reference genome, exiting..." >&2; exit 1; fi # Make sure it exists
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference directory
    if ! [[ $(find "${referenceDirectory}" -maxdepth 1 -name "${referenceName}.fai") ]]; then return 1; fi # Check to make sure we have the index file for our reference FASTA file
}

#   Export the function
export -f checkFaidx

#   A function to index the FASTA and exit
function faidxReference() {
    local reference="$1" # What is our reference FASTA file?
    echo "Indexing reference for SAM Processing, will quit upon completion..." >&2
    samtools faidx "${reference}" # Index our reference FASTA file
    echo "Please re-run sequence_handling to map reads" >&2
    exit 10 # Exit the script with a unique exit status
}

#   Export the function
export -f faidxReference

#   A function to make our outdirectories
function makeOutDirectories() {
    local outBase="$1"
    mkdir -p "${outBase}"/fixedHeader "${outBase}"/rawBAM/stats "${outBase}"/sorted/stats "${outBase}"/finished/stats
}

#   Export the function
export -f makeOutDirectories

#   Define a function to process SAM files using Picard
function SAMProcessingPicard() {
    #   Today's date
    YMD=`date +%Y-%m-%d`
    #   In order: SAM file being worked with, reference genome, outdirectory,
    #       sequencing platform, project name
    SAMFILE="$1"
    REF_SEQ="$3"
    OUTDIR="$4"
    TYPE="$5"
    PROJ="$6"
    #   Define path to Picard, either local or using MSI's Picard module
    /panfs/roc/itascasoft/picard/2.1.1
    if [[ -z "$2" ]]
    then
        #   Use MSI's Picard module, with changes to Java's heap memeory options
        module load picard
        PICARD_DIR=`echo $PTOOL | rev | cut -d " " -f 1 - | rev`
        PTOOL="java -Xmx15g -XX:MaxPermSize=10g -jar ${PICARD_DIR}"
    else
        #   Use local Picard installaiton, passed by arguments
        PICARD_DIR="$2"
        PTOOL="java -Xmx15g -jar ${PICARD_DIR}"
    fi
    #   Sample name, taken from full name of SAM file
    sampleName=`basename "${SAMFILE}" .sam`
#   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFile}" > "${out}"/fixedHeader/"${sampleName}"_FixedHeader.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${reference}" "${out}"/fixedHeader/"${sampleName}"_FixedHeader.sam > "${out}/rawBAM/${sampleName}_${YMD}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${out}/rawBAM/${sampleName}_${YMD}_raw.bam" > "${out}/rawBAM/stats/${sampleName}_${YMD}_raw_stats.out"
    #   Sort the raw BAM file
    "${PTOOL}"/picard.jar SortSam \
        INPUT="${OUTDIR}/raw/${SAMPLE_NAME}_${YMD}_raw.bam" \
        OUTPUT="${OUTDIR}/sorted/${SAMPLE_NAME}_${YMD}_sorted.bam" \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        TMP_DIR="${HOME}/Scratch/Picard_Tmp"
    #   Deduplicate the sorted BAM file
    "${PTOOL}"/picard.jar MarkDuplicates \
        INPUT="${OUTDIR}/sorted/${SAMPLE_NAME}_${YMD}_sorted.bam" \
        OUTPUT="${OUTDIR}/deduped/${SAMPLE_NAME}_${YMD}_deduplicated.bam" \
        METRICS_FILE="${OUTDIR}/deduped/${SAMPLE_NAME}_${YMD}_Metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        TMP_DIR="${HOME}/Scratch/Picard_Tmp" \
        MAX_RECORDS_IN_RAM=50000
    #   Then add read groups
    "${PTOOL}"/picard.jar AddOrReplaceReadGroups \
        INPUT="${OUTDIR}/deduped/${SAMPLE_NAME}_${YMD}_deduplicated.bam" \
        OUTPUT="${OUTDIR}/finished/${SAMPLE_NAME}_${YMD}_finished.bam" \
        RGID="${SAMPLE_NAME}" \
        RGPL="${TYPE}" \
        RGPU="${SAMPLE_NAME}" \
        RGSM="${SAMPLE_NAME}" \
        RGLB="${PROJ}_${SAMPLE_NAME}" \
        TMP_DIR="${HOME}/Scratch/Picard_Tmp" \
        CREATE_INDEX=True
    #   Create alignment statistics using SAMTools
    samtools flagstat "${out}/finished/${sampleName}_${YMD}_deduped.bam" > "${out}/finished/stats/${sampleName}_${YMD}_finished_stats.out"
}

#   Export the function
export -f SamProcessingPicard

#   Pass variables to and run the function in parallel
#cat ${SAMPLE_INFO} | parallel "dedup {} ${PICARD_DIR} ${REF_GEN} ${OUT} ${PLATFORM} ${PROJECT}"

#   A function to run the SAM processing
function SAM_Processing() {
    local SAMList="$1" # What is our list of samples?
    local outDirectory="$2"/SAM_Processing # Where are we storing our results?
    local referenceSequence="$3" # What is our reference sequence?
    local project="$4" # What do we call our results?
    makeOutDirectories "${outDirectory}" # Make our outdirectories
    parallel SAMToolsProcessing {} "${referenceSequence}" "${outDirectory}" :::: "${SAMList}" # Process our SAM files using SAMTools
    find "${outDirectory}/finished" -name "*_deduped.bam" | sort > "${outDirectory}"/"${project}"_Processed_BAM.txt # Create a list of finished files
    # samtools merge -r "${outDirectory}"/"${project}"_Merged.bam $(find "${outDirectory}/finished" -name "*_deduped.bam") # Merge the finished BAM files into one BAM file
}

#   Export the function
export -f SAM_Processing