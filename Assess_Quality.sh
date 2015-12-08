#!/bin/sh

#   This script runs FastQC on a series of samples
#   Please have FastQC installed before using this script

#PBS -l mem=4000mb,nodes=1:ppn=4,walltime=6:00:00
#PBS -m abe
#PBS -M user@example.com
#PBS -q lab

set -e
set -o pipefail

# CONFIG_FILE=$1
# UTILS=$2
#
# source ${CONFIG_FILE}
# source ${UTIlS}

# #   Check our sample list
# checkSamples ${RAW_SAMPLES}
#
# if [[ "$?" -eq 1 ]]; then echo "Error with samples"; exit 1; fi
#
# #   Check to see if FastQC was loaded or defined properly
# if ! `command -v fastqc > /dev/null 2> /dev/null`
# then
#     echo "Please make sure FastQC is either loaded or installed properly"
#     exit 1
# fi
#
# if ! `command -v parallel > dev/null 2> /dev/null`
# then
#     echo "Please make sure GNU Parallel is installed properly"
#     exit 1
# fi

#   What are the dependencies Assess Quality?
declare -a Assess_Quality_Dependencies=(fastqc parallel)

#   A function to run the quality assesment
function Assess_Quality() {
    local sampleList="$1" # What is our list of samples?
    local outDir="$2" # Where are we storing our results?
    local proj="$3" # What do we call our results?
    local out="${outDir}"/Quality_Assessment # Full path to output directory
    mkdir -p "${out}" # Make our output directory
    cat "${sampleList}" | parallel "fastqc --outdir ${out} {}" # Run FastQC in parallel
    find "${out}" -name "*.zip" | sort > "${outDir}"/"${proj}"_FastQC_ZipFiles.txt # Write a list of zipped files
}

export -f Assess_Quality

# echo "Assess_Quality ${RAW_SAMPLES} ${OUT_DIR} ${PROJECT}" | qsub -l "${QC_QSUB}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Assessment

# cat ${RAW_SAMPLES} | parallel "fastqc --outdir ${OUT} {}"
#
# # Create a list of ZIP files for use in counting read depth
# find ${OUT} -name "*.zip" | sort > ${OUT}/FastQC_zipfiles.txt
# echo "List of ZIP files can be found at"
# echo "${OUT}/FastQC_zipfiles.txt"
