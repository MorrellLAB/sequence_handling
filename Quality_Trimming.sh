#!/bin/env bash

#   This script performs quality and adapter trimming
#   on a series of FASTQ samples using Sickle, Seqqs,
#   and Scythe. Please install these before use.

#PBS -l mem=1gb,nodes=1:ppn=8,walltime=20:00:00
#PBS -m abe
#PBS -M user@example.com
#PBS -q lab

set -e
set -u
set -o pipefail

module load parallel

#   This script is a QSub submission for quality trimming a batch of files.
#   To use, on line 5, change the 'user@example.com' to your own email address
#       to get notifications on start and completion for this script
#   Place the full directory path to the Sequence Handling directory on line 66
#       This should look like:
#           SEQUENCE_HANLDING=${HOME}/sequence_handling
#       Use ${HOME}, as it is a link that the shell understands as your home directory
#   Add the full file path to list of samples on the 'SAMPLE_INFO' field on line 69
#       This should look like:
#           SAMPLE_INFO=${HOME}/Directory/list.txt
#   Specify the forward and reverse file extensions in the 'FORWARD_NAMING'
#       and 'REVERSE_NAMING' fields on lines 75 and 76
#       This should look like:
#           FORWARD_NAMING=_1_sequence.txt.gz
#           REVERSE_NAMING=_2_sequence.txt.gz
#   Name the project in the 'PROJECT' field on line 79
#       This should look lke:
#           PROJECT=Barley
#   Put the full directory path for the output in the 'SCRATCH' field on line 82
#       This should look like:
#           SCRATCH="${HOME}/Out_Directory"
#       Adjust for your own OUT directory.
#   List the adapters file on line 85
#       This should look like:
#           ADAPTERS=${HOME}/adapters.fa
#   Define the prior for Scythe on line 89
#       This should look like:
#           PRIOR=0.04
#   Set Sickle's quality threshold on line 94
#       This should look like:
#           THRESHOLD=20
#       Use 20 for normal trimming and 0 for no trimming
#   Define the platform usied for sequencing on line 97
#       This should look like:
#           PLATFORM=sanger
#   Define the R installation on line 101
#       If using MSI, leave the definition as is to use their installation of R
#       Otherwise, it should look like this:
#           R=${HOME}/software/R
#       Please be sure to comment out (put a '#' symbol in front of) the 'module load R' on line 100
#       And to uncomment (remove the '#' symbol) from lines 101 and 102
#   Run this script using the qsub command
#       qsub Quality_Trimming.sh
#   This script outputs gzipped FastQ files with the extension fq.qz
#   In the stats directory, there are text files with more details about the trimming
#       as well as a plots directory
#   In the plots directory, there are PDFs showing graphs of the quality before and after the trim
#   Finally, this script outputs a list of all trimmed FastQ files for use in the Read_Mapping.sh script
#       This is stored in ${SCRATCH}/${PROJECT}/Quality_Trimming, with the path dependent on how you name these fields.


#   Where is the Sequence Handling directory located? We're using some scripts located in here
SEQUENCE_HANDLING=

#   List of samples to be processed
SAMPLE_INFO=

#   Extension on forward and reverse read names to be trimmed by basename
#       Example:
#           _1_sequence.txt.gz  for forward
#           _2_sequence.txt.gz  for reverse
FORWARD_NAMING=
REVERSE_NAMING=

#   Project name
PROJECT=

#   Output directory
SCRATCH=

#   Adapters file
ADAPTERS=

#   What's the prior?
#   Defaults to 0.04
PRIOR=0.04

#   What's the threshold?
#   Use 20 for normal trimming based on quality
#   Use 0 for no trimming based on quality
THRESHOLD=

#   What is the platform used for sequencing?
PLATFORM=

#   Load the R Module
module load R
#R_DIR=
#export PATH='${PATH}':${R_DIR}

#   Test to see if there are equal numbers of forward and reverse reads
FORWARD_COUNT="`grep -cE "${FORWARD_NAMING}" ${SAMPLE_INFO}`"
REVERSE_COUNT="`grep -cE "${REVERSE_NAMING}" ${SAMPLE_INFO}`"

if [[ "${FORWARD_COUNT}" -eq "${REVERSE_COUNT}" ]]
then
    echo Equal numbers of forward and reverse samples
else
    exit 1
fi

#   Make the out directory
OUT=${SCRATCH}/${PROJECT}/Quality_Trimming
mkdir -p ${OUT}

#   Create arrays of forward and reverse samples
declare -a FORWARD=(`grep -E "${FORWARD_NAMING}" "${SAMPLE_INFO}"`)
declare -a REVERSE=(`grep -E "${REVERSE_NAMING}" "${SAMPLE_INFO}"`)

#   Create an array of sample names
declare -a SAMPLE_NAMES=()
counter=0
for s in "${FORWARD[@]}"
do
    SAMPLE_NAMES[`echo "$counter"`]="`basename $s ${FORWARD_NAMING}`"
    let "counter += 1"
done

#   A function to perform the trimming and plotting
#       Adapted from Tom Kono and Peter Morrell
function trimAutoplot() {
    #   Set the arguments for trimming
    sampleName="$1" #   Name of the sample
    forward="$2" #  Forward file
    reverse="$3" #  Reverse file
    out="$4"/"${sampleName}" #  Outdirectory
    adapters="$5" # Adapter file
    prior="$6" #    Prior informaiotn
    threshold="$7" #    Threshold Value
    encoding="$8" # Platform for sequencing
    seqHand="$9" #  The sequence_handling directory
    if [[ -d "${seqHand}"/Helper_Scripts ]] #  Check to see if helper scripts directory exists
    then
        helper="${seqHand}"/Helper_Scripts #    The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts!"
        exit 1
    fi
    #   Make the out directories
    stats="${out}"/stats
    plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Trim the sequences based on quality and adapters
    sickle pe -t "${encoding}" -q "${threshold}" \
        -f <(seqqs -e -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R1 "${forward}" | scythe -a "${adapters}" -p "${prior}" - 2> "${stats}"/"${sampleName}"_R1_scythe.stderr) \
        -r <(seqqs -e -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R2 "${reverse}" | scythe -a "${adapters}" -p "${prior}" - 2> "${stats}"/"${sampleName}"_R2_scythe.stderr) \
        -o >(seqqs -e -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R1 - | gzip > "${out}"/"${sampleName}"_R1_trimmed.fq.gz) \
        -p >(seqqs -e -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R2 - | gzip > "${out}"/"${sampleName}"_R2_trimmed.fq.gz) \
        -s >(seqqs -e -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_singles - | gzip > "${out}"/"${sampleName}"_singles_trimmed.fq.gz) > "${stats}"/"${sampleName}"_sickle.stderr
    #   Fix the quality scores
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R1_qual.txt 33
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R2_qual.txt 33
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R1_qual.txt 33
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R2_qual.txt 33
    #   Make the plots
    Rscript "${helper}"/plot_seqqs.R "${stats}" "${sample_name}"
}

#   Export the funciton to be used by Parallel
export -f trimAutoplot

#   Check to make sure Sickle, Scythe, and Seqqs are installed
if ! `command -v sickle > /dev/null 2> /dev/null` || ! `command -v seqqs > /dev/null 2> /dev/null` || ! `command -v scythe > /dev/null 2> /dev/null`
then
    echo "Cannot find Sickle, Scythe, and Seqqs!"
    echo "Please install and place in your PATH!"
    exit 1
fi

#   Run the job in parallel
parallel --xapply trimAutoplot {1} {2} {3} ${OUT} ${ADAPTERS} ${PRIOR} ${THRESHOLD} ${PLATFORM} ${SEQUENCE_HANDLING} ::: ${SAMPLE_NAMES[@]} ::: ${FORWARD[@]} ::: ${REVERSE[@]}

#   Create a list of outfiles to be used by read_mapping_start.sh
find ${OUT} -regex ".*_R[1-2]_trimmed.fq.gz" | sort > ${OUT}/"${PROJECT}"_samples_trimmed.txt
echo "List for read_mapping_start.sh can be found at"
echo "${OUT}"/"${PROJECT}"_samples_trimmed.txt
