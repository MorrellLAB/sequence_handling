#!/bin/bash

set -e
set -u
set -o pipefail

#   This script generates a series of QSub submissions for read mapping
#   The Burrows-Wheeler Aligner (BWA) and the Portable Batch System (PBS)
#   are required to use this script

usage() {
    echo -e "\
Usage: ./read_mapping_start.sh map scratch ref_gen sample_info project platform email
where:  scratch is the output directory for the read mapping \n\
\n\
        ref_gen is the reference genome for the read mapping \n\
\n\
        sample_info is the list of FASTQ files to be read mapped \n\
\n\
        project is the name of the project or capture facility \n\
\n\
        platform is the sequencing platform used for your samples \n\
\n\
        email is the email address at which you can be notified of the progress of the read mapping \n\
\n\
This uses the Burrows-Wheeler Aligner's (BWA) 'mem' algorithim with the following settings: \n\
    -t 8: Use 8 threads \n\
    -k 10: Seed length 25 \n\
    -r 1.0: Re-seed if match is greater than 1.0 * seed length \n\
    -M : mark split hits as secondary \n\
    -T 85: only output alignments greater than score 85 \n\
    -O 8: Gap open penalty \n\
    -E 1: Gap extension penalty \n\
\n\
By default, this script looks for forward and reverse samples starting ending with: \n\
    '_R1_trimmed.fq.gz' for forward and \n\
    '_R2_trimmed.fq.gz' for reverse \n\
\n\
This will generate a series of QSub submissions with the following settings: \n\
    8 GB memory \n\
    1 node \n\
    8 processor per node \n\
    16 hours walltime \n\
\n\
To change these default settings, open this file in your favourite text editor \n\
and edit lines 67, 69, 73, and 74 \n\
-------------------------------------------------------------------------------\n\
Read mapping requires an index of the reference genome to be made \n\
To do this, run: \n\
       ./read_mapping_start.sh index ref_gen email \n\
where:  ref_gen is the reference genome to be indexed \n
\n\
        email is the email address at which you can be notified of the progress of the read mapping \n\
\n\
This index process will be done in the same directory \n\
where the reference genome is stored, please make sure you have \n\
appropriate permissions for said directory\n\
" >&2
    exit 1
}

if [ "$#" -lt 1 ]; then
    usage;
fi

QUE_SETTINGS='-l mem=8gb,nodes=1:ppn=8,walltime=16:00:00'

case "$1" in
    "map" )
        #   Make sure there are enough arguments
        if [ "$#" -lt 7 ]; then
            usage;
        fi
        #   Assign variables for running BWA mem
        SCRATCH="$2"
        REF_GEN="$3"
        SAMPLE_INFO="$4"
        PROJ="$5"
        PLAT="$6"
        EMAIL="$7"
        #   BWA mem settings
        SETTINGS='-t 8 -k 10 -r 1.0 -M -T 85 -O 8 -E 1'
        #   Today's date
        YMD=`date +%Y-%m-%d`
        #   Create scratch directory if it doesn't exist
        mkdir -p "${SCRATCH}"/"${PROJ}"
        cd "${SCRATCH}"/"${PROJ}"
        #   Generate lists of forward and reverse reads that match
        FWD="_R1_trimmed.fq.gz"
        REV="_R2_trimmed.fq.gz"
        grep -E "$FWD" "${SAMPLE_INFO}" | sort > "${SCRATCH}"/fwd.txt
        FWD_FILE="${SCRATCH}"/fwd.txt
        grep -E "$REV" "${SAMPLE_INFO}" | sort > "${SCRATCH}"/rev.txt
        REV_FILE="${SCRATCH}"/rev.txt
        #   Check for equal numbers of forward and reverse reads
        if [ `wc -l < "$FWD_FILE"` = `wc -l < "$REV_FILE"` ]; then
            echo "Equal numbers of forwad and reverse reads"
        else
            exit 1
        fi
        #   Start a series of QSub submissions to run BWA
        for i in `seq $(wc -l < "$FWD_FILE")`
        do
            f=`head -"$i" "$FWD_FILE" | tail -1`
            r=`head -"$i" "$REV_FILE" | tail -1`
            s=`basename "$f" "$FWD"`
            RG="@RG\tID:$s\tLB:${PROJ}_$s\tPL:${PLAT}\tPU:$s\tSM:$s"
            echo -e "'"$RG"'" > RG.txt
            echo "module load bwa && bwa mem ${SETTINGS} -R `cat RG.txt` ${REF_GEN} ${f} ${r} > ${SCRATCH}/${s}_${YMD}.sam" | qsub "${QUE_SETTINGS}" -m abe -M "${EMAIL}" -N "$s"_Read_Mapping
        done
        ;;
    "index" )
        #   Make sure there are enough arguments
        if [ "$#" -lt 3]; then
            usage;
        fi
        #   Assign variables for running BWA index
        REF_GEN="$2"
        EMAIL="$3"
        s=`basename "${REF_GEN}"`
        #   Generate a QSub submission for BWA index
        echo "module load && bwa index ${REF_GEN}" | qsub "${QUE_SETTINGS}" -m abe -M "${EMAIL}" -N "$s"_Index
        ;;
    * )
        #   Anythin other than 'map' or 'index'
        usage
        ;;
esac
