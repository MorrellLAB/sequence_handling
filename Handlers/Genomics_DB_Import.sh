#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genomics_DB_Import_Dependencies=(java)

# Note with a large sample size intervals might need to be further subdivided to avoid memory problems
function GenomicsDBImport() {
    local sample_list="$1" # What is our sample list?
    local outFileName="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    local type="$4"  # 'WGS' or 'targeted'
    local intvlBED="$5"  # put NA for WGS, BED file name for interval of tageted sequencing
    local tmp="$6"    # temp directory
    # get the directory for the outputFile, and create it if it's missing
    local outDir=$(dirname "${outFileName}")
    local memory="$7" # How much memory can java use?
    mkdir -p "${outDir}"
    if ! [[ -s "${reference}" ]]; then echo "Cannot find readable reference genome, exiting..." >&2; exit 31; fi # Make sure it exists
    # create a file which list the chromosomes for -L option
    # In case of targeted sequences, GenomicDBImport doesn't work if hundreds of intervals are given
    # So, using the chromosomes which contain the targets are used.
    if [[ "${type}" == "targeted" ]]; then
        if ! [[ -s "${intvlBED}" ]]; then echo "Cannot find readable bed file for the target region, exiting..." >&2; exit 31; fi # Make sure it exists
        cut -f 1 "${intvlBED}" | sort | uniq > "${outDir}/intervals.list"
    else
        local dict="${reference%.*}.dict" # replace suffix (.fa or .fasta) with .dict
        # checkDict and creatDict must have been called in sequence_handling, but just checking again
        if ! [[ -s "${dict}" ]]; then echo "Cannot find readable reference dict genome (or bed file), exiting..." >&2; exit 31; fi # Make sure it exists
        chrom_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
        printf '%s\n' "${chrom_list[@]}" > "${outDir}/intervals.list"
    fi

    if [ $(cat "${outDir}/intervals.list" | wc -l) -gt 500 ]; then
        # When there are many small intervals (e.g exomes), following option increases performance.
        local mergeIntvl="--merge-input-intervals"
    else
        local mergeIntvl=""
    fi
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
	    GATK_IN+=(-V $s)
    done
    if [ -n "$tmp" ] ; then
	    tmp="--tmp-dir=${tmp}"
    fi
    # Check if mergeIntvl is an empty string, if so run without --merge-input-intervals flag
    if [ -z "${mergeIntvl}" ]; then
        echo "Interval list is <500."
        set -x
        gatk --java-options "-Xmx${memory}" \
            GenomicsDBImport \
            -R "${reference}" \
		    "${GATK_IN[@]}" \
		    -L "${outDir}/intervals.list" \
		    "${tmp}" \
		    --genomicsdb-workspace-path $outFileName
        set +x
    else
        echo "Interval list is >500, run with --merge-input-intervals flag."
        set -x
        gatk --java-options "-Xmx${memory}" \
            GenomicsDBImport \
            -R "${reference}" \
            "${GATK_IN[@]}" \
            -L "${outDir}/intervals.list" \
            "${mergeIntvl}" \
            "${tmp}" \
            --genomicsdb-workspace-path $outFileName
        set +x
    fi
    rm -f "${outDir}/intervals.list"
}

export -f GenomicsDBImport

#   A function to combine individual GVCFs into a single GVCF required for GATK4
# This is slow, so use GenomicsDBImport().
# keeping it here in case there is a situation when GenomicsDBImport doesn't work
function Combine_GVCFs() {
    local sample_list="$1" # What is our sample list?
    local outFileName="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
	    GATK_IN+=(-V $s)
    done
    # get the directory for the outputFile, and create it if it's missing
    local outDir=$(dirname "${outFileName}")
    mkdir -p "${outDir}"
    set -x; gatk CombineGVCFs \
		 -R "${reference}"\
		 "${GATK_IN[@]}" \
		 -O $outFileName
}

export -f Combine_GVCFs
