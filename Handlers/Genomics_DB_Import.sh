#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genomics_DB_Import_Dependencies=(java)

# Note with a large sample size intervals might need to be further subdivided to avoid memory problems
function GenomicsDBImport() {
    local sample_list="$1" # What is our sample list?
    local out_dir="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    local type="$4" # 'WGS' or 'targeted'
    local intvlFile="$5" # put NA for WGS, intervals for interval of targeted sequencing
    local tmp="$6" # temp directory
    local memory="$7" # How much memory can java use?
    local parallelize="$8" # Are we parallelizing across regions?
    # Get the directory for the outputFile, and create it if it's missing
    #local outDir=$(dirname "${out_dir}")
    # Check if out and temp dirs exists, if not make it
    mkdir -p "${out_dir}/Genotype_GVCFs" "${out_dir}/Genotype_GVCFs/combinedDB" "${tmp}"
    # Make sure reference exists
    if ! [[ -s "${reference}" ]]; then
        echo "Cannot find readable reference genome, exiting..." >&2
        exit 31
    fi
    # Note: -Xmx value should be less than total amount of physical memory by at least
    # a few GB. The native TileDB library requires additional memory on top of the Java memory.
    # So, we will subtract a few GB of mem from user provided memory
    mem_num=$(basename ${memory} G)
    new_mem_num=$[${mem_num} - 3]
    mem=$(printf ${new_mem_num}G)
    # Create a file which list the chromosomes for -L option
    # In case of targeted sequences, GenomicDBImport doesn't work if hundreds of intervals are given
    # So, each region in the custom intervals list will be submitted as its own job
    if [[ "${type}" == "targeted" ]]; then
        # Make sure file exists
        if ! [[ -s "${intvlFile}" ]]; then
            echo "Cannot find readable intervals file for the target region, exiting..." >&2
            exit 31
        fi
        cut -f 1 "${intvlFile}" | sort -V | uniq > "${out_dir}/Genotype_GVCFs/intervals.list"
    else
        local dict="${reference%.*}.dict" # replace suffix (.fa or .fasta) with .dict
        # checkDict and creatDict must have been called in sequence_handling, but just checking again
        if ! [[ -s "${dict}" ]]; then echo "Cannot find readable reference dict genome (or bed file), exiting..." >&2; exit 31; fi # Make sure it exists
        chrom_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
        printf '%s\n' "${chrom_list[@]}" > "${out_dir}/Genotype_GVCFs/intervals.list"
    fi
    # Check if we have > 500 intervals and are NOT parallelizing across regions
    if [ $(cat "${out_dir}/Genotype_GVCFs/intervals.list" | wc -l) -gt 500 ] && [ "${parallelize}" == "false" ]; then
        # When there are many small intervals (e.g exomes), following option increases performance.
        local mergeIntvl="--merge-input-intervals"
    else
        local mergeIntvl=""
    fi
    #   Put the sample list into array format
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}"))
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
        # Check if we are parallelizing across regions
        if [ "${parallelize}" == "true" ]; then
            echo "Parallelizing across regions."
            # Store list of custom intervals in an array
            intvl_arr=($(cat "${out_dir}/Genotype_GVCFs/intervals.list"))
            # What interval are we working on currently?
            local current_intvl="${intvl_arr[${PBS_ARRAYID}]}"
            # GATK 4 will throw an error when trying to make workspace if one already exists
            # So, check if directory exists, if so remove before running GenomicsDBImport
            # to make sure we are starting with a clean slate
            if [ -d "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${current_intvl}" ]; then
                echo "Directory for current interval exists, remove before proceeding." >&2
                rm -rf "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${current_intvl}"
            fi
            # Make output directory
            #mkdir -p "${out_dir}/gendb_wksp_${current_intvl}"
            set -x
            gatk --java-options "-Xmx${mem}" \
                GenomicsDBImport \
                -R "${reference}" \
                "${GATK_IN[@]}" \
                -L "${current_intvl}" \
                "${tmp}" \
                --genomicsdb-workspace-path "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${current_intvl}"
            set +x
        else
            echo "Interval list is <500 and we are NOT parallelizing across regions."
            set -x
            gatk --java-options "-Xmx${mem}" \
                GenomicsDBImport \
                -R "${reference}" \
                "${GATK_IN[@]}" \
                -L "${out_dir}/Genotype_GVCFs/intervals.list" \
                "${tmp}" \
                --genomicsdb-workspace-path "${out_dir}/Genotype_GVCFs/combinedDB"
            set +x
        fi
    else
        echo "Interval list is >500, run with --merge-input-intervals flag."
        set -x
        gatk --java-options "-Xmx${mem}" \
            GenomicsDBImport \
            -R "${reference}" \
            "${GATK_IN[@]}" \
            -L "${out_dir}/Genotype_GVCFs/intervals.list" \
            "${mergeIntvl}" \
            "${tmp}" \
            --genomicsdb-workspace-path "${out_dir}/Genotype_GVCFs/combinedDB"
        set +x
    fi
}

export -f GenomicsDBImport

#   A function to combine individual GVCFs into a single GVCF required for GATK4
# This is slow, so use GenomicsDBImport().
# keeping it here in case there is a situation when GenomicsDBImport doesn't work
function Combine_GVCFs() {
    local sample_list="$1" # What is our sample list?
    local out_dir="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
	    GATK_IN+=(-V $s)
    done
    # get the directory for the outputFile, and create it if it's missing
    #local out_dir=$(dirname "${out_dir}")
    mkdir -p "${out_dir}/Genotype_GVCFs" "${out_dir}/Genotype_GVCFs/combinedDB"
    set -x; gatk CombineGVCFs \
		 -R "${reference}"\
		 "${GATK_IN[@]}" \
		 -O "${out_dir}/Genotype_GVCFs/combinedDB"
}

export -f Combine_GVCFs
