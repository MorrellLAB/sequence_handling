#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

set -e
set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genomics_DB_Import_Dependencies=(java mkdir echo printf sort perl cut wc uniq)

# Note with a large sample size intervals might need to be further subdivided to avoid memory problems
function GenomicsDBImport() {
    local sample_list="$1" # What is our sample list?
    local out_dir="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    local type="$4" # 'WGS', 'targeted' or 'targetedHC'
    local intvlFile="$5" # put NA for WGS, intervals for interval of targeted sequencing
    local scaffolds="$6" # list of scaffolds or sequences not covered by chromosomes
    local memory="$7" # How much memory can java use?
    local parallelize="$8" # Are we parallelizing across regions?
    local tmp="$9" # temp directory
    # Check if out and temp dirs exists, if not make it
    mkdir -p "${out_dir}/Genotype_GVCFs" \
        "${out_dir}/Genotype_GVCFs/combinedDB"
    # Check if user has specified tmp dir
    if [ -z "${tmp}" ]; then
        echo "tmp veriable is empty, proceed without using tmp dir."
    else
        echo "Making tmp directory."
        mkdir -p "${tmp}"
    fi
    # Make sure reference exists
    if ! [[ -s "${reference}" ]]; then
        echo "Cannot find readable reference genome, exiting..." >&2
        exit 31
    fi
    # Note: -Xmx value should be less than total amount of physical memory by at least
    # a few GB. The native TileDB library requires additional memory on top of the Java memory.
    # So, we will subtract a few GB of mem from user provided memory
    mem_num=$(basename ${memory} g)
    # Check that we have requested the minimum required memory (>4g)
    if [ ${mem_num} -le 4 ]; then
        echo "You have 4g or less requested memory, please request more memory for this handler to work. Exiting..."
        exit 31
    fi
    new_mem_num="$((${mem_num}-4))"
    mem=$(printf "${new_mem_num}g")
    # Create a file which list the intervals/chromosomes for -L option
    # Each region in the custom intervals list will be submitted as its own job
    if [[ "${USE_PBS}" == "true" ]]; then
        # With PBS, multiple processes might write to this file, so
        # making unique file for each job.  Is PBS_JOBID better?
        local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals-${PBS_ARRAYID}.list")
    elif [[ "${USE_SLURM}" == true ]]; then
        local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals-${SLURM_ARRAY_TASK_ID}.list")
    else
	    local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals.list")
    fi
    if [[ "${type}" == "targeted" ]] || [[ "${type}" == "targeted-HC" ]] ; then
        # Make sure file exists
        if ! [[ -s "${intvlFile}" ]]; then
            echo "Cannot find readable intervals file for the target region, exiting..." >&2
            exit 31
        fi

        # Check if we have a .bed file or .intervals/.list file format
        if [[ "${intvlFile}" == *.bed ]]; then
            # If we have a .bed file, sort, convert to .list and store path in variable
            sort -k1,1 -k2,2n "${intvlFile}" | \
            perl -ne 'chomp; @a=split/\t/;print $a[0], ":", $a[1]+1,"-", $a[2], "\n"' > "${intervals_filepath}"
            # With the perl 1-liner
            # each line of *.bed becomes firstColumn:(2nd column+1)-(3rd column)
            # NOTE: https://software.broadinstitute.org/gatk/documentation/article?id=1319
            #   GATK interval uses 1-based (first position in the genome is position 1, not position 0).
            #   But bed file has the start position 0-based, and end position 1-based, so 1 is added to 2nd column

            # Naoki Comment: conversion of bed to GATK intervals list seems to make things simpler.
            # So I'm doing the conversion in Haplotype_Caller, Genomics_DB_Import and Genotype_GVCFs.
            #    With the previous implementation, it becomes weird when bed file contains
            #    more than 3 columns.
        else
	        if [ $(grep ":0-" "${intvlFile}" | wc -l) != 0 ]; then
                # For GATK-style .list or .intervals file in the form <chr>:<start>-<stop>
                # GATK 4 will throw an error if the first position starts at 0.
                echo "GATK-style .list or .intervals file uses 1-based coordinates. Please check your custom intervals file to make sure there are no positions that start at 0. If so, change so the starting position is 1. See docs: https://software.broadinstitute.org/gatk/documentation/article?id=11009. Exiting..." >&2
                exit 31
            fi
            cut -f 1 "${intvlFile}" | sort -V | uniq > "${intervals_filepath}"
        fi
    else
        # If analysisType="WGS", then use reference dict to create intervals list using chr
        local dict="${reference%.*}.dict" # replace suffix (.fa or .fasta) with .dict
        # checkDict and creatDict must have been called in sequence_handling, but just checking again
        if ! [[ -s "${dict}" ]]; then
            echo "Cannot find readable reference dict genome (or bed file), exiting..." >&2; exit 31
        fi # Make sure it exists
        # Check if we have scaffolds
        if [[ "${scaffolds}" != "false" ]]; then
            # We need to temporarily remove the scaffolds from our chromosome list (for compatibility with
            #   various cases, these will be added back at a later step)
            chrom_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':' | grep -vf ${scaffolds})) # Make an array of chromosome part names
        else
            # Assume dict doesn't have scaffolds
            chrom_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
        fi
        printf '%s\n' "${chrom_list[@]}" > "${intervals_filepath}"
    fi

    # Check if we have > 500 intervals and are NOT parallelizing across regions
    if [ $(cat "${intervals_filepath}" | wc -l) -gt 500 ] && [ "${parallelize}" == "false" ]; then
        # When there are many small intervals (e.g exomes), following option increases performance.
        local mergeIntvl="--merge-input-intervals"
	    echo "Interval list is >500, run with --merge-input-intervals flag."
    else
        local mergeIntvl=""
    fi

    #   Put the sample list into array format
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}"))
    #   Put the samples into a format that GATK can read, becomes a text string of "-V file1 -V file2 ..."
    local input_all_sample_vcf=$(printf -- "-V %s " "${sample_array[@]}")

    # Temp directory option
#    if [ -n "$tmp" ] ; then
#        # Length of string is non-zero
#        tmp="--tmp-dir ${tmp}"
#    fi

    # Check if we are parallelizing across regions
    if [ "${parallelize}" == "true" ]; then
        # The following will work with type="targeted", "targeted-HC", or "WGS"
        # since intervals_filepath is correctly set above.
        if [[ "${intvlFile}" == "NA" ]]; then
            echo "Parallelizing across chromosomes."
        else
            echo "Parallelizing across regions."
        fi

        # Store list of custom intervals in an array
        intvl_arr=($(cat "${intervals_filepath}"))
        # Prepare list of output names
        out_name_arr=("${intvl_arr[@]}")

        # If we have scaffolds or sequences not covered by the chromosomes.
        # append scaffolds list to array
        if [[ "${scaffolds}" != "false" ]]; then
            intvl_arr+=("${scaffolds}")
            out_name_arr+=("additional_intervals")
        fi
        # Input vcf
        #   For targeted-HC, we select the subset of vcf corresponding to the interval.
        #   For others, $input_all_sample_vcf works for any intervals.
        input_vcf_arr=()
        for ivl in "${out_name_arr[@]}" # going through each interval
        do
            if [[ "${type}" == "targeted-HC" ]]; then
                # the naming scheme of input vcf is sampleID_interval_RawGLs.g.
                # select subset of the input vcf files corresponding to this interval
                input_vcf_arr+=( "$(printf -- "-V %s\n" "${sample_array[@]}" | grep ${ivl} | tr '\n' ' ')" )
            else  # "targeted" or "WGS"
                input_vcf_arr+=( "${input_all_sample_vcf}" )
            fi
        done
        if [[ "${USE_PBS}" == "true" ]] || [[ "${USE_SLURM}" == true ]]; then
            # What interval are we working on currently?
            if [[ "$USE_PBS" == true ]]; then
                local current_intvl="${intvl_arr[${PBS_ARRAYID}]}"
                local current_output_dirname="${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${out_name_arr[${PBS_ARRAYID}]}"
                local current_input_vcf="${input_vcf_arr[${PBS_ARRAYID}]}"
            elif [[ "${USE_SLURM}" == true ]]; then
                local current_intvl="${intvl_arr[${SLURM_ARRAY_TASK_ID}]}"
                local current_output_dirname="${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${out_name_arr[${SLURM_ARRAY_TASK_ID}]}"
                local current_input_vcf="${input_vcf_arr[${SLURM_ARRAY_TASK_ID}]}"
            fi

            # GATK 4 will throw an error when trying to make workspace if one already exists
            # So, check if directory exists, if so remove before running GenomicsDBImport
            # to make sure we are starting with a clean slate
            if [ -d "${current_output_dirname}" ]; then
                echo "Directory for current interval exists, removing before proceeding." >&2
                rm -rf "${current_output_dirname}"
            fi
            set -x
            /panfs/jay/groups/9/morrellp/public/Software/gatk-4.1.8.0/gatk --java-options "-Xmx${mem}" \
                GenomicsDBImport \
                -R "${reference}" \
                $(printf -- '%s ' ${current_input_vcf}) \
                -L "${current_intvl}" \
                --tmp-dir "${tmp}" \
                --genomicsdb-workspace-path "${current_output_dirname}" \
                --overwrite-existing-genomicsdb-workspace true
            set +x
        else  # No PBS
            if [[ -z "${GDBI_THREADS}" ]]; then
                GDBI_THREADS='100%'   # Use all cores if not defined
            fi

            # clean up output workspace
            cleaned_flag=false
            for this_dir in "${out_name_arr[@]}"
            do
                local current_output_dirname="${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${this_dir}"
                if [ -d "${current_output_dirname}" ]; then
                    rm -rf "${current_output_dirname}"
                    cleaned_flag=true
                fi
            done
            if [[ "${cleaned_flag}" == true ]]; then
                echo "Directory for current interval exists, removing before proceeding." >&2
            fi
            # Using arrays to feed parallel can cause "Argument list too long", so making input file
            local temp_input_vcf_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_input_vcf.txt"
            local temp_intvl_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_intvl.txt"
            local temp_out_name_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_out_name.txt"
            printf '%s\n' "${input_vcf_arr[@]}" > "${temp_input_vcf_filepath}"
            printf '%s\n' "${intvl_arr[@]}" > "${temp_intvl_filepath}"
            printf '%s\n' "${out_name_arr[@]}" > "${temp_out_name_filepath}"
            
            set -x
            parallel --jobs ${GDBI_THREADS} /panfs/jay/groups/9/morrellp/public/Software/gatk-4.1.8.0/gatk --java-options "-Xmx${mem}" \
                    GenomicsDBImport \
                    -R "${reference}" \
                    '$(echo {1})' \
                    -L {2} \
                    --tmp-dir "${tmp}" \
                    --overwrite-existing-genomicsdb-workspace true \
                    --genomicsdb-workspace-path "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_{3}" \
                    :::: "${temp_input_vcf_filepath}" ::::+ "${temp_intvl_filepath}" ::::+ "${temp_out_name_filepath}"
                    # ::: "${input_vcf_arr[@]}" :::+ "${intvl_arr[@]}" :::+ "${out_name_arr[@]}"
            set +x

            rm -f "${temp_input_vcf_filepath}" "${temp_intvl_filepath}" "${temp_out_name_filepath}"
        fi
    else
        echo "NOT parallelizing across regions."
        # This would result in a single gendb workspace (NOT parallelizing)
        # Check if we have scaffolds in addition to custom intervals, if so append to list
        if [[ "${scaffolds}" != "false" ]]; then
            cat "${scaffolds}" >> "${intervals_filepath}"
        fi
    # clean up output workspace
        if [ -d "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp" ]; then
            echo "Directory for current interval exists, remove before proceeding." >&2
            rm -rf "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp"
        fi
        set -x
        /panfs/jay/groups/9/morrellp/public/Software/gatk-4.1.8.0/gatk --java-options "-Xmx${mem}" \
             GenomicsDBImport \
             -R "${reference}" \
             $(printf -- '%s ' ${input_all_sample_vcf}) \
             -L "${intervals_filepath}" \
            ${mergeIntvl} \
            --tmp-dir "${tmp}" \
            --genomicsdb-workspace-path "${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp" \
            --overwrite-existing-genomicsdb-workspace true
        set +x
    fi
    echo "Clean up ${intervals_filePath}" >&2
    rm -f "${intervals_filepath}"
}

export -f GenomicsDBImport
