#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_GenotypeGVCFs.job

set -e
set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genotype_GVCFs_Dependencies=(java realpath sort perl cut uniq printf)

#   A function to genotype the GVCFs without PBS
function Genotype_GVCFs() {
    local sample_list="$1" # Sample list (GATK 3). Note for GATK4 path to gendb workspace will be created automatically
    local out_dir="$2" # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local ploidy="$6" # What is the sample ploidy?
    local memory="$7" # How much memory can java use?
    local ref_dict="$8" # Where is the reference dictionary?
    local maxarray="$9" # What is the maximum array index?
    local gatkVer="${10}" # GATK version 3 or 4
    local intervals="${11}" # put NA for WGS, intervals for interval of targeted sequencing
    local parallelize="${12}" # Are we parallelizing across regions?
    local type="${13}" # 'WGS' or 'targeted'
    local scaffolds="${14}" # Any additional regions not covered by chromosomes
    local gg_combined_vcf="${15}" # Do we want to combine split vcfs into one vcf?
    local project="${16}"
    local tmp="${17}"
    # Making sure out_dir is absolute path
    # because of CL's note from Sept 23 2019 (gatk 4.1.2 can't use rel. or abs. path to gendb).
    # Naoki didn't have any problem with rel. path to gendb (gatk 4.1.7).
    # If this issue is resolved (i.e. no need to cd to the directory with DB),
    # no need for this conversion step with realpath (and nicer without it)
    out_dir="$(realpath "${out_dir}")"
    reference="$(realpath --no-symlinks "${reference}")"
    # temporary file to store intervals
    if [[ "${USE_PBS}" == "true" ]]; then
        ## With PBS, multiple processes might write to this file,
        ## so making unique file for each process.  Is PBS_JOBID better? Naoki
        local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals-${PBS_ARRAYID}.list")
    elif [[ "${USE_SLURM}" == "true" ]]; then
        local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals-${SLURM_ARRAY_TASK_ID}.list")
    else
	    local intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals.list")
    fi
    # There might be a left over from Genomics_DB_Import, so clean it up
    rm -f "${intervals_filepath}"
    
    # set up temporary intervals file
    if [[ "${type}" == "targeted" ]]; then
        if [[ "${intervals}" == *.bed ]]; then
	        # convert .bed to .list
	        sort -k1,1 -k2,2n "${intervals}" | \
		    perl -ne 'chomp; @a=split/\t/;print $a[0], ":", $a[1]+1,"-", $a[2], "\n"' > "${intervals_filepath}"
	    else
	        # .list/.intervals
            cut -f 1 "${intervals}" | sort -V | uniq > "${intervals_filepath}"
	    fi
    else
        # If analysisType="WGS", then use reference dict to create intervals list using chr
        # Make an array of chromosome part names
        intvl_arr=($(cut -f 2 ${ref_dict} | grep -E '^SN' | cut -f 2 -d ':'))
        printf '%s\n' "${chrom_list[@]}" > "${intervals_filepath}"
    fi

    # For GATK3, input files are bunch of vcf instead of gendb
    if [[ "$gatkVer" == 3 ]]; then
        # Put the sample list into array format
        declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}"))
	    # Put the samples into a format that GATK can read, becomes a text string of "-V file1.vcf -V file2.vcf ..."
	    local input_all_sample_vcf=$(printf -- '-V %s ' "${sample_array[@]}")
    fi

    # set up arrays:
    #    intvl_arr: stores values to -L option
    #    out_name_arr: values for -O/-o option
    #    input_opt_arr: values for input: e.g. "-V gendb://file" for GATK4 or "-V in-1.vcf -V in-2.vcf ..." for GATK3
    # Each element of the array is for different intervals
    if [[ "${parallelize}" == "true" ]]; then
        # Read in line by line and store in array called intvl_arr
        readarray -t intvl_arr < "${intervals_filepath}"
        # Prepare list of output names
	    #  example: "/out_dir_path/Genotype_GVCFs/vcf_split_regions/chr1:100-1000.vcf"
	    # about substitution: https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html
	    out_name_arr=( "${intvl_arr[@]/%/.vcf}" ) # append ".vcf"
	    # prepend with output directory
	    out_name_arr=( "${out_name_arr[@]/#/${out_dir}/Genotype_GVCFs/vcf_split_regions/}" )
	    # For parallelized WGS, output was in Genotype_GVCFs instead of in Genotype_GVCFs/vcf_split_regions/
	    # But it is in vcf_split_regions now.  I hope it's ok, Naoki
        if [[ "$gatkVer" == 3 ]]; then
            # Each element of this array becomes ${input_all_sample_vcf}
            input_opt_arr=( "${intvl_arr[@]/*/${input_all_sample_vcf}}" )
        else
            # prepend each interval with "-V gendb:..."
            input_opt_arr=( "${intvl_arr[@]/#/-V gendb://${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_}" )
        fi
    else
        # parallelize = false and we have a single gendb workspace
	    intvl_arr=("${intervals_filepath}")	 # needed for no PBS
        out_name_arr=($(echo "${out_dir}/Genotype_GVCFs/raw_variants.vcf"))
        if [[ "$gatkVer" == 3 ]]; then
            input_opt_arr=( "${input_all_sample_vcf}" )
        else
            input_opt_arr=( "-V gendb://${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp" )
        fi
    fi

    # If we have scaffolds or regions not covered by chromosomes,
    # append scaffolds list to array
    if [[ "${scaffolds}" != "false" ]]; then
	    scaffolds="$(realpath "${scaffolds}")"
        if [[ "${parallelize}" == "true" ]]; then
            intvl_arr+=("${scaffolds}")
            out_name_arr+=("${out_dir}/Genotype_GVCFs/vcf_split_regions/additional_intervals.vcf")
            if [[ "$gatkVer" == 3 ]]; then
                input_opt_arr+=("${input_all_sample_vcf}")
            else
                input_opt_arr+=("-V gendb://${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_additional_intervals")
            fi
        fi
	# When not parallelized, "-L ${scaffolds}" is attached later
    fi

    # Sanity check. Actually, if Config is changed between Haplotype_Caller,
    # Genomics_DB_Import and Genotype_GVCFs, this handler can fail.
    # It isn't currently checking this kind of error.
    if [ ${#intvl_arr[@]} -ne $[${maxarray} + 1] ]; then
        echo "Check NUM_CHR and CUSTOM_INTERVAL, it doesn't match with entries in the reference" >&2
        exit 1
    fi

    # Check GATK version and decide how to format flags and sample lists
    if [[ "$gatkVer" == 3 ]]; then
        analysisTypeOpt="-T GenotypeGVCFs"
        outFlag="-o"
        ploidyFlag="--sample_ploidy $ploidy"
    else
        analysisTypeOpt="GenotypeGVCFs"
        outFlag="-O"
        # GATK 4 can handle any ploidy (or mixed ploidies) intelligently, so not setting this option
        ploidyFlag=""
    fi

    #   Make sure the out directory exists
    if [[ "${parallelize}" == "true" ]]; then
        mkdir -p "${out_dir}/Genotype_GVCFs/vcf_split_regions"
    else
	    mkdir -p "${out_dir}/Genotype_GVCFs"
    fi

    if [[ "$USE_PBS" == "true" ]] || [[ "${USE_SLURM}" == "true" ]]; then
        if [[ "$USE_PBS" == "true" ]]; then
            if [[ "${parallelize}" == "true" ]]; then
                local arr_index=${PBS_ARRAYID}
            else
                local arr_index=0
            fi
        elif [[ "${USE_SLURM}" == true ]]; then
            if [[ "${parallelize}" == "true" ]]; then
                local arr_index=${SLURM_ARRAY_TASK_ID}
            else
                local arr_index=0
            fi
        fi
        #   What region of the genome are we working on currently?
        local current_intvl="${intvl_arr[${arr_index}]}"
        local current_input="${input_opt_arr[${arr_index}]}"
        local current_output="${out_name_arr[${arr_index}]}"
	    # -L argument(s)
	    local intvl_args=( '-L' "${current_intvl}" )
        if [[ "${parallelize}" != "true" ]] && [[ "${scaffolds}" != "false" ]]; then
            # not parallelized and there is a custom scaffold deal with scaffold
            intvl_args+=( '-L' "${scaffolds}" )
        fi

        if [[ "$gatkVer" == 3 ]]; then
            #   Run GATK using the parameters given
            java -Xmx"${memory}" -jar "${gatk}" \
                 "${analysisTypeOpt}" \
                 -R "${reference}" \
                 "${intvl_args[@]}" \
                 ${current_input} \
                 --heterozygosity "${heterozygosity}" \
                 $(printf -- "%s" "${ploidyFlag}") \
                 ${outFlag} "${current_output}"
        else  # GATK4
            #   As of Sep 23, 2019 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
            #   CL is unable to get it working with a relative or absolute filepath to the database
            #   Go into Genotype_GVCFs directory
	        # take basename after removing "-V gendb:" from the begining of the string
            if [[ "${USE_PBS}" == "true" ]]; then
                local gdb_workspace=$(basename "${input_opt_arr[${PBS_ARRAYID}]/#-V gendb:/}")
            elif [[ "${USE_SLURM}" == "true" ]]; then
                local gdb_workspace=$(basename "${input_opt_arr[${SLURM_ARRAY_TASK_ID}]/#-V gendb:/}")
            fi
            current_input="-V gendb://${gdb_workspace}"
            cd "${out_dir}/Genotype_GVCFs/combinedDB"
            if [[ -z "${tmp}" ]]; then
                # No tmp directory specified
                set -x
                gatk --java-options "-Xmx${memory}" \
                    "${analysisTypeOpt}" \
                    -R "${reference}" \
                    "${intvl_args[@]}" \
                    ${current_input} \
                    --heterozygosity "${heterozygosity}" \
                    $(printf -- "%s" "${ploidyFlag}") \
                    ${outFlag} "${current_output}"
                set +x
            else
                # tmp directory is specified
                set -x
                gatk --java-options "-Xmx${memory}" \
                    "${analysisTypeOpt}" \
                    -R "${reference}" \
                    "${intvl_args[@]}" \
                    ${current_input} \
                    --heterozygosity "${heterozygosity}" \
                    $(printf -- "%s" "${ploidyFlag}") \
                    ${outFlag} "${current_output}" \
                    --tmp-dir ${tmp}
                set +x
            fi
        fi
    else  # No PBS
	    if [[ -z "${GG_THREADS}" ]]; then
            GG_THREADS='100%'   # Use all cores if not defined
        fi
	    # Using arrays to feed parallel can cause "Argument list too long", so making input file
	    local temp_input_opt_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_input_opt.txt"
	    local temp_intvl_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_intvl.txt"
	    local temp_out_name_filepath="${out_dir}/Genotype_GVCFs/combinedDB/temp_out_name.txt"
	    printf '%s\n' "${input_opt_arr[@]}" > "${temp_input_opt_filepath}"
	    printf '%s\n' "${intvl_arr[@]}" > "${temp_intvl_filepath}"
	    printf '%s\n' "${out_name_arr[@]}" > "${temp_out_name_filepath}"

        parallel --jobs ${GG_THREADS} java -Xmx"${memory}" -jar "${gatk}" \
            "${analysisTypeOpt}" \
            -R "${reference}" \
            -L {1} \
            $( [[ "${parallelize}" != "true" ]] && [[ "${scaffolds}" != "false" ]] && printf -- '%s' "-L \"${scaffolds}\"" ) \
            '$(printf -- %s {2})' \
            --heterozygosity "${heterozygosity}" \
            $(printf -- "%s" "${ploidyFlag}") \
            ${outFlag} "{3}" \
            :::: "${temp_intvl_filepath}" ::::+ "${temp_input_opt_filepath}" ::::+ "${temp_out_name_filepath}"
            #::: "${intvl_arr[@]}" :::+ "${input_opt_arr[@]}" :::+ "${out_name_arr[@]}"
        rm -f "${temp_input_opt_filepath}" "${temp_intvl_filepath}" "${temp_out_name_filepath}"
    fi

    echo "Clean up ${intervals_filepath}" >&2
    rm -f "${intervals_filepath}"  # clean up the temp. intervals file
}

#   Export the function
export -f Genotype_GVCFs
