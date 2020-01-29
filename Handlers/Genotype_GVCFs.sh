#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_GenotypeGVCFs.job

set -e
set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genotype_GVCFs_Dependencies=(java)

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
    #set -x # For debugging purposes
    if [[ "${type}" == "targeted" ]]; then
        # Check if we parallelized across regions for GenomicsDBImport
        if [[ "${parallelize}" == "true" ]]; then
            if [[ "${intervals}" == *.bed ]]; then
                # If we have a .bed file, sort and store path in variable
                intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals.bed")
                # Read in line by line and store in array called intvl_arr
                readarray -t intvl_arr < ${intervals_filepath}
                # Prepare list of output names, assume tab delimited bed file
                out_name_arr=($(sed -e 's,\t,_,' -e 's,\t,-,' "${intervals_filepath}"))
            else
                intvl_arr=($(cat "${out_dir}/Genotype_GVCFs/intervals.list"))
                # Prepare list of output names
                out_name_arr=("${intvl_arr[@]}")
            fi
        else
            # parallelize = false and we have a single gendb workspace
            if [[ "${intervals}" == *.bed ]]; then
                # If we have a .bed file, sort and store path in variable
                intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals.bed")
                out_name_arr=($(echo "raw_variants"))
            else
                intervals_filepath=$(echo "${out_dir}/Genotype_GVCFs/intervals.list")
                out_name_arr=($(echo "raw_variants"))
            fi
        fi
    else
        # If analysisType="WGS", then use reference dict to create intervals list using chr
        if [[ "${maxarray}" -ne 0 ]]; then
            # Make an array of chromosome part names
            intvl_arr=($(cut -f 2 ${ref_dict} | grep -E '^SN' | cut -f 2 -d ':'))
        fi
        # A chromosome part name is used as the output name, or "custom_intervals" for the intervals in CUSTOM_INTERVAL
        out_name_arr=("${intvl_arr[@]}")
        if [[ ! -z "${intervals}" ]]; then # custom interval is appended
            intvl_arr+=("${intervals}")
            out_name_arr+=("custom_intervals")
        fi
        if [ ${#intvl_arr[@]} -ne $[${maxarray} + 1] ]; then
            echo "Check NUM_CHR and CUSTOM_INTERVAL, it doesn't match with entries in the reference" >&2
            exit 1
        fi
    fi
    # Check GATK version and decide how to format flags and sample lists
    if [[ "$gatkVer" == 3 ]]; then
        analysisTypeOpt="-T GenotypeGVCFs"
        outFlag="-o"
        ploidyFlag="--sample_ploidy"
        # Put the sample list into array format
        declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}"))
        #   Put the samples into a format that GATK can read
        GATK_IN=()
        for s in "${sample_array[@]}"
        do
            GATK_IN+=(-V $s)
        done
    else
        analysisTypeOpt="GenotypeGVCFs"
        outFlag="-O"
        # note that with GATK 4 the documentation says, no need to specify this
        ploidyFlag="--sample-ploidy"
        # Sample list format can remain as is
        GATK_IN=("${sample_list}")
    fi
    #   Make sure the out directory exists
    mkdir -p "${out_dir}/Genotype_GVCFs"
    if [[ "$USE_PBS" == "true" ]]; then
        if [[ "${type}" == "targeted" ]]; then
            #   What region of the genome are we working on currently?
            local current_intvl="${intvl_arr[${PBS_ARRAYID}]}"
            if [[ "${parallelize}" == "true" ]]; then
                local current_intvl_name="${out_name_arr[${PBS_ARRAYID}]}"
                local gdb_workspace=$(basename ${out_dir}/Genotype_GVCFs/combinedDB/gendb_wksp_${current_intvl_name})
                # Check if our out dir exists
                mkdir -p "${out_dir}/Genotype_GVCFs/vcf_split_regions"
                # Check if our input intervals was a .bed or .intervals/.list format
                if [[ "${intervals_filepath}" == *.bed ]]; then
                    # If input is .bed file, save current interval to temporary file
                    # and use as interval
                    printf "${current_intvl}\n" > "${out_dir}/Genotype_GVCFs/temp_intervals/temp_${current_intvl_name}.bed"
                    set -x
                    #   Assume we are running GATK 4
                    #   As of Sep 23, 2019 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
                    #   CL is unable to get it working with a relative or absolute filepath to the database
                    #   Go into Genotype_GVCFs directory
                    cd "${out_dir}/Genotype_GVCFs/combinedDB"
                    gatk --java-options "-Xmx${memory}" \
                        "${analysisTypeOpt}" \
                        -R "${reference}" \
                        -L "${out_dir}/Genotype_GVCFs/temp_intervals/temp_${current_intvl_name}.bed" \
                        -V "gendb://${gdb_workspace}" \
                        --heterozygosity "${heterozygosity}" \
                        "${ploidyFlag}" "${ploidy}" \
                        ${outFlag} "${out_dir}/Genotype_GVCFs/vcf_split_regions/${current_intvl_name}.vcf"
                    set +x
                else
                    set -x
                    #   Assume we are running GATK 4
                    #   As of Sep 23, 2019 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
                    #   CL is unable to get it working with a relative or absolute filepath to the database
                    #   Go into Genotype_GVCFs directory
                    cd "${out_dir}/Genotype_GVCFs/combinedDB"
                    gatk --java-options "-Xmx${memory}" \
                        "${analysisTypeOpt}" \
                        -R "${reference}" \
                        -L "${current_intvl}" \
                        -V "gendb://${gdb_workspace}" \
                        --heterozygosity "${heterozygosity}" \
                        "${ploidyFlag}" "${ploidy}" \
                        ${outFlag} "${out_dir}/Genotype_GVCFs/vcf_split_regions/${current_intvl_name}.vcf"
                    set +x
                fi
            else
                # parallelize = false and analysisType="WGS", we have a single gendb workspace
                local out_name="${out_name_arr[${PBS_ARRAYID}]}"
                if [[ "$gatkVer" == 3 ]]; then
                    set -x
                    #   Run GATK using the parameters given
                    java -Xmx"${memory}" -jar "${gatk}" \
                        "${analysisTypeOpt}" \
                        -R "${reference}" \
                        -L "${intervals_filepath}" \
                        "${GATK_IN[@]}" \
                        --heterozygosity "${heterozygosity}" \
                        "${ploidyFlag}" "${ploidy}" \
                        "${outFlag} ${out_dir}/Genotype_GVCFs/${out_name}.vcf"
                    set +x
                else
                    set -x
                    #   Assume we are running GATK 4
                    #   As of Sep 23, 2019 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
                    #   CL is unable to get it working with a relative or absolute filepath to the database
                    #   Go into Genotype_GVCFs directory
                    cd "${out_dir}/Genotype_GVCFs/combinedDB"
                    gatk --java-options "-Xmx${memory}" \
                        "${analysisTypeOpt}" \
                        -R "${reference}" \
                        -L "${intervals_filepath}" \
                        -V "gendb://gendb_wksp" \
                        --heterozygosity "${heterozygosity}" \
                        "${ploidyFlag}" "${ploidy}" \
                        "${outFlag} ${out_dir}/Genotype_GVCFs/${out_name}.vcf"
                    set +x
                fi
            fi
        else
            # analysisType="WGS", use reference dict to create intervals list using chr
            local current_chr="${intvl_arr[${PBS_ARRAYID}]}"
            local current_chr_name="${out_name_arr[${PBS_ARRAYID}]}"
            set -x
            #   Assume we are running GATK 4
            #   As of Sep 23, 2019 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
            #   CL is unable to get it working with a relative or absolute filepath to the database
            #   Go into Genotype_GVCFs directory
            cd "${out_dir}/Genotype_GVCFs/combinedDB"
            gatk --java-options "-Xmx${memory}" \
                "${analysisTypeOpt}" \
                -R "${reference}" \
                -L "${current_chr}" \
                -V "gendb://gendb_wksp_${current_chr_name}" \
                --heterozygosity "${heterozygosity}" \
                "${ploidyFlag}" "${ploidy}" \
                "${outFlag} ${out_dir}/Genotype_GVCFs/${current_chr_name}.vcf"
            set +x
        fi
    else
        set -x
        parallel java -Xmx"${memory}" -jar "${gatk}" \
            "${analysisTypeOpt}" \
            -R "${reference}" \
            -L {1} \
            "${GATK_IN[@]}" \
            --heterozygosity "${heterozygosity}" \
            "${ploidyFlag}" "${ploidy}" \
            "${outFlag} ${out_dir}/Genotype_GVCFs/{2}.vcf" ::: "${intvl_arr[@]}" :::+ "${out_name_arr[@]}"
        set +x
    fi
}

#   Export the function
export -f Genotype_GVCFs
