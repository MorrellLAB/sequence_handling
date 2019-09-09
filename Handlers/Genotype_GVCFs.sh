#!/bin/bash

#   This script creates VCF files for each
#   chromosome part using GVCFs as input.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_GenotypeGVCFs.job

set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genotype_GVCFs_Dependencies=(java)

#   A function to genotype the GVCFs without PBS
function Genotype_GVCFs() {
    local sample_list="$1" # What is our sample list? Or filename of combined GVCF with gatk4
    local out="$2"/Genotype_GVCFs # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local ploidy="$6" # What is the sample ploidy?
    local memory="$7" # How much memory can java use?
    local dict="$8" # Where is the reference dictionary?
    local maxarray="$9" # What is the maximum array index?
    local gatkVer="$10"  # GATK version  3 or 4
    local scaffolds="${11}" # Where is the scaffolds intervals file?
    local seqs_list=()
    if [[ "${maxarray}" -ne 0 ]]; then
    	seqs_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
    fi
    # a chromosome part name is used as the output name, or "custom_intervals" for the scaffolds in CUSTOM_INTERVAL
    out_name_arr=("${seqs_list[@]}")
    if [[ ! -z "${scaffolds}" ]]; then # custom interval is appended
        seqs_list+=("${scaffolds}")
        out_name_arr+=("custom_intervals")
    fi
    if [ ${#seqs_list[@]} -ne $[${maxarray} + 1] ]; then echo "Check NUM_CHR and CUSTOM_INTERVAL, it doesn't match with entries in the reference" >&2; exit 1; fi
    if [[ "$gatkVer" == 3 ]]; then
        analysisTypeOpt="-T GenotypeGVCFs"
        outFlag="-o"
        ploidyFlag="--sample_ploidy"
    else
        analysisTypeOpt="GenotypeGVCFs"
        outFlag="-O"
        # note that with GATK 4 the documentation says, no need to specify this
        ploidyFlag="--sample-ploidy"
    fi
    if [[ "$gatkVer" == 3 ]]; then
        declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
        #   Put the samples into a format that GATK can read
        GATK_IN=()
        for s in "${sample_array[@]}"
        do
            GATK_IN+=(-V $s)
        done
    else
        GATK_IN=("-V ${sample_list}")
    fi
    #   Make sure the out directory exists
    mkdir -p "${out}"
    if [[ "$USE_PBS" == "true" ]]; then
        #   What region of the genome are we working on currently?
        local current="${seqs_list[${PBS_ARRAYID}]}"
        local out_name="${out_name_arr[${PBS_ARRAYID}]}"
        #   Run GATK using the parameters given
        (set -x; java -Xmx"${memory}" -jar "${gatk}" \
            "${analysisTypeOpt}"  \
            -R "${reference}" \
            -L "${current}" \
            "${GATK_IN[@]}" \
            --heterozygosity "${heterozygosity}" \
            "${ploidyFlag}" "${ploidy}" \
            "${outFlag} ${out}/${out_name}.vcf")
    else
        set -x; parallel java -Xmx"${memory}" -jar "${gatk}" \
        "${analysisTypeOpt}"  \
            -R "${reference}" \
            -L {1} \
            "${GATK_IN[@]}" \
            --heterozygosity "${heterozygosity}" \
            "${ploidyFlag}" "${ploidy}" \
            "${outFlag} ${out}/{2}.vcf" ::: "${seqs_list[@]}" :::+ "${out_name_arr[@]}"
    fi
}

#   Export the function
export -f Genotype_GVCFs
