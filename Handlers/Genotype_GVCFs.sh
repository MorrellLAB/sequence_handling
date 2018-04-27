#!/bin/bash

#   This script creates VCF files for each 
#   chromosome part using GVCFs as input.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_GenotypeGVCFs.job

set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genotype_GVCFs_Dependencies=(java)

#   A function to genotype the GVCFs
function Genotype_GVCFs() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Genotype_GVCFs # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local ploidy="$6" # What is the sample ploidy?
    local memory="$7" # How much memory can java use?
    local dict="$8" # Where is the reference dictionary?
    local maxarray="$9" # What is the maximum array index?
    local scaffolds="${10}" # Where is the scaffolds intervals file?
    local seqs_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
    #   What region of the genome are we working on currently?
    if [[ "${PBS_ARRAYID}" = "${maxarray}" && ! -z "${scaffolds}" ]]
    then
        #   If this is the last array index AND we have a scaffolds file, use the scaffolds
        local current="${scaffolds}"
        local out_name="custom_intervals"
    else 
        #   If this isn't the last array index OR we don't have a scaffolds file, use a chromosome part name
        local current="${seqs_list[${PBS_ARRAYID}]}"
        local out_name="${current}"
    fi
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
		GATK_IN+=(-V $s)
	done
    #   Make sure the out directory exists
    mkdir -p "${out}"
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T GenotypeGVCFs \
	    -R "${reference}" \
        -L "${current}" \
	    "${GATK_IN[@]}" \
	    --heterozygosity "${heterozygosity}" \
	    --sample_ploidy "${ploidy}" \
	    -o "${out}/${out_name}.vcf")
}

#   Export the function
export -f Genotype_GVCFs