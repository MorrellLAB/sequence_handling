#!/bin/bash

#   This script creates VCF files for each 
#   chromosome part using GVCFs as input.

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
    local seqs_list="$8" # What is our array of chromosomes or chromosome parts?
    local current="${seqs_list[${PBS_ARRAYID}]}" # What is the current chromosome part we're working on?
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    # Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
		GATK_IN+=("-V $s")
	done
    # Make sure the out directory exists
    mkdir -p "${out}"
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T GenotypeGVCFs \
	    -R "${reference}" \
        -L "${current}" \
	    "${GATK_IN[@]}" \
	    --heterozygosity "${heterozygosity}" \
	    --sample_ploidy "${ploidy}" \
	    -o "${out}/${current}.vcf"
}

#   Export the function
export -f Genotype_GVCFs