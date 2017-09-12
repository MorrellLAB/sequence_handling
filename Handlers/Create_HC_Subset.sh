#!/bin/bash

#   This script creates a high-confidence
#   subset of variants to use in variant recalibration.

set -o pipefail

#   What are the dependencies for Create_HC_Subset?
declare -a Create_HC_Subset_Dependencies=(java vcftools vcflib python)

#   A function to remove indels
function filterIndels() {
    local sample="$1" # What is the vcf file we're filtering?
    local out="$2" # Where is the output directory?
    local name=$(basename ${sample} .vcf) # What is the name of the sample without the .vcf?
    vcftools --vcf "${sample}" --remove-indels --recode --recode-INFO-all --out "${out}/${sample}_no_indels.vcf" # Perform the filtering
    rm "${out}/${sample}_no_indels.log" # Remove the log file
}

#   Export the function
export -f filterIndels




#   A function to call each filtering step
function Create_HC_Subset() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Genotype_GVCFs # Where are we storing our results?

    declare -a sample_array=($(grep -E ".vcf" "${sample_list}")) # Put the sample list into array format
    echo "1. Filter indels with vcftools" >&2
    parallel --verbose filterIndels
    echo "2. "


    

    # Make sure the out directory exists
    mkdir -p "${out}"

}

#   Export the function
export -f Create_HC_Subset