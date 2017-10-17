#!/bin/bash

#   This script filters a final VCF file
#   based on user-specified parameters.

set -o pipefail

#   What are the dependencies for Variant_Filtering?
declare -a Variant_Filtering_Dependencies=(parallel vcftools R python2 vcfintersect python3)

#   A function to write a percentile table
function percentiles() {
    local vcf="$1"
    local out="$2"
    local project="$3"
    local filtered="$4"
    local seqhand="$5"
    #   Extract GQ information
    vcftools --vcf "${vcf}" --extract-FORMAT-info GQ --out "${out}/Intermediates/${project}_${filtered}"
    cut -f 3- "${out}/Intermediates/${project}_${filtered}.GQ.FORMAT" | sed 1,1d > "${out}/Intermediates/${project}_${filtered}.GQ.matrix"
    #   Extract DP per sample information
    vcftools --vcf "${vcf}" --geno-depth --out "${out}/Intermediates/${project}_${filtered}"
    sed 1d "${out}/Intermediates/${project}_${filtered}.gdepth" | cut -f 3- > "${out}/Intermediates/${project}_${filtered}.gdepth.matrix"
    sed -i -e 's/-1/NA/g' "${out}/Intermediates/${project}_${filtered}.gdepth.matrix" # Change -1 to NA
    #   Call the R script to write the percentiles table
    Rscript "${seqhand}/HelperScripts/percentiles.R" "${out}/Intermediates/${project}_${filtered}.GQ.matrix" "${out}/Intermediates/${project}_${filtered}.gdepth.matrix" "${out}/Percentile_Tables" "${project}" "${filtered}"
}

#   Export the function
export -f percentiles

#   A function to call each filtering step
function Variant_Filtering() {
    local vcf="$1" # What is our input VCF file?
    local out="$2"/Variant_Filtering # Where are we storing our results?
    local bed="$3" # Where is the capture regions bed file?
    local barley="$4" # Is this barley?
    local project="$5" # What is the name of this project?
    local seqhand="$6" # Where is sequence_handling located?
    
    local qual_cutoff="$7" # What is the quality cutoff?
    local gq_cutoff="$8" # What is the genotyping quality cutoff?
    local max_low_gq="$9" # What is the maximum number of low GQ samples?
    local dp_per_sample_cutoff="${10}" # What is the DP per sample cutoff?
    local max_low_dp="${11}" # What is the maximum number of samples with low DP?
    local max_het="${12}" # What is the maximum number of heterozygous samples?
    local max_missing="${13}" # What is the maximum number of missing samples?
    #   Make sure the out directories exist
    mkdir -p "${out}/Intermediates/Parts"
    mkdir -p "${out}/Percentile_Tables"
    #   1. Filter out indels using vcftools
    vcftools --vcf "${vcf}" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    #   2. If exome capture, filter out SNPs outside the exome capture region. If not, then do nothing
    if ! [[ "${bed}" == "NA" ]]
    then
        vcfintersect -b "${bed}" "${out}/Intermediates/${project}_no_indels.recode.vcf" > "${out}/Intermediates/${project}_capture_regions.vcf" # Perform the filtering
        local step4output="${out}/Intermediates/${project}_capture_regions.vcf"
    else
        local step4output="${out}/Intermediates/${project}_no_indels.recode.vcf"
    fi
    #   3. Create a percentile table for the unfiltered SNPs
    percentiles "${step4output}" "${out}" "${project}" "raw" "${seqhand}"
    #   4. Filter out unbalanced heterozygotes 
    python3 "${seqhand}/HelperScripts/filter_genotypes.py" "${step4output}" ... > "${out}/Intermediates/${project}_het_balanced.vcf"  
    #   5. Create a percentile table for the het balanced SNPs
    percentiles "${out}/Intermediates/${project}_het_balanced.vcf" "${out}" "${project}" "het_balanced" "${seqhand}"
    #   6. Filter out sites that are low quality
    python3 "${seqhand}/HelperScripts/filter_sites.py" "${step4output}" "${qual_cutoff}" "${max_het}" "${max_missing}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf"  
    #   7. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    vcftools --vcf "${out}/Intermediates/${project}_filtered.vcf" --mac 1 --recode --recode-INFO-all --out "${out}/${project}_final"
    mv "${out}/${project}_final.recode.vcf" "${out}/${project}_final.vcf" # Rename the output file
    #   8. Create a percentile table for the final SNPs
    percentiles "${out}/Intermediates/${project}_final.vcf" "${out}" "${project}" "final" "${seqhand}"
    #   9. Remove intermediates to clear space
    #rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Filtering