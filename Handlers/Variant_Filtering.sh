#!/bin/bash

#   This script filters a final VCF file
#   based on user-specified parameters.

set -o pipefail

#   What are the dependencies for Variant_Filtering?
declare -a Variant_Filtering_Dependencies=(vcftools R vcfintersect python3)

#   A function to call each filtering step
function Variant_Filtering() {
    local vcf="$1" # What is our input VCF file?
    local out="$2"/Variant_Filtering # Where are we storing our results?
    local bed="$3" # Where is the capture regions bed file?
    local barley="$4" # Is this barley?
    local project="$5" # What is the name of this project?
    local seqhand="$6" # Where is sequence_handling located?
    local mindp="$7" # What is the minimum number of reads needed to support a genotype?
    local maxdp="$8" # What is the maximum number of reads allowed to support a genotype?
    local maxdev="$9" # What is the maximum percent deviation from 50/50 reference/alternative reads allowed in heterozygotes?
    local dp_per_sample_cutoff="${10}" # What is the DP per sample cutoff?
    local gq_cutoff="${11}" # What is the genotyping quality cutoff?
    local max_het="${12}" # What is the maximum number of heterozygous samples?
    local max_bad="${13}" # What is the maximum number of bad samples?
    local qual_cutoff="${14}" # What is the quality cutoff?
    #   Make sure the out directories exist
    mkdir -p "${out}/Intermediates"
    mkdir -p "${out}/Percentile_Tables"
    #   0. Filter out SNPs that did not pass Variant_Recalibrator
    #   According to the GATK docs, SNPs tagged with "VQSRTrancheSNP99.90to100.00" in the FILTER field are considered to be false positives by the model
    #   https://software.broadinstitute.org/gatk/documentation/article.php?id=39
    (set -x; grep -v "VQSRTrancheSNP99.90to100.00" "${vcf}" > "${out}/Intermediates/${project}_recal_filtered.vcf")
    #   1. If exome capture, filter out SNPs outside the exome capture region. If not, then do nothing
    if ! [[ "${bed}" == "NA" ]]
    then
        (set -x; vcfintersect -b "${bed}" "${out}/Intermediates/${project}_recal_filtered.vcf" > "${out}/Intermediates/${project}_capture_regions.vcf") # Perform the filtering
        local step1output="${out}/Intermediates/${project}_capture_regions.vcf"
    else
        local step1output="${out}/Intermediates/${project}_recal_filtered.vcf"
    fi
    #   2. Filter out indels using vcftools
    vcftools --vcf "${step1output}" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    #   3. Create a percentile table for the unfiltered SNPs
    source "${seqhand}/HelperScripts/percentiles.sh"
    percentiles "${out}/Intermediates/${project}_no_indels.recode.vcf" "${out}" "${project}" "raw" "${seqhand}"
    #   4. Filter out unbalanced heterozygotes 
    (set -x; python3 "${seqhand}/HelperScripts/filter_genotypes.py" "${out}/Intermediates/${project}_no_indels.recode.vcf" "${mindp}" "${maxdp}" "${maxdev}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_het_balanced.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_genotypes.py, exiting..." >&2; exit 24; fi # If something went wrong with the python script, exit
    #   5. Remove any sites that aren't polymorphic (minor allele count of 0). This is because filter_genotypes.py has the potential to create monomorphic sites via filtering
    vcftools --vcf "${out}/Intermediates/${project}_het_balanced.vcf" --mac 1 --recode --recode-INFO-all --out "${out}/Intermediates/${project}_het_balanced_poly"
    local num_sites=$(grep -v "#" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" | wc -l) # Get the number of sites left after filtering out unbalanced heterozygotes
    if [[ num_sites == 0 ]]; then echo "No sites left after filtering out unbalanced heterozygotes! Try using less stringent criteria. Exiting..." >&2; exit 8; fi # If no sites left, error out with message
    #   6. Create a percentile table for the het balanced SNPs
    percentiles "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${out}" "${project}" "het_balanced" "${seqhand}"
    #   7. Filter out sites that are low quality
    (set -x; python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 25; fi # If something went wrong with the python script, exit
    #   8. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    vcftools --vcf "${out}/Intermediates/${project}_filtered.vcf" --mac 1 --recode --recode-INFO-all --out "${out}/${project}_final"
    mv "${out}/${project}_final.recode.vcf" "${out}/${project}_final.vcf" # Rename the output file
    local num_sites_final=$(grep -v "#" "${out}/${project}_final.vcf" | wc -l) # Get the number of sites left after filtering
    if [[ num_sites_final == 0 ]]; then echo "No sites left after filtering out low quality sites! Try using less stringent criteria. Exiting..." >&2; exit 9; fi # If no sites left, error out with message
    #   9. Create a percentile table for the final SNPs
    percentiles "${out}/${project}_final.vcf" "${out}" "${project}" "final" "${seqhand}"
    #   10. Remove intermediates to clear space
    #rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Filtering
