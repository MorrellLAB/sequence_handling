#!/bin/bash

#   This script creates a high-confidence
#   subset of variants to use in variant recalibration.

set -o pipefail

#   What are the dependencies for Create_HC_Subset?
declare -a Create_HC_Subset_Dependencies=(parallel vcftools R vcfintersect python3)

#   A function to call each filtering step
function Create_HC_Subset() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Create_HC_Subset # Where are we storing our results?
    local bed="$3" # Where is the capture regions bed file?
    local barley="$4" # Is this barley?
    local project="$5" # What is the name of this project?
    local seqhand="$6" # Where is sequence_handling located?
    local qual_cutoff="$7" # What is the quality cutoff?
    local gq_cutoff="$8" # What is the genotyping quality cutoff?
    local dp_per_sample_cutoff="$9" # What is the DP per sample cutoff?
    local max_het="${10}" # What is the maximum number of heterozygous samples?
    local max_bad="${11}" # What is the maximum number of bad samples?
    #   Make sure the out directories exist
    mkdir -p "${out}/Intermediates/Parts"
    mkdir -p "${out}/Percentile_Tables"
    #   1. Gzip all the chromosome part VCF files
    source "${seqhand}/HelperScripts/gzip_parts.sh"
    parallel -v gzip_parts {} "${out}/Intermediates/Parts" :::: "${sample_list}" # Do the gzipping in parallel, preserve original files
    "${seqhand}/HelperScripts/sample_list_generator.sh" .vcf.gz "${out}/Intermediates/Parts" gzipped_parts.list # Make a list of the gzipped files for the next step
    #   2. Use vcftools to concatenate all the gzipped VCF files
    vcf-concat -f "${out}/Intermediates/Parts/gzipped_parts.list" > "${out}/Intermediates/${project}_concat.vcf"
    #   3. If exome capture, filter out SNPs outside the exome capture region. If not, then do nothing
    if ! [[ "${bed}" == "NA" ]]
    then
        (set -x; vcfintersect -b "${bed}" "${out}/Intermediates/${project}_concat.vcf" > "${out}/Intermediates/${project}_capture_regions.vcf") # Perform the filtering
        local step3output="${out}/Intermediates/${project}_capture_regions.vcf"
    else
        local step3output="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   4. Filter out indels using vcftools
    vcftools --vcf "${step3output}" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    #   5. Create a percentile table for the unfiltered SNPs
    source "${seqhand}/HelperScripts/percentiles.sh"
    percentiles "${out}/Intermediates/${project}_no_indels.recode.vcf" "${out}" "${project}" "unfiltered" "${seqhand}"
    if [[ "$?" -ne 0 ]]; then echo "Error creating raw percentile tables, exiting..." >&2; exit 32; fi # If something went wrong with the R script, exit
    #   6. Filter out sites that are low quality
    (set -x; python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Intermediates/${project}_no_indels.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 22; fi # If something went wrong with the python script, exit
    local num_sites=$(grep -v "#" "${out}/Intermediates/${project}_filtered.vcf" | wc -l) # Get the number of sites left after filtering
    if [[ num_sites == 0 ]]; then echo "No sites left after filtering! Try using less stringent criteria. Exiting..." >&2; exit 23; fi # If no sites left, error out with message
    #   7. Create a percentile table for the filtered SNPs
    percentiles "${out}/Intermediates/${project}_filtered.vcf" "${out}" "${project}" "filtered" "${seqhand}"
    if [[ "$?" -ne 0 ]]; then echo "Error creating filtered percentile tables, exiting..." >&2; exit 33; fi # If something went wrong with the R script, exit
    #   8. If barley, convert the parts positions into pseudomolecular positions. If not, then do nothing
    if [[ "${barley}" == true ]]
    then
        (set -x; python3 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" "${out}/Intermediates/${project}_filtered.vcf" > "${out}/Intermediates/${project}_pseudo.vcf")
        local step8output="${out}/Intermediates/${project}_pseudo.vcf"
    else
        local step8output="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   9. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    vcftools --vcf "${step8output}" --mac 1 --recode --recode-INFO-all --out "${out}/${project}_high_confidence_subset"
    mv "${out}/${project}_high_confidence_subset.recode.vcf" "${out}/${project}_high_confidence_subset.vcf" # Rename the output file
    #   10. Remove intermediates to clear space
    #rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Create_HC_Subset
