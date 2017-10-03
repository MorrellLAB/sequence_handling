#!/bin/bash

#   This script creates a high-confidence
#   subset of variants to use in variant recalibration.

set -o pipefail

#   What are the dependencies for Create_HC_Subset?
declare -a Create_HC_Subset_Dependencies=(parallel vcftools R python2 vcfintersect python3)

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

#   A function to do gzipping while preserving the original file
function gzip_parts() {
    local sample="$1" # The vcf file to be gzipped
    local out="$2" # Where to put the gzipped file
    local name=$(basename ${sample} .vcf) # The name of the sample
    gzip -c "${sample}" > "${out}/${name}.vcf.gz" # Perform the gzipping without altering the original file
}

#   Export the function
export -f gzip_parts

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
    local max_low_gq="$9" # What is the maximum number of low GQ samples?
    local dp_per_sample_cutoff="${10}" # What is the DP per sample cutoff?
    local max_low_dp="${11}" # What is the maximum number of samples with low DP?
    local max_het="${12}" # What is the maximum number of heterozygous samples?
    local max_missing="${13}" # What is the maximum number of missing samples?
    #   Make sure the out directories exist
    mkdir -p "${out}/Intermediates/Parts"
    mkdir -p "${out}/Percentile_Tables"
    #   1. Gzip all the chromosome part VCF files
    parallel -v gzip_parts {} "${out}/Intermediates/Parts" :::: "${sample_list}" # Do the gzipping in parallel, preserve original files
    "${seqhand}/HelperScripts/sample_list_generator.sh" .vcf.gz "${out}/Intermediates/Parts" gzipped_parts.list # Make a list of the gzipped files for the next step
    #   2. Use vcftools to concatenate all the gzipped VCF files
    vcf-concat -f "${out}/Intermediates/Parts/gzipped_parts.list" > "${out}/Intermediates/${project}_concat.vcf"
    #   3. Filter out indels using vcftools
    vcftools --vcf "${out}/Intermediates/${project}_concat.vcf" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    #   4. If exome capture, filter out SNPs outside the exome capture region. If not, then do nothing
    if ! [[ "${bed}" == "NA" ]]
    then
        vcfintersect -b "${bed}" "${out}/Intermediates/${project}_no_indels.recode.vcf" > "${out}/Intermediates/${project}_capture_regions.vcf" # Perform the filtering
        local step4output="${out}/Intermediates/${project}_capture_regions.vcf"
    else
        local step4output="${out}/Intermediates/${project}_no_indels.recode.vcf"
    fi
    #   5. Create a percentile table for the unfiltered SNPs
    percentiles "${step4output}" "${out}" "${project}" "unfiltered" "${seqhand}"
    #   6. Filter out sites that are low quality
    python3 "${seqhand}/HelperScripts/filter_sites_stringent.py" "${step4output}" "${qual_cutoff}" "${max_het}" "${max_missing}" "${gq_cutoff}" "${max_low_gq}" "${dp_per_sample_cutoff}" "${max_low_dp}" > "${out}/Intermediates/${project}_filtered.vcf"
    #   7. Create a percentile table for the filtered SNPs
    percentiles "${out}/Intermediates/${project}_filtered.vcf" "${out}" "${project}" "filtered" "${seqhand}"
    #   8. If barley, convert the parts positions into pseudomolecular positions. If not, then do nothing
    if [[ "${barley}" == true ]]
    then
        python2 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" "${out}/Intermediates/${project}_filtered.vcf" > "${out}/Intermediates/${project}_pseudo.vcf"
        local step8output="${out}/Intermediates/${project}_pseudo.vcf"
    else
        local step8output="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   9. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    vcftools --vcf "${step8output}" --mac 1 --recode --recode-INFO-all --out "${out}/${project}_high_confidence_subset"
    mv "${out}/${project}_high_confidence_subset.recode.vcf" "${out}/${project}_high_confidence_subset.vcf" # Rename the output file
    #   10. Remove intermediates to clear space
    rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Create_HC_Subset
