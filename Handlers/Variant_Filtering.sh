#!/bin/bash

#   This script filters a final VCF file
#   based on user-specified parameters.

set -e
set -o pipefail

#   What are the dependencies for Variant_Filtering?
declare -a Variant_Filtering_Dependencies=(vcftools R vcfintersect python3)

function Variant_Filtering_GATK4() {
    local vcf="$1" # Where is our recalibrated VCF file?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference genome?
    local project="$4" # What is the name of this project?
    local seqhand="$5" # Where is sequence_handling located?
    local mindp="$6" # What is the minimum number of reads needed to support a genotype?
    local maxdp_percent="$7" # What is the maximum number of reads allowed to support a genotype?
    local maxdev="$8" # What is the maximum percent deviation from 50/50 reference/alternative reads allowed in heterozygotes?
    local dp_per_sample_cutoff="$9" # What is the DP per sample cutoff?
    local gq_cutoff_percent="${10}" # What is the genotyping quality cutoff?
    local max_het="${11}" # What is the maximum number of heterozygous samples?
    local max_bad="${12}" # What is the maximum number of bad samples?
    local qual_cutoff="${13}" # What is the quality cutoff?
    local temp_dir="${14}" # Where is our temporary directory?
    # Check if out directories exist, if not make them
    mkdir -p ${out_dir} \
        ${out_dir}/Variant_Filtering \
        ${out_dir}/Variant_Filtering/Intermediates \
        ${out_dir}/Variant_Filtering/Percentile_Tables \
        ${temp_dir}

    # 1. Remove SNPs that did not pass Variant_Recalibrator
    #   According to GATK docs: https://gatk.broadinstitute.org/hc/en-us/articles/360035531612
    #   Variants that are above the threshold pass the filter, so the FILTER field will contain PASS.
    #   Variants that are below the threshold will be filtered out; they will be written to the output
    #       file, but in the FILTER field they will have the name of the tranche they belonged to.
    #       So VQSRTrancheSNP99.90to100.00 means that the variant was in the range of VQSLODs
    #       corresponding to the remaining 0.1% of the truth set, which are considered false positives.
    #   Here, we will pull out only "PASS" variants.
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        # We are working with a VCF produced from Variant_Recalibrator
        out_prefix=$(basename ${vcf} .recalibrated.vcf.gz)
    else
        # We can't assume the file suffix, use the project to name the file
        out_prefix=${project}
    fi
    # Select PASS variants only
    # For large VCF files, this can take a long time, so we will check if the index has been
    #   created to indicate completion of the process. This allows for the handler to not
    #   re-run this time consuming step upon resubmission of partially completed job.
    if [ -f ${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz.tbi ]; then
        echo "Pass sites have been selected and the VCF has been indexed. Proceeding with existing file: ${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz"
    else
        echo "Selecting pass sites only from VCF file..."
        if [[ -z "${temp_dir}" ]]; then
            # No tmp directory specified
            gatk SelectVariants \
                -R ${ref} \
                -V ${vcf} \
                --exclude-filtered true \
                --create-output-variant-index true \
                -O ${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz
        else
            # tmp directory specified
            gatk SelectVariants \
                -R ${ref} \
                -V ${vcf} \
                --exclude-filtered true \
                --create-output-variant-index true \
                --tmp-dir ${temp_dir} \
                -O ${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz
        fi
        echo "Done selecting pass sites only."
    fi

    # 2. Create a percentile table for the unfiltered SNPs (pass recalibration sites in this case)
    source "${seqhand}/HelperScripts/percentiles.sh"
    percentiles "${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz" \
        "${out_dir}/Variant_Filtering" \
        "${project}" \
        "recal_pass_sites" \
        "${seqhand}"
    # Set maxdp based on percentiles table
    echo "Max percentile of reads allowed is currently set to: ${maxdp_percent}"
    local maxdp=$(grep -w ${maxdp_percent} ${out_dir}/Variant_Filtering/Percentile_Tables/${project}_recal_pass_sites_DP_per_sample.txt | cut -f 2)
    echo "Max number of reads allowed ${maxdp_percent} corresponds to the value ${maxdp}"
    # Set gq_cutoff based on percentiles table
    echo "GQ cutoff percentile is currently set to: ${gq_cutoff_percent}"
    local gq_cutoff=$(grep -w ${gq_cutoff_percent} ${out_dir}/Variant_Filtering/Percentile_Tables/${project}_recal_pass_sites_GQ.txt | cut -f 2)
    echo "GQ cutoff ${gq_cutoff_percent} corresponds to the value ${gq_cutoff}"

    # 3. Filter out unbalanced heterozygotes
    python3 "${seqhand}/HelperScripts/filter_genotypes.py" \
        "${out_dir}/Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz" \
        "${mindp}" \
        "${maxdp}" \
        "${maxdev}" \
        "${gq_cutoff}" \
        "${dp_per_sample_cutoff}" > "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced.vcf"
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_genotypes.py, exiting..." >&2; exit 24; fi # If something went wrong with the python script, exit

    # 4. Remove sites that aren't polymorphic (minor allele count of 0).
    #   This is because filter_genotypes.py has the potential to create monomorphic sites via filtering.
    #   It does still keep the monomorphic sites for the ALT allele.
    #   For large VCF files, this can take a long time, so we will check if the index has been
    #   created to indicate completion of the process. This allows for the handler to not
    #   re-run this time consuming step upon resubmission of partially completed job.
    if [ -f ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz.tbi ]; then
        echo "Sites that aren't polymorphic have already been removed, proceeding with existing file: ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz"
    else
        echo "Removing sites that aren't polymorphic..."
        if [[ -z "${temp_dir}" ]]; then
            # No tmp directory specified
            gatk SelectVariants \
                -V "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced.vcf" \
                --exclude-non-variants \
                -O "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz"
        else
            # tmp directory specified
            gatk SelectVariants \
                -V "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced.vcf" \
                --exclude-non-variants \
                -O "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz" \
                --tmp-dir ${temp_dir}
        fi
        echo "Done removing sites that aren't polymorphic."
    fi
    # Get the number of sites left after filtering out unbalanced heterozygotes
    hbp_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz"
    if [[ "${hbp_vcf}" == *".gz"* ]]; then
        local num_sites=$(zgrep -v "#" "${hbp_vcf}" | wc -l)
    else
        local num_sites=$(grep -v "#" "${hbp_vcf}" | wc -l)
    fi
    if [[ "${num_sites}" == 0 ]]; then echo "No sites left after filtering out unbalanced heterozygotes! Try using less stringent criteria. Exiting..." >&2; exit 8; fi # If no sites left, error out with message

    # 5. Create a percentile table for the het balanced SNPs
    percentiles ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz \
        "${out_dir}/Variant_Filtering" \
        "${project}" \
        "het_balanced" \
        "${seqhand}"

    # 6. Filter out sites that are low quality
    python3 "${seqhand}/HelperScripts/filter_sites.py" \
        "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_het_balanced_poly.vcf.gz" \
        "${qual_cutoff}" \
        "${max_het}" \
        "${max_bad}" \
        "${gq_cutoff}" \
        "${dp_per_sample_cutoff}" > "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered.vcf"
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 25; fi # If something went wrong with the python script, exit

    # 7. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    if [[ -z "${temp_dir}" ]]; then
        # No tmp directory specified
        gatk SelectVariants \
            -V "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered.vcf" \
            --exclude-non-variants \
            -O "${out_dir}/Variant_Filtering/${out_prefix}_final.vcf.gz"
    else
        # tmp directory specified
        gatk SelectVariants \
            -V "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered.vcf" \
            --exclude-non-variants \
            -O "${out_dir}/Variant_Filtering/${out_prefix}_final.vcf.gz" \
            --tmp-dir ${temp_dir}
    fi
    # Get the number of sites left after filtering
    final_vcf="${out_dir}/Variant_Filtering/${out_prefix}_final.vcf.gz"
    if [[ "${final_vcf}" == *".gz"* ]]; then
        local num_sites_final=$(zgrep -v "#" "${final_vcf}" | wc -l)
    else
        local num_sites_final=$(grep -v "#" "${final_vcf}" | wc -l)
    fi
    if [[ "${num_sites_final}" == 0 ]]; then echo "No sites left after filtering out low quality sites! Try using less stringent criteria. Exiting..." >&2; exit 9; fi

    # 8. Create a percentile table for the final SNPs
    percentiles "${out_dir}/Variant_Filtering/${out_prefix}_final.vcf.gz" \
        "${out_dir}/Variant_Filtering" \
        "${project}" \
        "final" \
        "${seqhand}"

    # 9. Remove intermediates to clear space
    #rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

export -f Variant_Filtering_GATK4

#   A function to call each filtering step
function Variant_Filtering_GATK3() {
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
    # vcftools 0.1.6 has a bug and convert '.|.' to '.'
    #    vcftools --vcf "${step1output}" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    gatk SelectVariants -V "${step1output}" --select-type-to-include SNP -O "${out}/Intermediates/${project}_no_indels.recode.vcf"
    #   3. Create a percentile table for the unfiltered SNPs
    source "${seqhand}/HelperScripts/percentiles.sh"
    percentiles "${out}/Intermediates/${project}_no_indels.recode.vcf" "${out}" "${project}" "raw" "${seqhand}"
    #   4. Filter out unbalanced heterozygotes 
    (set -x; python3 "${seqhand}/HelperScripts/filter_genotypes.py" "${out}/Intermediates/${project}_no_indels.recode.vcf" "${mindp}" "${maxdp}" "${maxdev}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_het_balanced.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_genotypes.py, exiting..." >&2; exit 24; fi # If something went wrong with the python script, exit
    #   5. Remove any sites that aren't polymorphic (minor allele count of 0). This is because filter_genotypes.py has the potential to create monomorphic sites via filtering. It does still keep the monomorphic sites for the ALT allele. 
    #    vcftools --vcf "${out}/Intermediates/${project}_het_balanced.vcf" --non-ref-ac 1 --recode --recode-INFO-all --out "${out}/Intermediates/${project}_het_balanced_poly"
    gatk SelectVariants -V "${out}/Intermediates/${project}_het_balanced.vcf" --exclude-non-variants -O "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf"
    local num_sites=$(grep -v "#" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" | wc -l) # Get the number of sites left after filtering out unbalanced heterozygotes
    if [[ "${num_sites}" == 0 ]]; then echo "No sites left after filtering out unbalanced heterozygotes! Try using less stringent criteria. Exiting..." >&2; exit 8; fi # If no sites left, error out with message
    #   6. Create a percentile table for the het balanced SNPs
    percentiles "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${out}" "${project}" "het_balanced" "${seqhand}"
    #   7. Filter out sites that are low quality
    (set -x; python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 25; fi # If something went wrong with the python script, exit
    #   8. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    # vcftools --vcf "${out}/Intermediates/${project}_filtered.vcf" --non-ref-ac 1 --recode --recode-INFO-all --out "${out}/${project}_final"
    gatk SelectVariants -V "${out}/Intermediates/${project}_filtered.vcf" --exclude-non-variants -O  "${out}/${project}_final.recode.vcf"
    mv "${out}/${project}_final.recode.vcf" "${out}/${project}_final.vcf" # Rename the output file
    local num_sites_final=$(grep -v "#" "${out}/${project}_final.vcf" | wc -l) # Get the number of sites left after filtering
    if [[ "${num_sites_final}" == 0 ]]; then echo "No sites left after filtering out low quality sites! Try using less stringent criteria. Exiting..." >&2; exit 9; fi # If no sites left, error out with message
    #   9. Create a percentile table for the final SNPs
    percentiles "${out}/${project}_final.vcf" "${out}" "${project}" "final" "${seqhand}"
    #   10. Remove intermediates to clear space
    rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Filtering_GATK3
