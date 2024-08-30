#!/bin/bash

#   This script filters a final VCF file
#   based on user-specified parameters.

set -e
set -o pipefail

#   What are the dependencies for Variant_Filtering?
declare -a Variant_Filtering_Dependencies=(vcftools R vcfintersect python3 java bcftools)
source HelperScripts/graph_annotations.sh

function check_filter_stringency() {
    local vcf_file="$1"
    local num_sites="$2"
    if [[ ${num_sites} == 0 ]]; then
        echo "No sites left after filtering. Try using less stringent criteria. File with no sites is: ${vcf_file}. Exiting..." >&2
        exit 8 # If not sites left, error out with message
    fi
}

export -f check_filter_stringency

function count_sites() {
    local vcf="$1"
    local log_file="$2"
    # Get the number of sites left after filtering
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        num_sites=$(zgrep -v "#" ${vcf} | wc -l)
    else
        # We are working with uncompressed vcf
        num_sites=$(grep -v "#" ${vcf} | wc -l)
    fi
    # Append the number of sites remaining to file
    printf "${vcf}\t${num_sites}\n" >> ${log_file}
    check_filter_stringency ${vcf} ${num_sites}
}

export -f count_sites

function Variant_Filtering_GATK4() {
    local vcf="$1" # Where is our recalibrated VCF file?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference genome?
    local project="$4" # What is the name of this project?
    local seqhand="$5" # Where is sequence_handling located?
    local qual_cutoff="$6"
    local mindp="$7"
    local maxdp="$8"
    local gq_cutoff="$9"
    local qd_cutoff="${10}"
    local sor_cutoff="${11}"
    local fs_cutoff="${12}"
    local mq_cutoff="${13}"
    local prop_missing="${14}"
    local prop_het="${15}"
    local gen_num="${16}"
    local gen_len="${17}"
    local vf_diploid="${18}"
    local excess_het_cutoff="${19}"
    local temp_dir="${20}"
    # 0. Prepare directories and out file prefix
    # Check if out directories exist, if not make them
    mkdir -p ${out_dir} \
        ${out_dir}/Variant_Filtering \
        ${out_dir}/Variant_Filtering/Intermediates \
        ${temp_dir}
    
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        # We are working with a VCF produced from Variant_Recalibrator and
        #   Variant_Filtering handlers
        if [[ "${vcf}" == *".gz"* ]]; then
            out_prefix=$(basename ${vcf} .recalibrated.pass_sites.vcf.gz)
        else
            out_prefix=$(basename ${vcf} .recalibrated.pass_sites.vcf)
        fi
    else
        # Working with VCF not processed through sequence_handling
        out_prefix=${project}
    fi

    # At each step in the VCF filtering process, output the number of variants in the VCF to a file
    # Variant count in input VCF
    num_sites_log_file="${out_dir}/Variant_Filtering/num_sites_in_vcfs.txt"
    # Prepare header
    printf "Filename\tNum_Sites\n" > ${num_sites_log_file}
    # Count number of sites and append to file
    count_sites ${vcf} ${num_sites_log_file}

    # 1. Filter out lower quality sites
    echo "Filtering out low quality sites..."
    bcftools filter -e "QUAL < ${qual_cutoff}" ${vcf} > ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_qual_filtered.vcf
    # Store path in variable
    step1_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_qual_filtered.vcf"
    # Count number of sites and append to file
    count_sites ${step1_out_vcf} ${num_sites_log_file}

    # 2. Filter out low depth and extremely high depth sites
    echo "Filtering out low depth and extremely high depth sites..."
    bcftools filter -e "INFO/DP < ${mindp} | INFO/DP > ${maxdp}" ${step1_out_vcf} > ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_dp_filtered.vcf
    # Store path in variable
    step2_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_dp_filtered.vcf"
    # Count number of sites and append to file
    count_sites ${step2_out_vcf} ${num_sites_log_file}

    # 3. Filter out low GQ sites
    echo "Filter out low GQ sites..."
    bcftools filter -e "FORMAT/GQ < ${gq_cutoff}" ${step2_out_vcf} > ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_gq_filtered.vcf
    # Store path in variable
    step3_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_gq_filtered.vcf"
    # Count number of sites and append to file
    count_sites ${step3_out_vcf} ${num_sites_log_file}

    # 4. Filter sites with high proportions of missingness
    echo "Filter on missingness..."
    bcftools view -e "F_MISSING > ${prop_missing}" ${step3_out_vcf} > ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missing_filtered.vcf
    # Store path in variable
    step4_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missing_filtered.vcf"
    # Count number of sites and append to file
    count_sites ${step4_out_vcf} ${num_sites_log_file}

    # 5a. Filter out highly heterozygous sites
    echo "Filter highly heterozygous sites..."
    bcftools filter -i "COUNT(GT='het')/(N_SAMPLES-N_MISSING) < ${prop_het}" ${step4_out_vcf} > ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered_het.vcf
    # Store path in variable
    step5a_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered_het.vcf"
    # Count number of sites and append to file
    count_sites ${step5a_out_vcf} ${num_sites_log_file}

    # 5b. If conditions are met, filter on Excess Heterozygosity (see config for linked documentation)
    if [ ${vf_diploid} == "true" ] && [ ${excess_het_cutoff} != "NA" ]; then
        echo "Filtering on excess heterozygosity..."
        if [[ -z "${temp_dir}" ]]; then
            # No tmp directory specified
            gatk SelectVariants \
                -R ${ref} \
                -V ${step5a_out_vcf} \
                --selectExpressions "ExcessHet < ${excess_het_cutoff}" \
                -O ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered_excesshet.vcf.gz
        else
            # tmp directory specified
            gatk SelectVariants \
                -R ${ref} \
                -V ${step5a_out_vcf} \
                --selectExpressions "ExcessHet < ${excess_het_cutoff}" \
                -O ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered_excesshet.vcf.gz \
                --tmp-dir ${temp_dir}
        fi
        # Store vcf in variable
        step5b_out_vcf="${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_filtered_excesshet.vcf.gz"
        # Count number of sites and append to file
        count_sites ${step5b_out_vcf} ${num_sites_log_file}
    else
        # Conditions are not met, set vcf path to "NA"
        step5b_out_vcf="NA"
    fi

    if [ ${step5b_out_vcf} == "NA" ]; then
        # Set step5 out vcf to step5a_out_vcf
        step5_out_vcf=${step5a_out_vcf}
    else
        # We filtered on excess heterozygosity
        step5_out_vcf=${step5b_out_vcf}
    fi

    # 6. Remove sites that aren't polymorphic (minor allele count of 0) and unused alternate alleles.
    #   For large VCF files, this can take a long time, so we will check if the index has been
    #   created to indicate completion of the process. This allows for the handler to not
    #   re-run this time consuming step upon resubmission of partially completed job.
    if [ -f ${out_dir}/Variant_Filtering/${out_prefix}_polymorphic.vcf.gz.tbi ]; then
        echo "Sites that aren't polymorphic have already been removed, proceeding with existing file: ${out_dir}/Variant_Filtering/${out_prefix}_polymorphic.vcf.gz"
    else
        echo "Removing sites that aren't polymorphic and unused alternate alleles..."
        if [[ -z "${temp_dir}" ]]; then
            # No tmp directory specified
            gatk SelectVariants \
                -V "${step5_out_vcf}" \
                --exclude-non-variants true \
                --remove-unused-alternates true \
                -O "${out_dir}/Variant_Filtering/${out_prefix}_polymorphic.vcf.gz"
        else
            # tmp directory specified
            gatk SelectVariants \
                -V "${step5_out_vcf}" \
                --exclude-non-variants true \
                --remove-unused-alternates true \
                -O "${out_dir}/Variant_Filtering/${out_prefix}_polymorphic.vcf.gz" \
                --tmp-dir ${temp_dir}
        fi
        echo "Done removing sites that aren't polymorphic and unused alternate alleles."
    fi
    # Get the number of sites left after filtering out unbalanced heterozygotes
    step6_out_vcf="${out_dir}/Variant_Filtering/${out_prefix}_polymorphic.vcf.gz"
    # Count number of sites and append to file
    count_sites ${step6_out_vcf} ${num_sites_log_file}

    # 7. Filter to only biallelic sites
    echo "Filter to only biallelic sites..."
    bcftools view -m2 -M2 "${step6_out_vcf}" -O z -o "${out_dir}/Variant_Filtering/${project}_final.vcf.gz"
    # Store vcf in variable
    step7_out_vcf="${out_dir}/Variant_Filtering/${project}_final.vcf.gz"
    # Index VCF
    tabix -p vcf ${step7_out_vcf}
    # Count number of sites and append to file
    count_sites ${step7_out_vcf} ${num_sites_log_file}

    # 9. Generate plots following filtering
    final_vcf=${step7_out_vcf}
    # Generate graphs showing distributions of variant annotations
    echo "Graphing annotations..."
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        vcf_prefix="PassRecalSites" # Raw pass sites VCF file
    else
        vcf_prefix="Raw" # Raw VCF file
    fi
    filtered_vcf_prefix="Final_Filtered"
    graph_annotations \
        "${vcf}" \
        "${final_vcf}" \
        "${out_dir}/Variant_Filtering" \
        "${project}" \
        "${ref}" \
        "${seqhand}" \
        "${gen_num}" \
        "${gen_len}" \
        "${vcf_prefix}" \
        "${filtered_vcf_prefix}" \
        "pair"

    # Visualize missingness present in VCF after filtering
    echo "Visualizing missingness in data..."
    #   Use vcftools to identify amount of missingness in VCF
    if [[ ${final_vcf} == *".gz"* ]]; then
        # Use gzip flag
        vcftools --gzvcf ${final_vcf} \
            --missing-indv \
            --out ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness
        # Missing data per site
        vcftools --gzvcf ${final_vcf} \
            --missing-site \
            --out ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness
    else
        # Uncompressed VCF
        vcftools --vcf ${final_vcf} \
            --missing-indv \
            --out ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness
        # Missing data per site
        vcftools --vcf ${final_vcf} \
            --missing-site \
            --out ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness
    fi

    # Visualize missingness
    "${seqhand}/HelperScripts/graph_missingness.R" \
        ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness.imiss \
        ${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_missingness.lmiss \
        ${out_dir}/Variant_Filtering

    # Visualize Heterozygosity, Excess Heterozygosity, and Inbreeding Coefficients
    echo "Visualizing heterozygosity, excess heterozygosity, and inbreeding coefficients..."
    # Generate variants table for visualization
    gatk VariantsToTable \
        -V "${final_vcf}" \
        -F CHROM -F POS -F TYPE -F DP \
        -F ExcessHet -F InbreedingCoeff \
        -F HET -F HOM-REF -F HOM-VAR -F NCALLED \
        -O "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_Variants_HetInfo.table"
    # Generate the plots
    ${seqhand}/HelperScripts/graph_heterozygotes.R \
        "${out_dir}/Variant_Filtering/Intermediates/${out_prefix}_Variants_HetInfo.table" \
        "${out_dir}/Variant_Filtering"

    # 10. Remove intermediates to clear space
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
    gatk SelectVariants \
        -V "${out}/Intermediates/${project}_het_balanced.vcf" \
        --exclude-non-variants \
        -O "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf"
    local num_sites=$(grep -v "#" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" | wc -l) # Get the number of sites left after filtering out unbalanced heterozygotes
    if [[ "${num_sites}" == 0 ]]; then
        echo "No sites left after filtering out unbalanced heterozygotes! Try using less stringent criteria. Exiting..." >&2; exit 8 # If no sites left, error out with message
    else
        echo "Number of sites left after filtering out unbalanced heterozygotes: ${num_sites}"
    fi
    #   6. Create a percentile table for the het balanced SNPs
    percentiles "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${out}" "${project}" "het_balanced" "${seqhand}"
    #   7. Filter out sites that are low quality
    (set -x; python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Intermediates/${project}_het_balanced_poly.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 25; fi # If something went wrong with the python script, exit
    #   8. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    # vcftools --vcf "${out}/Intermediates/${project}_filtered.vcf" --non-ref-ac 1 --recode --recode-INFO-all --out "${out}/${project}_final"
    gatk SelectVariants \
        -V "${out}/Intermediates/${project}_filtered.vcf" \
        --exclude-non-variants \
        -O "${out}/${project}_final.recode.vcf"
    mv "${out}/${project}_final.recode.vcf" "${out}/${project}_final.vcf" # Rename the output file
    local num_sites_final=$(grep -v "#" "${out}/${project}_final.vcf" | wc -l) # Get the number of sites left after filtering
    if [[ "${num_sites_final}" == 0 ]]; then
        echo "No sites left after final filtering! Try using less stringent criteria. Exiting..." >&2; exit 9 # If no sites left, error out with message
    else
        echo "Number of sites left for *_final.vcf: ${num_sites_final}"
    fi
    #   9. Create a percentile table for the final SNPs
    percentiles "${out}/${project}_final.vcf" "${out}" "${project}" "final" "${seqhand}"
    #   10. Remove intermediates to clear space
    rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Filtering_GATK3
