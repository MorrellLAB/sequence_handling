#!/bin/bash

#   This script generates graphs of vcf annotations and percentile tables for the VCF prior 
#       to running the Variant_Filtering handler.
#       If given a GATK4 recalibrated VCF, this script selects only PASS sites from 
#       the Variant_Recalibrator output VCF

set -e
set -o pipefail

#   What are the dependencies for Variant_Filtering?
declare -a Pre_Variant_Filtering_Dependencies=(vcftools R java)

function Pre_Variant_Filtering_GATK4() {
    local vcf="$1" # Where is our input VCF file?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference genome?
    local project="$4" # What is the name of this project?
    local seqhand="$5" # Where is sequence_handling located?
    local gen_num="$6"
    local gen_len="$7"
    local temp_dir="$8" # Where is our temporary directory?
    # Check if out directories exist, if not make them
    mkdir -p ${out_dir} \
        ${out_dir}/Pre_Variant_Filtering \
        ${out_dir}/Pre_Variant_Filtering/Intermediates \
        ${out_dir}/Pre_Variant_Filtering/Percentile_Tables \
        ${temp_dir}

    # 1. Preparation steps
    #   For GATK4 recalibrated VCF, we need to select PASS sites only first
    # Check if we are working with GATK 4 recalibrated vcf or non-GATK vcf
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        # We are working with a VCF produced from Variant_Recalibrator
        out_prefix=$(basename ${vcf} .recalibrated.vcf.gz)
        # Remove SNPs that did not pass Variant_Recalibrator (if this hasn't been done already)
        #   According to GATK docs: https://gatk.broadinstitute.org/hc/en-us/articles/360035531612
        #   Variants that are above the threshold pass the filter, so the FILTER field will contain PASS.
        #   Variants that are below the threshold will be filtered out; they will be written to the output
        #       file, but in the FILTER field they will have the name of the tranche they belonged to.
        #       So VQSRTrancheSNP99.90to100.00 means that the variant was in the range of VQSLODs
        #       corresponding to the remaining 0.1% of the truth set, which are considered false positives.
        #   Here, we will select PASS variants only
        # For large VCF files, this can take a long time, so we will check if the index has been
        #   created to indicate completion of the process. This allows for the handler to not
        #   re-run this time consuming step upon resubmission of partially completed job.
        if [ -f ${out_dir}/Variant_Recalibrator/${out_prefix}.recalibrated.pass_sites.vcf.gz.tbi ] || [ -f ${out_dir}/Pre_Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz.tbi ]; then
            echo "Pass sites have been selected and the VCF has been indexed. Proceeding with existing file."
        else
            echo "Selecting pass sites only from VCF file..."
            if [[ -z "${temp_dir}" ]]; then
                # No tmp directory specified
                gatk SelectVariants \
                    -R ${ref} \
                    -V ${vcf} \
                    --exclude-filtered true \
                    --create-output-variant-index true \
                    -O ${out_dir}/Pre_Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz
            else
                # tmp directory specified
                gatk SelectVariants \
                    -R ${ref} \
                    -V ${vcf} \
                    --exclude-filtered true \
                    --create-output-variant-index true \
                    --tmp-dir ${temp_dir} \
                    -O ${out_dir}/Pre_Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz
            fi
            echo "Done selecting pass sites only."
        fi
        vcf_file="${out_dir}/Pre_Variant_Filtering/${out_prefix}.recalibrated.pass_sites.vcf.gz"
        # Set variant table prefix for graphing annotations
        vcf_prefix="PassRecalSites_${out_prefix}"
    else
        # Working with non-GATK generated VCF file
        # We can't assume the file suffix, use the project to name the file
        out_prefix=${project}
        vcf_file="${vcf}"
        # Set variant table prefix for graphing annotations
        vcf_prefix="${project}"
    fi

    # 2. Generate graphs showing distributions of variant annotations
    echo "Generating graphs of variant annotation distributions..."
    source "${seqhand}/HelperScripts/graph_annotations.sh"
    # Check if we are working with GATK 4 recalibrated vcf or non-GATK vcf
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        # Graph recal sites and recal pass sites side by side
        raw_vcf_prefix="RawRecalSites"
        graph_annotations \
            "${vcf}" \
            "${vcf_file}" \
            "${out_dir}/Pre_Variant_Filtering" \
            "${project}" \
            "${ref}" \
            "${seqhand}" \
            "${gen_num}" \
            "${gen_len}" \
            "${raw_vcf_prefix}" \
            "${vcf_prefix}" \
            "pair"
    else
        # Graph just input vcf
        graph_annotations \
            "${vcf_file}" \
            "NA" \
            "${out_dir}/Pre_Variant_Filtering" \
            "${project}" \
            "${ref}" \
            "${seqhand}" \
            "${gen_num}" \
            "${gen_len}" \
            "${vcf_prefix}" \
            "NA" \
            "single"
    fi

    # 3. Visualize missingness present in VCF
    echo "Visualizing missingness in data..."
    #   Use vcftools to identify amount of missingness in VCF
    if [[ ${vcf_file} == *".gz"* ]]; then
        # Use gzip flag
        vcftools --vcf ${vcf_file} \
            --gzvcf \
            --missing-indv \
            --out ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness
        # Missing data per site
        vcftools --vcf ${vcf_file} \
            --gzvcf \
            --missing-site \
            --out ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness
    else
        # Uncompressed VCF
        vcftools --vcf ${vcf_file} \
            --missing-indv \
            --out ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness
        # Missing data per site
        vcftools --vcf ${vcf_file} \
            --missing-site \
            --out ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness
    fi

    # Visualize missingness
    "${seqhand}/HelperScripts/graph_missingness.R" \
        ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness.imiss \
        ${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_missingness.lmiss \
        ${out_dir}/Pre_Variant_Filtering

    # 4. Visualize Heterozygosity, Excess Heterozygosity, and Inbreeding Coefficients
    echo "Visualizing heterozygosity, excess heterozygosity, and inbreeding coefficients..."
    # Generate variants table for visualization
    gatk VariantsToTable \
        -V "${vcf_file}" \
        -F CHROM -F POS -F TYPE -F DP \
        -F ExcessHet -F InbreedingCoeff \
        -F HET -F HOM-REF -F HOM-VAR -F NCALLED \
        -O "${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_Variants_HetInfo.table"
    # Generate the plots
    ${seqhand}/HelperScripts/graph_heterozygotes.R \
        "${out_dir}/Pre_Variant_Filtering/Intermediates/${out_prefix}_Variants_HetInfo.table" \
        "${out_dir}/Pre_Variant_Filtering"

    # 5. Create a percentile table and plot for the unfiltered SNPs (pass recalibration sites in this case)
    # Generate graphs showing distributions of variant annotations
    echo "Generating percentiles tables..."
    source "${seqhand}/HelperScripts/percentiles.sh"
    # Check if we are working with GATK 4 recalibrated vcf or non-GATK vcf
    if [[ "${vcf}" == *"recalibrated"* ]]; then
        filtered_status="recal_pass_sites"
    else
        filtered_status="unfiltered_vcf"
    fi
    # Percentiles table
    percentiles "${vcf_file}" \
            "${out_dir}/Pre_Variant_Filtering" \
            "${project}" \
            "${filtered_status}" \
            "${seqhand}"
    echo "Done running Pre_Variant_Filtering handler."
}

export -f Pre_Variant_Filtering_GATK4
