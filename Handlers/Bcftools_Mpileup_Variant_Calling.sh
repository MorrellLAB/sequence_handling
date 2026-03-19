#!/bin/bash

#   This script performs SNP-only variant calling using bcftools mpileup
#   as an alternative lightweight path after BAM files are sorted and indexed.
#   This handler generates only SNP variants (no indels) for short-read data.

set -e
set -o pipefail

#   What are the dependencies for Bcftools_Mpileup_Variant_Calling?
declare -a Bcftools_Mpileup_Variant_Calling_Dependencies=(bcftools)

function Bcftools_Mpileup_Variant_Calling() {
    local sample_list="$1" # What is our sample list (BAM files)?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference sequence?
    local threads="$4" # How many threads can we use?
    local tmp="${5:-}" # Temporary directory (optional)
    
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
    
    # Create output directories
    mkdir -p "${out_dir}/Bcftools_Mpileup_Variant_Calling"
    mkdir -p "${out_dir}/Bcftools_Mpileup_Variant_Calling/Intermediates"
    
    # Check if temp variable is provided
    if [[ -z "${tmp}" ]]; then
        echo "Not using temp directory, proceeding..."
    else
        echo "Making temp directory: ${tmp}"
        mkdir -p "${tmp}"
    fi
    
    # Check if reference genome is indexed for bcftools
    if ! [[ -f "${ref}.fai" ]]; then
        echo "Indexing reference genome for bcftools..."
        bcftools faidx "${ref}"
    fi
    
    # Verify all BAM files are indexed
    for bam in "${sample_array[@]}"; do
        if ! [[ -f "${bam}.bai" ]]; then
            echo "Indexing BAM file: ${bam}"
            samtools index "${bam}"
        fi
    done
    
    # Perform mpileup variant calling
    local mpileup_vcf="${out_dir}/Bcftools_Mpileup_Variant_Calling/Intermediates/mpileup_raw.vcf.gz"
    local filtered_vcf="${out_dir}/Bcftools_Mpileup_Variant_Calling/mpileup_snps_only.vcf.gz"
    local sample_log="${out_dir}/Bcftools_Mpileup_Variant_Calling/variant_calling.log"
    
    echo "[Bcftools_Mpileup] Starting SNP-only variant calling with bcftools mpileup" | tee -a "${sample_log}"
    echo "[Bcftools_Mpileup] Input BAM files:" | tee -a "${sample_log}"
    for bam in "${sample_array[@]}"; do
        echo "[Bcftools_Mpileup]   ${bam}" | tee -a "${sample_log}"
    done
    echo "[Bcftools_Mpileup] Reference genome: ${ref}" | tee -a "${sample_log}"
    
    # Run bcftools mpileup
    if bcftools mpileup \
        --fasta-ref "${ref}" \
        --output "${mpileup_vcf}" \
        --output-type z \
        --threads "${threads}" \
        --annotate FORMAT/AD,FORMAT/DP \
        "${sample_array[@]}" >> "${sample_log}" 2>&1; then
        
        echo "[Bcftools_Mpileup] Successfully completed mpileup" | tee -a "${sample_log}"
        
        # Filter to SNPs only (exclude indels)
        echo "[Bcftools_Mpileup] Filtering to SNPs only..." | tee -a "${sample_log}"
        
        if bcftools view \
            "${mpileup_vcf}" \
            --types snps \
            --output "${filtered_vcf}" \
            --output-type z >> "${sample_log}" 2>&1; then
            
            echo "[Bcftools_Mpileup] Successfully filtered to SNPs only" | tee -a "${sample_log}"
            echo "[Bcftools_Mpileup] Output VCF: ${filtered_vcf}" | tee -a "${sample_log}"
            
            # Index the final VCF
            if bcftools index "${filtered_vcf}" >> "${sample_log}" 2>&1; then
                echo "[Bcftools_Mpileup] Successfully indexed VCF" | tee -a "${sample_log}"
            else
                echo "[WARN] Failed to index VCF file" >&2
            fi
        else
            echo "[ERROR] Failed to filter VCF to SNPs only" >&2
            echo "[ERROR] Check log file: ${sample_log}" >&2
            return 1
        fi
    else
        echo "[ERROR] Bcftools mpileup failed" >&2
        echo "[ERROR] Check log file: ${sample_log}" >&2
        return 1
    fi
    
    echo "[Bcftools_Mpileup] SNP-only variant calling complete"
}

export -f Bcftools_Mpileup_Variant_Calling
