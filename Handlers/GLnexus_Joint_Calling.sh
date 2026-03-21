#!/bin/bash

#   This script performs joint genotyping on per-sample GVCFs
#   using GLnexus, consolidating calls from Clair3, DeepVariant,
#   or other callers that emit reference-confidence blocks.

set -euo pipefail

#   What are the dependencies for GLnexus_Joint_Calling?
declare -a GLnexus_Joint_Calling_Dependencies=(glnexus bcftools samtools)

function GLnexus_Joint_Calling() {
    local vcf_list="$1" # What is our list of single-sample GVCFs?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference sequence?
    local glnexus_config="${4:-DeepVariant}" # Which GLnexus config to use (default: DeepVariant)
    local threads="${5:-1}" # How many threads can we use?
    local tmp="${6:-}" # Temporary directory (optional)
    
    mapfile -t vcf_array < <(grep -E '\.(g\.vcf|gvcf|vcf)(\.gz)?$|\.bcf$' "${vcf_list}")
    
    # Validate input
    if [[ ${#vcf_array[@]} -eq 0 ]]; then
        echo "[ERROR] No GVCF files found in list: ${vcf_list}" >&2
        return 1
    fi
    
    # Create output directories
    mkdir -p "${out_dir}/GLnexus_Joint_Calling"
    mkdir -p "${out_dir}/GLnexus_Joint_Calling/Intermediates/bcf_indices"
    
    local glnexus_log="${out_dir}/GLnexus_Joint_Calling/glnexus.log"
    
    echo "[GLnexus] Starting joint genotyping" | tee -a "${glnexus_log}"
    echo "[GLnexus] Input GVCF count: ${#vcf_array[@]}" | tee -a "${glnexus_log}"
    echo "[GLnexus] GVCF list file: ${vcf_list}" | tee -a "${glnexus_log}"
    echo "[GLnexus] Reference genome: ${ref}" | tee -a "${glnexus_log}"
    echo "[GLnexus] GLnexus config: ${glnexus_config}" | tee -a "${glnexus_log}"
    echo "[GLnexus] Using ${threads} thread(s)" | tee -a "${glnexus_log}"
    
    # Check if temp variable is provided
    if [[ -z "${tmp}" ]]; then
        echo "[GLnexus] Not using temp directory, proceeding..." | tee -a "${glnexus_log}"
    else
        echo "[GLnexus] Making temp directory: ${tmp}" | tee -a "${glnexus_log}"
        mkdir -p "${tmp}"
    fi
    
    # Index reference genome if needed
    if ! [[ -f "${ref}.fai" ]]; then
        echo "[GLnexus] Indexing reference genome..." | tee -a "${glnexus_log}"
        samtools faidx "${ref}"
    fi
    
    # Prepare BCF index directory for GLnexus.
    # GLnexus consumes gVCF content so that non-variant reference blocks are preserved.
    echo "[GLnexus] Preparing GVCF-derived BCF files and indices..." | tee -a "${glnexus_log}"
    
    local bcf_index_dir="${out_dir}/GLnexus_Joint_Calling/Intermediates/bcf_indices"
    declare -a bcf_file_array=()
    
    for vcf in "${vcf_array[@]}"; do
        if [[ ! -f "${vcf}" ]]; then
            echo "[ERROR] GVCF file not found: ${vcf}" >&2
            return 1
        fi

        if [[ "${vcf}" != *.bcf ]]; then
            if ! bcftools view -h "${vcf}" | grep -q 'NON_REF'; then
                echo "[ERROR] Input does not appear to be a GVCF with reference-confidence blocks: ${vcf}" >&2
                echo "[ERROR] GLnexus joint genotyping requires per-sample GVCFs to avoid missing reference states." >&2
                return 1
            fi
        fi
        
        local vcf_name
        vcf_name=$(basename "${vcf}" .g.vcf.gz)
        vcf_name=$(basename "${vcf_name}" .gvcf.gz)
        vcf_name=$(basename "${vcf_name}" .vcf.gz)
        vcf_name=$(basename "${vcf_name}" .g.vcf)
        vcf_name=$(basename "${vcf_name}" .gvcf)
        vcf_name=$(basename "${vcf_name}" .vcf)
        vcf_name=$(basename "${vcf_name}" .bcf)
        local bcf_file="${bcf_index_dir}/${vcf_name}.bcf"
        
        echo "[GLnexus] Converting GVCF to BCF: ${vcf_name}" | tee -a "${glnexus_log}"
        
        if [[ "${vcf}" == *.bcf ]]; then
            cp "${vcf}" "${bcf_file}"
        elif bcftools view "${vcf}" \
            --output-type b \
            --output "${bcf_file}" \
            --threads "${threads}" \
            >> "${glnexus_log}" 2>&1; then
            true
        else
            echo "[ERROR] Failed to convert GVCF to BCF: ${vcf}" >&2
            return 1
        fi

        if bcftools index -f "${bcf_file}" >> "${glnexus_log}" 2>&1; then
            echo "[GLnexus] Successfully indexed BCF: ${bcf_file}" | tee -a "${glnexus_log}"
            bcf_file_array+=("${bcf_file}")
        else
            echo "[WARN] Failed to index BCF file: ${bcf_file}" >&2
            return 1
        fi
    done
    
    # Run GLnexus
    local glnexus_work_dir="${out_dir}/GLnexus_Joint_Calling/Intermediates/glnexus_workspace"
    echo "[GLnexus] Running GLnexus joint genotyping..." | tee -a "${glnexus_log}"
    echo "[GLnexus] Work directory: ${glnexus_work_dir}" | tee -a "${glnexus_log}"
    
    if glnexus \
        --dir "${glnexus_work_dir}" \
        --config "${glnexus_config}" \
        "${bcf_file_array[@]}" \
        >> "${glnexus_log}" 2>&1; then
        
        echo "[GLnexus] Successfully completed joint genotyping" | tee -a "${glnexus_log}"
        
        # Extract GLnexus output to VCF format
        local glnexus_joint_vcf="${out_dir}/GLnexus_Joint_Calling/glnexus_joint.vcf.gz"
        
        echo "[GLnexus] Extracting joint variants to VCF..." | tee -a "${glnexus_log}"
        if bcftools view \
            "${glnexus_work_dir}/variants" \
            --output-type z \
            --output "${glnexus_joint_vcf}" \
            --threads "${threads}" \
            >> "${glnexus_log}" 2>&1; then
            
            echo "[GLnexus] Successfully extracted joint VCF: ${glnexus_joint_vcf}" | tee -a "${glnexus_log}"
            
            # Index the output VCF
            if bcftools index "${glnexus_joint_vcf}" >> "${glnexus_log}" 2>&1; then
                echo "[GLnexus] Successfully indexed joint VCF" | tee -a "${glnexus_log}"
            else
                echo "[WARN] Failed to index joint VCF" >&2
            fi
            
            echo "[GLnexus] Joint genotyping pipeline complete" | tee -a "${glnexus_log}"
            echo "[GLnexus] Output VCF: ${glnexus_joint_vcf}" | tee -a "${glnexus_log}"
            
        else
            echo "[ERROR] Failed to extract GLnexus results to VCF" >&2
            return 1
        fi
        
    else
        echo "[ERROR] GLnexus joint calling failed" >&2
        echo "[ERROR] Check log file: ${glnexus_log}" >&2
        return 1
    fi
}

export -f GLnexus_Joint_Calling
