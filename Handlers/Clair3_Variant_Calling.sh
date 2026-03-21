#!/bin/bash

#   This script performs SNP/indel calling on long-read BAM samples
#   using Clair3, which bypasses the GATK germline workflow.
#   Clair3 should emit per-sample GVCFs so downstream joint genotyping
#   via GLnexus retains homozygous-reference states.

set -euo pipefail

#   What are the dependencies for Clair3_Variant_Calling?
declare -a Clair3_Variant_Calling_Dependencies=(clair3 samtools bcftools)

function Clair3_Variant_Calling() {
    local sample_list="$1" # What is our sample list (BAM files)?
    local out_dir="$2" # Where are we storing our results?
    local ref="$3" # Where is the reference sequence?
    local seq_platform="$4" # What is the sequencing platform? (ONT, PACBIO, PACBIO_HIFI)
    local basecaller="$5" # What basecaller was used? (e.g., guppy, dorado) - optional for ONT
    local threads="$6" # How many threads can we use?
    local memory="$7" # How much memory can we use?
    local tmp="${8:-}" # Temporary directory (optional)
    local generate_all_sites="${9:-false}" # Generate AllSites VCF for pixy? (default: false)
    
    # Clair3 supports GVCF output and GLnexus joint genotyping depends on it
    # to retain reference-confidence blocks instead of introducing missing data.
    mapfile -t sample_array < <(grep -E '\.bam$' "${sample_list}")

    if [[ ${#sample_array[@]} -eq 0 ]]; then
        echo "[ERROR] No BAM files found in list: ${sample_list}" >&2
        return 1
    fi
    
    # Create output directories
    mkdir -p "${out_dir}/Clair3_Variant_Calling"
    mkdir -p "${out_dir}/Clair3_Variant_Calling/Intermediates"
    
    # Check if temp variable is provided
    if [[ -z "${tmp}" ]]; then
        echo "Not using temp directory, proceeding..."
    else
        echo "Making temp directory: ${tmp}"
        mkdir -p "${tmp}"
    fi
    
    # Index reference genome if needed
    if ! [[ -f "${ref}.fai" ]]; then
        echo "Indexing reference genome..."
        samtools faidx "${ref}"
    fi
    
    # Process each BAM file
    for sample in "${sample_array[@]}"; do
        local sample_name
        sample_name=$(basename "${sample}" .bam)
        local sample_gvcf="${out_dir}/Clair3_Variant_Calling/${sample_name}.gvcf.gz"
        local sample_log="${out_dir}/Clair3_Variant_Calling/${sample_name}.log"
        local clair3_output_dir="${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}"
        local clair3_output_gvcf="${clair3_output_dir}/merge_output.gvcf.gz"
        local clair3_output_vcf="${clair3_output_dir}/merge_output.vcf.gz"
        
        echo "[Clair3] Starting variant calling for sample: ${sample_name}" | tee -a "${sample_log}"
        echo "[Clair3] Input BAM: ${sample}" | tee -a "${sample_log}"
        echo "[Clair3] Sequencing platform: ${seq_platform}" | tee -a "${sample_log}"
        echo "[Clair3] Memory budget: ${memory}" | tee -a "${sample_log}"
        
        # Determine Clair3 model based on sequencing platform
        local model=""
        case "${seq_platform}" in
            ONT)
                if [[ -n "${basecaller}" ]]; then
                    # ONT with basecaller information
                    case "${basecaller}" in
                        dorado) model="ont_small" ;;
                        guppy) model="ont_small" ;;
                        *) model="ont_small" ;;
                    esac
                else
                    model="ont_small"
                fi
                ;;
            ONT_Q20|ONT-Q20)
                model="ont_small"
                ;;
            PACBIO|PACBIO_HIFI|PACBIO-HiFi)
                model="hifi_small"
                ;;
            *)
                echo "[ERROR] Unknown sequencing platform: ${seq_platform}. Expected: ONT, ONT_Q20, PACBIO, PACBIO_HIFI" >&2
                echo "[ERROR] Failed for sample: ${sample_name}" >&2
                return 1
                ;;
        esac
        
        echo "[Clair3] Using model: ${model}" | tee -a "${sample_log}"
        
        # Run Clair3 with GVCF output enabled for downstream GLnexus joint genotyping.
        if clair3 \
            --bam_fn="${sample}" \
            --ref_fn="${ref}" \
            --output="${clair3_output_dir}" \
            --model_path="${model}" \
            --num_workers="${threads}" \
            --chunk_size=5000 \
            --chunk_num=5 \
            --sample_name="${sample_name}" \
            --gvcf \
            --enable_phasing \
            >> "${sample_log}" 2>&1; then
            
            echo "[Clair3] Successfully completed variant calling for sample: ${sample_name}" | tee -a "${sample_log}"
            
            # Clair3 may expose the final merged output as either merge_output.gvcf.gz
            # or merge_output.vcf.gz depending on version, but the content must be GVCF here.
            if [[ -f "${clair3_output_gvcf}" ]]; then
                mv "${clair3_output_gvcf}" "${sample_gvcf}"
                if [[ -f "${clair3_output_gvcf}.tbi" ]]; then
                    mv "${clair3_output_gvcf}.tbi" "${sample_gvcf}.tbi"
                fi
            elif [[ -f "${clair3_output_vcf}" ]]; then
                mv "${clair3_output_vcf}" "${sample_gvcf}"
                if [[ -f "${clair3_output_vcf}.tbi" ]]; then
                    mv "${clair3_output_vcf}.tbi" "${sample_gvcf}.tbi"
                fi
            fi

            if [[ -f "${sample_gvcf}" ]]; then
                echo "[Clair3] Output GVCF: ${sample_gvcf}" | tee -a "${sample_log}"
                echo "[Clair3] This GVCF is ready for joint genotyping via GLnexus" | tee -a "${sample_log}"
                echo "[Clair3] See: https://github.com/dnanexus-rnd/GLnexus" | tee -a "${sample_log}"

                if bcftools index -f "${sample_gvcf}" >> "${sample_log}" 2>&1; then
                    echo "[Clair3] Successfully indexed GVCF" | tee -a "${sample_log}"
                else
                    echo "[WARN] Failed to index GVCF file" >&2
                fi
                
                # If AllSites VCF for pixy is requested
                if [[ "${generate_all_sites}" == "true" ]]; then
                    local sample_allsites="${out_dir}/Clair3_Variant_Calling/${sample_name}_allsites.vcf.gz"
                    echo "[Clair3] Generating AllSites VCF for pixy analysis: ${sample_allsites}" | tee -a "${sample_log}"
                    
                    # Note: Full AllSites VCF generation requires bcftools mpileup on the entire genome
                    # This is a placeholder comment for now; see pixy documentation at:
                    # https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
                    echo "[Clair3] AllSites VCF generation requires additional steps - see pixy documentation" | tee -a "${sample_log}"
                fi
            else
                echo "[ERROR] Clair3 output GVCF not found for sample: ${sample_name}" >&2
                return 1
            fi
        else
            echo "[ERROR] Clair3 failed for sample: ${sample_name}" >&2
            echo "[ERROR] Check log file: ${sample_log}" >&2
            return 1
        fi
    done
    
    echo "[Clair3] Variant calling complete for all samples"
}

export -f Clair3_Variant_Calling
