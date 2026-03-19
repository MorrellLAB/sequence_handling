#!/bin/bash

#   This script performs SNP/indel calling on long-read BAM samples
#   using Clair3, which bypasses the GATK germline workflow.
#   Clair3 outputs standard VCFs suitable for joint calling via GLnexus
#   or other non-GATK joint callers.

set -e
set -o pipefail

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
    
    # Note: VCF output is the standard primary output, suitable for joint calling with GLnexus
    # GVCF format (with hom-ref calls) is NOT required for GLnexus joint calling
    
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
    
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
        local sample_name=$(basename "${sample}" .bam)
        local sample_vcf="${out_dir}/Clair3_Variant_Calling/${sample_name}.vcf.gz"
        local sample_log="${out_dir}/Clair3_Variant_Calling/${sample_name}.log"
        
        echo "[Clair3] Starting variant calling for sample: ${sample_name}" | tee -a "${sample_log}"
        echo "[Clair3] Input BAM: ${sample}" | tee -a "${sample_log}"
        echo "[Clair3] Sequencing platform: ${seq_platform}" | tee -a "${sample_log}"
        
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
        
        # Run Clair3 with GVCF output enabled
        if clair3 \
            --bam_fn="${sample}" \
            --ref_fn="${ref}" \
            --output="${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}" \
            --model_path="${model}" \
            --num_workers="${threads}" \
            --chunk_size=5000 \
            --chunk_num=5 \
            --sample_name="${sample_name}" \
            --enable_phasing \
            >> "${sample_log}" 2>&1; then
            
            echo "[Clair3] Successfully completed variant calling for sample: ${sample_name}" | tee -a "${sample_log}"
            
            # Check if the output VCF was created and is not empty
            if [[ -f "${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}/merge_output.vcf.gz" ]]; then
                # Move to final location
                mv "${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}/merge_output.vcf.gz" "${sample_vcf}"
                if [[ -f "${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}/merge_output.vcf.gz.tbi" ]]; then
                    mv "${out_dir}/Clair3_Variant_Calling/Intermediates/${sample_name}/merge_output.vcf.gz.tbi" "${sample_vcf}.tbi"
                fi
                echo "[Clair3] Output VCF: ${sample_vcf}" | tee -a "${sample_log}"
                
                # Note: Standard VCF output from Clair3 is suitable for joint calling with GLnexus
                # GLnexus is the recommended joint caller for long-read variant consolidation
                echo "[Clair3] This VCF is ready for joint calling via GLnexus (recommended for long-read data)" | tee -a "${sample_log}"
                echo "[Clair3] See: https://github.com/dnanexus-rnd/GLnexus" | tee -a "${sample_log}"
                
                # Optional: Index the VCF for use by downstream tools
                if bcftools index "${sample_vcf}" >> "${sample_log}" 2>&1; then
                    echo "[Clair3] Successfully indexed VCF" | tee -a "${sample_log}"
                else
                    echo "[WARN] Failed to index VCF file" >&2
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
                echo "[ERROR] Clair3 output VCF not found for sample: ${sample_name}" >&2
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
