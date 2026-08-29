#!/bin/bash

# Run Octopus individual calling as an independent short-read branch.

set -euo pipefail

declare -a Octopus_Variant_Calling_Dependencies=(octopus bcftools)

function Octopus_Variant_Calling() {
    local sample_list="$1"
    local out_dir="$2"
    local ref="$3"
    local error_model="$4"
    local calling_model="$5"
    local ploidy="$6"
    local threads="$7"

    if [[ ! -f "${sample_list}" ]]; then
        echo "[ERROR] BAM list not found: ${sample_list}" >&2
        return 1
    fi
    if [[ ! -f "${ref}" ]]; then
        echo "[ERROR] Reference genome not found: ${ref}" >&2
        return 1
    fi
    if [[ -z "${error_model}" ]]; then
        echo "[ERROR] OCTOPUS_ERROR_MODEL must specify a library preparation and sequencer, for example PCR.NOVASEQ." >&2
        return 1
    fi
    if [[ "${calling_model}" != "individual" && "${calling_model}" != "population" ]]; then
        echo "[ERROR] OCTOPUS_CALLING_MODEL must be individual or population." >&2
        return 1
    fi

    mapfile -t sample_array < <(grep -E '\.bam$' "${sample_list}")
    if [[ ${#sample_array[@]} -eq 0 ]]; then
        echo "[ERROR] No BAM files found in list: ${sample_list}" >&2
        return 1
    fi

    mkdir -p "${out_dir}/Octopus_Variant_Calling"

    if [[ "${calling_model}" == "population" ]]; then
        local population_vcf="${out_dir}/Octopus_Variant_Calling/population.vcf.gz"
        local population_log="${out_dir}/Octopus_Variant_Calling/population.log"

        echo "[Octopus] Starting population variant calling" | tee -a "${population_log}"
        echo "[Octopus] Version: $(octopus --version)" | tee -a "${population_log}"
        echo "[Octopus] Input BAM list: ${sample_list}" | tee -a "${population_log}"
        echo "[Octopus] Sequence error model: ${error_model}" | tee -a "${population_log}"

        octopus \
            --reference "${ref}" \
            --reads-file "${sample_list}" \
            --output "${population_vcf}" \
            --caller population \
            --sequence-error-model "${error_model}" \
            --organism-ploidy "${ploidy}" \
            --threads "${threads}" \
            >> "${population_log}" 2>&1

        if [[ ! -f "${population_vcf}" ]]; then
            echo "[ERROR] Octopus population output VCF not found." >&2
            return 1
        fi

        bcftools index --force "${population_vcf}" >> "${population_log}" 2>&1
        echo "[Octopus] Output VCF: ${population_vcf}" | tee -a "${population_log}"
        return 0
    fi

    for sample in "${sample_array[@]}"; do
        local sample_name
        sample_name=$(basename "${sample}" .bam)
        local sample_vcf="${out_dir}/Octopus_Variant_Calling/${sample_name}.vcf.gz"
        local sample_log="${out_dir}/Octopus_Variant_Calling/${sample_name}.log"

        if [[ ! -f "${sample}" ]]; then
            echo "[ERROR] BAM file not found: ${sample}" >&2
            return 1
        fi

        echo "[Octopus] Starting individual variant calling for sample: ${sample_name}" | tee -a "${sample_log}"
        echo "[Octopus] Version: $(octopus --version)" | tee -a "${sample_log}"
        echo "[Octopus] Input BAM: ${sample}" | tee -a "${sample_log}"
        echo "[Octopus] Reference genome: ${ref}" | tee -a "${sample_log}"
        echo "[Octopus] Sequence error model: ${error_model}" | tee -a "${sample_log}"

        octopus \
            --reference "${ref}" \
            --reads "${sample}" \
            --output "${sample_vcf}" \
            --caller individual \
            --sequence-error-model "${error_model}" \
            --organism-ploidy "${ploidy}" \
            --threads "${threads}" \
            >> "${sample_log}" 2>&1

        if [[ ! -f "${sample_vcf}" ]]; then
            echo "[ERROR] Octopus output VCF not found for sample: ${sample_name}" >&2
            return 1
        fi

        bcftools index --force "${sample_vcf}" >> "${sample_log}" 2>&1
        echo "[Octopus] Output VCF: ${sample_vcf}" | tee -a "${sample_log}"
    done
}

export -f Octopus_Variant_Calling