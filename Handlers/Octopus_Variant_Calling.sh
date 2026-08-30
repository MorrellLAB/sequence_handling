#!/bin/bash

# Run Octopus variant calling as an independent short-read branch.
# Supports individual, population, and trio calling models.
#
# For "trio" mode, sample_list stays a flat list of BAMs (same format
# as individual/population mode, and the same file you'd point
# OCTOPUS_BAM_LIST at). Two of the sample names in that list are
# designated as the maternal and paternal samples via
# maternal_sample/paternal_sample (OCTOPUS_MATERNAL_SAMPLE /
# OCTOPUS_PATERNAL_SAMPLE in Config_fastp). Every other BAM in the
# list is treated as a progeny and gets its own trio call against
# those same two parent BAMs.
#
# Sample names are derived from each BAM's basename (minus .bam),
# matching the convention used elsewhere in this handler.

set -euo pipefail

declare -a Octopus_Variant_Calling_Dependencies=(octopus bcftools)

function Octopus_Variant_Calling() {
    local sample_list="$1"
    local out_dir="$2"
    local cohort_name="$3"
    local ref="$4"
    local error_model="$5"
    local calling_model="$6"
    local ploidy="$7"
    local threads="$8"
    local maternal_sample="${9:-}"
    local paternal_sample="${10:-}"

    if [[ ! -f "${sample_list}" ]]; then
        echo "[ERROR] BAM list not found: ${sample_list}" >&2
        return 1
    fi

    if [[ ! -f "${ref}" ]]; then
        echo "[ERROR] Reference genome not found: ${ref}" >&2
        return 1
    fi

    if [[ ! "${cohort_name}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "[ERROR] OCTOPUS_COHORT_NAME must contain only letters, numbers, periods, underscores, or hyphens." >&2
        return 1
    fi

    if [[ -z "${error_model}" ]]; then
        echo "[ERROR] OCTOPUS_ERROR_MODEL must specify a library preparation and sequencer, for example PCR.NOVASEQ." >&2
        return 1
    fi

    if [[ "${calling_model}" != "individual" && "${calling_model}" != "population" && "${calling_model}" != "trio" ]]; then
        echo "[ERROR] OCTOPUS_CALLING_MODEL must be individual, population, or trio." >&2
        return 1
    fi

    mapfile -t sample_array < <(grep -E '\.bam$' "${sample_list}")

    if [[ ${#sample_array[@]} -eq 0 ]]; then
        echo "[ERROR] No BAM files found in list: ${sample_list}" >&2
        return 1
    fi

    for sample in "${sample_array[@]}"; do
        if [[ ! -f "${sample}" ]]; then
            echo "[ERROR] BAM file not found: ${sample}" >&2
            return 1
        fi
    done

    local cohort_out_dir="${out_dir}/Octopus_Variant_Calling/${cohort_name}"
    mkdir -p "${cohort_out_dir}"

    # ------------------------------------------------------------------
    # Trio mode
    # ------------------------------------------------------------------
    if [[ "${calling_model}" == "trio" ]]; then
        if [[ -z "${maternal_sample}" || -z "${paternal_sample}" ]]; then
            echo "[ERROR] OCTOPUS_MATERNAL_SAMPLE and OCTOPUS_PATERNAL_SAMPLE must both be set for trio calling." >&2
            return 1
        fi

        local maternal_bam="" paternal_bam=""
        local -a progeny_array=()

        for sample in "${sample_array[@]}"; do
            local sample_name
            sample_name=$(basename "${sample}" .bam)
            if [[ "${sample_name}" == "${maternal_sample}" ]]; then
                maternal_bam="${sample}"
            elif [[ "${sample_name}" == "${paternal_sample}" ]]; then
                paternal_bam="${sample}"
            else
                progeny_array+=("${sample}")
            fi
        done

        if [[ -z "${maternal_bam}" ]]; then
            echo "[ERROR] OCTOPUS_MATERNAL_SAMPLE (${maternal_sample}) not found among BAMs in ${sample_list}." >&2
            return 1
        fi

        if [[ -z "${paternal_bam}" ]]; then
            echo "[ERROR] OCTOPUS_PATERNAL_SAMPLE (${paternal_sample}) not found among BAMs in ${sample_list}." >&2
            return 1
        fi

        if [[ ${#progeny_array[@]} -eq 0 ]]; then
            echo "[ERROR] No progeny BAMs found in ${sample_list} (only the maternal/paternal samples were present)." >&2
            return 1
        fi

        for sample in "${progeny_array[@]}"; do
            local progeny_name
            progeny_name=$(basename "${sample}" .bam)

            local trio_vcf="${cohort_out_dir}/${progeny_name}.trio.vcf.gz"
            local trio_log="${cohort_out_dir}/${progeny_name}.trio.log"

            echo "[Octopus] Starting trio variant calling for progeny: ${progeny_name} (mother=${maternal_sample}, father=${paternal_sample})" | tee -a "${trio_log}"
            echo "[Octopus] Version: $(octopus --version)" | tee -a "${trio_log}"
            echo "[Octopus] Sequence error model: ${error_model}" | tee -a "${trio_log}"

            octopus \
                --reference "${ref}" \
                --reads "${maternal_bam}" "${paternal_bam}" "${sample}" \
                --maternal-sample "${maternal_sample}" \
                --paternal-sample "${paternal_sample}" \
                --output "${trio_vcf}" \
                --sequence-error-model "${error_model}" \
                --organism-ploidy "${ploidy}" \
                --threads "${threads}" \
                >> "${trio_log}" 2>&1

            if [[ ! -f "${trio_vcf}" ]]; then
                echo "[ERROR] Octopus trio output VCF not found for progeny: ${progeny_name}" >&2
                return 1
            fi

            bcftools index --force "${trio_vcf}" >> "${trio_log}" 2>&1
            echo "[Octopus] Output VCF: ${trio_vcf}" | tee -a "${trio_log}"
        done
        return 0
    fi

    # ------------------------------------------------------------------
    # Population mode
    # ------------------------------------------------------------------
    if [[ "${calling_model}" == "population" ]]; then
        local population_vcf="${cohort_out_dir}/population.vcf.gz"
        local population_log="${cohort_out_dir}/population.log"

        echo "[Octopus] Starting population variant calling for cohort: ${cohort_name}" | tee -a "${population_log}"
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

    # ------------------------------------------------------------------
    # Individual mode
    # ------------------------------------------------------------------
    for sample in "${sample_array[@]}"; do
        local sample_name
        sample_name=$(basename "${sample}" .bam)
        local sample_vcf="${cohort_out_dir}/${sample_name}.vcf.gz"
        local sample_log="${cohort_out_dir}/${sample_name}.log"

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
