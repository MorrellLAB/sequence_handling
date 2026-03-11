#!/bin/bash

set -o pipefail

# Legacy compatibility wrapper. New long-read mapping is centralized in Long_Read_Mapping.
declare -a NP_Read_Mapping_Dependencies=(minimap2 samtools)

function NP_Read_Mapping() {
    local sample_list="$1"
    local out_dir="$2"
    local reference="$4"
    local threads="${13}"

    if [[ -z "${threads}" ]]; then
        threads=8
    fi

    # Preserve historical NP behavior by defaulting to ONT preset unless user overrides.
    local preset="${MINIMAP2_PRESET:-map-ont}"
    if [[ -n "${MINIMAP2_READ_TYPE:-}" ]]; then
        case "${MINIMAP2_READ_TYPE}" in
            ONT_Q20|ONT-Q20) preset="lr:hq" ;;
            ONT) preset="map-ont" ;;
        esac
    fi

    source "${SEQUENCE_HANDLING}/Handlers/Long_Read_Mapping.sh"
    Long_Read_Mapping "${sample_list}" "${preset}" "${PROJECT}" "ONT" "${out_dir}" "${reference}" "${threads}"
}

export -f NP_Read_Mapping
