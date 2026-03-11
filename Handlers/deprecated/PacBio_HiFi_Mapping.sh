#!/bin/bash

set -o pipefail

# Legacy compatibility wrapper. Mapping implementation is centralized in Long_Read_Mapping.
declare -a PacBioHiFiMapping_DEPENDENCIES=(minimap2 samtools)

function PacBio_HiFi_Mapping() {
    local project_dir="$1"
    local handler="$2"
    local config="$3"
    source "${config}"

    if [[ -z "${THREADS:-}" ]]; then
        THREADS=8
    fi

    source "${SEQUENCE_HANDLING}/Handlers/Long_Read_Mapping.sh"
    Long_Read_Mapping "${PACBIO_SAMPLES}" "map-hifi" "${PROJECT}" "PACBIO" "${OUT_DIR}" "${REF_GEN}" "${THREADS}"
}

export -f PacBio_HiFi_Mapping