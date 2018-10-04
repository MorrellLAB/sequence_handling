#!/bin/sh

set -o pipefail

#   What are the dependencies for Read_Mapping?
declare -a NP_Read_Mapping_Dependencies=(minimap2)

#   A function to run read mapping using minimap2
function NP_Read_Mapping() {
    local sample_list="$1" # The list of paths to samples to be mapped
    local out="$2/Read_Mapping" # The location to output mapped files
    local output_format="$3" # What file type should be outputted?
    local reference="$4" # The reference genome to map to
    local reference_index="$5" # The index to use for the reference genome
    local bandwidth="$6" # Mapping parameter
    local matching_score="$7" # Mapping parameter
    local mismatch_penalty="$8" # Mapping parameter
    local gap_open_penalty="${9}" # Mapping parameter
    local gap_extension_penalty="${10}" # Mapping parameter
    local z_drop_score="${11}" # Mapping parameter
    local minimal_peak_DP_score="${12}" # Mapping parameter
    local threads="${13}" # Mapping parameter
    #   Make sure the output directory exists
    mkdir -p "${out}" 
    #   Make an array of samples
    declare -a sample_array=($(grep -E ".fastq|.fastq.gz|.fasta|.fasta.gz|.fa|.fq|.fa.gz|.fq.gz" "${sample_list}"))
    #   Get which sample in the list we are working on
    local sample="${sample_array[${PBS_ARRAYID}]}"
    #   Get the name of the sample without the path in front
    #   If it is gzipped, this will strip off the .gz
    local sample_name=$(basename "${sample}" .gz)
    #   Get the name of the sample without the file extension, whatever it is
    local base_name="${sample_name%.*}"
    #   Convert the output format into a setting that minimap2 can understand
    if [[ "${output_format}" == "SAM" ]]
    then
        local minimap_output=(-a)
    else
        local minimap_output=()
    fi
    #   Run minimap2
    (set -x; minimap2 \
        -r "${bandwidth}" \
        -A "${matching_score}" \
        -B "${mismatch_penalty}" \
        -O "${gap_open_penalty}" \
        -E "${gap_extension_penalty}" \
        -z "${z_drop_score}" \
        -s "${minimal_peak_DP_score}" \
        -t "${threads}" \
        -L \
        "${minimap_output}" \
        "${reference_index}" \
        "${sample}" > "${out}/${base_name}.sam" )
}

export -f NP_Read_Mapping
