#!/bin/sh

set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a NP_Adapter_Trimming_Dependencies=(python porechop)

#   A function to run adapter trimming
function NP_Adapter_Trimming() {
    local sample_list="$1" # The list of paths to samples to be trimmed
    local out="$2/Adapter_Trimming" # The location to output trimmed files
    local project="$3" # The name of the project that the samples belong to
    #   Make sure the output directory exists
    mkdir -p "${out}" 
    #   Get which sample in the list we are working on
    local sample="${sample_list[${PBS_ARRAYID}]}"
    #   Get the name of the sample without the path in front
    #   If it is gzipped, this will strip off the .gz
    local sample_name=$(basename "${sample}" .gz)
    #   Get the name of the sample without the file extension, whatever it is
    local base_name="${sample_name%.*}"
    #   Get the file extension only (or nothing for files without an extension)
    local extension=$([[ "${sample_name}" = *.* ]] && echo ".${sample_name##*.}" || echo '')
    #   Use the file extension to figure out how to name the output file
    case "${extension}" in
        .fasta | .fa )
            (set -x; porechop -i "${sample}" -o "${out}/${base_name}.fasta")
            ;;
        .fastq | .fq )
            (set -x; porechop -i "${sample}" -o "${out}/${base_name}.fastq")
            ;;
        * )
            # If the file isn't fastq or fasta, then we have a problem
            echo "ERROR: Unexpected file type encountered." >&2
            exit 1
            ;;
    esac
}

export -f NP_Adapter_Trimming