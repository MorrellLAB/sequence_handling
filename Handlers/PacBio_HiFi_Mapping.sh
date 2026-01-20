#!/bin/bash

#   This script maps PacBio HiFi reads to a reference genome
#   using minimap2 with optimal settings for HiFi data.
#   It is designed to process multiple samples efficiently.

set -e
set -o pipefail

#   What are the dependencies for this handler?
declare -a PacBioHiFiMapping_DEPENDENCIES=(minimap2 samtools)

#   A function to check the dependencies
function check_PacBioHiFiMapping_Dependencies() {
    module_Dependencies PacBioHiFiMapping_DEPENDENCIES
}

#   A function to run the handler
function PacBio_HiFi_Mapping() {
    #   Get the dependencies
    check_PacBioHiFiMapping_Dependencies

    #   User specified variables
    local PROJECT_DIR="$1"
    local HANDLER="$2"
    local CONFIG="$3"

    #   Source the config file
    source "${CONFIG}"

    #   Set up our I/O variables
    local READS_DIR="${OUT_DIR}/SAM_Processing"
    local OUT_DIR="${OUT_DIR}/PacBio_HiFi_Mapping"
    mkdir -p "${OUT_DIR}"
    mkdir -p "${OUT_DIR}/logs"

    #   Get the reference genome path
    local REFERENCE="${REF_GEN}"
    if [[ ! -f "${REFERENCE}" ]]; then
        echo "Reference genome not found at ${REFERENCE}"
        exit 1
    fi

    #   Create a samples array from the FASTQ list
    declare -a SAMPLES=()
    while read -r FASTQ_FILE; do
        if [[ -f "${FASTQ_FILE}" ]]; then
            SAMPLE_NAME=$(basename "${FASTQ_FILE}" | sed 's/\.fastq\.gz$//')
            SAMPLES+=("${SAMPLE_NAME}:${FASTQ_FILE}")
        fi
    done < "${PACBIO_SAMPLES}"

    echo "Found ${#SAMPLES[@]} samples to process"

    #   Create the job submission scripts
    for SAMPLE_PAIR in "${SAMPLES[@]}"; do
        IFS=':' read -r SAMPLE_NAME FASTQ_FILE <<< "${SAMPLE_PAIR}"
        
        #   Create job script
        echo "#!/bin/bash
#PBS -N PBMap_${SAMPLE_NAME}
#PBS -q ${QUEUE}
#PBS -m abe
#PBS -M ${EMAIL}
#PBS -l walltime=${WALLTIME}
#PBS -l nodes=1:ppn=${MAPPING_THREADS}
#PBS -l mem=${MAPPING_MEM}gb
set -e
set -o pipefail

module load minimap2
module load samtools

SAMPLE_DIR=\"${OUT_DIR}/${SAMPLE_NAME}\"
mkdir -p \"\${SAMPLE_DIR}\"

echo \"Mapping PacBio HiFi reads for sample ${SAMPLE_NAME}\"
minimap2 -ax map-hifi \\
         -t ${MAPPING_THREADS} \\
         -R \"@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\" \\
         ${REFERENCE} \\
         ${FASTQ_FILE} | \\
samtools sort -@ ${MAPPING_THREADS} \\
             -m 4G \\
             -o \"\${SAMPLE_DIR}/${SAMPLE_NAME}.sorted.bam\"

# Index BAM file
samtools index \"\${SAMPLE_DIR}/${SAMPLE_NAME}.sorted.bam\"

# Generate mapping statistics
samtools flagstat \"\${SAMPLE_DIR}/${SAMPLE_NAME}.sorted.bam\" > \"\${SAMPLE_DIR}/${SAMPLE_NAME}.flagstat\"
samtools stats \"\${SAMPLE_DIR}/${SAMPLE_NAME}.sorted.bam\" > \"\${SAMPLE_DIR}/${SAMPLE_NAME}.stats\"

echo \"Completed processing sample: ${SAMPLE_NAME}\"
" > "${OUT_DIR}/logs/PBMap_${SAMPLE_NAME}.sh"

        #   Make the job script executable
        chmod +x "${OUT_DIR}/logs/PBMap_${SAMPLE_NAME}.sh"
        
        #   Submit the job
        if [[ "${SCHEDULER}" == "PBS" ]]; then
            qsub "${OUT_DIR}/logs/PBMap_${SAMPLE_NAME}.sh"
        elif [[ "${SCHEDULER}" == "SLURM" ]]; then
            sbatch "${OUT_DIR}/logs/PBMap_${SAMPLE_NAME}.sh"
        else
            echo "Unknown scheduler: ${SCHEDULER}"
            exit 1
        fi
    done
}