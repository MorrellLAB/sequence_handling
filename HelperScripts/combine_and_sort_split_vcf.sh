#!/bin/bash

# This script uses GATK 4 to combine and sort split VCF files output from the
# Genotype_GVCFs handler. IMPORTANT: This is only necessary if the user doesn't plan
# to run the Create_HC_Subset handler and the Variant_Recalibrator handler.
# In other words, only run this if you are planning to hard filter your VCF file.

set -e
set -o pipefail

# Dependencies
module load java/openjdk-8_202
module load gatk/4.1.2

# User provided arguments are pulled from the Config file
# Provide the full filepath to the config file
CONFIG=
# Requested resources, modify as appropriate
QSUB="mem=22gb,nodes=1:ppn=4,walltime=24:00:00"
# Which MSI queue are we using?
QUEUE="small"

#--------------------------------------

function Combine_and_Sort_Split_VCFs() {
    local out=$1
    local project=$2
    out_subdir="${out}/Genotype_GVCFs/vcf_split_regions"
    ls ${out_subdir}/*.vcf > ${out_subdir}/temp-FileList.list
    # Since we parallelized across regions, combine split vcfs into one vcf file here
    # Check if we have specified a TMP directory (Note: pulled from global variable specified in config)
    if [ -z "${TMP}" ]; then
        echo "No temp directory specified. Proceeding without using temp directory."
        # This works for GATK 4, but not sure about GATK 3
        gatk SortVcf \
             -I ${out_subdir}/temp-FileList.list \
             -O ${out}/Genotype_GVCFs/${project}_raw_variants.vcf
    else
        echo "Proceed using temp directory: ${TMP}"
        # This works for GATK 4, but not sure about GATK 3
        gatk SortVcf \
                --TMP_DIR ${TMP} \
                -I ${out_subdir}/temp-FileList.list \
                -O ${out}/Genotype_GVCFs/${project}_raw_variants.vcf
    fi
    # Cleanup
    rm -f ${out_subdir}/temp-FileList.list
}

export -f Combine_and_Sort_Split_VCFs

# Submit process as a job
echo "source ${CONFIG} && Combine_and_Sort_Split_VCFs ${OUT_DIR} ${PROJECT}" | qsub -q "${QUEUE}" -l "${QSUB}" -m abe -M "${EMAIL}"
