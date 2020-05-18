#!/bin/bash
#PBS -l mem=16gb,nodes=1:ppn=4,walltime=24:00:00
#PBS -m abe
#PBS -M username@umn.edu
#PBS -q mesabi

set -e
set -o pipefail

# This helper script takes a list of VCF files and converts them to hdf5 format.
# To run, fill out filepaths under the section "User provided arguments" and submit as a job
# on MSI. This script outputs .h5 files in the same directory as the vcf file(s).

# Dependencies
module load python3/3.6.3_anaconda5.0.1
module load parallel/20180922

# User provided arguments
# Full filepath to list of vcf files to convert to hdf5 format
VCF_LIST=
# Full filepath to our sequence_handling directory
SEQ_HANDLING_DIR=

#-----------------------

# Source directory where vcf_to_h5.py script is located
# Directory where custom script is located
export PATH=${PATH}:${SEQ_HANDLING_DIR}/HelperScripts

function vcf2h5() {
    local vcf=$1
    # Convert VCF to hdf5
    vcf_to_h5.py ${vcf}
}

export -f vcf2h5

# Store vcf list in array
VCF_ARR=($(cat ${VCF_LIST}))

# Run program
parallel vcf2h5 {} ::: ${VCF_ARR[@]}
