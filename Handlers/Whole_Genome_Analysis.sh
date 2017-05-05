#!/bin/bash

#   This script analyzes whole genome sequencing data for 10x Genomics data.
#   Additional documentation can be found at:
#       https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/wgs

set -o pipefail

#   What are the dependencies of Whole_Genome_Analysis?
declare -a Whole_Genome_Analysis_Dependencies=(longranger java)

#   A function to run Whole Genome Whole_Genome_Analysis
function WGS_Analysis() {
    local prefix="$1" # What are our sample prefixes?
    local fastqDir="$2" # Where is our FASTQ directory?
    local referenceDir="$3" # Where is our compatible reference?
    local variantCaller="$4" # What variant caller are we using?
    local sex="$5" # What sex is our organism?
    local outDir="$6" # Where are we storing our results?
    local out="${outDir}"/Whole_Genome_Analysis # Full path to output directory
    # local vcf="$5" # Are we using a pre-called VCF file?

    #   Create out directory name
    sampleID=`basename ${out}`

    #   Go into output directory
    #   Long Ranger outputs results in current working directory
    cd ${out}
    #   Run longranger wgs
    longranger wgs --id=${sampleID} \
                   --sample=${prefix} \
                   --fastqs=${fastqDir} \
                   --reference=${referenceDir} \
                   --vcmode=${variantCaller} \
                   --sex=${sex}
}

export -f WGS_Analysis
