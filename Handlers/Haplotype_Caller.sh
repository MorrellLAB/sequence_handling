#!/bin/bash

#   This script performs preliminary SNP calling
#   on a single BAM sample using GATK.

set -o pipefail

#   What are the dependencies for Haplotype_Caller?
declare -a Haplotype_Caller_Dependencies=(java)

#   A function to run the SNP calling
function Haplotype_Caller() {
    local sample_list="$1" # What is our sample list?
    local sample="${sample_list[${PBS_ARRAYID}]}" # Which sample are we working on currently?
    local out="$2"/Haplotype_Caller # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local memory="$6" # How much memory can java use?
    local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
    mkdir -p "${out}" # Make sure the out directory exists
    #   Run GATK using the parameters given
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
	    -T HaplotypeCaller \
	    -R "${reference}" \
	    -I "${sample}" \
	    -o "${out}/${sample_name}_RawGLs.g.vcf" \
	    -nct 4 \
	    --genotyping_mode DISCOVERY \
	    --heterozygosity "${heterozygosity}" \
	    --emitRefConfidence GVCF \
	    -variant_index_type LINEAR \
	    -variant_index_parameter 128000)
}

#   Export the function
export -f Haplotype_Caller