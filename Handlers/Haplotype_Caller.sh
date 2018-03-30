#!/bin/bash

#   This script performs preliminary SNP calling
#   on a single BAM sample using GATK.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_HaplotypeCaller.job

set -o pipefail

#   What are the dependencies for Haplotype_Caller?
declare -a Haplotype_Caller_Dependencies=(java)

#   A function to run the SNP calling
function Haplotype_Caller() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Haplotype_Caller # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local memory="$6" # How much memory can java use?
    local qscores="$7" # Do we fix quality scores?
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
    local sample="${sample_array[${PBS_ARRAYID}]}" # Which sample are we working on currently?
    local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
    mkdir -p "${out}" # Make sure the out directory exists
    if [[ "${qscores}" == true ]]
    then
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
        --fix_misencoded_quality_scores \
	    -variant_index_parameter 128000)
    else
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
    fi
}

#   Export the function
export -f Haplotype_Caller