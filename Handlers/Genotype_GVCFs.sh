#!/bin/bash

#   This script creates VCF files for each 
#   chromosome part using GVCFs as input.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_GenotypeGVCFs.job

set -o pipefail

#   What are the dependencies for Genotype_GVCFs?
declare -a Genotype_GVCFs_Dependencies=(java)

#   A function to genotype the GVCFs without PBS
function Genotype_GVCFs() {
    local sample_list="$1" # What is our sample list? Or filename of combined GVCF with gatk4
    local out="$2"/Genotype_GVCFs # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local ploidy="$6" # What is the sample ploidy?
    local memory="$7" # How much memory can java use?
    local dict="$8" # Where is the reference dictionary?
    local maxarray="$9" # What is the maximum array index?
    local gatkVer="$10"  # GATK version  3 or 4
    local scaffolds="${11}" # Where is the scaffolds intervals file?
    local seqs_list=()
    if [[ "${maxarray}" -ne 0 ]]; then
    	seqs_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names
    fi
    # a chromosome part name is used as the output name, or "custom_intervals" for the scaffolds in CUSTOM_INTERVAL
    out_name_arr=("${seqs_list[@]}")    
    if [[ ! -z "${scaffolds}" ]]; then # custom interval is appended
        seqs_list+=("${scaffolds}")
        out_name_arr+=("custom_intervals")
    fi
    if [ ${#seqs_list[@]} -ne $[${maxarray} + 1] ]; then echo "Check NUM_CHR and CUSTOM_INTERVAL, it doesn't match with entries in the reference" >&2; exit 1; fi
    if [[ "$gatkVer" == 3 ]]; then
        analysisTypeOpt="-T GenotypeGVCFs"
        outFlag="-o"
        ploidyFlag="--sample_ploidy"
    else
        analysisTypeOpt="GenotypeGVCFs"
        outFlag="-O"
        # note that with GATK 4 the documentation says, no need to specify this
        ploidyFlag="--sample-ploidy"
    fi
    if [[ "$gatkVer" == 3 ]]; then
        declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format	
        #   Put the samples into a format that GATK can read	
        GATK_IN=()	
        for s in "${sample_array[@]}"		 
        do	    
            GATK_IN+=(-V $s)	    
        done
    else
        GATK_IN=("-V ${sample_list}")
    fi
    #   Make sure the out directory exists
    mkdir -p "${out}"
    if [[ "$USE_PBS" == "true" ]]; then
        #   What region of the genome are we working on currently?	
        local current="${seqs_list[${PBS_ARRAYID}]}"	
        local out_name="${out_name_arr[${PBS_ARRAYID}]}"	
        #   Run GATK using the parameters given
        (set -x; java -Xmx"${memory}" -jar "${gatk}" \
            "${analysisTypeOpt}"  \
            -R "${reference}" \
            -L "${current}" \
            "${GATK_IN[@]}" \
            --heterozygosity "${heterozygosity}" \
            "${ploidyFlag}" "${ploidy}" \
            "${outFlag} ${out}/${out_name}.vcf")
    else
        set -x; parallel java -Xmx"${memory}" -jar "${gatk}" \
        "${analysisTypeOpt}"  \
            -R "${reference}" \
            -L {1} \
            "${GATK_IN[@]}" \
            --heterozygosity "${heterozygosity}" \
            "${ploidyFlag}" "${ploidy}" \
            "${outFlag} ${out}/{2}.vcf" ::: "${seqs_list[@]}" :::+ "${out_name_arr[@]}"	
    fi   
}

#   Export the function
export -f Genotype_GVCFs

#   A function to combine individual GVCFs into a single GVCF required for GATK4
# This is slow, so use GenomicsDBImport().
# keeping it here in case there is a situation when GenomicsDBImport doesn't work
function Combine_GVCFs() {
    local sample_list="$1" # What is our sample list?
    local outFileName="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
	    GATK_IN+=(-V $s)
    done    
    # get the directory for the outputFile, and create it if it's missing
    local outDir=$(dirname "${outFileName}")
    mkdir -p "${outDir}"
    set -x; gatk CombineGVCFs \
		 -R "${reference}"\
		 "${GATK_IN[@]}" \
		 -O $outFileName  
}

export -f Combine_GVCFs

# Note with a large sample size intervals might need to be further subdivided to avoid memory problems
function GenomicsDBImport() {
    local sample_list="$1" # What is our sample list?
    local outFileName="$2" # Where are we storing our results?
    local reference="$3" # Where is the reference sequence?
    local type="$4"  # 'WGS' or 'targeted'
    local intvlBED="$5"  # put NA for WGS, BED file name for interval of tageted sequencing
    local tmp="$6"    # temp directory
    # get the directory for the outputFile, and create it if it's missing
    local outDir=$(dirname "${outFileName}")
    mkdir -p "${outDir}"
    if ! [[ -s "${reference}" ]]; then echo "Cannot find readable reference genome, exiting..." >&2; exit 31; fi # Make sure it exists
    # create a file which list the chromosomes for -L option
    # In case of targeted sequences, GenomicDBImport doesn't work if hundreds of intervals are given
    # So, using the chromosomes which contain the targets are used.
    if [[ "${type}" == "targeted" ]]; then
        if ! [[ -s "${intvlBED}" ]]; then echo "Cannot find readable bed file for the target region, exiting..." >&2; exit 31; fi # Make sure it exists
        cut -f 1 "${intvlBED}" | sort | uniq > "${outDir}/intervals.list"
    else	
        local dict="${reference%.*}.dict" # replace suffix (.fa or .fasta) with .dict	
        # checkDict and creatDict must have been called in sequence_handling, but just checking again	
        if ! [[ -s "${dict}" ]]; then echo "Cannot find readable reference dict genome (or bed file), exiting..." >&2; exit 31; fi # Make sure it exists	
        chrom_list=($(cut -f 2 ${dict} | grep -E '^SN' | cut -f 2 -d ':')) # Make an array of chromosome part names	
        printf '%s\n' "${chrom_list[@]}" > "${outDir}/intervals.list"
    fi
    local mergeIntvl=""
    if [ $(cat "${outDir}/intervals.list"|wc -l) -gt  500 ]; then
        # When there are many small intervals (e.g exomes), following option increases performance.
        mergeIntvl="--merge-input-intervals"
    fi
    declare -a sample_array=($(grep -E ".g.vcf" "${sample_list}")) # Put the sample list into array format
    #   Put the samples into a format that GATK can read
    GATK_IN=()
    for s in "${sample_array[@]}"
    do
	    GATK_IN+=(-V $s)
    done    
    if [ -n "$tmp" ] ; then
	    tmp="--tmp-dir=${tmp}"
    fi    
    set -x; gatk GenomicsDBImport \
		 -R "${reference}" \
		 "${GATK_IN[@]}" \
		 -L "${outDir}/intervals.list" \
		 ${mergeIntvl} \
		 "$tmp" \
		 --genomicsdb-workspace-path $outFileName
    rm -f "${outDir}/intervals.list"
}
