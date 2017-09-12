#!/bin/bash

#   This script generates summary statistics from a VCF file.

set -o pipefail

#   What are the dependencies for Variant_Analysis?
declare -a Variant_Analysis_Dependencies=(python parallel)

#   This function generates population statistics for 18 barley loci to be compared with Morrell et al. 2006
function Loci_Analysis() {
    local vcf="$1" # What is our vcf file?
    local bed="$2" # Which BED file are we looking at?
    local seqhand="$3" # Where is the sequence_handling directory located?
    local out="$4" # Where is the output directory?
    local project="$5" # What is the name of this project?
    #   Get the name of the loci from the bed file
    local loci=$(basename ${bed} .bed)
    #   Create a new VCF with only variants located in the loci of interest
    vcfintersect -b "${seqhand}/HelperScripts/18_barley_loci/${loci}.bed" "${vcf}" > "${out}/${loci}.vcf"
    #   Calculate the number of variants
    local n_vars=$(grep -v "#" "${out}/${loci}.vcf" | wc -l)
    #   If we actually have variants, do more analysis
    if [[ "${n_vars}" -ne 0 ]]
    then
        #   Convert the new VCF to Hudson table format
        python2 "${seqhand}/HelperScripts/VCF_To_Htable.py" "${out}/${loci}.vcf" > "${out}/${loci}.htable"
        #   Use molpopgen compute to calculate summary statistics for this loci
        compute -h "${out}/${loci}.htable" > "${out}/${loci}.stats"
        #   Cut out the statistics we care about
        local variants=$(grep "${loci}" "${out}/${loci}.stats" | tail -n 1 | cut -f 5)
        local theta_w=$(grep "${loci}" "${out}/${loci}.stats" | tail -n 1 | cut -f 12)
        local theta_pi=$(grep "${loci}" "${out}/${loci}.stats" | tail -n 1 | cut -f 13)
        local taj_d=$(grep "${loci}" "${out}/${loci}.stats" | tail -n 1 | cut -f 14)
        #   Remove the intermediate files
        rm "${out}/${loci}.htable"
        rm "${out}/${loci}.stats"
    else    #   Otherwise, set the parameters to "NA"
        local variants=0
        local theta_w="NA"
        local theta_pi="NA"
        local taj_d="NA"
    fi
    #   Print the statistics to the summary file, don't use n_vars because sometimes molpopgen doesn't load all the variants correctly
    echo -e "${loci}"'\t'"${variants}"'\t'"${theta_w}"'\t'"${theta_pi}"'\t'"${taj_d}" >> "${out}/${project}_variant_summary_unfinished.txt"
    #   Remove the intermediate files
    rm "${out}/${loci}.vcf"
}

#   Export the function
export -f Loci_Analysis

#   A function to run each analysis
function Variant_Analysis() {
    local vcf="$1" # What is our VCF file?
    local out="$2"/Variant_Analysis # Where are we storing our results?
    local seqhand="$3" # Where is the sequence_handling directory located?
    local project="$4" # What is the name of this project?
    local barley="$5" # Is this barley?
    #   Make sure the out directory exists
    mkdir -p "${out}" 
    #   If we have barley, do some additional analysis
    if [[ "${barley}" == true ]]
    then
        #   Generate a list of full file paths to the bed files included in sequence_handling
        "${seqhand}/HelperScripts/sample_list_generator.sh" ".bed" "${seqhand}/HelperScripts/18_barley_loci" "18_loci_beds.list"
        #   Call the loci analysis function in parallel
        parallel -v Loci_Analysis "${vcf}" {} "${seqhand}" "${out}" "${project}" :::: "${seqhand}/HelperScripts/18_barley_loci/18_loci_beds.list"
        #   Create the header for the summary file
        echo -e "Loci\tVariants\tTheta_W\tTheta_Pi\tTajimas_D" > "${out}/${project}_variant_summary.txt"
        #   Sort the summary file 
        sort "${out}/${project}_variant_summary_unfinished.txt" >> "${out}/${project}_variant_summary.txt"
        #   Remove the intermediate file
        rm "${out}/${project}_variant_summary_unfinished.txt"
    fi
}

#   Export the function
export -f Variant_Analysis
