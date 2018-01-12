#!/bin/bash

#   This script generates summary statistics from a VCF file.

set -o pipefail

#   What are the dependencies for Variant_Analysis?
declare -a Variant_Analysis_Dependencies=(python3 parallel compute vcfintersect bcftools R vcftools)

#   This function generates population statistics for 18 barley loci to be compared with Morrell et al. 2006
function Loci_Analysis() {
    local vcf="$1" # What is our vcf file?
    local bed="$2" # Which BED file are we looking at?
    local seqhand="$3" # Where is the sequence_handling directory located?
    local out="$4" # Where is the output directory?
    local name="$5" # What is the name of the file?
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
        python3 "${seqhand}/HelperScripts/VCF_To_Htable.py" "${out}/${loci}.vcf" > "${out}/${loci}.htable"
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
    echo -e "${loci}"'\t'"${variants}"'\t'"${theta_w}"'\t'"${theta_pi}"'\t'"${taj_d}" >> "${out}/${name}_loci_summary_unfinished.txt"
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
    local barley="$4" # Is this barley?
    #   Make sure the out directory exists
    mkdir -p "${out}" 
    #   What's the name of the vcf file?
    local name=$(basename ${vcf} .vcf)
    #   Generate some pdf plots using bcftools, python-epd, and texlive
    cd "${out}"
    plot-vcfstats <(bcftools stats ${vcf}) -p plots
    #   Rename the combined pdf
    mv "${out}/plots-summary.pdf" "${out}/${name}_summary.pdf"
    #   Remove the pdf parts and junk files
    rm plots-counts_by_af.snps.pdf plots-counts_by_af.snps.png plots-tstv_by_af.0.pdf plots-tstv_by_af.0.png plots-tstv_by_qual.0.dat plots-tstv_by_qual.0.pdf plots-tstv_by_qual.0.png plots-plot.py plots-plot-vcfstats.log plots-substitutions.0.pdf plots-substitutions.0.png plots-summary.aux plots-summary.log plots-summary.tex
    #   Create a file with minor allele frequencies (MAF)
    python3 "${seqhand}/HelperScripts/VCF_MAF.py" "${vcf}" > "${out}/${name}_MAF.txt"
    #   Use R to plot the MAF file as a histogram
    Rscript "${seqhand}/HelperScripts/plot_maf.R" "${out}/${name}_MAF.txt" "${name}" "${out}/${name}_MAF.pdf"
    #   Calculate inbreeding coefficients and heterozygosity
    vcftools --vcf "${vcf}" --het --out "${out}/${name}_unsorted"
    echo -e "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF" > "${out}/${name}_heterozygosity.txt"
    tail -n +2 "${out}/${name}_unsorted.het" | sort -k 5 >> "${out}/${name}_heterozygosity.txt"
    rm "${out}/${name}_unsorted.het"
    #   Calculate missingness per individual
    vcftools --vcf "${vcf}" --missing-indv --out "${out}/${name}_unsorted"
    sort -rk 5 "${out}/${name}_unsorted.imiss" > "${out}/${name}_missingness.txt"
    rm "${out}/${name}_unsorted.imiss"
    #   If we have barley, do some additional analysis
    if [[ "${barley}" == true ]]
    then
        #   Check to see if the VCF is pseudomolecular positions or parts positions
        local positions=$(grep -v "#" ${vcf} | head -n 1 | grep "part")
        #   If it's parts positions, convert the VCF to pseudomolecular positions
        if ! [[ -z "${positions}" ]]
        then
            python3 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" "${vcf}" > "${out}/${name}_pseudo.vcf"
            local to_analyze="${out}/${name}_pseudo.vcf"
        else
            local to_analyze="${vcf}"
        fi
        #   Generate a list of full file paths to the bed files included in sequence_handling
        "${seqhand}/HelperScripts/sample_list_generator.sh" ".bed" "${seqhand}/HelperScripts/18_barley_loci" "18_loci_beds.list"
        #   Call the loci analysis function in parallel
        parallel -v Loci_Analysis "${to_analyze}" {} "${seqhand}" "${out}" "${name}" :::: "${seqhand}/HelperScripts/18_barley_loci/18_loci_beds.list"
        #   Create the header for the summary file
        echo -e "Loci\tVariants\tTheta_W\tTheta_Pi\tTajimas_D" > "${out}/${name}_loci_summary.txt"
        #   Sort the summary file 
        sort "${out}/${name}_loci_summary_unfinished.txt" >> "${out}/${name}_loci_summary.txt"
        #   Remove the intermediate file
        rm "${out}/${name}_loci_summary_unfinished.txt"
    fi
}

#   Export the function
export -f Variant_Analysis