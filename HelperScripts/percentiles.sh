#!/bin/bash

#   A function to write a percentile table
function percentiles() {
    local vcf="$1"
    local out="$2"
    local project="$3"
    local filtered="$4"
    local seqhand="$5"
    #   For very large VCF files, any step with vcftools will take a very long time
    #   (ex. for ~1.8TB vcf file, vcftools takes ~24 hours just to pull out GQ info)
    #   So, to manage large vcf file we will add a check to see if the file exists, if so we will continue to the next step
    #   Extract GQ information
    if [ -n "$(ls -A ${out}/Intermediates/${project}_${filtered}.GQ.FORMAT 2>/dev/null)" ]; then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/${project}_${filtered}.GQ.FORMAT"
    else
        vcftools --vcf "${vcf}" --extract-FORMAT-info GQ --out "${out}/Intermediates/${project}_${filtered}"
        echo "Finished generating ${out}/Intermediates/${project}_${filtered}.GQ.FORMAT"
    fi

    if [ -n "$(ls -A ${out}/Intermediates/${project}_${filtered}.GQ.matrix 2>/dev/null)" ]; then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/${project}_${filtered}.GQ.matrix"
    else
        cut -f 3- "${out}/Intermediates/${project}_${filtered}.GQ.FORMAT" | sed 1,1d > "${out}/Intermediates/${project}_${filtered}.GQ.matrix"
        echo "Finished generating ${out}/Intermediates/${project}_${filtered}.GQ.matrix"
    fi

    #   Extract DP per sample information
    if [ -n "$(ls -A ${out}/Intermediates/${project}_${filtered}.gdepth 2>/dev/null)" ]; then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/${project}_${filtered}.gdepth"
    else
        vcftools --vcf "${vcf}" --geno-depth --out "${out}/Intermediates/${project}_${filtered}"
        echo "Finished generating ${out}/Intermediates/${project}_${filtered}.gdepth"
    fi

    if [ -n "$(ls -A ${out}/Intermediates/${project}_${filtered}.gdepth.matrix 2>/dev/null)" ]; then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/${project}_${filtered}.gdepth.matrix"
    else
        sed 1d "${out}/Intermediates/${project}_${filtered}.gdepth" | cut -f 3- > "${out}/Intermediates/${project}_${filtered}.gdepth.matrix"
        sed -i -e 's/-1/NA/g' "${out}/Intermediates/${project}_${filtered}.gdepth.matrix" # Change -1 to NA
        echo "Finished generating ${out}/Intermediates/${project}_${filtered}.gdepth.matrix"
    fi

    #   Call the R script to write the percentiles table
    Rscript "${seqhand}/HelperScripts/percentiles.R" "${out}/Intermediates/${project}_${filtered}.GQ.matrix" "${out}/Intermediates/${project}_${filtered}.gdepth.matrix" "${out}/Percentile_Tables" "${project}" "${filtered}"
    echo "Finished generating the files:"
    echo "${out}/Percentile_Tables/${project}_${filtered}_GQ.txt"
    echo "${out}/Percentile_Tables/${project}_${filtered}_DP_per_sample.txt"
}

#   Export the function
export -f percentiles