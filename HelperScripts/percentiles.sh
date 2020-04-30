#!/bin/bash

#   A function to write a percentile table
function percentiles() {
    local vcf="$1"
    local out="$2"
    local project="$3"
    local filtered="$4"
    local seqhand="$5"
    local temp_dir="$6"
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

    #   Generate the percentiles tables
    #   Prepare flattened matrix and remove missing values denoted by "."
    #   For .GQ.matrix
    #   Check if file exists, if so proceed to next step
    if [ -n "$(ls -A ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt 2>/dev/null)" ]
    then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt"
    else
        # Missing "." values will be removed
        cat ${out}/Intermediates/${project}_${filtered}.GQ.matrix | tr '\t' '\n' | grep -vw '\.' > ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt
        echo "Finished generating temp_flattened_${project}_${filtered}.GQ.matrix.txt"
    fi
    #   For .gdepth.matrix, missing values denoted by "NA"
    if [ -n "$(ls -A ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt 2>/dev/null)" ]
    then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt"
    else
        # Missing "NA" values will be removed
        cat ${out}/Intermediates/${project}_${filtered}.gdepth.matrix | tr '\t' '\n' | grep -vw "NA" > ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt
        echo "Finished generating temp_flattened_${project}_${filtered}.gdepth.matrix.txt"
    fi
    #   Use the flattened matrix to save on memory usage and reduce processing done inside R
    Rscript "${seqhand}/HelperScripts/percentiles.R" "${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt" "${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt" "${out}/Percentile_Tables" "${project}" "${filtered}"
    echo "Finished calculating percentiles, see .log file in the same directory for additional info."
}

#   Export the function
export -f percentiles
