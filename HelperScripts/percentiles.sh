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

    #   Call the R script to write the percentiles table
    # Rscript "${seqhand}/HelperScripts/percentiles.R" "${out}/Intermediates/${project}_${filtered}.GQ.matrix" "${out}/Intermediates/${project}_${filtered}.gdepth.matrix" "${out}/Percentile_Tables" "${project}" "${filtered}"
    #   Prepare flattened matrix and remove missing values denoted by "."
    #   For .GQ.matrix
    #   Check if file exists, if so proceed to next step
    if [ -n "$(ls -A ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt 2>/dev/null)" ]
    then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt"
    else
        cat ${out}/Intermediates/${project}_${filtered}.GQ.matrix | tr '\t' '\n' | grep -vw '\.' > ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt
        echo "Finished generating temp_flattened_${project}_${filtered}.GQ.matrix.txt"
    fi
    #   For .gdepth.matrix, missing values denoted by "NA"
    if [ -n "$(ls -A ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt 2>/dev/null)" ]
    then
        echo "Proceeding to next step, this file exists: ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt"
    else
        cat ${out}/Intermediates/${project}_${filtered}.gdepth.matrix | tr '\t' '\n' | grep -vw "NA" > ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt
        echo "Finished generating temp_flattened_${project}_${filtered}.gdepth.matrix.txt"
    fi
    #   Use the flattened matrix to save on memory usage and reduce processing done inside R
    Rscript "${seqhand}/HelperScripts/percentiles.R" "${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt" "${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt" "${out}/Percentile_Tables" "${project}" "${filtered}"
    # #   Use datamash to calculate percentiles instead, should be more efficient than R or Python
    # percentile_arr=($(seq 1 100))
    # #   Temporary directory to save outputs prior to concatenating
    # mkdir -p ${out}/Intermediates/temp_percentiles
    # #   Set up function so we can run this in parallel
    # function datamash_percentile() {
    #     local p_num=$1
    #     local flattened_matrix=$2
    #     local suffix=$3
    #     local out_dir=$4
    #     # Use datamash to calculate percentile
    #     temp_percentile=$(datamash perc:${p_num} 1 < ${flattened_matrix})
    #     # Save to a file
    #     echo ${p_num}$'%' ${temp_percentile} | tr ' ' '\t' > ${out_dir}/Intermediates/temp_percentiles/temp_${p_num}_${suffix}.txt
    # }

    # export -f datamash_percentile

    # # Calculate percentiles in parallel for .GQ.matrix
    # #parallel --compress --tmpdir ${temp_dir} datamash_percentile {} ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt ${filtered}_GQ ${out} ::: ${percentile_arr[@]}
    # for i in ${percentile_arr[@]}
    # do
    #     datamash_percentile ${i} ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt ${filtered}_GQ ${out}
    # done
    # # Check that we have the expected number of output files and concatenate
    # if [ "$(ls ${out}/Intermediates/temp_percentiles/temp_*${filtered}_GQ.txt | wc -l)" -eq "100" ]; then
    #     echo "We have calculated percentiles 1-100 for *GQ.txt, concatenating files..."
    #     cat ${out}/Intermediates/temp_percentiles/temp_*${filtered}_GQ.txt | sort -V -k1,1 > ${out}/Percentile_Tables/${project}_${filtered}_GQ.txt
    #     # Cleanup
    #     rm ${out}/Intermediates/temp_percentiles/temp_*${filtered}_GQ.txt
    # else
    #     echo "Error when calculating percentiles in parallel for file ${out}/Intermediates/temp_flattened_${project}_${filtered}.GQ.matrix.txt"
    #     exit 32 # If something went wrong, exit
    # fi

    # # Calculate percentiles in parallel for .gdepth.matrix
    # #parallel --compress --tmpdir ${temp_dir} datamash_percentile {} ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt ${filtered}_DP_per_sample ${out} ::: ${percentile_arr[@]}
    # for i in ${percentile_arr[@]}
    # do
    #     datamash_percentile ${i} ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt ${filtered}_DP_per_sample ${out}
    # done
    # # Check that we have the expected number of output files and concatenate
    # if [ "$(ls ${out}/Intermediates/temp_percentiles/temp_*${filtered}_DP_per_sample.txt | wc -l)" -eq "100" ]; then
    #     echo "We have calculated percentiles 1-100 for *_DP_per_sample.txt, concatenating files..."
    #     cat ${out}/Intermediates/temp_percentiles/temp_*${filtered}_DP_per_sample.txt | sort -V -k1,1 > ${out}/Percentile_Tables/${project}_${filtered}_DP_per_sample.txt
    #     # Cleanup
    #     rm ${out}/Intermediates/temp_percentiles/temp_*${filtered}_DP_per_sample.txt
    # else
    #     echo "Error when calculating percentiles in parallel for file ${out}/Intermediates/temp_flattened_${project}_${filtered}.gdepth.matrix.txt"
    #     exit 32 # If something went wrong, exit
    # fi
    echo "Finished generating the files:"
    echo "${out}/Percentile_Tables/${project}_${filtered}_GQ.txt"
    echo "${out}/Percentile_Tables/${project}_${filtered}_DP_per_sample.txt"
}

#   Export the function
export -f percentiles
