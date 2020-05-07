#!/bin/bash

#   This script creates a high-confidence
#   subset of variants to use in variant recalibration.

set -e
set -o pipefail

#   What are the dependencies for Create_HC_Subset?
declare -a Create_HC_Subset_Dependencies=(parallel vcftools R vcfintersect python3 bcftools)

function Create_HC_Subset_GATK4() {
    local vcf_list="$1" # What is our sample list?
    local out="$2" # Where are we storing our results?
    local barley="$3" # Is this barley?
    local project="$4" # What is the name of this project?
    local seqhand="$5" # Where is sequence_handling located?
    local qual_cutoff="$6" # What is the quality cutoff?
    local gq_cutoff="$7" # What is the genotyping quality cutoff?
    local dp_per_sample_cutoff="$8" # What is the DP per sample cutoff?
    local max_het="$9" # What is the maximum number of heterozygous samples?
    local max_bad="${10}" # What is the maximum number of bad samples?
    local temp_dir="${11}" # Where can we store temporary files while running parallel?
    mkdir -p ${out}/Create_HC_Subset \
             ${out}/Create_HC_Subset/Intermediates \
             ${out}/Create_HC_Subset/Percentile_Tables \
             ${temp_dir}
    # Note: gzipping large files takes a long time and bcftools concat (based on a quick time comparison) runs faster with uncompressed vcf files. So, we will not gzip the VCF file.
    # 1. Use bcftools to concatenate all the split VCF files into a raw VCF file
    # Many users will want to back this file up since it is the most raw form of the SNP calls
    # Check if files are already concatenated, if so skip this time consuming step
    if [ -n "$(ls -A ${out}/Create_HC_Subset/${project}_concat_raw.vcf 2>/dev/null)" ]; then
        echo "File exists, proceed with next step using file: ${out}/Create_HC_Subset/${project}_concat_raw.vcf"
    else
        echo "File doesn't exist, concatenate vcf."
        # File doesn't exist, concatenate vcf
        bcftools concat -f ${vcf_list} > "${out}/Create_HC_Subset/${project}_concat_raw.vcf"
    fi

    # 2. Filter out indels using vcftools
    # Check if our concatenated file is larger than 200GB (equivalent to 214,748,364,800 bytes) in size, if
    # so we will work with the split vcf files and remove indels in parallel to speed up the processing time.
    # First, get the size of our concatenated VCF
    vcf_size=$(stat -c %s ${out}/Create_HC_Subset/${project}_concat_raw.vcf)
    if [ ${vcf_size} -gt "214748364800" ]; then
        echo "File is larger than 200GB, we will run remove indels on the split vcf files and concatenate after."
        # Check if we have already filtered out indels, if so skip and proceed to next step
        if [ -n "$(ls -A ${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf 2>/dev/null)" ]; then
            echo "Already filtered out indels, proceed with next step using file: ${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf"
        else
            echo "Removing indels from split vcf files output from Genotype_GVCFs handler."
            # Make a temporary directory to store these intermediate files
            mkdir -p "${out}/Create_HC_Subset/Intermediates/temp_split_vcfs_no_indels"
            # Store split vcf list in array
            vcf_arr=($(cat ${vcf_list}))
            # Create a small function that makes it easier to filter out indels in parallel
            function filter_indels() {
                local vcf_file=$1
                local out_dir=$2
                # Prefix for current vcf file
                prefix=$(basename ${vcf_file} .vcf)
                # Filter out indels
                vcftools --vcf ${vcf_file} --remove-indels --recode --recode-INFO-all --out "${out_dir}/Create_HC_Subset/Intermediates/temp_split_vcfs_no_indels/${prefix}_no_indels"
            }

            export -f filter_indels # export function
            # Filter out indels in parallel
            parallel filter_indels {} ${out} ::: ${vcf_arr[@]}
            # Concatenate no_indels vcfs
            cd "${out}/Create_HC_Subset/Intermediates/temp_split_vcfs_no_indels"
            find $(pwd -P) -name "*_no_indels.recode.vcf" | sort -V > "${out}/Create_HC_Subset/Intermediates/all_chr_no_indels_vcf_list.txt"
            # Check that this list is the same length as the original VCF list, if not exit out and ask
            # user to re-run this handler
            if [ "$(wc -l < ${out}/Create_HC_Subset/Intermediates/all_chr_no_indels_vcf_list.txt)" -eq "$(wc -l < ${vcf_list})" ]
            then
                # Line counts are equal
                echo "All split vcfs have indels removed, concatenating files."
                bcftools concat -f "${out}/Create_HC_Subset/Intermediates/all_chr_no_indels_vcf_list.txt" > "${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf"
                # Cleanup to save space
                rm -rf "${out}/Create_HC_Subset/Intermediates/temp_split_vcfs_no_indels" # Comment out this line if you need to debug this handler
            else
                # Line counts are NOT equal
                echo "Not all split vcfs and indels removed successfully, please manually check the files in ${out}/Create_HC_Subset/Intermediates/temp_split_vcfs_no_indels and re-run handler."
                exit 32 # If there are mismatches in number of files, exit
            fi
        fi
    else
        # Concatenated vcf is smaller than 500GB, we will proceed with the concatenated VCF file
        # Check if we have already filtered out indels, if so skip and proceed to next step
        if [ -n "$(ls -A ${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf 2>/dev/null)" ]; then
            echo "Already filtered out indels, proceed to next step."
        else
            # We need to filter out indels
            vcftools --vcf "${out}/Create_HC_Subset/${project}_concat_raw.vcf" --remove-indels --recode --recode-INFO-all --out "${out}/Create_HC_Subset/Intermediates/${project}_no_indels"
        fi
    fi

    # 3. Create a percentile table for the unfiltered SNPs
    source "${seqhand}/HelperScripts/percentiles.sh"
    # Check if we have already created the percentiles table for the unfiltered SNPs
    # The *_unfiltered_DP_per_sample.txt should be the last file that gets created
    if [ -n "$(ls -A ${out}/Create_HC_Subset/Percentile_Tables/${project}_unfiltered_DP_per_sample.txt 2>/dev/null)" ]
    then
        echo "Already generated percentiles table for unfiltered SNPs, proceed to next step."
    else
        percentiles "${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf" "${out}/Create_HC_Subset" "${project}" "unfiltered" "${seqhand}" "${temp_dir}"
    fi
    if [[ "$?" -ne 0 ]]; then
        echo "Error creating raw percentile tables, exiting..." >&2
        exit 32 # If something went wrong with the R script, exit
    fi

    # 4. Filter out sites that are low quality
    # We don't want to re-run this time consuming step if all of our cutoffs are the same as the previous run
    #   (in the case that we are re-running this handler due to exceeding walltime, etc.).
    # We only want to run this step if we have new cutoffs provided
    # If file exists, check header line for current VCF file's cutoffs used
    if [ -n "$(ls -A ${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf 2>/dev/null)" ]
    then
        echo "Filtered vcf file exists, checking if cutoffs have been changed..."
        # Filtered vcf file exists, check if cutoffs have been changed
        set -x
        # First check if we have a header line that starts with ##Create_HC_Subset_filter_cutoffs
        # This is used to check if our current cutoffs have been updated
        if grep -q "##Create_HC_Subset_filter_cutoffs" ${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf
        then
            # Expected header exists, check if cutoffs have been updated
            cutoffs_in_config=($(echo "Quality:"${qual_cutoff} "Het:"${max_het} "Max_bad:"${max_bad} "Genotype_Quality:"${gq_cutoff} "DP_per_sample:"${dp_per_sample_cutoff} | tr ' ' '\n'))
            cutoffs_in_vcf=($(grep "##Create_HC_Subset_filter_cutoffs" ${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf | cut -d'=' -f 2 | tr ',' '\n'))
        else
            # Expected header doesn't exist
            echo "VCF header line starting with ##Create_HC_Subset_filter_cutoffs doesn't exist."
            echo "If you have updated your cutoffs, the easiest solution is to delete your filtered VCF file (${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf) and re-run this handler."
            echo "If you know you have NOT updated your cutoffs, please read the following:"
            echo "Please check your VCF and make sure it has the header line starting with ##Create_HC_Subset_filter_cutoffs containing your previous run's cutoff values. If not, please add it."
            echo "Here is the format using the current cutoffs specified in your config:"
            echo "##Create_HC_Subset_filter_cutoffs=""Quality:"${qual_cutoff}",Het:"${max_het}",Max_bad:"${max_bad}",Genotype_Quality:"${gq_cutoff}",DP_per_sample:"${dp_per_sample_cutoff}
            exit 22 # Exit
        fi

        for i in "${!cutoffs_in_vcf[@]}"
        do
            if [ ${cutoffs_in_vcf[$i]} != ${cutoffs_in_config[$i]} ]; then
                echo "We have one or more updated cutoffs, re-run filtering with new cutoffs."
                echo "Cutoff in VCF" ${cutoffs_in_vcf[$i]}
                echo "Updated cutoff in Config" ${cutoffs_in_config[$i]}
                # Cutoffs have been changed, re-run filtering with updated cutoffs
                python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf"
                if [[ "$?" -ne 0 ]]; then
                    echo "Error with filter_sites.py, exiting..." >&2
                    exit 22 # If something went wrong with the python script, exit
                fi
                break # Break out of loop
            else
                echo ${cutoffs_in_vcf[$i]} "Cutoff in VCF"
                echo ${cutoffs_in_config[$i]} "Cutoff in Config"
            fi
        done
        set +x
    else
        # This is our first time filtering the vcf file
        python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Create_HC_Subset/Intermediates/${project}_no_indels.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf"
        if [[ "$?" -ne 0 ]]; then
            echo "Error with filter_sites.py, exiting..." >&2
            exit 22 # If something went wrong with the python script, exit
        fi
    fi

    # Get the number of sites left after filtering
    local num_sites=$(grep -v "#" "${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf" | wc -l)
    if [[ "${num_sites}" == 0 ]]; then
        echo "No sites left after filtering! Try using less stringent criteria. Exiting..." >&2
        exit 23 # If no sites left, error out with message
    fi

    # 5. Create a percentile table for the filtered SNPs
    percentiles "${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf" "${out}/Create_HC_Subset" "${project}" "filtered" "${seqhand}" "${temp_dir}"
    if [[ "$?" -ne 0 ]]; then
        echo "Error creating filtered percentile tables, exiting..." >&2
        exit 33 # If something went wrong with the R script, exit
    fi

    # 6. If barley, convert the parts positions into pseudomolecular positions. If not, then do nothing
    if [[ "${barley}" == true ]]
    then
        # Make sure the dictionary of positions in convert_parts_to_pseudomolecules.py match the current reference genome version
        echo "Converting barley parts positions to pseudomolecular positions..."
        python3 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" --vcf "${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf" > "${out}/Create_HC_Subset/Intermediates/${project}_filtered_pseudo.vcf"
        echo "Finished converting barley parts positions to pseudomolecular positions."
        local vcfoutput="${out}/Create_HC_Subset/Intermediates/${project}_filtered_pseudo.vcf"
    else
        local vcfoutput="${out}/Create_HC_Subset/Intermediates/${project}_filtered.vcf"
    fi

    # 7. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    echo "Removing sites that aren't polymorphic."
    vcftools --vcf "${vcfoutput}" --non-ref-ac 1 --recode --recode-INFO-all --out "${out}/Create_HC_Subset/${project}_high_confidence_subset"
    mv "${out}/Create_HC_Subset/${project}_high_confidence_subset.recode.vcf" "${out}/Create_HC_Subset/${project}_high_confidence_subset.vcf" # Rename the output file
    echo "Finished removing sites that aren't polymorphic."

    # 8. Remove intermediates to clear space
    # Since the user likely will run this handler multiple times, don't remove intermediate files
    # that don't need to be re-generated (i.e., filtered indels vcf) when handler is re-run.
    # Let the user decide what to remove when they are done.
    # rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

export -f Create_HC_Subset_GATK4

#   A function to call each filtering step
function Create_HC_Subset_GATK3() {
    local sample_list="$1" # What is our sample list?
    local out="$2"/Create_HC_Subset # Where are we storing our results?
    local bed="$3" # Where is the capture regions bed file?
    local barley="$4" # Is this barley?
    local project="$5" # What is the name of this project?
    local seqhand="$6" # Where is sequence_handling located?
    local qual_cutoff="$7" # What is the quality cutoff?
    local gq_cutoff="$8" # What is the genotyping quality cutoff?
    local dp_per_sample_cutoff="$9" # What is the DP per sample cutoff?
    local max_het="${10}" # What is the maximum number of heterozygous samples?
    local max_bad="${11}" # What is the maximum number of bad samples?
    local temp_dir="${12}" # Where can we store temporary files while running parallel?
    #   Make sure the out directories exist
    mkdir -p "${out}/Intermediates/Parts"
    mkdir -p "${out}/Percentile_Tables"
    #   1. Gzip all the chromosome part VCF files
    source "${seqhand}/HelperScripts/gzip_parts.sh"
    # Note: gzipping large files takes a long time and bcftools concat (based on a quick time comparison) runs faster
    # with uncompressed vcf files. So, for large VCF files, you may not want to gzip the files.
    parallel -v gzip_parts {} "${out}/Intermediates/Parts" :::: "${sample_list}" # Do the gzipping in parallel, preserve original files
    "${seqhand}/HelperScripts/sample_list_generator.sh" .vcf.gz "${out}/Intermediates/Parts" gzipped_parts.list # Make a list of the gzipped files for the next step
    #   2. Use bcftools to concatenate all the gzipped VCF files
    bcftools concat -f "${out}/Intermediates/Parts/gzipped_parts.list" > "${out}/Intermediates/${project}_concat.vcf"
    #   3. If exome capture, filter out SNPs outside the exome capture region. If not, then do nothing (This is necessary for GATK v3 since we would have called SNPs in all regions, but not always necessary for GATK 4 if we provided GATK4 with regions)
    if ! [[ "${bed}" == "NA" ]]
    then
        (set -x; vcfintersect -b "${bed}" "${out}/Intermediates/${project}_concat.vcf" > "${out}/Intermediates/${project}_capture_regions.vcf") # Perform the filtering
        local step3output="${out}/Intermediates/${project}_capture_regions.vcf"
    else
        local step3output="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   4. Filter out indels using vcftools
    vcftools --vcf "${step3output}" --remove-indels --recode --recode-INFO-all --out "${out}/Intermediates/${project}_no_indels" # Perform the filtering
    #   5. Create a percentile table for the unfiltered SNPs
    source "${seqhand}/HelperScripts/percentiles.sh"
    percentiles "${out}/Intermediates/${project}_no_indels.recode.vcf" "${out}" "${project}" "unfiltered" "${seqhand}" "${temp_dir}"
    if [[ "$?" -ne 0 ]]; then echo "Error creating raw percentile tables, exiting..." >&2; exit 32; fi # If something went wrong with the R script, exit
    #   6. Filter out sites that are low quality
    (set -x; python3 "${seqhand}/HelperScripts/filter_sites.py" "${out}/Intermediates/${project}_no_indels.recode.vcf" "${qual_cutoff}" "${max_het}" "${max_bad}" "${gq_cutoff}" "${dp_per_sample_cutoff}" > "${out}/Intermediates/${project}_filtered.vcf")
    if [[ "$?" -ne 0 ]]; then echo "Error with filter_sites.py, exiting..." >&2; exit 22; fi # If something went wrong with the python script, exit
    local num_sites=$(grep -v "#" "${out}/Intermediates/${project}_filtered.vcf" | wc -l) # Get the number of sites left after filtering
    if [[ "${num_sites}" == 0 ]]; then echo "No sites left after filtering! Try using less stringent criteria. Exiting..." >&2; exit 23; fi # If no sites left, error out with message
    #   7. Create a percentile table for the filtered SNPs
    percentiles "${out}/Intermediates/${project}_filtered.vcf" "${out}" "${project}" "filtered" "${seqhand}" "${temp_dir}"
    if [[ "$?" -ne 0 ]]; then echo "Error creating filtered percentile tables, exiting..." >&2; exit 33; fi # If something went wrong with the R script, exit
    #   8. If barley, convert the parts positions into pseudomolecular positions. If not, then do nothing
    if [[ "${barley}" == true ]]
    then
        (set -x; python3 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" "${out}/Intermediates/${project}_filtered.vcf" > "${out}/Intermediates/${project}_pseudo.vcf")
        local step8output="${out}/Intermediates/${project}_pseudo.vcf"
    else
        local step8output="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   9. Remove any sites that aren't polymorphic (minor allele count of 0). This is just a safety precaution
    vcftools --vcf "${step8output}" --non-ref-ac 1 --recode --recode-INFO-all --out "${out}/${project}_high_confidence_subset"
    mv "${out}/${project}_high_confidence_subset.recode.vcf" "${out}/${project}_high_confidence_subset.vcf" # Rename the output file
    #   10. Remove intermediates to clear space
    # rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Create_HC_Subset_GATK3
