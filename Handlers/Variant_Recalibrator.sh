#!/bin/bash

#   This script trains a model to recognize good SNPs
#   and then applies that model on a VCF file.

set -e
set -o pipefail

#   What are the dependencies for Variant_Recalibrator?
declare -a Variant_Recalibrator_Dependencies=(java parallel vcftools python3)

#   A function to parse the resources and determine appropriate settings
function ParseResources() {
    # Additional resources
    local res1="$1" # Where is the first VCF file to use as a training reference?
    local res2="$2" # Where is the second VCF file to use as a training reference?
    local res3="$3" # Where is the third VCF file to use as a training reference?
    local res4="$4" # Where is the fourth VCF file to use as a training reference?
    # Priors
    local p1="$5" # What is the prior for the first resource?
    local p2="$6" # What is the prior for the second resource?
    local p3="$7" # What is the prior for the third resource?
    local p4="$8" # What is the prior for the fourth resource?
    # Are the following resources known variants?
    local known1="$9"
    local known2="${10}"
    local known3="${11}"
    local known4="${12}"
    # Can the resources be used for training the model?
    local train1="${13}"
    local train2="${14}"
    local train3="${15}"
    local train4="${16}"
    # Are any of the resources truth sets?
    local truth1="${17}"
    local truth2="${18}"
    local truth3="${19}"
    local truth4="${20}"
    local gatk_version="${21}"
    # Create an array of resource arguments for GATK
    arguments=() # set as global variable, can be accessed after running function
    # The syntax for GATK 3 and GATK 4 differ slightly, check for version
    if [ ${gatk_version} == "4" ]; then
        # We are running GATK 4
        # If res1 exists, add it to the arguments
        if ! [[ "${res1}" == "NA" ]]; then
            arguments+=(--resource:first,known=${known1},training=${train1},truth=${truth1},prior=${p1} ${res1})
            # Make sure resource VCF file is indexed
            if [ -n "$(ls -A ${res1}.idx 2>/dev/null)" ]; then
                echo "Resource 1 VCF file is already indexed, add to arguments."
            else
                echo "Indexing Resource 1 VCF file..."
                gatk IndexFeatureFile -F ${res1}
                echo "Finished indexing Resource 1 vcf file"
            fi
        fi
        # If res2 exists, add it to the arguments
        if ! [[ "${res2}" == "NA" ]]; then
            arguments+=(--resource:second,known=${known2},training=${train2},truth=${truth2},prior=${p2} ${res2})
            # Make sure resource VCF file is indexed
            if [ -n "$(ls -A ${res2}.idx 2>/dev/null)" ]; then
                echo "Resource 2 VCF file is already indexed, add to arguments."
            else
                echo "Indexing Resource 2 VCF file..."
                gatk IndexFeatureFile -F ${res2}
                echo "Finished indexing Resource 2 vcf file"
            fi
        fi
        # If res3 exists, add it to the arguments
        if ! [[ "${res3}" == "NA" ]]; then
            arguments+=(--resource:third,known=${known3},training=${train3},truth=${truth3},prior=${p3} ${res3})
            # Make sure resource VCF file is indexed
            if [ -n "$(ls -A ${res3}.idx 2>/dev/null)" ]; then
                echo "Resource 3 VCF file is already indexed, add to arguments."
            else
                echo "Indexing Resource 3 VCF file..."
                gatk IndexFeatureFile -F ${res3}
                echo "Finished indexing Resource 3 vcf file"
            fi
        fi
        # If res4 exists, add it to the arguments
        if ! [[ "${res4}" == "NA" ]]; then
            arguments+=(--resource:fourth,known=${known4},training=${train4},truth=${truth4},prior=${p4} ${res4})
            # Make sure resource VCF file is indexed
            if [ -n "$(ls -A ${res4}.idx 2>/dev/null)" ]; then
                echo "Resource 4 VCF file is already indexed, add to arguments."
            else
                echo "Indexing Resource 4 VCF file..."
                gatk IndexFeatureFile -F ${res4}
                echo "Finished indexing Resource 4 vcf file"
            fi
        fi
    else
        # Assume we are running GATK 3
        # If res1 exists, add it to the arguments
        if ! [[ "${res1}" == "NA" ]]; then
            arguments+=(-resource:first,known=${known1},training=${train1},truth=${truth1},prior=${p1} ${res1})
        fi
        # If res2 exists, add it to the arguments
        if ! [[ "${res2}" == "NA" ]]; then
            arguments+=(-resource:second,known=${known2},training=${train2},truth=${truth2},prior=${p2} ${res2})
        fi
        # If res3 exists, add it to the arguments
        if ! [[ "${res3}" == "NA" ]]; then arguments+=(-resource:third,known=${known3},training=${train3},truth=${truth3},prior=${p3} ${res3})
        fi
        # If res4 exists, add it to the arguments
        if ! [[ "${res4}" == "NA" ]]; then
            arguments+=(-resource:fourth,known=${known4},training=${train4},truth=${truth4},prior=${p4} ${res4})
        fi
    fi
}

#   Export the function
export -f ParseResources

#   A function to run GATK Variant Recalibrator
function Variant_Recalibrator_GATK4() {
    local vcf_raw_concat="$1"
    local vcf_list="$2" # What is our VCF list?
    local out="$3" # Where are we storing our results?
    local gatk="$4" # Where is the GATK jar?
    local reference="$5" # Where is the reference sequence?
    local memory="$6" # How much memory can java use?
    local project="$7" # What is the name of the project?
    local seqhand="$8" # Where is sequence_handling located?
    # Resources used for recalibration
    local hc_subset="$9" # Where are the high-confidence variants?
    local res1="${10}" # Where is the first VCF file to use as a training reference?
    local res2="${11}" # Where is the second VCF file to use as a training reference?
    local res3="${12}" # Where is the third VCF file to use as a training reference?
    local res4="${13}" # Where is the fourth VCF file to use as a training reference?
    # Priors
    local hc_prior="${14}" # What is the prior for the high-confidence variants?
    local p1="${15}" # What is the prior for the first resource?
    local p2="${16}" # What is the prior for the second resource?
    local p3="${17}" # What is the prior for the third resource?
    local p4="${18}" # What is the prior for the fourth resource?
    # Are the following resources known variants?
    local hc_known="${19}"
    local known1="${20}"
    local known2="${21}"
    local known3="${22}"
    local known4="${23}"
    # Can the resources be used for training the model?
    local hc_train="${24}"
    local train1="${25}"
    local train2="${26}"
    local train3="${27}"
    local train4="${28}"
    # Are any of the resources truth sets?
    local hc_truth="${29}"
    local truth1="${30}"
    local truth2="${31}"
    local truth3="${32}"
    local truth4="${33}"
    # Are we working with barley?
    local barley="${34}"
    local gatk_version="${35}"
    local ts_filter_level="${36}"
    #   Check if diretory exists, if not make it
    mkdir -p ${out}/Variant_Recalibrator \
             ${out}/Variant_Recalibrator/Intermediates

    set -x # for testing, remove after done
    #   Check if we need to concatenate our split VCF files into a single raw VCF file
    #   If we don't have a split VCF file, assume that we have a concatentated vcf file
    if [ ${vcf_list} == "NA" ]; then
        echo "Using the concatenated raw VCF file provided: ${vcf_raw_concat}"
        to_recal_vcf=${vcf_raw_concat}
    else
        # Assume we have split VCF files that we need to concatenate
        echo "We need to concatenate our split VCF files, concatenating..."
        bcftools concat -f ${vcf_list} > "${out}/Variant_Recalibrator/${project}_concat_raw.vcf"
        to_recal_vcf="${out}/Variant_Recalibrator/${project}_concat_raw.vcf"
        echo "Done concatenating."
    fi

    # Check if high confidence subset vcf file is indexed, if not index it
    echo ${hc_subset} # for testing
    if [ -n "$(ls -A ${hc_subset}.idx 2>/dev/null)" ]; then
        echo "High confidence subset VCF file is already indexed, proceeding to recalibration step..."
    else
        echo "Indexing high confidence subset VCF file..."
        gatk IndexFeatureFile -F ${hc_subset}
        echo "Finished indexing high confidence subset VCF file."
    fi

    #   Get the GATK settings for the additional resources
    #   Function returns global variable: arguments
    ParseResources ${res1} ${res2} ${res3} ${res4} ${p1} ${p2} ${p3} ${p4} ${known1} ${known2} ${known3} ${known4} ${train1} ${train2} ${train3} ${train4} ${truth1} ${truth2} ${truth3} ${truth4} ${gatk_version}
    local settings=$(echo -n ${arguments[@]}) # Strip trailing newline
    #   Build the recalibration model for SNPs
    #   For GATK 4, SNPs and indels must be recalibrated in separate runs, but
    #       it is not necessary to separate them into different files
    #   Recalibrate indels
    #   Note: Removed -an MQ because: "For filtering indels, most annotations related to mapping quality have been removed since there is a conflation with the length of an indel in a read and the degradation in mapping quality that is assigned to the read by the aligner. This covariation is not necessarily indicative of being an error in the same way that it is for SNPs."
    echo "Starting indel recalibration..."
    gatk --java-options "-Xmx${memory}" VariantRecalibrator \
        -R "${reference}" \
        -V "${to_recal_vcf}" \
        -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
        -mode INDEL \
        -O "${out}/Variant_Recalibrator/Intermediates/${project}_recal_indels.txt" \
        --resource:highconfidence,known=${hc_known},training=${hc_train},truth=${hc_truth},prior=${hc_prior} ${hc_subset} \
        ${settings} \
        --tranches-file ${out}/Variant_Recalibrator/Intermediates/${project}_tranches_indels.txt \
        --rscript-file ${out}/Variant_Recalibrator/Intermediates/${project}_indels.plots.R
    echo "Finished indel recalibration."

    #   Recalibrate SNPs
    echo "Starting snp recalibration..."
    gatk --java-options "-Xmx${memory}" VariantRecalibrator \
        -R "${reference}" \
        -V "${to_recal_vcf}" \
        -an FS -an ReadPosRankSum -an MQ -an MQRankSum -an QD -an SOR -an DP \
        -mode SNP \
        -O "${out}/Variant_Recalibrator/Intermediates/${project}_recal_snps.txt" \
        --resource:highconfidence,known=${hc_known},training=${hc_train},truth=${hc_truth},prior=${hc_prior} ${hc_subset} \
        ${settings} \
        --tranches-file ${out}/Variant_Recalibrator/Intermediates/${project}_tranches_snps.txt \
        --rscript-file ${out}/Variant_Recalibrator/Intermediates/${project}_snps.plots.R
    echo "Finished snp recalibration."

    #   Add Rscripts to environment PATH
    export PATH=${out}/Variant_Recalibrator/Intermediates/${project}_indels.plots.R:${PATH}
    export PATH=${out}/Variant_Recalibrator/Intermediates/${project}_snps.plots.R:${PATH}

    #   Now, successively apply the indel and SNP recalibrations to the full callset to
    #       produce a final filtered callset
    #   We use --ts_filter 99.9 to take 99.9% of true positives from the model, which is recommended in the GATK docs
    #   Filter indels on VQSLOD using ApplyVQSR, outputs an indel filtered callset
    echo "Filtering indels on VQSLOD using ApplyVQSR..."
    gatk --java-options "-Xmx${memory}" ApplyVQSR \
        -R "${reference}" \
        -V "${to_recal_vcf}" \
        -mode INDEL \
        --truth-sensitivity-filter-level ${ts_filter_level} \
        --recal-file "${out}/Variant_Recalibrator/Intermediates/${project}_recal_indels.txt" \
        --tranches-file "${out}/Variant_Recalibrator/Intermediates/${project}_tranches_indels.txt" \
        --create-output-variant-index true \
        -O "${out}/Variant_Recalibrator/Intermediates/${project}_indel.recalibrated.vcf"
    echo "Finished filtering indels on VQSLOD. This outputs an indel filtered callset: ${out}/Variant_Recalibrator/Intermediates/${project}_indel.recalibrated.vcf"

    #   Now, filter SNP variants
    echo "Filtering snps on VQSLOD using ApplyVQSR..."
    gatk --java-options "-Xmx${memory}" ApplyVQSR \
        -R "${reference}" \
        -V "${out}/Variant_Recalibrator/Intermediates/${project}_indel.recalibrated.vcf" \
        -mode SNP \
        --truth-sensitivity-filter-level ${ts_filter_level} \
        --recal-file "${out}/Variant_Recalibrator/Intermediates/${project}_recal_snps.txt" \
        --tranches-file "${out}/Variant_Recalibrator/Intermediates/${project}_tranches_snps.txt" \
        --create-output-variant-index true \
        -O "${out}/Variant_Recalibrator/${project}_snps.recalibrated.vcf"
    echo "Finished filtering snps on VQSLOD. This outputs a SNP filtered callset, which is also the final filtered callset since we gave it the indel filtered callset: ${out}/Variant_Recalibrator/${project}_snps.recalibrated.vcf"
    echo "Note: filtered means that variants failing the requested tranche cutoff are marked as filtered in the output VCF, these are NOT discarded yet."
    set +x # for testing, remove after done
}

export -f Variant_Recalibrator_GATK4

#   A function to run GATK Variant Recalibrator
function Variant_Recalibrator_GATK3() {
    local vcf_raw_concat="$1"
    local vcf_list="$2" # What is our VCF list?
    local out="$3" # Where are we storing our results?
    local gatk="$4" # Where is the GATK jar?
    local reference="$5" # Where is the reference sequence?
    local memory="$6" # How much memory can java use?
    local project="$7" # What is the name of the project?
    local seqhand="$8" # Where is sequence_handling located?
    # Resources used for recalibration
    local hc_subset="$9" # Where are the high-confidence variants?
    local res1="${10}" # Where is the first VCF file to use as a training reference?
    local res2="${11}" # Where is the second VCF file to use as a training reference?
    local res3="${12}" # Where is the third VCF file to use as a training reference?
    local res4="${13}" # Where is the fourth VCF file to use as a training reference?
    # Priors
    local hc_prior="${14}" # What is the prior for the high-confidence variants?
    local p1="${15}" # What is the prior for the first resource?
    local p2="${16}" # What is the prior for the second resource?
    local p3="${17}" # What is the prior for the third resource?
    local p4="${18}" # What is the prior for the fourth resource?
    # Are the following resources known variants?
    local hc_known="${19}"
    local known1="${20}"
    local known2="${21}"
    local known3="${22}"
    local known4="${23}"
    # Can the resources be used for training the model?
    local hc_train="${24}"
    local train1="${25}"
    local train2="${26}"
    local train3="${27}"
    local train4="${28}"
    # Are any of the resources truth sets?
    local hc_truth="${29}"
    local truth1="${30}"
    local truth2="${31}"
    local truth3="${32}"
    local truth4="${33}"
    # Are we working with barley?
    local barley="${34}"
    local gatk_version="${35}"
    local ts_filter_level="${36}"
    mkdir -p "${out}/Intermediates/Parts" # Make sure the out directory exists
    #   Gzip all the chromosome part VCF files, because they must be gzipped to combine
    source "${seqhand}/HelperScripts/gzip_parts.sh"
    parallel -v gzip_parts {} "${out}/Intermediates/Parts" :::: "${vcf_list}" # Do the gzipping in parallel, preserve original files
    "${seqhand}/HelperScripts/sample_list_generator.sh" .vcf.gz "${out}/Intermediates/Parts" gzipped_parts.list # Make a list of the gzipped files for the next step
    #   Use vcftools to concatenate all the gzipped VCF files
    vcf-concat -f "${out}/Intermediates/Parts/gzipped_parts.list" > "${out}/Intermediates/${project}_concat.vcf"
    #   Change the concatenated VCF to pseudomolecular positions if barley. If not barley, do nothing.
    if [[ "${barley}" == true ]]
    then
        python3 "${seqhand}/HelperScripts/convert_parts_to_pseudomolecules.py" "${out}/Intermediates/${project}_concat.vcf" > "${out}/Intermediates/${project}_pseudo.vcf"
        local to_recal="${out}/Intermediates/${project}_pseudo.vcf"
    else
        local to_recal="${out}/Intermediates/${project}_concat.vcf"
    fi
    #   Get the GATK settings for the resources
    local settings=$(ParseResources ${res1} ${res2} ${res3} ${res4} ${p1} ${p2} ${p3} ${p4})
    #   Build the recalibration model for SNPs
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
        -T VariantRecalibrator \
        -an MQ \
        -an MQRankSum \
        -an DP \
        -an ReadPosRankSum \
        -mode SNP \
        -input "${to_recal}" \
        -R "${reference}" \
        -recalFile "${out}/Intermediates/${project}_recal_file.txt" \
        -tranchesFile "${out}/Intermediates/${project}_tranches_file.txt" \
        -resource:highconfidence,known=false,training=true,truth=false,prior="${hc_prior}" "${hc_subset}" \
        ${settings})
    #   Now, actually apply it
    #   We use --ts_filter 99.9 to take 99.9% of true positives from the model, which is recommended in the GATK docs
    (set -x; java -Xmx"${memory}" -jar "${gatk}" \
        -T ApplyRecalibration \
        -R "${reference}" \
        -input "${to_recal}" \
        -mode SNP \
        --ts_filter_level ${ts_filter_level} \
        -recalFile "${out}/Intermediates/${project}_recal_file.txt" \
        -tranchesFile "${out}/Intermediates/${project}_tranches_file.txt" \
        -o "${out}/${project}_recalibrated.vcf")
    #   Remove the intermediates
    rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Recalibrator_GATK3
