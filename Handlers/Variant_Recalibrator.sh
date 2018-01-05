#!/bin/bash

#   This script trains a model to recognize good SNPs
#   and then applies that model on a VCF file.

set -o pipefail

#   What are the dependencies for Variant_Recalibrator?
declare -a Variant_Recalibrator_Dependencies=(java parallel vcftools python3)

#   A function to parse the resources and determine appropriate settings
function ParseResources() {
    local res1="$1" # Where is the first VCF file to use as a training reference?
    local res2="$2" # Where is the second VCF file to use as a training reference?
    local res3="$3" # Where is the third VCF file to use as a training reference?
    local res4="$4" # Where is the fourth VCF file to use as a training reference?
    local p1="$5" # What is the prior for the first resource?
    local p2="$6" # What is the prior for the second resource?
    local p3="$7" # What is the prior for the third resource?
    local p4="$8" # What is the prior for the fourth resource?
    arguments=() # Create an array of resource arguments for GATK
    if ! [[ "${res1}" == "NA" ]]; then arguments+=(-resource:first,known=true,training=true,truth=true,prior=${p1} ${res1}); fi # If res1 exists, add it to the arguments
    if ! [[ "${res2}" == "NA" ]]; then arguments+=(-resource:second,known=true,training=true,truth=true,prior=${p2} ${res2}); fi # If res2 exists, add it to the arguments
    if ! [[ "${res3}" == "NA" ]]; then arguments+=(-resource:third,known=true,training=true,truth=true,prior=${p3} ${res3}); fi # If res3 exists, add it to the arguments
    if ! [[ "${res4}" == "NA" ]]; then arguments+=(-resource:fourth,known=true,training=true,truth=true,prior=${p4} ${res4}); fi # If res4 exists, add it to the arguments
    echo -n ${arguments[@]} # Return our settings
}

#   Export the function
export -f ParseResources

#   A function to run GATK Variant Recalibrator
function Variant_Recalibrator() {
    local vcf_list="$1" # What is our VCF list?
    local out="$2"/Variant_Recalibrator # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local hc_subset="$5" # Where are the high-confidence variants?
    local memory="$6" # How much memory can java use?
    local project="$7" # What is the name of the project?
    local res1="$8" # Where is the first VCF file to use as a training reference?
    local res2="$9" # Where is the second VCF file to use as a training reference?
    local res3="${10}" # Where is the third VCF file to use as a training reference?
    local res4="${11}" # Where is the fourth VCF file to use as a training reference?
    local seqhand="${12}" # Where is sequence_handling located?
    local p1="${13}" # What is the prior for the first resource?
    local p2="${14}" # What is the prior for the second resource?
    local p3="${15}" # What is the prior for the third resource?
    local p4="${16}" # What is the prior for the fourth resource?
    local hc_prior="${17}" # What is the prior for the high-confidence variants?
    local barley="${18}" # Is this barley?
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
        --ts_filter_level 99.9 \
        -recalFile "${out}/Intermediates/${project}_recal_file.txt" \
        -tranchesFile "${out}/Intermediates/${project}_tranches_file.txt" \
        -o "${out}/${project}_recalibrated.vcf")
    #   Remove the intermediates
    #rm -Rf "${out}/Intermediates" # Comment out this line if you need to debug this handler
}

#   Export the function
export -f Variant_Recalibrator