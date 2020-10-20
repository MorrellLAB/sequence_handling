#!/bin/bash

# Generate variant annotation files for graphing distributions
function graph_annotations() {
    local raw_vcf="$1"
    local out="$2"
    local project="$3"
    local ref_gen="$4"
    local seqhand="$5"
    local gen_num="$6"
    local gen_len="$7"

# Assign values if the gen_num and gen_len variables are unset (Bedtools default values):
if [ -z "${gen_num}" ]; then
    echo "Number of Genomic Regions to subset not specified. Will be set to 1000000"
    gen_num=1000000
fi
if [ -z "${gen_len}" ]; then
    echo "Length of Genomic Regions to subset not specified. Will be set to 100"
    gen_len=100
fi
# Create .bed file with random genome intervals (to speed up computational time)
echo "Creating intervals file to subset genome randomly at ${gen_num} x ${gen_len}bp regions"
awk -v OFS='\t' {'print $1,$2'} "${ref_gen}.fai" > "${out}/Intermediates/GenomeFile.txt"

bedtools random -l ${gen_len} -n ${gen_num} -seed 65 \
-g "${out}/Intermediates/GenomeFile.txt" | \
sort -k 1,1 -k2,2n > "${out}/Intermediates/Genome_Random_Intervals.bed"

# Obtain annotation information for raw variants in intervals:
echo "Creating a table of annotation scores for raw variants in intervals"

gatk VariantsToTable \
     -V ${raw_vcf} \
     -L "${out}/Intermediates/Genome_Random_Intervals.bed" \
     -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O "${out}/Intermediates/RawVariants.table"

hc_subset="${out}/${project}_high_confidence_subset.vcf"

# high-confidence variants need to be indexed for VariantsToTable function
# check if it has already been indexed
if [ -f ${hc_subset}.idx ]; then
    echo "high-confidence VCF file is already indexed."
else
    echo "Indexing high-confidence VCF file..."
    gatk IndexFeatureFile -F ${hc_subset}
    echo "Finished indexing high-confidence VCF file"
fi

# Obtain annotation information for high-confidence variants in intervals:
echo "Creating a table of annotation scores for the variants that passed filtering in intervals"

gatk VariantsToTable \
     -V ${hc_subset} \
     -L "${out}/Intermediates/Genome_Random_Intervals.bed" \
     -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O "${out}/Intermediates/HCVariants.table"

# Make graphs of annotation distributions
echo "Calculating annotation distributions for variant sets"

num_name=$(( ${gen_num}/1000000 ))
suffix="${num_name}Mx${gen_len}bp"

Rscript "${seqhand}/HelperScripts/graph_annotations.R" \
    "${out}" \
    "${out}/Intermediates/RawVariants.table" \
    "${out}/Intermediates/HCVariants.table" \
    "${suffix}"

#rm "${out}/Intermediates/GenomeFile.txt"
#rm "${out}/Intermediates/Genome_Random_intervals.bed"
#rm "${out}/Intermediates/RawVariants.table"
#rm "${out}/Intermediates/HCVariants.table"
}

export -f graph_annotations