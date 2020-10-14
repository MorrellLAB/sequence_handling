#!/bin/bash

# Generate variant annotation files for graphing distributions
function graph_annotations() {
	local raw_vcf="$1"
    local out="$2"
    local project="$3"
    local ref_gen="$4"
    local seqhand="$5"

# Create .bed file with random genome intervals (to speed up computational time)
echo "Creating intervals file to subset genome randomly"
awk -v OFS='\t' {'print $1,$2'} "${ref_gen}.fai" > "${out}/Intermediates/GenomeFile.txt"

bedtools random -l 100 -n 1000000 -seed 65 \
-g "${out}/Intermediates/GenomeFile.txt" | \
sort -k 1,1 -k2,2n > "${out}/Intermediates/Genome_Random_Intervals.bed"

# Obtain annotation information for raw variants in intervals:
echo "Creating a table of annotation scores for raw variants in intervals"

gatk VariantsToTable \
     -V ${raw_vcf} \
     -L "${out}/Intermediates/Genome_Random_Intervals.bed" \
     -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O "${out}/Intermediates/RawVariants.table"

# Obtain annotation information for high-confidence variants in intervals:
echo "Creating a table of annotation scores for the variants that passed filtering in intervals"

hc_subset="${out}/${project}_high_confidence_subset.vcf"

# high-confidence variants need to be indexed for VariantsToTable function:
gatk IndexFeatureFile \
     -F ${hc_subset}

gatk VariantsToTable \
     -V ${hc_subset} \
     -L "${out}/Intermediates/Genome_Random_Intervals.bed" \
     -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O "${out}/Intermediates/HCVariants.table"

# Make graphs of annotation distributions
echo "Calculating annotation distributions for variant sets"

Rscript "${seqhand}/HelperScripts/graph_annotations.R" "${out}" "Intermediates/RawVariants.table" "Intermediates/HCVariants.table"

#rm "${out}/Intermediates/GenomeFile.txt"
#rm "${out}/Intermediates/Genome_Random_intervals.bed"
#rm "${out}/Intermediates/RawVariants.table"
#rm "${out}/Intermediates/HCVariants.table"
}

export -f graph_annotations