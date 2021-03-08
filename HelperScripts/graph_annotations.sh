#!/bin/bash

# Generate variant annotation files for graphing distributions
function index_vcf() {
    local vcf="$1"
    if [ -f ${vcf}.tbi ]; then
        echo "VCF is already indexed, proceed to obtain annotation info."
    else
        echo "Indexing VCF file..."
        gatk IndexFeatureFile -F ${vcf}
        echo "Finished indexing VCF file."
    fi
}

export -f index_vcf

function graph_annotations() {
    local vcf1="$1"
    local vcf2="$2"
    local out="$3"
    local project="$4"
    local ref_gen="$5"
    local seqhand="$6"
    local gen_num="$7"
    local gen_len="$8"
    local vcf1_out_prefix="$9"
    local vcf2_out_prefix="${10}"
    local plot_type="${11}" # Valid options: "single" or "pair"
    # Check if out subdirectory exists, if not make it
    mkdir -p ${out}/Plots
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
    # First, need genome file
    awk -v OFS='\t' {'print $1,$2'} "${ref_gen}.fai" > "${out}/Intermediates/GenomeFile.txt"

    # Assign variables for filenames
    num_name=$(echo "scale=2; $gen_num/1000000" | bc)
    if [[ "$(echo $num_name | cut -c 2-)" == ".00" ]]; then
        num_name=$(echo "scale=0; $gen_num/1000000" | bc)
    fi
    suffix="${num_name}Mx${gen_len}bp"

    echo "Creating intervals file to subset genome randomly at ${gen_num} ${gen_len}bp regions"
    bedtools random -l ${gen_len} -n ${gen_num} -seed 65 \
    -g "${out}/Intermediates/GenomeFile.txt" | \
    sort -k 1,1 -k2,2n > "${out}/Intermediates/Genome_Random_Intervals_${suffix}.bed"

    # Check if we are plotting one vcf or two vcf files (unfiltered and filtered)
    # This is setup this way so this script can be called in multiple handlers
    #   (i.e., Create_HC_Subset, Variant_Recalibrator, and Variant_Filtering)
    if [ ${plot_type} == "single" ];then
        # Plotting one VCF
        echo "Plotting one VCF file..."
        # Check if VCF file is index, if not index file
        index_vcf ${vcf1}
        
        # Obtain annotation information for raw variants in intervals:
        echo "Creating a table of annotation scores for first VCF in intervals"
        gatk VariantsToTable \
            -V ${vcf1} \
            -L "${out}/Intermediates/Genome_Random_Intervals_${suffix}.bed" \
            -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
            -O "${out}/Intermediates/${vcf1_out_prefix}_Variants_in${suffix}.table"
        
        # Make graph of annotation distributions
        # Single VCF
        echo "Calculating annotation distributions for variant sets"
        vcf2_out_prefix="NA" # Since we are only plotting a single plot
        Rscript "${seqhand}/HelperScripts/graph_annotations.R" \
            "${out}" \
            "${out}/Intermediates/${vcf1_out_prefix}_Variants_in${suffix}.table" \
            "NA" \
            "${suffix}" \
            "${plot_type}" \
            "${vcf1_out_prefix}" \
            "${vcf2_out_prefix}"
    elif [ ${plot_type} == "pair" ]; then
        # Plotting two VCF files
        echo "Plotting pair of VCF file..."
        # Check if VCF files are indexed, if not index them
        index_vcf ${vcf1}
        index_vcf ${vcf2}

        # Obtain annotation information for raw variants in intervals:
        echo "Creating a table of annotation scores for first VCF in intervals"
        gatk VariantsToTable \
            -V ${vcf1} \
            -L "${out}/Intermediates/Genome_Random_Intervals_${suffix}.bed" \
            -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
            -O "${out}/Intermediates/${vcf1_out_prefix}_Variants_in${suffix}.table"
        echo "Creating a table of annotation scores for second VCF in intervals"
        gatk VariantsToTable \
            -V ${vcf2} \
            -L "${out}/Intermediates/Genome_Random_Intervals_${suffix}.bed" \
            -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
            -O "${out}/Intermediates/${vcf2_out_prefix}_Variants_in${suffix}.table"

        # Make graph of annotation distributions
        # Two VCF files
        echo "Calculating annotation distributions for variant sets"
        Rscript "${seqhand}/HelperScripts/graph_annotations.R" \
            "${out}" \
            "${out}/Intermediates/${vcf1_out_prefix}_Variants_in${suffix}.table" \
            "${out}/Intermediates/${vcf2_out_prefix}_Variants_in${suffix}.table" \
            "${suffix}" \
            "${plot_type}" \
            "${vcf1_out_prefix}" \
            "${vcf2_out_prefix}"
    fi

    #rm "${out}/Intermediates/GenomeFile.txt"
    #rm "${out}/Intermediates/Genome_Random_Intervals_${suffix}.bed"
    #rm "${out}/Intermediates/RawVariants_in${suffix}.table"
    #rm "${out}/Intermediates/${vcf2_out_prefix}_Variants_in${suffix}.table"
}

export -f graph_annotations
