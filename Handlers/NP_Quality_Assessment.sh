#!/bin/sh

set -o pipefail

#   What are the dependencies for Quality_Assessment?
declare -a NP_Quality_Assessment_Dependencies=(python NanoPlot)

#   A function to run quality assessment on fast5 files using Poretools
function NP_Quality_Assessment_Poretools() {
    local fast5_directory="$1" # Where is the directory of FAST5 files?
    local out="$2/Quality_Assessment" # Where are the results stored?
    local project="$3" # What do we call the results?
    # Make the output directory
    mkdir -p "${out}" 
    # Generate read statistics (total #, base pairs, averages)
    poretools stats "${fast5_directory}" > "${out}/read_statistics.txt"
    # Generate nucleotide distribution stats
    poretools nucdist "${fast5_directory}" > "${out}/nucleotide_statistics.txt"
    # Generate quality distributions stats
    poretools qualdist "${fast5_directory}" > "${out}/quality_statistics.txt"
    # Generate yield plots for reads over time
    poretools yield_plot --plot-type reads --saveas "${out}/read_yield_plot.pdf" "${fast5_directory}"
    # Generate yield plots for base pairs over time
    poretools yield_plot --plot-type basepairs --saveas "${out}/basepair_yield_plot.pdf" "${fast5_directory}"
    # Generate a histogram of read lengths
    poretools hist --min-length 50 --saveas "${out}/read_lengths_histogram.pdf" "${fast5_directory}"
    # Generate a boxplot of quality scores over position
    poretools qualpos --saveas "${out}/quality_versus_position.pdf" "${fast5_directory}"
    # Generate a diagram of flowcell occupancy
    poretools occupancy --saveas "${out}/flowcell_occupancy.pdf" "${fast5_directory}"
    # Merge all pdf files into one
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${out}/${project}_quality_summary.pdf" "${out}/read_lengths_histogram.pdf" "${out}/read_yield_plot.pdf" "${out}/basepair_yield_plot.pdf" "${out}/quality_versus_position.pdf" "${out}/flowcell_occupancy.pdf"
    # Delete the individual pdfs
    #rm "${out}/read_lengths_histogram.pdf" "${out}/read_yield_plot.pdf" "${out}/basepair_yield_plot.pdf" "${out}/quality_versus_position.pdf" "${out}/flowcell_occupancy.pdf"
}

export -f NP_Quality_Assessment_Poretools

#   A function to run quality assessment on fastq files using NanoPlot
function NP_Quality_Assessment() {
    local sample_list="$1" # Where is the list of samples?
    local out="$2/Quality_Assessment" # Where are the results stored?
    local project="$3" # What do we call the results?
    # Make an array of samples
    declare -a sample_array=($(grep -E ".fastq|.fastq.gz" "${sample_list}"))
    # Which sample are we working on now?
    local sample="${sample_array[${PBS_ARRAYID}]}"
    # Get the name of the sample without the path in front
    # If it is gzipped, this will strip off the .gz
    local base_name=$(basename ${sample} .gz)
    # Get the name of the sample without the file extension
    local sample_name=$(basename ${base_name} .fastq)
    # Make the output directory
    mkdir -p "${out}/${sample_name}"
    # Generate some plots
    NanoPlot --outdir "${out}/${sample_name}" --prefix "${sample_name}_" --format pdf --plots dot --plots kde --plots hex --plots pauvre --N50 --title "${project} ${sample_name} Quality Assessment" --fastq_rich "${sample}"
}

export -f NP_Quality_Assessment