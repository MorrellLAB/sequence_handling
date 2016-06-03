#!/bin/bash

#   This script is designed to generate a coverage histogram
#   from BAM files using BEDTools. Please have BEDTools installed
#   before using this script

set -e
set -o pipefail

#   What are the dependencies for Coverage_Mapping?
declare -a Coverage_Mapping_Dependencies=(bedtools Rscript)

#   A function to make our outdirectories
function makeOutDirectories() {
    local outPrefix="$1"
    mkdir -p "${outPrefix}"/CoverageMaps "${outPrefix}"/CoveragePlots
}

#   Export the function
export -f makeOutDirectories

#   A function to map the coverage
function mapCoverage() {
    local BAMFile="$1"
    local name="$2"
    local out="$3"
    local referenceAnnotation="$4"
    bedtools coverage -hist -abam "${BAMFile}" -b "${referenceAnnotation}" > "${out}"/"${name}".coverage.hist.txt
}

#   Export the function
export -f mapCoverage

#   A function to plot the coverage
function plotCoverage() {
    local sample="$1" # Figure out what this sample is
    local out="$2" # Where do we store our output files?
    local helperScripts="$3" # Where is our helper script located?
    local name="$(basename $sample .coverage.hist.txt)" # Get the name of the sample
    rm -f "${out}"/"${name}"_* # Remove any pipes
    mkfifo "${out}"/"${name}"_genome.PIPE # Make a pipe for the genome map
    mkfifo "${out}"/"${name}"_exon.PIPE # Make a pipe for the exon map
    mkfifo "${out}"/"${name}"_gene.PIPE # Make a pipe for the gene map
    grep 'all' "$sample" > "${out}"/"${name}"_genome.PIPE & # Make a map for the genome
    grep 'exon' "$sample" > "${out}"/"${name}"_exon.PIPE & # Make a map for exons
    grep 'gene' "$sample" > "${out}"/"${name}"_gene.PIPE & # Make a map for genes
    Rscript "${helperScripts}"/plot_cov.R "${out}"/"${name}"_genome.PIPE "${out}"/"${name}"_exon.PIPE "$sample" > "${out}"/"${name}"_gene.PIPE "${out}" "${name}"
}

#   Export the function
export -f plotCoverage

function Coverage_Mapping() {
    local sampleList="$1" # What is our list of samples?
    local outDirectory="$2"/Coverage_Mapping # Where do we store our results?
    local referenceAnnotation="$3" # What is our reference annotation file?
    local helperScripts="$4"/HelperScripts # Where do we keep our helper scripts?
    local -a sampleNames=($(parallel basename {} .bam :::: "${sampleList}")) # Create an array of names
    makeOutDirectories "${outDirectory}" # Make our output directories
    parallel --jobs 2 --xapply mapCoverage {1} {2} "${outDirectory}/CoverageMaps" "${referenceAnnotation}" :::: "${sampleList}" ::: "${sampleNames}" # Generate our coverage maps in text histogram form
    local -a histograms=($(find ${outDirectory}/CoverageMaps -name "*.coverage.hist.txt" | sort)) # Get a list of our coverage histograms
    parallel plotCoverage {} "${outDirectory}/CoveragePlots" "${helperScripts}" ::: "${histograms[@]}" # Generate our coverage plots
    find "${outDirectory}/CoveragePlots" -type p -exec rm {} \; # Remove any excess pipes
}

#   Export the function
export -f Coverage_Mapping
