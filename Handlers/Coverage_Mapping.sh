#!/bin/bash

#   This script is designed to generate a coverage histogram
#   from BAM files using BEDTools. Please have BEDTools installed
#   before using this script

set -e
set -o pipefail

#   What are the dependencies for Coverage_Mapping?
declare -a Coverage_Mapping_Dependencies=(bedtools Rscript parallel)

#   A function to find BAM files
function findBAM() {
    local bamDirectory="$1" # Where are the BAM files stored?
    local project="$2" # What is our project called?
    local BAMList="${bamDirectory}"/"${project}"_BAM.txt # Create a name for our sample list
    find -L "${bamDirectory}" -name "*.bam" | sort > "${BAMList}" # Create our list
    echo "${BAMList}" # Return the name of our sample list
}

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
    local name="$(basename ${sample} .coverage.hist.txt)" # Get the name of the sample
    Rscript "${helperScripts}"/plot_cov.R "${sample}" "${out}" "${name}"
}

#   Export the function
export -f plotCoverage

function Coverage_Mapping() {
    local sampleList="$1" # What is our list of samples?
    local outDirectory="$2"/Coverage_Mapping # Where do we store our results?
    local referenceAnnotation="$3" # What is our reference annotation file?
    local helperScripts="$4"/HelperScripts # Where do we keep our helper scripts?
    echo "Collecting sample names..." >&2
    local -a sampleNames=($(xargs --arg-file="${sampleList}" -I @ --max-args=1 basename @ .bam)) # Create an array of names
    makeOutDirectories "${outDirectory}" # Make our output directories
    echo "Mapping coverage..." >&2
    parallel --jobs 2 --xapply mapCoverage {1} {2} "${outDirectory}/CoverageMaps" "${referenceAnnotation}" :::: "${sampleList}" ::: "${sampleNames[@]}" # Generate our coverage maps in text histogram form
    echo "Finding coverage histograms..." >&2
    local -a histograms=($(find ${outDirectory}/CoverageMaps -name "*.coverage.hist.txt" | sort)) # Get a list of our coverage histograms
    echo "Plotting coverage..." >&2
    parallel plotCoverage {} "${outDirectory}/CoveragePlots" "${helperScripts}" ::: "${histograms[@]}" # Generate our coverage plots
}

#   Export the function
export -f Coverage_Mapping
