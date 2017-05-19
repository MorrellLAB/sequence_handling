#!/bin/bash

#   This script is designed to generate a coverage histogram
#   from BAM files using BEDTools. Please have BEDTools installed
#   before using this script

set -o pipefail

#   What are the dependencies for Coverage_Mapping?
declare -a Coverage_Mapping_Dependencies=(bedtools Rscript parallel)

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
    bedtools coverage -hist -abam "${BAMFile}" -b "${referenceAnnotation}" > "${out}/${name}.coverage.hist.txt"
}

#   Export the function
export -f mapCoverage

#   A function to plot the coverage
function plotCoverage() {
    local sample="$1" # Figure out what this sample is
    local out="$2" # Where do we store our output files?
    local sequenceHandling="$3" # Where is sequence_handling?
    local helperScripts="${sequenceHandling}"/HelperScripts
    local name="$(basename ${sample} .coverage.hist.txt)" # Get the name of the sample
    Rscript "${helperScripts}"/plot_cov.R "${sample}" "${out}" "${name}" "${sequenceHandling}"
}

#   Export the function
export -f plotCoverage

function sampleCoverage() {
    local sample="$1" # What sample are we working with?
    local outRoot="$2" # Root of our output directory
    #   Create an output name and directory
    local sampleExtension="$(echo ${sample} | rev | cut -f 1 -d '.' | rev)"
    local sampleName="$(basename ${sample} ${sampleExtension})"
    local outDirectory="${outRoot}/${sampleName}"
    (set -x; mkdir -p "${outdirectories}")
    #   What are the chromosomes that we can calculate coverage for?
    local chromosomes=($(cut -f 1 "${sample}" | sort -u))
}

#   A function to summarize the coverage
function summarizeCoverage() {
    return -8
}

function Coverage_Mapping() {
    local sampleList="$1" # What is our list of samples?
    local outDirectory="$2"/Coverage_Mapping # Where do we store our results?
    local referenceAnnotation="$3" # What is our reference annotation file?
    local sequenceHandling="$4" # Where is sequence_handling
    echo "Collecting sample names..." >&2
    local -a sampleNames=($(xargs --arg-file="${sampleList}" -I @ --max-args=1 basename @ .bam)) # Create an array of names
    makeOutDirectories "${outDirectory}" # Make our output directories
    echo "Mapping coverage..." >&2
    parallel --jobs 2 --xapply mapCoverage {1} {2} "${outDirectory}/CoverageMaps" "${referenceAnnotation}" :::: "${sampleList}" ::: "${sampleNames[@]}" # Generate our coverage maps in text histogram form
    echo "Finding coverage histograms..." >&2
    local -a histograms=($(find ${outDirectory}/CoverageMaps -name "*.coverage.hist.txt" | sort)) # Get a list of our coverage histograms
    echo "Plotting coverage..." >&2
    Rscript "${sequenceHandling}"/HelperScripts/install_dependencies.R "${sequenceHandling}" 'Hmisc'
    parallel plotCoverage {} "${outDirectory}/CoveragePlots" "${sequenceHandling}" ::: "${histograms[@]}" # Generate our coverage plots
}

#   Export the function
export -f Coverage_Mapping
