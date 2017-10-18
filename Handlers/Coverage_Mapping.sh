#!/bin/bash

set -o pipefail

#   What are the dependencies for Coverage_Mapping?
declare -a Coverage_Mapping_Dependencies=(bedtools parallel)

#   Makes the outdirectories
function makeOutDirectories() {
    local outPrefix="$1"
    mkdir -p "${outPrefix}"/Histograms
    mkdir -p "${outPrefix}"/Plots
}

export -f makeOutDirectories

#   A function to plot the coverage - DEPRICATED, might be added back later
# function plotCoverage() {
#     local sample="$1" # Figure out what this sample is
#     local out="$2" # Where do we store our output files?
#     local sequenceHandling="$3" # Where is sequence_handling?
#     local helperScripts="${sequenceHandling}"/HelperScripts
#     local name="$(basename ${sample} .coverage.hist.txt)" # Get the name of the sample
#     Rscript "${helperScripts}"/plot_cov.R "${sample}" "${out}" "${name}" "${sequenceHandling}"
# }

# export -f plotCoverage

#   This is to calculate histograms and summary statistics for exome capture data
function EC_Coverage() {
    local bam_file="$1"
    local bam_dir=$(dirname "${bam_file}")
    local region_file="$2"
    local out_dir="$3"
    local project="$4"
    #   Get the sample name without the .bam
    local sampleName=$(basename "${bam_file}" .bam)
    #   Generate coverage histograms as text files
    bedtools coverage -hist -abam "${bam_file}" -b "${region_file}" > ${out_dir}/Histograms/${sampleName}.hist
    #   Begin calculating statistics per bp
    #   The minimum is the coverage on the first line of the "all" fields since they're already sorted
    local min=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n 1 | awk -F "\t" '{print $2}')
    #   The maximum is the coverage on the last line of the "all" fields
    local max=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | tail -n 1 | awk -F "\t" '{print $2}')
    #   The mean is the sum of (each coverage * the percent of the genome at that coverage)
    local mean=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | awk '{ sum += $2*$5 } END { print sum }')
    #   The mode is the coverage that has the highest percent of the genome at that coverge (excludes zero coverage)
    local mode=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | tail -n +2 | sort -grk5,5 | head -1 | awk -F "\t" '{print $2}')
    #   The quantiles are a bit tricky...
    #   row_count will count how many rows down the "all" fields we are
    local row_count="0"
    #   freq_sum will be the sum of the frequency fields (column 5) from row 0 to row_count
    local freq_sum="0"
    #   While freq_sum < 0.25
    while [ $(echo "if (${freq_sum} < 0.25) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        #   freq is the value of the frequency field (column 5) on the row corresponding to row_count
        local freq=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        #   Add freq to freq_sum until the while loop exits
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    #   The first quantile is the coverage on the row at which the cumulative frequency hits 0.25 or greater
    local Q1=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Repeat for Q2 (median)
    while [ $(echo "if (${freq_sum} < 0.5) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        local freq=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    local Q2=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Repeat for Q3
    while [ $(echo "if (${freq_sum} < 0.75) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        local freq=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    local Q3=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Append the statistics to the summary file
    echo -e "${sampleName}"'\t'"${min}"'\t'"${Q1}"'\t'"${mode}"'\t'"${Q2}"'\t'"${mean}"'\t'"${Q3}"'\t'"${max}" >> "${out_dir}/${project}_coverage_summary_unfinished.txt"
    #   Put a call to plotCoverage here
}

export -f EC_Coverage

#   This is to calculate histograms and summary statistics for whole genome data
function WG_Coverage() {
    local bam_file="$1"
    local bam_dir=$(dirname "${bam_file}")
    local out_dir="$2"
    local project="$3"
    #   Get the sample name without the .bam
    local sampleName=$(basename "${bam_file}" .bam)
    #   Generate coverage histograms as text files
    bedtools genomecov -ibam "${bam_file}" > ${out_dir}/Histograms/${sampleName}.hist
    #   Begin calculating statistics per bp
    #   The minimum is the coverage on the first line of the "genome" fields since they're already sorted
    local min=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n 1 | awk -F "\t" '{print $2}')
    #   The maximum is the coverage on the last line of the "all" fields
    local max=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | tail -n 1 | awk -F "\t" '{print $2}')
    #   The mean is the sum of (each coverage * the percent of the genome at that coverage)
    local mean=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | awk '{ sum += $2*$5 } END { print sum }')
    #   The mode is the coverage that has the highest percent of the genome at that coverge (excludes zero coverage)
    local mode=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | tail -n +2 | sort -grk5,5 | head -1 | awk -F "\t" '{print $2}')
    #   The quantiles are a bit tricky...
    #   row_count will count how many rows down the "all" fields we are
    local row_count="0"
    #   freq_sum will be the sum of the frequency fields (column 5) from row 0 to row_count
    local freq_sum="0"
    #   While freq_sum < 0.25
    while [ $(echo "if (${freq_sum} < 0.25) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        #   freq is the value of the frequency field (column 5) on the row corresponding to row_count
        local freq=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        #   Add freq to freq_sum until the while loop exits
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    #   The first quantile is the coverage on the row at which the cumulative frequency hits 0.25 or greater
    local Q1=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Repeat for Q2 (median)
    while [ $(echo "if (${freq_sum} < 0.5) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        local freq=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    local Q2=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Repeat for Q3
    while [ $(echo "if (${freq_sum} < 0.75) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        local freq=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $5}')
        local freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done
    local Q3=$(grep "genome" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')
    #   Append the statistics to the summary file
    echo -e "${sampleName}"'\t'"${min}"'\t'"${Q1}"'\t'"${mode}"'\t'"${Q2}"'\t'"${mean}"'\t'"${Q3}"'\t'"${max}" >> "${out_dir}/${project}_coverage_summary_unfinished.txt"
    #   Put a call to plotCoverage here
}

export -f WG_Coverage

#   The main function that sets up and calls the various others
function Coverage_Mapping() {
    local sampleList="$1" # What is our list of samples?
    local outDirectory="$2"/Coverage_Mapping # Where do we store our results?
    local regions="$3" # What is our regions file?
    local proj="$4" # What is the name of the project?
    makeOutDirectories "${outDirectory}" # Make our output directories
    if ! [[ -f "${REGIONS_FILE}" ]]
    then # Whole-genome sequencing
        proj="${regions}" # Because regions was empty, the code read the project variable into the regions slot. This fixes it.
        #   Make the header for the summary file
        echo -e "Sample name\tMin\t1st Q\tMode\tMedian\tMean\t3rd Q\tMax" >> "${outDirectory}/${proj}_coverage_summary_unfinished.txt"
        parallel --jobs 4 --xapply WG_Coverage {1} "${outDirectory}" "${proj}" :::: "${sampleList}"
    else # Exome capture
        #   Make the header for the summary file
        echo -e "Sample name\tMin\t1st Q\tMode\tMedian\tMean\t3rd Q\tMax" >> "${outDirectory}/${proj}_coverage_summary_unfinished.txt"
        parallel --jobs 4 --xapply EC_Coverage {1} "${regions}" "${outDirectory}" "${proj}" :::: "${sampleList}"
    fi
    #   Make the header for the sorted summary file
    echo -e "Sample name\tMin\t1st Q\tMode\tMedian\tMean\t3rd Q\tMax" >> "${outDirectory}/${proj}_coverage_summary.txt"
    #   Sort the summary file based on sample name
    tail -n +2 "${outDirectory}/${proj}_coverage_summary_unfinished.txt" | sort >> "${outDirectory}/${proj}_coverage_summary.txt"
    #   Remove the unsorted file
    rm "${outDirectory}/${proj}_coverage_summary_unfinished.txt"
}

export -f Coverage_Mapping
