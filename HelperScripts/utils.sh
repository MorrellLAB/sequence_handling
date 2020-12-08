#!/bin/bash

#   Check to make sure our samples exist
function checkSamples() {
    local sample_list="$1" # Sample list
    if [[ -f "${sample_list}" ]] # If the sample list exists
    then
        if [[ -r "${sample_list}" ]] # If the sample list is readable
        then
            for sample in $(cat "${sample_list}") # For each sample in the sample list
            do
                if ! [[ -f "${sample}" ]] # If the sample doesn't exist
                then
                    echo "The sample ${sample} does not exist, exiting..." >&2 # Exit out with error
                    return 1
                else
                    if ! [[ -r "${sample}" ]] # If the sample isn't readable
                    then
                        echo "The sample ${sample} does not have read permissions, exiting..." >&2
                        return 1
                    fi
                fi
            done
        else # If the sample isn't readable
            echo "The sample list ${sample_list} does not have read permissions, exiting..." >&2
            return 1
        fi
    else # If the sample list doesn't exist
        echo "The sample list $sample_list does not exist, exiting..." >&2 # Exit out with error
        return 1
    fi
}

#   Export the function to be used elsewhere
export -f checkSamples

#   Check only custom job array samples exist for Read Mapping handler
function checkSamplesCustomJobArrRM() {
    local sample_list="$1" # Sample list
    local custom_job_arr="$2" # Custom job arrays following -t flag
    local singles_trimmed="$3"
    local forward_trimmed="$4"
    local reverse_trimmed="$5"
    # Prep custom job arrays, expand numbers in a range (i.e., 3-6)
    temp=($(echo ${custom_job_arr} | tr ',' '\n'))
    new_custom_arr=()
    for i in ${temp[@]}
    do
        if [[ ${i} == *"-"* ]]; then
            temp_range=($(echo ${i} | tr '-' '\n'))
            # Add to array
            new_custom_arr+=($(seq ${temp_range[0]} ${temp_range[1]}))
        else
            # Add to array
            new_custom_arr+=(${i})
        fi
    done
    if [[ -f "${sample_list}" ]] # If the sample list exists
    then
        if [[ -r "${sample_list}" ]] # If the sample list exists
        then
            single_samples=($(grep -E "${singles_trimmed}" "${sample_list}")) # Get the single-end samples
            forward_samples=($(grep -E "${forward_trimmed}" "${sample_list}")) # Get the forward samples
            reverse_samples=($(grep -E "${reverse_trimmed}" "${sample_list}")) # Get the reverse samples
            # If we have paired end samples
            if [[ ! -z "${forward_samples[@]}" && ! -z "${reverse_samples[@]}" ]] # If we have paired-end samples
            then
                for i in ${new_custom_arr[@]}
                do
                    # Check if samples exist
                    # Forward samples
                    if ! [[ -f "${forward_samples[$i]}" ]] # If the sample doesn't exist
                    then
                        echo "The sample ${forward_samples[$i]} does not exist, exiting..." >&2
                        return 1 # Exit out with error
                    else
                        if ! [[ -r "${forward_samples[$i]}" ]] # If the sample isn't readable
                        then
                            echo "The sample ${forward_samples[$i]} does not have read permissions, exiting..." >&2
                            return 1
                        fi
                    fi
                    # Reverse samples
                    if ! [[ -f "${reverse_samples[$i]}" ]] # If the sample doesn't exist
                    then
                        echo "The sample ${reverse_samples[$i]} does not exist, exiting..." >&2
                        return 1 # Exit out with error
                    else
                        if ! [[ -r "${reverse_samples[$i]}" ]] # If the sample isn't readable
                        then
                            echo "The sample ${reverse_samples[$i]} does not have read permissions, exiting..." >&2
                            return 1
                        fi
                    fi
                done
            else # We have single end samples
                for i in ${new_custom_arr[@]}
                do
                    # Check if samples exist
                    if ! [[ -f "${single_samples[$i]}" ]] # If the sample doesn't exist
                    then
                        echo "The sample ${single_samples[$i]} does not exist, exiting..." >&2
                        return 1
                    else
                        if ! [[ -r "${single_samples[$i]}" ]] # If the sample isn't readable
                        then
                            echo "The sample ${single_samples[$i]} does not have read permissions, exiting..." >&2
                            return 1
                        fi
                    fi
                done
            fi
        else # If the sample isn't readable
            echo "The sample list ${sample_list} does not have read permissions, exiting..." >&2
            return 1 # Exit out with error
        fi
    else # If the sample list doesn't exist
        echo "The sample list ${sample_list} does not exist, exiting..." >&2
        return 1
    fi
}

export -f checkSamplesCustomJobArrRM

#   Check to make sure our dependencies are installed
function checkDependencies() {
    local dependencies=("${!1}") # BASH array to hold dependencies
    for dep in "${dependencies[@]}" # For each dependency
    do
        if ! `command -v "$dep" > /dev/null 2> /dev/null` # If it's not installed
        then
            echo "Failed to find $dep installation, exiting..." >&2 # Write error message
            return 1 # Exit out with error
        fi
    done
}

#   Export the function to be used elsewhere
export -f checkDependencies

# Compares two tuple-based, dot-delimited version numbers a and b (possibly
# with arbitrary string suffixes). Returns:
# 1 if a<b
# 2 if equal
# 3 if a>b
# Everything after the first character not in [0-9.] is compared
# lexicographically using ASCII ordering if the tuple-based versions are equal.
function compare-versions() {
    if [[ $1 == $2 ]]; then
        return 2
    fi
    local IFS=.
    local i a=(${1%%[^0-9.]*}) b=(${2%%[^0-9.]*})
    local arem=${1#${1%%[^0-9.]*}} brem=${2#${2%%[^0-9.]*}}
    for ((i=0; i<${#a[@]} || i<${#b[@]}; i++)); do
        if ((10#${a[i]:-0} < 10#${b[i]:-0})); then
            return 1
        elif ((10#${a[i]:-0} > 10#${b[i]:-0})); then
            return 3
        fi
    done
    if [ "$arem" '<' "$brem" ]; then
        return 1
    elif [ "$arem" '>' "$brem" ]; then
        return 3
    fi
    return 2
}
export -f compare-versions

#   Check versions of tools to see if it meets the minimum version
function checkVersion() {
    local tool="$1"
    local minVersion="$2"
    # Make sure ${tool} is installed
    if [[ "${tool}" -eq "gatk" ]]; then
        GATK_JAR=$(checkGATK ${GATK_JAR})
        "${tool}" --version > /dev/null 2>&1
        retVal=$?
    else
        "${tool}" --version > /dev/null 2>&1
        retVal=$?
    fi
    if [ $retVal -ne 0 ]; then
	    echo "Please make sure ${tool} are installed and under your PATH"
	    return 1
    fi
    # Get the installed version
    local installedVers=""
    if [[ "${tool}" == "samtools" ]]; then
	    installedVer=$(samtools --version-only | perl -pe 's/\+htslib-[\d\.]*\s*$//')
    elif [[ "${tool}" == "bedtools" ]]; then
	    installedVer=$(bedtools --version | perl -pe 's/^bedtools\s+v//')
    elif [[ "${tool}" == "gatk" ]]; then
        # GATK3 uses "GenomeAnalysisTK.jar" file while in GATK4, this jar file naming no longer exists
        if [[ ${GATK_JAR} == *"GenomeAnalysisTK.jar"* ]]
        then
            installedVer=$(java -jar ${GATK_JAR} --version | cut -d'-' -f 1)
	elif [[ ${GATK_JAR} == *.jar ]]
	then
		installedVer=$(java -jar ${GATK_JAR} --version | grep "Genome Analysis Toolkit" | grep -Eo '[0-9]+[\.0-9]*')
        else # if GATK_JAR is pointing to the wrapper script
	        installedVer=$(${GATK_JAR} --version | grep "Genome Analysis Toolkit" | grep -Eo '[0-9]+[\.0-9]*')
        fi
    else
	    echo "ERROR: checkVersion() in utils.sh doesn't know how to check the version of ${tool}"
	    return 1
    fi
    compare-versions $minVersion $installedVer
    if [ $? -gt 2 ]; then
	    return 1  # fail
    else
	    return 0  # min version met
    fi
}

#   Export the function to be used elsewhere
export -f checkVersion

#   Figure out memory requirements based on Qsub or sbatch settings
#   This code written by Paul Hoffman for the RNA version of sequence handling at https://github.com/LappalainenLab/sequence_handling/
#   Function modified to work with Slurm
function getMemory() {
    local settings="$1"
    local wkld_manager="$2" # Are we using PBS or Slurm?
    # Pull out memory from sbatch settings
    if [[ ${wkld_manager} == "pbs" ]]; then
        # qsub settings
        mem_raw=$(echo "${settings}" | grep -oE 'mem=[[:alnum:]]+' | cut -f 2 -d '=')
    elif [[ ${wkld_manager} == "slurm" ]]; then
        # sbatch settings
        mem_raw=$(echo ${settings} | tr ' ' '\n' | grep "mem" | cut -d'=' -f 2)
    fi
    # Get just the number
    mem_digits=$(echo "${mem_raw}" | grep -oE '[[:digit:]]+')
    # Format memory for use in handler
    if $(echo "${mem_raw}" | grep -i 'g' > /dev/null 2> /dev/null)
    then
        max_mem="${mem_digits}g"
    elif $(echo "${mem_raw}" | grep -i 'm' > /dev/null 2> /dev/null)
    then
        max_mem="${mem_digits}M"
    elif $(echo "${mem_raw}" | grep -i 'k' > /dev/null 2> /dev/null)
    then
        max_mem="${mem_digits}K"
    else
        max_mem="${mem_digits}"
    fi
    echo "${max_mem}" # Return just the memory setting
}

export -f getMemory

function getThreads() {
    local settings="$1"
    local wkld_manager="$2" # Are we using PBS or Slurm?
    if [[ ${wkld_manager} == "pbs" ]]; then
        # qsub settings
        threads=$(echo "${settings}" | grep -oE 'ppn=[[:alnum:]]+' | cut -d '=' -f 2)
    elif [[ ${wkld_manager} == "slurm" ]]; then
        # sbatch settings
        threads=$(echo ${settings} | tr ' ' '\n' | grep "ntasks-per-node" | cut -d '=' -f 2)
    fi
    # Return the number of threads
    echo "${threads}"
}

export -f getThreads

#   A function to check to make sure Picard is where it actually is
function checkPicard() {
    local Picard="$1" # Where is Picard?
    if ! [[ -f "${Picard}" ]]; then echo "Failed to find Picard, exiting..." >&2; return 1; fi # If we can't find Picard, exit with error
}

#   Export the function
export -f checkPicard

#   A function to check to make sure Picard is where it actually is
#   If the variable isn't set, it will check GATK_LOCAL_JAR, and return the path
function checkGATK() {
    local GATK="$1" # Where is GATK?
    if ! [[ -f "${GATK}" ]]; then
        if [[ -f "${GATK_LOCAL_JAR}" ]]; then
            echo "${GATK_LOCAL_JAR}"  # with gatk4, the env var should be set
            return 0
        else
            echo "Failed to find GATK, exiting..." >&2
            echo 1 # If we can't find GATK, exit with error
            return 1
        fi
    else
	    echo "${GATK}"
	    return 0
    fi
}

#   Export the function to be used elsewhere
export -f checkGATK

#   A function to check the first line of a VCF file
#       If the first line isn't ##fileformat= then Variant_Recalibrator will not be able to read it
function checkVCF() {
    local vcf="$1"
    if ! [[ -r "${vcf}" ]]; then echo "${vcf} does not have read permissions, exiting..." >&2; return 11; fi # If the vcf isn't readable, exit
    local firstline=$(head -n 1 "${vcf}" | grep "##fileformat=")
    if [[ -z "${firstline}" ]]; then echo "${vcf} is not parseable by Variant_Recalibrator. Make sure that the first line is ##fileformat. Exiting..." >&2; return 12; fi
}

#   Export the function to be used elsewhere
export -f checkVCF

#   A function to see if our referenced FASTA has a .dict file
function checkDict() {
    local reference="$1" # What is our reference FASTA file?
    if ! [[ -f "${reference}" ]]; then echo "Cannot find reference genome, exiting..." >&2; exit 31; fi # Make sure it exists
    if ! [[ -r "${reference}" ]]; then echo "Reference genome does not have read permissions, exiting..." >&2; exit 30; fi # Make sure we can read it
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference
    local referenceBase="${referenceName%.*}" # Get the basename of the reference without extension since it could be either .fa or .fasta
    if [[ ! $(ls "${referenceDirectory}" | grep "${referenceBase}.dict" ) ]]; then return 1; fi # Check that we have the dict file, if we don't then return 1 (not exit 1) so that we can make it
    if [[ ! -r "${referenceDirectory}"/"${referenceBase}.dict" ]]; then echo "Reference dictionary file does not have read permissions, exiting..." >&2; exit 29; fi # Make sure we can read the dict file if it exists
}

#   Export the function
export -f checkDict

#   A function to generate a dictionary file for a reference
function createDict() {
    local reference="$1" # What is our reference FASTA file?
    local memory="$2" # How much memory can Picard use?
    local picard="$3" # Where is the Picard jar?
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference
    local referenceBase="${referenceName%.*}" # Get the basename of the reference without extension since it could be either .fa or .fasta
    if ! [[ -w "${referenceDirectory}" ]]; then echo "Cannot create reference dictionary file because you do not have write permissions for ${referenceDirectory}, exiting..." >&2; exit 28; fi # Make sure we can create the dict file
    # Don't need to check for Java because createDict is only used with GATK handlers, which already check for it
    checkPicard "${picard}"
    if [[ "$?" -ne 0 ]]; then exit 32; fi
    # Make the dict file
    (set -x; java -Xmx"${memory}" -jar "${picard}" CreateSequenceDictionary \
        R="${reference}" \
        O="${referenceDirectory}/${referenceBase}.dict")
    if [[ "$?" -ne 0 ]]; then echo "Error creating reference dictionary, exiting..."; exit 27; fi
}

#   Export the function
export -f createDict

#   Check to make sure our BAM files are indexed
function checkBaiIndex() {
    local sample_list="$1" # Sample list of BAM files, already checked by checkSamples (above)
    for sample in $(cat "${sample_list}") # For each sample in the sample list
    do
        local basename=$(basename "${sample}" .bam)
        local dirname=$(dirname "${sample}")
        if [[ ! -f "${dirname}/${basename}.bai" && ! -f "${dirname}/${basename}.bam.bai" ]] # If the sample doesn't have a .bai index file of either naming convention
        then
            echo "The sample ${sample} does not have a .bai index, exiting..." >&2
            exit 32
        fi
    done
}

#   Export the function
export -f checkBaiIndex

#   A function to see if our reference FASTA has a .mmi file
function checkMinimap2Index() {
    local reference="$1" # What is our reference FASTA file?
    if ! [[ -f "${reference}" ]]; then echo "Cannot find reference genome, exiting..." >&2; exit 31; fi # Make sure it exists
    if ! [[ -r "${reference}" ]]; then echo "Reference genome does not have read permissions, exiting..." >&2; exit 30; fi # Make sure we can read it
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference
    local referenceBase="${referenceName%.*}" # Get the basename of the reference without extension since it could be either .fa or .fasta
    if [[ ! $(ls "${referenceDirectory}" | grep "${referenceBase}.mmi" ) ]]; then return 1; fi # Check that we have the .mmi file, if we don't then return 1 (not exit 1) so that we can make it
    if [[ ! -r "${referenceDirectory}"/"${referenceBase}.mmi" ]]; then echo "Reference index file does not have read permissions, exiting..." >&2; exit 29; fi # Make sure we can read the .mmi file if it exists
}

#   Export the function
export -f checkMinimap2Index

# Check to make sure our .g.vcf files (output from Haplotype Caller) are indexed
function checkGvcfIndex() {
    local sample_list="$1" # List full filepaths to .g.vcf files
    local out_dir="$2" # OUT_DIR listed in config file
    # Remove existing temp list if it exists so we have a clean list to start appending to
    if [[ -f ${out_dir}/temp_re-index_gvcf_list.txt ]]
    then
        rm ${out_dir}/temp_re-index_gvcf_list.txt
    fi
    # Check if index exists, if not, add it to temporary list to re-run
    for sample in $(cat "${sample_list}")
    do
        if [[ ! -f "${sample}.idx" ]]
        then
            echo "The sample ${sample} does not have a .idx index file." >&2
            # Save sample to re-run and index in temporary file
            printf "${sample}\n" >> ${out_dir}/temp_re-index_gvcf_list.txt
        fi
    done
    # Check if temporary list to re-index exists
    if [[ -f ${out_dir}/temp_re-index_gvcf_list.txt ]]
    then
        echo "List of all .g.vcf files to index can be found at: ${out_dir}/temp_re-index_gvcf_list.txt" >&2
        echo "Please re-run Haplotype Caller to index these .g.vcf files before proceeding with Genomics_DB_Import, exiting..." >&2
        exit 32
    fi
}

export -f checkGvcfIndex

