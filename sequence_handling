#!/bin/bash

set -o pipefail

function Usage() {
    echo -e "\
Usage:  `basename $0` <handler> <config> [-t custom_array_indices]
Where:  <handler> is one of:
            Recommended workflow:
                1 | Quality_Assessment
                2 | Adapter_Trimming
                3 | Read_Mapping
                4 | SAM_Processing
                5 | Coverage_Mapping
                6 | Haplotype_Caller
                7 | Genomics_DB_Import
                8 | Genotype_GVCFs
                9 | Create_HC_Subset
                10 | Variant_Recalibrator
                11 | Pre_Variant_Filtering
                12 | Variant_Filtering
                13 | Variant_Analysis
            Other handlers:
                14 | GBS_Demultiplex (in progress)
                15 | Quality_Trimming
                16 | Realigner_Target_Creator
                17 | Indel_Realigner
                18 | Freebayes_Variant_Calling
            Nanopore workflow:
                1NP | NP_Quality_Assessment (in progress)
                2NP | NP_Adapter_Trimming (in progress)
                3NP | NP_Read_Mapping (in progress)
                4NP | NP_SAM_Processing (in progress)
And:    <config> is the full file path to the configuration file

Optional arguments:
        [-t custom_array_indices] is a range of arrays and/or comma separated list of specific arrays to run.
        Important: No spaces allowed when providing range or list of arrays.
            and -t flag must be provided as the 3rd argument on the command line (see example below)
        If left blank, the DEFAULT runs all samples.
        This is helpful if only some of your jobs arrays fail and you need to re-run only those.

        Example: ./sequence_handling SAM_Processing /path/to/config -t 1-5,10,12
" >&2
    exit 1
}

export -f Usage

#   Where is 'sequence_handling' located?
#   Mod by Naoki
SCRIPT=$(realpath $0)  # Track back the symlink to find the actual location
SEQUENCE_HANDLING=$(dirname "${SCRIPT}")

#   A list of valid sequencing platforms
VALID_SEQ_PLATFORMS=('CAPILLARY' 'LS454' 'ILLUMINA' 'SOLID' 'HELICOS' 'IONTORRENT' 'ONT' 'PACBIO')

#   If we have less than two arguments
if [[ "$#" -lt 1 ]]; then Usage; fi # Display the usage message and exit

ROUTINE="$1" # What routine are we running?
CONFIG=$(realpath $2) # Where is our config file?

#   If the specified config exists
if [[ -f "${CONFIG}" ]]
then
    source "${CONFIG}" # Source it, providing parameters and software
    bash "${CONFIG}" > /dev/null 2> /dev/null # Load any modules
    source "${SEQUENCE_HANDLING}"/HelperScripts/utils.sh # And the utils script
else # If it doesn't
    echo "Please specify a valid config file." >&2 # Print error message
    exit 1 # Exit with non-zero exit status
fi

#   After loading Config, make sure that we will be able to properly output files
mkdir -p "${OUT_DIR}"
if ! [[ -w "${OUT_DIR}" ]]; then echo "You don't have write permissions for the output directory ${OUT_DIR}, exiting..." >&2; exit 1; fi

#   Custom list of job arrays to run using a -t flag
#   Currently, -t flag has to be provided as 3rd argument
case "$3" in
    -t)
        echo "Using custom list of job arrays provided from -t flag on the command line..." >&2
        #   Do we have a list of custom job arrays to run?
        #   Acceptable formats include range separated by hyphen, comma separated list, or mix of both
        CUSTOM_JOB_ARR="$4"
    ;;
esac

#   Where do we output the standard error and standard output files?
ERROR="${OUT_DIR}"/Error_Files
mkdir -p "${ERROR}"

#   Run sequence_handling
case "${ROUTINE}" in
    1 | Quality_Assessment)
        echo "$(basename $0): Assessing quality..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Quality_Assessment.sh
        checkDependencies Quality_Assessment_Dependencies[@] # Check to see if dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${QA_SAMPLES}" # Check that the samples and sample list exist
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == "true" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Quality_Assessment.sh && Quality_Assessment ${QA_SAMPLES} ${OUT_DIR} ${PROJECT} ${TARGET} ${MSI}" | qsub -q "${QA_QUEUE}" -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Assessment
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            sbatch --job-name=${PROJECT}_Quality_Assessment ${QA_SBATCH} -o ${ERROR}/${PROJECT}_Quality_Assessment.%j.out -e ${ERROR}/${PROJECT}_Quality_Assessment.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Quality_Assessment.job ${CONFIG} ${SEQUENCE_HANDLING}
        else
            (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Quality_Assessment.sh" && Quality_Assessment "${QA_SAMPLES}" "${OUT_DIR}" "${PROJECT}" "${TARGET} ${MSI}") > "${ERROR}/Quality_Assessment.log" 2>&1
        fi
        ;;
    2 | Adapter_Trimming )
        echo "$(basename $0): Trimming adapters..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Adapter_Trimming.sh
        checkDependencies Adapter_Trimming_Dependencies[@] # Check to see if dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${RAW_SAMPLES}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        if ! [[ -f "${ADAPTERS}" ]]; then echo "Please specify a valid adapters file" >&2; exit 1; fi # Check for a valid adapters file
        if [[ -z "${QUAL_ENCODING}" ]]; then echo "Please specify the QUAL_ENCODING in the config file" >&2; exit 1; fi # Make sure the platform is filled out
        #   Run Adapter_Trimming using a task array
        declare -a AT_LIST=($(grep -E ".fastq|.fastq.gz" "${RAW_SAMPLES}"))
        SINGLE_ARRAY_LIMIT=$[${#AT_LIST[@]} - 1] # Get the maximum number of Torque tasks we're doing for our samples
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == "true" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Adapter_Trimming.sh && Adapter_Trimming ${RAW_SAMPLES} ${OUT_DIR} ${PROJECT} ${FORWARD_NAMING} ${REVERSE_NAMING} ${ADAPTERS} ${PRIOR} ${QUAL_ENCODING}" | qsub -t "${job_array_index_range}" -q "${AT_QUEUE}" -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Adapter_Trimming
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            sbatch --job-name=${PROJECT}_Adapter_Trimming --array=${job_array_index_range} ${AT_SBATCH} -o ${ERROR}/${PROJECT}_Adapter_Trimming.%A_%a.out -e ${ERROR}/${PROJECT}_Adapter_Trimming.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Adapter_Trimming.job ${CONFIG} ${SEQUENCE_HANDLING}
	    else
	        #   Non PBS processing
	        (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Adapter_Trimming.sh && Adapter_Trimming ${RAW_SAMPLES} ${OUT_DIR} ${PROJECT} ${FORWARD_NAMING} ${REVERSE_NAMING} ${ADAPTERS} ${PRIOR} ${QUAL_ENCODING}) > ${ERROR}/Adapter_Trimming.log 2>&1
            # I'm adding this here since default Config has ADAPTED_LIST pointing to this file, Naoki June 5, 2020
            # CL Note: Moved within non-PBS processing because this only works when using gnu parallel and not PBS job arrays.
	        ${SEQUENCE_HANDLING}/HelperScripts/sample_list_generator.sh .fastq.gz ${OUT_DIR}/Adapter_Trimming ${PROJECT}_trimmed_adapters.txt
	    fi
        ;;
    3 | Read_Mapping )
        echo "`basename $0`: Mapping reads...">&2
        source "${SEQUENCE_HANDLING}"/Handlers/Read_Mapping.sh
        checkDependencies Read_Mapping_Dependencies[@] # Check to see if dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        # Check samples
        if [[ -z ${CUSTOM_JOB_ARR+x} ]] # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
        then
            # Check to see if all samples and sample list exists
            checkSamples "${TRIMMED_LIST}" # Check to see if samples and sample list exists
        else
            # Only check that samples in CUSTOM_JOB_ARR exists
            # Purpose: if we need to re-run some samples and do not want to re-download all the samples
            # (due to space limitation), this allows checking of only those samples to re-run without having
            # to change our sample list. The sample list can still include all the samples.
            echo "Using custom job arrays following -t flag. Only check if custom job array samples exist."
            checkSamplesCustomJobArrRM "${TRIMMED_LIST}" "${CUSTOM_JOB_ARR}" "${SINGLES_TRIMMED}" "${FORWARD_TRIMMED}" "${REVERSE_TRIMMED}"
        fi
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkIndex "${REF_GEN}" # Check to make sure our reference genome is indexed
        if [[ "$?" -ne 0 ]]; then echo "Reference genome is not indexed for BWA mem..." >&2; indexReference "${REF_GEN}"; fi # If not, index it and exit
        if [[ -z "${SEQ_PLATFORM}" ]]; then echo "Please specify the SEQ_PLATFORM in the config file" >&2; exit 1; fi # Make sure the platform is filled out
        [[ "${VALID_SEQ_PLATFORMS[@]}" =~ "${SEQ_PLATFORM}" ]] || (echo -e "'${SEQ_PLATFORM}' is not a valid platform\nPlease choose from:" >&2; for PLAT in ${VALID_SEQ_PLATFORMS[@]}; do echo -e "\t${PLAT}"; done; exit 1)
        declare -a SINGLE_SAMPLES=($(grep -E "${SINGLES_TRIMMED}" "${TRIMMED_LIST}")) # Get the single-end samples
        declare -a FORWARD_SAMPLES=($(grep -E "${FORWARD_TRIMMED}" "${TRIMMED_LIST}")) # Get the forward samples
        declare -a REVERSE_SAMPLES=($(grep -E "${REVERSE_TRIMMED}" "${TRIMMED_LIST}")) # Get the reverse samples
        if [[ ! -z "${FORWARD_SAMPLES[@]}" && ! -z "${REVERSE_SAMPLES[@]}" ]] # If we have paired-end samples
        then
            declare -a PAIRED_NAMES=($(parallel basename {} "${FORWARD_TRIMMED}" ::: "${FORWARD_SAMPLES[@]}")) # Create an array of paired-end sample names
        fi
	    if ! [[ -z "${SINGLE_SAMPLES[@]}" ]]; then # If we have single-end samples
            MODE="single";
            SINGLE_ARRAY_LIMIT=$[${#SINGLE_SAMPLES[@]} - 1]
        elif [[ ! -z "${FORWARD_SAMPLES[@]}" && ! -z "${REVERSE_SAMPLES[@]}" ]]; then # If we have paired-end samples
            MODE="paired";
            SINGLE_ARRAY_LIMIT=$[${#FORWARD_SAMPLES[@]} - 1]
            # Make sure we have equal numbers of forward and reverse samples
            if [[ "${#FORWARD_SAMPLES[@]}" -ne "${#REVERSE_SAMPLES[@]}" ]]; then echo "Unequal numbers of forward and reverse reads, exiting..." >&2; exit 1; fi
        else
            echo "ERROR: No samples match the provided suffix \"${FORWARD_TRIMMED}\", exiting..." >&2;
            exit 17;
        fi
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
	    elif [[ "$USE_PBS" == "true" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Read_Mapping.sh && Read_Mapping ${TRIMMED_LIST} ${SINGLES_TRIMMED} ${FORWARD_TRIMMED} ${REVERSE_TRIMMED} ${MODE} ${PROJECT} ${SEQ_PLATFORM} ${OUT_DIR} ${REF_GEN}" | qsub -t "${job_array_index_range}" -q "${RM_QUEUE}" -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Read_Mapping
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            sbatch --job-name=${PROJECT}_Read_Mapping --array=${job_array_index_range} ${RM_SBATCH} -o ${ERROR}/${PROJECT}_Read_Mapping.%A_%a.out -e ${ERROR}/${PROJECT}_Read_Mapping.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Read_Mapping.job ${CONFIG} ${SEQUENCE_HANDLING} ${MODE}
	    else # without PBS
	        (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Read_Mapping.sh && Read_Mapping ${TRIMMED_LIST} ${SINGLES_TRIMMED} ${FORWARD_TRIMMED} ${REVERSE_TRIMMED} ${MODE} ${PROJECT} ${SEQ_PLATFORM} ${OUT_DIR} ${REF_GEN}) > ${ERROR}/Read_Mapping.log 2>&1
	    fi
        ;;
    4 | SAM_Processing )
        case "${METHOD}" in
            samtools )
                echo "$(basename $0): Processing SAM files using SAMtools..." >&2
                source "${SEQUENCE_HANDLING}"/Handlers/SAM_Processing_SAMtools.sh
                checkDependencies SAM_Processing_Dependencies[@] # Check to see if the dependencies are installed
                if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
                checkSamples "${MAPPED_LIST}" # Check to see if samples and sample list exists
                if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
                checkFaidx "${REF_GEN}" # Check to see if reference genome is indexed
                if [[ "$?" -ne 0 ]]; then echo "Reference genome is not indexed for SAM Processing...">&2; fadixReference "${REF_GEN}"; fi # If not, index and exit
                checkVersion 'samtools' '1.3' # Check SAMtools version 1.3 or higher
                if [[ "$?" -ne 0 ]]; then echo "Please use SAMtools version 1.3 or higher" >&2; exit 1; fi
                #   Run SAM_Processing using a task array
                declare -a SAM_LIST=($(grep -E ".sam" "${MAPPED_LIST}"))
                SINGLE_ARRAY_LIMIT=$[${#SAM_LIST[@]} - 1] # Get the maximum number of Torque tasks we're doing for SAM samples
                echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
                #   Run SAM_Processing Samtools
                if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                    echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                    exit 1
		        elif [[ ${USE_PBS} == true ]]; then
                    echo "PBS is our workload manager/job scheduler."
                    # If we have enough samples for a task array
                    if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                        # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                        if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                        then
                            echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                            job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                        else
                            # Use job arrays following -t flag
                            echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                            job_array_index_range="${CUSTOM_JOB_ARR}"
                        fi
                    else # If we only have one sample
                        echo "Using a single job"
                        job_array_index_range="0"
                    fi
                    echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/SAM_Processing_SAMtools.sh && SAM_Processing ${MAPPED_LIST} ${OUT_DIR} ${REF_GEN} ${PROJECT}" | qsub -t "${job_array_index_range}" -l "${SP_QSUB}" -q "${SP_QUEUE}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_SAM_Processing
                elif [[ "${USE_SLURM}" == true ]]; then
                    echo "Slurm is our workload manager/job scheduler."
                    # If we have enough samples for a task array
                    if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                        # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                        if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                        then
                            echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                            job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                        else
                            # Use job arrays following -t flag
                            echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                            job_array_index_range="${CUSTOM_JOB_ARR}"
                        fi
                    else # If we only have one sample
                        echo "Using a single job"
                        job_array_index_range="0"
                    fi
                    sbatch --job-name=${PROJECT}_SAM_Processing --array=${job_array_index_range} ${SP_SBATCH} -o ${ERROR}/${PROJECT}_SAM_Processing.%A_%a.out -e ${ERROR}/${PROJECT}_SAM_Processing.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/SAM_Processing_SAMtools.job ${CONFIG} ${SEQUENCE_HANDLING}
		        else
		            (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/SAM_Processing_SAMtools.sh && SAM_Processing ${MAPPED_LIST} ${OUT_DIR} ${REF_GEN} ${PROJECT}) > ${ERROR}/SAM_Processing_SAMtools.log 2>&1
		        fi
                ;;
            picard )
                echo "$(basename $0): Processing SAM files using Picard..." >&2
                source "${SEQUENCE_HANDLING}"/Handlers/SAM_Processing_Picard.sh
                checkDependencies SAM_Processing_Dependencies[@] # Check to see if the dependencies are installed
                if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
                checkSamples "${MAPPED_LIST}" # Check to see if samples and sample list exists
                if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
                checkPicard "${PICARD_JAR}" # Check to make sure Picard is installed
                if [[ "$?" -ne 0 ]]; then exit 1; fi # If we don't have Picard, exit with error
                checkVersion 'samtools' '1.3' # Check SAMtools version 1.3 or higher
                if [[ "$?" -ne 0 ]]; then echo "Please use SAMtools version 1.3 or higher" >&2; exit 1; fi
                if [[ -z "${SEQ_PLATFORM}" ]]; then echo "Please specify the SEQ_PLATFORM in the config file" >&2; exit 1; fi # Make sure the platform is filled out
                [[ "${VALID_SEQ_PLATFORMS[@]}" =~ "${SEQ_PLATFORM}" ]] || (echo -e "'${SEQ_PLATFORM}' is not a valid platform\nPlease choose from:" >&2; for PLAT in ${VALID_SEQ_PLATFORMS[@]}; do echo -e "\t${PLAT}"; done; exit 1)
                # Figure out memory requirements based on the Qsub or Slurm settings
                if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                    echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                    exit 1
                elif [[ ${USE_PBS} == true ]]; then
                    SP_MEM=$(getMemory "${SP_QSUB}" "pbs")
                elif [[ ${USE_SLURM} == true ]]; then
                    SP_MEM=$(getMemory "${SP_SBATCH}" "slurm")
                fi
                #   Create the header for the mapping stats summary file
                mkdir -p "${OUT_DIR}/SAM_Processing/Picard/Statistics"
                echo -e "Sample name\tTotal reads\tPercent mapped\tPercent paired\tPercent singletons\tFraction with mate mapped to different chr" > "${OUT_DIR}/SAM_Processing/Picard/Statistics/${PROJECT}_mapping_summary.tsv"
                #   Run SAM_Processing using a task array
                declare -a SAM_LIST=($(grep -E ".sam" "${MAPPED_LIST}"))
                SINGLE_ARRAY_LIMIT=$[${#SAM_LIST[@]} - 1] # Get the maximum number of Torque tasks we're doing for SAM samples
                echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
                # Run the handler
                if [[ "${USE_PBS}" == "true" ]]; then
                    echo "PBS is our workload manager/job scheduler."
                    # If we have enough samples for a task array
                    if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                        # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                        if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                        then
                            echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                            job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                        else
                            # Use job arrays following -t flag
                            echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                            job_array_index_range="${CUSTOM_JOB_ARR}"
                        fi
                    else # If we only have one sample
                        echo "Using a single job"
                        job_array_index_range="0"
                    fi
                    echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/SAM_Processing_Picard.sh && SAM_Processing ${MAPPED_LIST} ${OUT_DIR} ${PICARD_JAR} ${SEQ_PLATFORM} ${SP_MEM} ${MAX_FILES} ${PROJECT} ${PICARD_MAX_REC_IN_RAM} ${SORTING_COLL_SIZE_RATIO} ${TMP}" | qsub -t "${job_array_index_range}" -l "${SP_QSUB}" -q "${SP_QUEUE}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_SAM_Processing
                elif [[ "${USE_SLURM}" == true ]]; then
                    echo "Slurm is our workload manager/job scheduler."
                    # If we have enough samples for a task array
                    if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                        # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                        if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                        then
                            echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                            job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                        else
                            # Use job arrays following -t flag
                            echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                            job_array_index_range="${CUSTOM_JOB_ARR}"
                        fi
                    else # If we only have one sample
                        echo "Using a single job"
                        job_array_index_range="0"
                    fi
                    sbatch --job-name=${PROJECT}_SAM_Processing --array=${job_array_index_range} ${SP_SBATCH} -o ${ERROR}/${PROJECT}_SAM_Processing.%A_%a.out -e ${ERROR}/${PROJECT}_SAM_Processing.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/SAM_Processing_Picard.job ${CONFIG} ${SEQUENCE_HANDLING} ${SP_MEM}
		        else   # NO PBS
		            if [[ -z "${SAM_PROCESSING_THREADS}" ]]; then
			            SAM_PROCESSING_THREADS='100%'   # Use all cores if not defined
		            fi
		            (set -eo pipefail; source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/SAM_Processing_Picard.sh && printf '%s\n' "${SAM_LIST[@]}" | parallel --jobs ${SAM_PROCESSING_THREADS} "SAM_Processing {} ${OUT_DIR} ${PICARD_JAR} ${SEQ_PLATFORM} ${SP_MEM} ${MAX_FILES} ${PROJECT} ${TMP} ${PICARD_MAX_REC_IN_RAM} ${SORTING_COLL_SIZE_RATIO}") > ${ERROR}/SAM_Processing_Picard.log 2>&1
		        fi
                ;;
            * )
                echo "Invalid method"
                exit 1
                ;;
        esac
        ;;
    5 | Coverage_Mapping )
        echo "$(basename $0): Mapping coverage..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Coverage_Mapping.sh
        checkDependencies Coverage_Mapping_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkVersion bedtools 2.17.0 # Check to see that we have bedtools version 2.17.0 or newer
        if [[ "$?" -ne 0 ]]; then echo "Please use Bedtools version 2.17.0" >&2; exit 1; fi # If not, exit out with error
	    # With version 2.24.0 or newer, the behavior of bedtools coverage was changed.
        checkVersion bedtools 2.24.0 # Check to see that we have bedtools version 2.24.0
	    if [[ "$?" -eq 0 ]]; then
	        bedtoolsPre2_24_0="false"
	    else
	        bedtoolsPre2_24_0="true"
	    fi
        #checkVersion R 3.3 # Check to see that we have R version 3.3.X
        #if [[ "$?" -ne 0 ]]; then echo "Please use R version 3.3.X" >&2; exit 1; fi # If not, exit out with error
        checkSamples "${BAM_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        if ! [[ -f "${REGIONS_FILE}" ]]; then
            echo "No regions file found, assuming whole-genome sequencing data..." >&2
        else
            echo "Using regions file at ${REGIONS_FILE} for exome capture data..." >&2
        fi
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
	    elif [[ "$USE_PBS" == "true" ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Coverage_Mapping.sh && Coverage_Mapping ${BAM_LIST} ${OUT_DIR} ${PROJECT} ${bedtoolsPre2_24_0} ${REGIONS_FILE}" | qsub -q "${CM_QUEUE}" -l "${CM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Coverage_Mapping
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            sbatch --job-name=${PROJECT}_Coverage_Mapping ${CM_SBATCH} -o ${ERROR}/${PROJECT}_Coverage_Mapping.%j.out -e ${ERROR}/${PROJECT}_Coverage_Mapping.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Coverage_Mapping.job ${CONFIG} ${SEQUENCE_HANDLING} ${bedtoolsPre2_24_0}
	    else
	        (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Coverage_Mapping.sh && Coverage_Mapping ${BAM_LIST} ${OUT_DIR} ${PROJECT} ${bedtoolsPre2_24_0} ${REGIONS_FILE}) > ${ERROR}/Coverage_Mapping.log 2>&1
	    fi
        ;;
    6 | Haplotype_Caller )
        echo "$(basename $0): Calling SNPs with GATK Haplotype_Caller..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Haplotype_Caller.sh
        checkDependencies Haplotype_Caller_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${FINISHED_BAM_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkBaiIndex "${FINISHED_BAM_LIST}" # Check to see if our samples are indexed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If they're not indexed, exit out with error
	    # if GATK_JAR is not set in Config, it will check GATK_LOCAL_JAR
        GATK_JAR=$(checkGATK "${GATK_JAR}") # Check to make sure GATK is installed
        if [[ "${GATK_JAR}" == 1 ]]; then exit 1; fi # If we don't have GATK, exit with error
	    if checkVersion gatk 4.0 ; then
	        gatkVer=4
	    elif checkVersion gatk 3.8.0; then
	        gatkVer=3
	    else
	        echo "Please install GATK 3.8.0 or newer" >&2
	        exit 1
	    fi
        echo "GATK version: ${gatkVer}"
        if [[ -z "${DO_NOT_TRIM_ACTIVE_REGIONS}}" ]]; then echo "Please put 'true' or 'false' for DO_NOT_TRIM_ACTIVE_REGIONS in the config file" >&2; exit 1; fi # Make sure barley is filled out
        if [[ -z "${FORCE_ACTIVE}" ]]; then echo "Please put 'true' or 'false' for FORCE_ACTIVE in the config file" >&2; exit 1; fi # Make sure barley is filled out
        # Figure out memory requirements and number of threads we have based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            HC_MEM=$(getMemory "${HC_QSUB}" "pbs")
            HC_THREADS=$(getThreads "${HC_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            HC_MEM=$(getMemory "${HC_SBATCH}" "slurm")
            HC_THREADS=$(getThreads "${HC_SBATCH}" "slurm")
        fi
        echo "Checking to make sure our reference genome has a dict file..."
        checkDict "${REF_GEN}" "${PICARD_JAR}" # Check to make sure our reference genome has a dict file
        declare -a HC_LIST=($(grep -E ".bam" "${FINISHED_BAM_LIST}")) # Create an array of the BAM files

        # Check if we need to create *-fixed.bam for each input file
        # On the cluster, this function will be called within the Haplotyp_Caller_GATK4 function
        # Non-cluster, it will be called here
        if [[ "${FIX_QUALITY_SCORES}" == "true" ]] && [[ "${gatkVer}" == 4 ]] ; then
            if [[ "${USE_PBS}" == "true" ]] || [[ "${USE_SLURM}" == "true" ]]; then
                # CL Note: This needs to be moved to within the Haplotype_Caller.sh script for running on the cluster since there is no easy way to proceed to the next step, Haplotype Caller, with job arrays
                echo "Fix quality scores will be included in the Haplotype Caller run."
            else # No PBS
                (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Haplotype_Caller.sh" && Fix_Qscores "${FINISHED_BAM_LIST}" "${FIX_QUALITY_SCORES}") > "${ERROR}/Haplotype_Caller-Fix_Qscores.log" 2>&1
            fi
        fi

        # If we are parallelizing across custom regions for every sample, we will have
        # a total array limit of num_intervals * num_samples
        if [ "${HC_PARALLELIZE}" == "true" ] && [ "${HC_CUSTOM_INTERVALS}" != false ]; then
            # Get the maximum number of Torque tasks
            if [[ "${HC_SCAFFOLDS}" != "false" ]] && [ "${HC_PARALLELIZE}" == "true" ]; then
                # with additional scaffolds: (num_intervals+1) * num_samples - 1
                SINGLE_ARRAY_LIMIT=$(( ($(wc -l < "${HC_CUSTOM_INTERVALS}") + 1) * $(wc -l < "${FINISHED_BAM_LIST}") - 1 ))
            else
                # no additional scaffolds: num_intervals * num_samples - 1
                SINGLE_ARRAY_LIMIT=$(( $(wc -l < "${HC_CUSTOM_INTERVALS}") * $(wc -l < "${FINISHED_BAM_LIST}") - 1))
            fi
        elif [[ "${HC_PARALLELIZE}" == "false" ]] && [[ "${HC_SCAFFOLDS}" != "false" ]]; then
                # We are not parallelizing across regions so scaffolds list will just be appended to the run
                #   instead of run as a separate array index
                SINGLE_ARRAY_LIMIT=$[${#HC_LIST[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        else
            # Not parallelizing, total array limit is equal to the number of samples
            SINGLE_ARRAY_LIMIT=$[${#HC_LIST[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
            # If we also have scaffolds, add one to the array limit
            if [[ "${HC_SCAFFOLDS}" != "false" ]]; then
                SINGLE_ARRAY_LIMIT=${#HC_LIST[@]}
            fi
        fi
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == "true" ]]; then
            echo "PBS is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
            then # If we have enough samples for a task array
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we have only one sample
                job_array_index_range=0
            fi
            # Check if GATK 3 or GATK 4
            if checkVersion gatk 4.0; then
                echo "Using GATK 4..."
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Haplotype_Caller.sh && Haplotype_Caller_GATK4 ${FINISHED_BAM_LIST} ${OUT_DIR} ${GATK_JAR} ${REF_GEN} ${THETA} ${HC_MEM} ${HC_THREADS} ${FIX_QUALITY_SCORES} ${DO_NOT_TRIM_ACTIVE_REGIONS} ${FORCE_ACTIVE} ${HC_PARALLELIZE} ${HC_CUSTOM_INTERVALS} ${HC_SCAFFOLDS} ${TMP}" | qsub -t "${job_array_index_range}" -q "${HC_QUEUE}" -l "${HC_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Haplotype_Caller
            else
                # Assume we are using GATK 3
                echo "Using GATK 3..."
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Haplotype_Caller.sh && Haplotype_Caller_GATK3 ${FINISHED_BAM_LIST} ${OUT_DIR} ${GATK_JAR} ${REF_GEN} ${THETA} ${HC_MEM} ${HC_THREADS} ${FIX_QUALITY_SCORES} ${DO_NOT_TRIM_ACTIVE_REGIONS} ${FORCE_ACTIVE} ${HC_PARALLELIZE} ${HC_CUSTOM_INTERVALS} ${HC_SCAFFOLDS} ${TMP}" | qsub -t "${job_array_index_range}" -q "${HC_QUEUE}" -l "${HC_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Haplotype_Caller
            fi
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
            then # If we have enough samples for a task array
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we have only one sample
                job_array_index_range=0
            fi
            # Check if GATK 3 or GATK 4
            if checkVersion gatk 4.0; then
                echo "Using GATK 4..."
                sbatch --job-name=${PROJECT}_Haplotype_Caller --array=${job_array_index_range} ${HC_SBATCH} -o ${ERROR}/${PROJECT}_Haplotype_Caller.%A_%a.out -e ${ERROR}/${PROJECT}_Haplotype_Caller.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Haplotype_Caller_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING} ${HC_MEM} ${HC_THREADS}
            else
                # Assume we are using GATK 3
                echo "Using GATK 3..."
                sbatch --job-name=${PROJECT}_Haplotype_Caller --array=${job_array_index_range} ${HC_SBATCH} -o ${ERROR}/${PROJECT}_Haplotype_Caller.%A_%a.out -e ${ERROR}/${PROJECT}_Haplotype_Caller.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Haplotype_Caller_GATK3.job ${CONFIG} ${SEQUENCE_HANDLING} ${HC_MEM} ${HC_THREADS}
            fi
        else # without PBS
            if checkVersion gatk 4.0; then
                echo "Using GATK 4..."
                (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Haplotype_Caller.sh" && Haplotype_Caller_GATK4 "${FINISHED_BAM_LIST}" "${OUT_DIR}" "${GATK_JAR}" "${REF_GEN}" ${THETA} ${HC_MEM} ${HC_THREADS} ${FIX_QUALITY_SCORES} ${DO_NOT_TRIM_ACTIVE_REGIONS} ${FORCE_ACTIVE} ${gatkVer} ${HC_PARALLELIZE} "${HC_CUSTOM_INTERVALS}" "${HC_SCAFFOLDS}" "${TMP}") > ${ERROR}/Haplotype_Caller.log 2>&1
            else
                # Assume we are using GATK 3
                echo "Using GATK 3..."
                (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Haplotype_Caller.sh" && Haplotype_Caller_GATK3 "${FINISHED_BAM_LIST}" "${OUT_DIR}" "${GATK_JAR}" "${REF_GEN}" ${THETA} ${HC_MEM} ${HC_THREADS} ${FIX_QUALITY_SCORES} ${DO_NOT_TRIM_ACTIVE_REGIONS} ${FORCE_ACTIVE} ${gatkVer} ${HC_PARALLELIZE} "${HC_CUSTOM_INTERVALS}" "${HC_SCAFFOLDS}" "${TMP}") > ${ERROR}/Haplotype_Caller.log 2>&1
            fi
            # Currently this only works for non-PBS runs due to inability to track when all array indices have completed running
            if [[ "${FIX_QUALITY_SCORES}" == "true" ]] && [[ "${gatkVer}" == 4 ]] ; then
                # clean up *-fixed.bam
                for index in "${!HC_LIST[@]}"; do rm -f ${HC_LIST[${index}]}-fixed.bam; done
            fi
        fi
        ;;
    7 | Genomics_DB_Import )
        echo "$(basename $0): Genomics DB Import..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Genomics_DB_Import.sh
        checkDependencies Genomics_DB_Import_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${GVCF_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkGvcfIndex "${GVCF_LIST}" "${OUT_DIR}" # Check if .g.vcf files are indexed
	    # if GATK_JAR is not set in Config, it will check GATK_LOCAL_JAR
	    GATK_JAR=$(checkGATK "${GATK_JAR}") # Check to make sure GATK is installed
        if [[ "${GATK_JAR}" == 1 ]]; then exit 1; fi # If we don't have GATK, exit with error
	    rm -f "${ERROR}/Genomics_DB_Import.log"; touch "${ERROR}/Genomics_DB_Import.log"
        GDBI_MEM=$(getMemory "${GDBI_QSUB}")
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            GDBI_MEM=$(getMemory "${GDBI_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            GDBI_MEM=$(getMemory "${GDBI_SBATCH}" "slurm")
        fi
        # If we used custom intervals for the Haplotype_Caller step, we will need
        # some special processing of the data to get the arrays setup correctly
        #   Note: Check if HC_PARALLELIZE variable is set in case users did not use
        #   sequence_handling to run the Haplotype_Caller step
        if checkVersion gatk 4.0; then
            # With gatk4, individual vcf files have to be combined into a DB
            gatkVer=4
            # GenomicsDBImport() requires reference dict
            checkDict "${REF_GEN}" "${PICARD_JAR}" # Check to make sure our reference genome has a dict file

            ## Naoki Comments:
            ## Previously, I put the code to reuse pre-existing Genomics_DB because
            ## it was combined with Genotype_GVCFs handler.  So I thought that there
            ## could be a time when one might want to run Genotype_GVCFs multiple times.
            ## However, with the new structure with a separate Genomics_DB_Import,
            ## it is unlikely that we want to reuse the pre-exisiting Genomics_DB,
            ## So I removed the corresponding code.
            ## I also reorganized the logics a little bit to enable Parallelization
            ## at the HaplotypeCaller for non-PBS

            # When HC_PARALLIZE was used, these settings need to be overridden
            parallelizeQ=${PARALLELIZE}
            this_custom_intervals=${CUSTOM_INTERVALS}
            this_scaffolds=${SCAFFOLDS}
            # If variable is not unset and hc parallelize is true
            if [ ! -z ${HC_PARALLELIZE+x} ] && [ ${HC_PARALLELIZE} == "true" ] && [ ! -z "${HC_CUSTOM_INTERVALS}" ] && [ "${HC_CUSTOM_INTERVALS}" != "false" ] ; then
                analysisType="targeted-HC"
                # Parallelization occurs at Haplotype Caller phase
                echo "INFO: HC_PARALLELIZE is true, so PARRALELIZE is forced to be true."
                echo "Ignoring CUSTOM_INTERVALS & SCAFFOLDS and using HC_CUSTOM_INTERVALS & HC_SCAFFOLDS, instead."
                parallelizeQ=true
                this_custom_intervals="${HC_CUSTOM_INTERVALS}"
                if [ ! -z ${HC_SCAFFOLDS+x} ]; then
                    this_scaffolds="${HC_SCAFFOLDS}"
                else
                    this_scaffolds="false"
                fi
            fi
            # Otherwise, if we are using GATK 4.0 and are either NOT providing intervals or
            # providing custom intervals STARTING at the Genomics_DB_Import step, run the following
            # If we have custom intervals set, check if we are parallelizing across regions
            if [[ "${this_custom_intervals}" != "false" ]]; then
                analysisType="${analysisType:-targeted}" # if $analysisType isn't set, this var becomes "targeted"
                intvlFile="${this_custom_intervals}"
                # If we are parallelizing across regions, we will have one gendb
                # workspace for every custom interval
                if [ "${parallelizeQ}" == "true" ]; then
                    echo "We will parallelize across the provided custom intervals (i.e., regions)."
                    # Use CUSTOM_INTERVALS to determine array limit
                    NUM_INTVLS=$(wc -l "${this_custom_intervals}" | cut -d ' ' -f 1)
                    SINGLE_ARRAY_LIMIT=$[${NUM_INTVLS} - 1]
                    # If we also have scaffolds in addition to custom intervals
                    # and we are parallelizing across regions, add one to the array limit
                    if [[ "${SCAFFOLDS}" != "false" ]]; then
                        SINGLE_ARRAY_LIMIT=${NUM_INTVLS}
                    fi
                else
                    # We are NOT parallelizing across regions
                    # List of custom intervals is given to GATK, but only a single
                    # gendb workspace is created.
                    echo "Custom intervals list provided, but we are NOT parallelizing across regions or chromosomes/scaffolds."
                    SINGLE_ARRAY_LIMIT="0"
                fi
            else
                analysisType="WGS"
                intvlFile="NA"
                if [ "${parallelizeQ}" == "true" ]; then
                    echo "Analysis type is: WGS. We will parallelize across chromosomes and scaffolds if available."
                    # Use the number of chromosomes to set array limit
                    SINGLE_ARRAY_LIMIT=$[${NUM_CHR} - 1] # Get the maximum number of Torque tasks (# in array - 1)
                    # If we also have scaffolds in addition to the chromosomes, add one to the array limit
                    if [[ "${SCAFFOLDS}" != "false" ]]; then
                        SINGLE_ARRAY_LIMIT=${NUM_CHR}
                    fi
                else
                    echo "Not parallelizing across regions or chromosomes/scaffolds."
                    SINGLE_ARRAY_LIMIT="0"
                fi
            fi
            # Prepare arrays
            echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
            # Check if we have enough samples for a task
            # Run handler
            if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                exit 1
            elif [[ "${USE_PBS}" == true ]]; then
                echo "PBS is our workload manager/job scheduler."
                if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                    # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                    if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                    then
                        echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                        job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                    else
                        # Use job arrays following -t flag
                        echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                        job_array_index_range="${CUSTOM_JOB_ARR}"
                    fi
                else # If we are creating a single gendb workspace
                    echo "Using a single job"
                    job_array_index_range="0"
                fi
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Genomics_DB_Import.sh && GenomicsDBImport ${GVCF_LIST} ${OUT_DIR} ${REF_GEN} ${analysisType} ${intvlFile} ${this_scaffolds} ${GDBI_MEM} ${parallelizeQ} ${TMP}" | qsub -t "${job_array_index_range}" -l "${GDBI_QSUB}" -q "${GDBI_QUEUE}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_DBImport
            elif [[ "${USE_SLURM}" == true ]]; then
                echo "Slurm is our workload manager/job scheduler."
                if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                    # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                    if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                    then
                        echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                        job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                    else
                        # Use job arrays following -t flag
                        echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                        job_array_index_range="${CUSTOM_JOB_ARR}"
                    fi
                else # If we are creating a single gendb workspace
                    echo "Using a single job"
                    job_array_index_range="0"
                fi
                sbatch --job-name=${PROJECT}_DBImport --array=${job_array_index_range} ${GDBI_SBATCH} -o ${ERROR}/${PROJECT}_DBImport.%A_%a.out -e ${ERROR}/${PROJECT}_DBImport.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Genomics_DB_Import_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING} ${GDBI_MEM} ${analysisType} ${intvlFile} ${this_scaffolds} ${parallelizeQ}
            else # NO PBS
                (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Genotype_GVCFs.sh" && GenomicsDBImport "${GVCF_LIST}" "${OUT_DIR}" "${REF_GEN}" ${analysisType} "${intvlFile}" "${this_scaffolds}" "${GDBI_MEM}" "${parallelizeQ}" "${TMP}") >> "${ERROR}/Genomics_DB_Import.log" 2>&1
            fi
            echo "Making database. Please run Genotype_GVCFs after the database is complete."
            exit 1
        elif checkVersion gatk 3.8.0 ; then
            gatkVer=3
            echo "Don't need to make database for GATK version 3.8.0 or earlier. Please proceed with Genotype_GVCFs."
        else
            echo "Please install GATK 3.8.0 or newer" >&2
            exit 1
        fi
        ;;
    8 | Genotype_GVCFs )
        echo "$(basename $0): Genotyping GVCFs..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Genotype_GVCFs.sh
        checkDependencies Genotype_GVCFs_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${GVCF_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        # if GATK_JAR is not set in Config, it will check GATK_LOCAL_JAR
        GATK_JAR=$(checkGATK "${GATK_JAR}") # Check to make sure GATK is installed
        if [[ "${GATK_JAR}" == 1 ]]; then exit 1; fi # If we don't have GATK, exit with error
	    rm -f "${ERROR}/Genotype_GVCFs.log"; touch "${ERROR}/Genotype_GVCFs.log"
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            GG_MEM=$(getMemory "${GG_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            GG_MEM=$(getMemory "${GG_SBATCH}" "slurm")
        fi
        # GenomicsDBImport() for GATK4 and parallelized WGS requires reference dict
        checkDict "${REF_GEN}" "${PICARD_JAR}" # Check to make sure our reference genome has a dict file
        # Check if THETA variable is defined in config since this is taken from the Haplotype_Caller section of the config
        if [ -z ${THETA+x} ]; then
            echo "THETA variable in the config is unset. Please fill out THETA under the Haplotype_Caller section of the config. Exiting..."
            exit 1
        fi

        if checkVersion gatk 4.0 ; then   # check gatk version
            gatkVer=4
            # With gatk4, individual vcf files have to be combined into a DB
            # Check if DB exists, if so, proceed with running Genotype_GVCFs
            if [ -d "${OUT_DIR}"/Genotype_GVCFs/combinedDB ] && [ $(ls "${OUT_DIR}/Genotype_GVCFs/combinedDB" | wc -l) -gt 0 ]; then
                echo "INFO: *** Using pre-existing ${OUT_DIR}/Genotype_GVCFs/combinedDB ***" >> "${ERROR}/Genotype_GVCFs.log"
            else
                echo "Please run Genomics_DB_Import first to make the database. Exiting..."
                exit 1
            fi
            # Note the gendb:// prefix to the database input directory path is needed for GATK4
            # As of Sep 23, 2019 with 4.1.2 it seems like we need to be in Genotype_GVCFs for GATK4 to find database
            # CL is unable to get it working with a relative or absolute filepath to the database
            # Naoki: June 4, 2020 relative path with gendb:// works with 4.1.7, but
            # The following variable: input_gvcf is ignored
            # in Handlers/GenotypeGVCFs.sh, so I'm putting "DUMMY" here.
            input_gvcf="DUMMY"
        elif checkVersion gatk 3.8.0 ; then
                gatkVer=3
                input_gvcf="${GVCF_LIST}"
        else
            echo "Please install GATK 3.8.0 or newer" >&2
            exit
        fi

        SINGLE_ARRAY_LIMIT=$[${NUM_CHR} - 1] # Get the maximum number of Torque tasks (# in array - 1)

        parallelizeQ="${PARALLELIZE}"
        scaffoldsFile="${SCAFFOLDS}"
        # If we are parallelizing across regions, we will have one gendb
            # workspace for every custom interval
        # With the parallelization at Haplotype_Caller stage,
        # we need to use HC_CUSTOM_INTERVALS instead of CUSTOM_INTERVALS
        if [ ! -z ${HC_PARALLELIZE+x} ] && [ ${HC_PARALLELIZE} == "true" ] && [ ! -z "${HC_CUSTOM_INTERVALS}" ] && [ "${HC_CUSTOM_INTERVALS}" != "false" ]; then
            analysisType="targeted"
            intvlFile="${HC_CUSTOM_INTERVALS}"
            parallelizeQ="true"
            if [ ! -z ${HC_SCAFFOLDS+x} ] && [ "${HC_SCAFFOLDS}" != "false" ]; then
                scaffoldsFile="${HC_SCAFFOLDS}"
            fi
        elif [[ "${CUSTOM_INTERVALS}" != "false" ]] ; then
            analysisType="targeted"
            intvlFile="${CUSTOM_INTERVALS}"
        else
            analysisType="WGS"
            intvlFile="NA"
        fi

        # set SINGLE_ARRAY_LIMIT
        if [[ "${parallelizeQ}" == "true" ]]; then
            if [[ "${analysisType}" == "targeted" ]]; then
                # Use CUSTOM_INTERVALS to determine array limit
                NUM_INTVLS=$(wc -l "${intvlFile}" | cut -d ' ' -f 1)
                SINGLE_ARRAY_LIMIT=$[${NUM_INTVLS} - 1]
            else  # WGS
                # Use the number of chromosomes to set array limit
                SINGLE_ARRAY_LIMIT=$[${NUM_CHR} - 1] # Get the maximum number of Torque tasks (# in array - 1)
            fi
            # If we also have scaffolds in addition to the chromosomes, add one to the array limit
            if [[ "${scaffoldsFile}" != "false" ]]; then
                SINGLE_ARRAY_LIMIT=${NUM_INTVLS}
            fi
        else
            # We are NOT parallelizing across regions
            # List of custom intervals is given to GATK, but only a single
            # gendb workspace is created.
            # If we have scaffolds, add one to the maximum array index
            SINGLE_ARRAY_LIMIT="0"
        fi

        if [ "${SINGLE_ARRAY_LIMIT}" -lt 0 ]; then
            echo "Please set NUM_CHR > 0. Or if NUM_CHR=0, set CUSTOM_INTERVALS" >&2
            exit 1
        fi

        # Prepare job arrays
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        # Run the Handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
            then # If we have enough samples for a task array
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we have only one sample
                job_array_index_range="0"
            fi
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Genotype_GVCFs.sh && Genotype_GVCFs ${input_gvcf} ${OUT_DIR} ${GATK_JAR} ${REF_GEN} ${THETA} ${PLOIDY} ${GG_MEM} ${REF_DICT} ${SINGLE_ARRAY_LIMIT} ${gatkVer} ${intvlFile} ${parallelizeQ} ${analysisType} ${scaffoldsFile} ${GG_COMBINED_VCF} ${PROJECT} ${TMP}" | qsub -t "${job_array_index_range}" -q "${GG_QUEUE}" -l "${GG_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Genotype_GVCFs
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
            then # If we have enough samples for a task array
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we have only one sample
                job_array_index_range="0"
            fi
            sbatch --job-name=${PROJECT}_Genotype_GVCFs --array=${job_array_index_range} ${GG_SBATCH} -o ${ERROR}/${PROJECT}_Genotype_GVCFs.%A_%a.out -e ${ERROR}/${PROJECT}_Genotype_GVCFs.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Genotype_GVCFs.job ${CONFIG} ${SEQUENCE_HANDLING} ${GG_MEM} ${input_gvcf} ${SINGLE_ARRAY_LIMIT} ${gatkVer} ${intvlFile} ${parallelizeQ} ${analysisType} ${scaffoldsFile}
        else # Without PBS
            source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Genotype_GVCFs.sh && Genotype_GVCFs "${input_gvcf}" "${OUT_DIR}" "${GATK_JAR}" "${REF_GEN}" ${THETA} ${PLOIDY} ${GG_MEM} "${REF_DICT}" ${SINGLE_ARRAY_LIMIT} ${gatkVer} "${intvlFile}" ${parallelizeQ} ${analysisType} "${scaffoldsFile} ${GG_COMBINED_VCF} ${PROJECT} ${TMP}" >> "${ERROR}/Genotype_GVCFs.log" 2>&1
            # Combine and sort split vcfs
            if [[ "${parallelize}" == true && "${gg_combined_vcf}" == true ]]; then
                out_subdir="${out_dir}/Genotype_GVCFs/vcf_split_regions"
                ls ${out_subdir}/*.vcf > ${out_subdir}/temp-FileList.list # note sufix has to be .list
                gatk --java-options "-Xmx${memory}" SortVcf \
                     -I ${out_subdir}/temp-FileList.list \
                     -O ${out_dir}/Genotype_GVCFs/${project}_raw_variants.vcf >> "${ERROR}/Genotype_GVCFs.log" 2>&1
                rm -f ${out_subdir}/temp-FileList.list # Cleanup
            fi
        fi
        # CL Note: Combining split vcfs after the process above only works for non-PBS.
        # On PBS, different regions/parts will finish running at different times due to how the
        # job scheduler is set up to run task arrays.
        # For those using PBS AND don't plan to run the Create_HC_Subset handler and GATK's variant recalibrator,
        # the script HelperScripts/combine_and_sort_split_vcf.sh will combine and sort the split VCF files.
        # This is only necessary if user is NOT running Create_HC_Subset.
        ;;
    9 | Create_HC_Subset )
        echo "$(basename $0): Creating a high-confidence subset of variants..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Create_HC_Subset.sh
        checkDependencies Create_HC_Subset_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        if [[ ${CHS_VCF_LIST} != "NA" ]]; then
            # Check sample list
            checkSamples "${CHS_VCF_LIST}" # Check to see if samples and sample list exists
        else
            # Check single VCF file
            if ! [[ -f "${CHS_RAW_VCF}" ]]; then # If the sample doesn't exist
                echo "The sample ${CHS_RAW_VCF} does not exist, exiting..." >&2 # Exit out with error
                return 1
            else
                if ! [[ -r "${CHS_RAW_VCF}" ]]; then # If the sample isn't readable
                    echo "The sample ${CHS_RAW_VCF} does not have read permissions, exiting..." >&2
                    return 1
                fi
            fi
        fi
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            CHS_MEM=$(getMemory "${CHS_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            CHS_MEM=$(getMemory "${CHS_SBATCH}" "slurm")
        fi
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        if [[ -z "${BARLEY}" ]]; then echo "Please specify whether or not the organism is barley in the config file" >&2; exit 1; fi # Make sure barley is filled out
        if [[ -z "${CAPTURE_REGIONS}" ]]; then echo "Please either specify the capture regions bed file or put NA for CAPTURE_REGIONS" >&2; exit 1; fi # Make sure CAPTURE_REGIONS is filled out
        # NOTE: Currently (Aug 17, 2020), we need both CHS_RAW_VCF and CHS_VCF_LIST variables in the Config file to make this handler backwards compatible with GATK 3 and the Create_HC_Subset_GATK3 function
        # Check if both CHS_RAW_VCF and CHS_VCF_LIST are filled out appropriately
        if test -z "${CHS_RAW_VCF}"; then
            echo "CHS_RAW_VCF variable is empty, please fill out with either NA or filepath to raw variants."
            exit 1
        fi
        if test -z "${CHS_VCF_LIST}"; then
            echo "CHS_VCF_LIST variable is empty, please fill out with either NA or filepath to list of split vcf files."
            exit 1
        fi
        # Check if we are using a job scheduler or not
	    if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            # Check if GATK 3 or GATK 4
            if checkVersion gatk 4.0; then
                # Run GATK 4 code
                echo "Using GATK 4..."
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Create_HC_Subset.sh && Create_HC_Subset_GATK4 ${CHS_RAW_VCF} ${CHS_VCF_LIST} ${OUT_DIR} ${BARLEY} ${PROJECT} ${SEQUENCE_HANDLING} ${CHS_QUAL_CUTOFF} ${CHS_GQ_CUTOFF} ${CHS_MAX_LOWGQ} ${CHS_DP_PER_SAMPLE_CUTOFF} ${CHS_MAX_HET} ${CHS_MAX_MISS} ${CHS_MEM} ${REF_GEN} ${CHS_SUBSET_N} ${CHS_SUBSET_LEN}" | qsub -q "${CHS_QUEUE}" -l "${CHS_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Create_HC_Subset
                # Add new line to make it easier to see job ID
                printf "\nIMPORTANT REMINDER: To save processing time, this handler checks if the ${OUT_DIR}/Create_HC_Subset/${PROJECT}_raw_variants.vcf.gz file already exist and skips the first step if the file exists. If you know your file is truncated (due to running out of walltime, etc.), please delete the file(s) that are incomplete and re-run this handler. If everything looks ok, you can ignore this message."
                echo "The same checks apply for each step within the ${SEQUENCE_HANDLING}/HelperScripts/percentiles.sh script. If you suspect any files are truncated (due to exceeding walltime, etc.), please delete relevant files and re-run this handler."
            else
                # Assume we are using GATK 3
                echo "Using GATK 3..."
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Create_HC_Subset.sh && Create_HC_Subset_GATK3 ${CHS_VCF_LIST} ${OUT_DIR} ${CAPTURE_REGIONS} ${BARLEY} ${PROJECT} ${SEQUENCE_HANDLING} ${CHS_QUAL_CUTOFF} ${CHS_GQ_CUTOFF} ${CHS_DP_PER_SAMPLE_CUTOFF} ${CHS_MAX_HET} ${CHS_MAX_MISS} ${TMP} ${REF_GEN} ${CHS_SUBSET_N} ${CHS_SUBSET_LEN}" | qsub -q "${CHS_QUEUE}" -l "${CHS_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Create_HC_Subset
            fi
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            # Check if GATK 3 or GATK 4
            if checkVersion gatk 4.0; then
                # Run GATK 4 code
                echo "Using GATK 4..."
                sbatch --job-name=${PROJECT}_Create_HC_Subset ${CHS_SBATCH} -o ${ERROR}/${PROJECT}_Create_HC_Subset.%j.out -e ${ERROR}/${PROJECT}_Create_HC_Subset.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Create_HC_Subset_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING} ${CHS_MEM}
                # Add new line to make it easier to see job ID
                printf "\nIMPORTANT REMINDER: To save processing time, this handler checks if the ${OUT_DIR}/Create_HC_Subset/${PROJECT}_raw_variants.vcf.gz file already exist and skips the first step if the file exists. If you know your file is truncated (due to running out of walltime, etc.), please delete the file(s) that are incomplete and re-run this handler. If everything looks ok, you can ignore this message."
                echo "The same checks apply for each step within the ${SEQUENCE_HANDLING}/HelperScripts/percentiles.sh script. If you suspect any files are truncated (due to exceeding walltime, etc.), please delete relevant files and re-run this handler."
            else
                # Assume we are using GATK 3
                echo "Using GATK 3..."
                sbatch --job-name=${PROJECT}_Create_HC_Subset ${CHS_SBATCH} -o ${ERROR}/${PROJECT}_Create_HC_Subset.%j.out -e ${ERROR}/${PROJECT}_Create_HC_Subset.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Create_HC_Subset_GATK3.job ${CONFIG} ${SEQUENCE_HANDLING}
            fi
	    else
            # Non-PBS
            if checkVersion gatk 4.0; then
                # Run GATK 4 code
                (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Create_HC_Subset.sh && Create_HC_Subset_GATK4 ${CHS_RAW_VCF} ${CHS_VCF_LIST} ${OUT_DIR} ${BARLEY} ${PROJECT} ${SEQUENCE_HANDLING} ${CHS_QUAL_CUTOFF} ${CHS_GQ_CUTOFF} ${CHS_MAX_LOWGQ} ${CHS_DP_PER_SAMPLE_CUTOFF} ${CHS_MAX_HET} ${CHS_MAX_MISS} ${CHS_MEM} ${REF_GEN} ${CHS_SUBSET_N} ${CHS_SUBSET_LEN}) > "${ERROR}/Create_HC_Subset.log" 2>&1
                echo "IMPORTANT REMINDER: To save processing time, this handler checks if the ${OUT_DIR}/Create_HC_Subset/${PROJECT}_raw_variants.vcf.gz file and filtered ${OUT_DIR}/Create_HC_Subset/Intermediates/${PROJECT}_no_indels.recode.vcf file already exist and skips the first two steps if the files exist. If you know your file is truncated (due to running out of walltime, etc.), please delete the file(s) that are incomplete and re-run this handler. If everything looks ok, you can ignore this message."
                echo "The same checks apply for each step within the ${SEQUENCE_HANDLING}/HelperScripts/percentiles.sh script. If you suspect any files are truncated (due to exceeding walltime, etc.), please delete relevant files and re-run this handler."
            else
	            (source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Create_HC_Subset.sh && Create_HC_Subset_GATK3 ${CHS_VCF_LIST} ${OUT_DIR} ${CAPTURE_REGIONS} ${BARLEY} ${PROJECT} ${SEQUENCE_HANDLING} ${CHS_QUAL_CUTOFF} ${CHS_GQ_CUTOFF} ${CHS_DP_PER_SAMPLE_CUTOFF} ${CHS_MAX_HET} ${CHS_MAX_MISS} ${TMP} ${REF_GEN} ${CHS_SUBSET_N} ${CHS_SUBSET_LEN}) > "${ERROR}/Create_HC_Subset.log" 2>&1
            fi
	    fi
        ;;
    10 | Variant_Recalibrator )
        echo "$(basename $0): Training model and recalibrating quality scores of variants..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Variant_Recalibrator.sh
        checkDependencies Variant_Recalibrator_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        if [ ${VR_VCF_LIST} == "NA" ]; then
            if ! [[ -f "${VR_RAW_VCF}" ]]; then # if the sample doesn't exist
                echo "The file ${VR_RAW_VCF} does not exist, exiting..." >&2
                # Exit out with error
                return 1
            fi
        else
            checkSamples "${VR_VCF_LIST}" # Check to see if samples and sample list exists
        fi
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkGATK "${GATK_JAR}" # Check to make sure GATK is installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we don't have GATK, exit with error
        if [[ -z "${BARLEY}" ]]; then echo "Please specify whether or not the organism is barley in the config file" >&2; exit 1; fi # Make sure barley is filled out
        checkVCF "${HC_SUBSET}" # Make sure the VCF is formatted properly
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If it's not formatted properly, exit with error
        # Check if resource variables are specified correctly in config
        # Currently, if resource variable is not used, user should put "NA" and not leave it blank
        if [[ -z "${RESOURCE_1}" ]]; then echo "Variable RESOURCE_1 in Config is empty, please put NA if you are not planning to use that variable and re-run this handler. Exiting..."; exit 1; fi
        if [[ -z "${RESOURCE_2}" ]]; then echo "Variable RESOURCE_2 in Config is empty, please put NA if you are not planning to use that variable and re-run this handler. Exiting..."; exit 1; fi
        if [[ -z "${RESOURCE_3}" ]]; then echo "Variable RESOURCE_3 in Config is empty, please put NA if you are not planning to use that variable and re-run this handler. Exiting..."; exit 1; fi
        if [[ -z "${RESOURCE_4}" ]]; then echo "Variable RESOURCE_4 in Config is empty, please put NA if you are not planning to use that variable and re-run this handler. Exiting..."; exit 1; fi
        # Check VCF formatting of resource files
        if ! [[ "${RESOURCE_1}" == "NA" ]]; then checkVCF "${RESOURCE_1}"; fi # Make sure the VCF is formatted properly
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If it's not formatted properly, exit with error
        if ! [[ "${RESOURCE_2}" == "NA" ]]; then checkVCF "${RESOURCE_2}"; fi # Make sure the VCF is formatted properly
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If it's not formatted properly, exit with error
        if ! [[ "${RESOURCE_3}" == "NA" ]]; then checkVCF "${RESOURCE_3}"; fi # Make sure the VCF is formatted properly
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If it's not formatted properly, exit with error
        if ! [[ "${RESOURCE_4}" == "NA" ]]; then checkVCF "${RESOURCE_4}"; fi # Make sure the VCF is formatted properly
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If it's not formatted properly, exit with error
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            VR_MEM=$(getMemory "${VR_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            VR_MEM=$(getMemory "${VR_SBATCH}" "slurm")
        fi
        checkDict "${VR_REF}" ${PICARD_JAR} # Check to make sure our reference genome has a dict file
        if checkVersion gatk 4.0; then
            GATK_VERSION="4"
            # Check if we have a valid RECAL_MODE specified, error out if it is not one of the valid values
            if [ "${RECAL_MODE}" == "BOTH" ] || [ "${RECAL_MODE}" == "INDELS_ONLY" ] || [ "${RECAL_MODE}" == "SNPS_ONLY" ]; then
                printf "\nRECAL_MODE value is valid and is currently set to: ${RECAL_MODE}\n\n"
                # Generate expected output filename for file checks
                if [ "${VR_VCF_LIST}" == "NA" ]; then
                    # Use raw concatenated VCF filename
                    if [[ ${VR_RAW_VCF} == *.vcf.gz ]]; then
                        temp_vcf_filename=$(basename ${VR_RAW_VCF} .vcf.gz)
                    else
                        # Asssume vcf files ends in .vcf extension
                        temp_vcf_filename=$(basename ${VR_RAW_VCF} .vcf)
                    fi
                else
                    # Generate expected filename after concatenating split vcf files
                    temp_vcf_filename=$(basename "${OUT_DIR}/Variant_Recalibrator/${PROJECT}_raw_variants.vcf.gz" .vcf)
                fi
                # Check if we have already separated indels or snps only using GATK SelectVariants
                # For very large raw VCF files (>1 TB), this can be a time consuming step (>12 hours)
                # so, this comes in handy if we ran out of walltime but the select variants step completed successfully
                if [ "${RECAL_MODE}" == "INDELS_ONLY" ]; then
                    if [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/${temp_vcf_filename}_indels.vcf 2>/dev/null)" ]; then
                        echo "Indels only VCF file exists: ${OUT_DIR}/Variant_Recalibrator/${temp_vcf_filename}_indels.vcf"
                        printf "Proceeding with this file. If you suspect your vcf file is truncated (e.g., due to walltime issues, etc.), please delete this file and re-run the handler.\n\n"
                    else
                        printf "Indels only VCF file doesn't exist, we will pull out the indels first before recalibrating.\n\n"
                    fi
                elif [ "${RECAL_MODE}" == "SNPS_ONLY" ]; then
                    if [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/${temp_vcf_filename}_snps.vcf 2>/dev/null)" ]; then
                        echo "SNPs only VCF file exists: ${OUT_DIR}/Variant_Recalibrator/${temp_vcf_filename}_snps.vcf"
                        printf "Proceeding with this file. If you suspect your vcf file is truncated (e.g., due to walltime issues, etc.), please delete this file and re-run the handler.\n\n"
                    else
                        printf "SNPs only VCF file doesn't exist, we will pull out the SNPs first before recalibrating.\n\n"
                    fi
                else
                    printf "Recalibrating both indels and SNPs.\n\n"
                fi
            else
                printf "\nInvalid RECAL_MODE: ${RECAL_MODE}. Please set to one of the following - BOTH, INDELS_ONLY, SNPS_ONLY"
                exit 1 # Exiting out of program
            fi
            # Check if we are running on MSI or locally (Note: running locally will be a feature added in the future)
            if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                exit 1
            elif [[ "${USE_PBS}" == true ]]; then
                echo "PBS is our workload manager/job scheduler."
                # Run GATK 4 code
                echo "GATK Version ${GATK_VERSION}"
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Variant_Recalibrator.sh && Variant_Recalibrator_GATK4 ${VR_RAW_VCF} ${VR_VCF_LIST} ${OUT_DIR} ${GATK_JAR} ${VR_REF} ${VR_MEM} ${PROJECT} ${SEQUENCE_HANDLING} ${HC_SUBSET} ${RESOURCE_1} ${RESOURCE_2} ${RESOURCE_3} ${RESOURCE_4} ${HC_PRIOR} ${PRIOR_1} ${PRIOR_2} ${PRIOR_3} ${PRIOR_4} ${HC_KNOWN} ${KNOWN_1} ${KNOWN_2} ${KNOWN_3} ${KNOWN_4} ${HC_TRAIN} ${TRAINING_1} ${TRAINING_2} ${TRAINING_3} ${TRAINING_4} ${HC_TRUTH} ${TRUTH_1} ${TRUTH_2} ${TRUTH_3} ${TRUTH_4} ${BARLEY} ${GATK_VERSION} ${TS_FILTER_LEVEL} ${RECAL_MODE} ${TMP}" | qsub -q "${VR_QUEUE}" -l "${VR_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Variant_Recalibrator
                # Add new line to make it easier to see job ID
                if [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_sitesonly.vcf.gz 2>/dev/null)" ] || [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_indels_sitesonly.vcf.gz 2>/dev/null)" ] || [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_snps_sitesonly.vcf.gz 2>/dev/null)" ]; then
                    printf "\nIMPORTANT REMINDER: \n"
                    echo "To save processing time, this handler checks if the" ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}*_sitesonly.vcf.gz "file already exist and skips the MakeSitesOnlyVcf step if the file exists. If you know your file is truncated (due to running out of walltime, etc.), please delete the file(s) that are incomplete and re-run this handler. If everything looks ok, you can ignore this message."
                fi
            elif [[ "${USE_SLURM}" == true ]]; then
                echo "Slurm is our workload manager/job scheduler."
                # Run GATK 4 code
                echo "GATK Version ${GATK_VERSION}"
                source ${CONFIG} # to get sbatch settings loaded
                echo "Error files are located in directory: ${ERROR}"
                sbatch --job-name=${PROJECT}_Variant_Recalibrator ${VR_SBATCH} -o ${ERROR}/${PROJECT}_Variant_Recalibrator.%j.out -e ${ERROR}/${PROJECT}_Variant_Recalibrator.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Variant_Recalibrator_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING} ${VR_MEM} ${GATK_VERSION}
                # Add new line to make it easier to see job ID
                if [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_sitesonly.vcf.gz 2>/dev/null)" ] || [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_indels_sitesonly.vcf.gz 2>/dev/null)" ] || [ -n "$(ls -A ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}_snps_sitesonly.vcf.gz 2>/dev/null)" ]; then
                    printf "\nIMPORTANT REMINDER: \n"
                    echo "To save processing time, this handler checks if the" ${OUT_DIR}/Variant_Recalibrator/Intermediates/${PROJECT}*_sitesonly.vcf.gz "file already exist and skips the MakeSitesOnlyVcf step if the file exists. If you know your file is truncated (due to running out of walltime, etc.), please delete the file(s) that are incomplete and re-run this handler. If everything looks ok, you can ignore this message."
                fi
            fi
        else
            # Assume we are running GATK 3
            GATK_VERSION="3"
            if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                exit 1
            elif [[ "${USE_PBS}" == true ]]; then
                echo "PBS is our workload manager/job scheduler."
                echo "GATK Version ${GATK_VERSION}"
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Variant_Recalibrator.sh && Variant_Recalibrator_GATK3 {VR_RAW_VCF} ${VR_VCF_LIST} ${OUT_DIR} ${GATK_JAR} ${VR_REF} ${VR_MEM} ${PROJECT} ${SEQUENCE_HANDLING} ${HC_SUBSET} ${RESOURCE_1} ${RESOURCE_2} ${RESOURCE_3} ${RESOURCE_4} ${HC_PRIOR} ${PRIOR_1} ${PRIOR_2} ${PRIOR_3} ${PRIOR_4} ${HC_KNOWN} ${KNOWN_1} ${KNOWN_2} ${KNOWN_3} ${KNOWN_4} ${HC_TRAIN} ${TRAINING_1} ${TRAINING_2} ${TRAINING_3} ${TRAINING_4} ${HC_TRUTH} ${TRUTH_1} ${TRUTH_2} ${TRUTH_3} ${TRUTH_4} ${BARLEY} ${GATK_VERSION} ${TS_FILTER_LEVEL} ${TMP}" | qsub -q "${VR_QUEUE}" -l "${VR_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Variant_Recalibrator
            elif [[ "${USE_SLURM}" == true ]]; then
                echo "Slurm is our workload manager/job scheduler."
                # Assume we are using GATK 3
                echo "GATK Version ${GATK_VERSION}"
                source ${CONFIG} # to get sbatch settings loaded
                echo "Error files are located in directory: ${ERROR}"
                sbatch --job-name=${PROJECT}_Variant_Recalibrator ${VR_SBATCH} -o ${ERROR}/${PROJECT}_Variant_Recalibrator.%j.out -e ${ERROR}/${PROJECT}_Variant_Recalibrator.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Variant_Recalibrator_GATK3.job ${CONFIG} ${SEQUENCE_HANDLING} ${VR_MEM} ${GATK_VERSION}
            fi
        fi
        ;;
    11 | Pre_Variant_Filtering )
        echo "$(basename $0): Preparing to filter variants, generating annotation graphs and percentile tables..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Pre_Variant_Filtering.sh
        checkDependencies Pre_Variant_Filtering_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Pre_Variant_Filtering.sh && Pre_Variant_Filtering_GATK4 ${PVF_VCF} ${OUT_DIR} ${REF_GEN} ${PROJECT} ${SEQUENCE_HANDLING} ${PVF_GEN_NUM} ${PVF_GEN_LEN} ${TMP}" | qsub -q "${PVF_QUEUE}" -l "${PVF_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Pre_Variant_Filtering
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            echo "Error files are located in directory: ${ERROR}"
            sbatch --job-name=${PROJECT}_Pre_Variant_Filtering ${PVF_SBATCH} -o ${ERROR}/${PROJECT}_Pre_Variant_Filtering.%j.out -e ${ERROR}/${PROJECT}_Pre_Variant_Filtering.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Pre_Variant_Filtering_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING}
        fi
        ;;
    12 | Variant_Filtering )
        echo "$(basename $0): Filtering variants..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Variant_Filtering.sh
        checkDependencies Variant_Filtering_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        if [[ -z "${BARLEY}" ]]; then echo "Please specify whether or not the organism is barley in the config file" >&2; exit 1; fi # Make sure barley is filled out
	    if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "$USE_PBS" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            # Check if GATK 4 or GATK 3
            if checkVersion gatk 4.0; then
                GATK_VERSION="4"
                echo "GATK Version ${GATK_VERSION}"
                # Check if we have a TMP directory set
                if [ -z "${TMP}" ]; then
                    echo "TMP variable is empty in the config, please set the TMP directory before running this handler. Exiting..."
                    exit 1
                fi
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Variant_Filtering.sh && Variant_Filtering_GATK4 ${VF_VCF} ${OUT_DIR} ${REF_GEN} ${PROJECT} ${SEQUENCE_HANDLING} ${VF_MIN_DP} ${VF_MAX_DP} ${VF_MAX_DEV} ${VF_DP_PER_SAMPLE_CUTOFF} ${VF_GQ_CUTOFF} ${VF_MAX_HET} ${VF_MAX_BAD} ${VF_QUAL_CUTOFF} ${TMP}" | qsub -q "${VF_QUEUE}" -l "${VF_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Variant_Filtering
            else
                GATK_VERSION="3"
                echo "GATK Version ${GATK_VERSION}"
                if [[ -z "${VF_CAPTURE_REGIONS}" ]]; then echo "Please either specify the capture regions bed file or put NA for VF_CAPTURE_REGIONS" >&2; exit 1; fi # Make sure CAPTURE_REGIONS is filled out
                echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Variant_Filtering.sh && Variant_Filtering_GATK3 ${VF_VCF} ${OUT_DIR} ${VF_CAPTURE_REGIONS} ${BARLEY} ${PROJECT} ${SEQUENCE_HANDLING} ${VF_MIN_DP} ${VF_MAX_DP} ${VF_MAX_DEV} ${VF_DP_PER_SAMPLE_CUTOFF} ${VF_GQ_CUTOFF} ${VF_MAX_HET} ${VF_MAX_BAD} ${VF_QUAL_CUTOFF}" | qsub -q "${VF_QUEUE}" -l "${VF_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Variant_Filtering
            fi
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            # Check if GATK 4 or GATK 3
            if checkVersion gatk 4.0; then
                GATK_VERSION="4"
                echo "GATK Version ${GATK_VERSION}"
                # Check if we have a TMP directory set
                if [ -z "${TMP}" ]; then
                    echo "TMP variable is empty in the config, please set the TMP directory before running this handler. Exiting..."
                    exit 1
                fi
                echo "Error files are located in directory: ${ERROR}"
                sbatch --job-name=${PROJECT}_Variant_Filtering ${VF_SBATCH} -o ${ERROR}/${PROJECT}_Variant_Filtering.%j.out -e ${ERROR}/${PROJECT}_Variant_Filtering.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Variant_Filtering_GATK4.job ${CONFIG} ${SEQUENCE_HANDLING}
            else
                GATK_VERSION="3"
                echo "GATK Version ${GATK_VERSION}"
                if [[ -z "${VF_CAPTURE_REGIONS}" ]]; then echo "Please either specify the capture regions bed file or put NA for VF_CAPTURE_REGIONS" >&2; exit 1; fi # Make sure CAPTURE_REGIONS is filled out
                echo "Error files are located in directory: ${ERROR}"
                sbatch --job-name=${PROJECT}_Variant_Filtering ${VF_SBATCH} -o ${ERROR}/${PROJECT}_Variant_Filtering.%j.out -e ${ERROR}/${PROJECT}_Variant_Filtering.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Variant_Filtering_GATK3.job ${CONFIG} ${SEQUENCE_HANDLING}
            fi
	    else
	        (source ${CONFIG} && source "${SEQUENCE_HANDLING}/Handlers/Variant_Filtering.sh" && Variant_Filtering "${VF_VCF}" "${OUT_DIR}" "${VF_CAPTURE_REGIONS}" ${BARLEY} ${PROJECT} "${SEQUENCE_HANDLING}" ${MIN_DP} ${MAX_DP} ${MAX_DEV} ${DP_PER_SAMPLE_CUTOFF} ${GQ_CUTOFF} ${MAX_HET} ${MAX_BAD} ${QUAL_CUTOFF}) > "${ERROR}/Variant_Filtering.log" 2>&1
	    fi
        ;;
    13 | Variant_Analysis )
        echo "$(basename $0): Generating summary statistics..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Variant_Analysis.sh
        checkDependencies Variant_Analysis_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        if [[ -z "${VA_VCF}" ]]; then echo "Please specify the VCF file to be analyzed in the config file" >&2; exit 1; fi # Make sure the VCF file is filled out
        if [[ -z "${BARLEY}" ]]; then echo "Please specify whether or not the organism is barley in the config file" >&2; exit 1; fi # Make sure barley is filled out
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
	    elif [[ "$USE_PBS" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Variant_Analysis.sh && Variant_Analysis ${VA_VCF} ${OUT_DIR} ${SEQUENCE_HANDLING} ${BARLEY}" | qsub -l "${VA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Variant_Analysis
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            echo "Error files are located in directory: ${ERROR}"
            sbatch --job-name=${PROJECT}_Variant_Analysis ${VA_SBATCH} -o ${ERROR}/${PROJECT}_Variant_Analysis.%j.out -e ${ERROR}/${PROJECT}_Variant_Analysis.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Variant_Analysis.job ${CONFIG} ${SEQUENCE_HANDLING}
	    else
	        (source "${CONFIG}" && source "${SEQUENCE_HANDLING}/Handlers/Variant_Analysis.sh" && Variant_Analysis "${VA_VCF}" "${OUT_DIR}" "${SEQUENCE_HANDLING}" ${BARLEY}) > "${ERROR}/Variant_Analysis.log" 2>&1
	    fi
        ;;
    14 | GBS_Demultiplex )
        echo "$(basename $0): Splitting files based on barcodes..." >&2
        echo "GBS_Demultiplex is not yet functional, exiting..." >&2
        exit 1
        source "${SEQUENCE_HANDLING}"/Handlers/GBS_Demultiplex.sh
        checkDependencies GBS_Demultiplex_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${GBS_SAMPLES}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        barcodeGenerator "${SEQUENCE_HANDLING}" "${KEY_FILE}" "${OUT_DIR}" "${PROJECT}"
        declare -a BARCODE_LIST=($(grep -E ".barcode" "${OUT_DIR}/GBS_Demultiplex/barcodes/${PROJECT}_barcode_list.txt"))
        #   Run GBS_Demultiplexer
        SINGLE_ARRAY_LIMIT=$[${#BARCODE_LIST[@]} - 1] # Get the maximum number of Torque tasks we're doing
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        echo -e "#!/bin/bash\n#PBS -l ${GD_QSUB}\n#PBS -e ${ERROR}\n#PBS -o ${ERROR}\n#PBS -m abe\n#PBS -M ${EMAIL}\nset -e\nset -o pipefail\nsource ${CONFIG}\nsource ${SEQUENCE_HANDLING}/Handlers/GBS_Demultiplex.sh\ndeclare -a BARCODE_LIST=($(grep -E ".barcode" "${OUT_DIR}/GBS_Demultiplex/barcodes/${PROJECT}_barcode_list.txt"))\nSINGLE_ARRAY_LIMIT=\$[\${#BARCODE_LIST[@]} - 1]\nGBS_Demultiplex \${BARCODE_LIST[\${PBS_ARRAYID}]} ${GBS_SAMPLES} ${OUT_DIR} ${LINE_ENDING} ${MISMATCH_TOL} ${PARTIAL} ${FILE_TYPE} ${PROJECT}" > ${PROJECT}_GBS_Demultiplex
        if [[ -z ${CUSTOM_JOB_ARR+x} ]]
        then
            # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
            echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
            qsub -t 0-"${SINGLE_ARRAY_LIMIT}" "${PROJECT}"_GBS_Demultiplex
        else
            # Use job arrays following -t flag
            echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
            qsub -t "${CUSTOM_JOB_ARR}" "${PROJECT}"_GBS_Demultiplex
        fi
        ;;
    15 | Quality_Trimming )
        echo "$(basename $0): Trimming based on quality..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Quality_Trimming.sh
        checkDependencies Quality_Trimming_Dependencies[@] # Check to see if dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${ADAPTED_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        # Run the handler
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "${USE_PBS}" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            #   Run Quality_Trimming
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Quality_Trimming.sh && Quality_Trimming ${ADAPTED_LIST} ${FORWARD_ADAPTED} ${REVERSE_ADAPTED} ${SINGLES_ADAPTED} ${OUT_DIR} ${QT_THRESHOLD} ${QUAL_ENCODING} ${SEQUENCE_HANDLING} ${PROJECT}"| qsub -q "${QT_QUEUE}" -l "${QT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Trimming
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            sbatch --job-name=${PROJECT}_Quality_Trimming ${QT_SBATCH} -o ${ERROR}/${PROJECT}_Quality_Trimming.%j.out -e ${ERROR}/${PROJECT}_Quality_Trimming.%j.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Quality_Trimming.job ${CONFIG} ${SEQUENCE_HANDLING}
        fi
        ;;
    16 | Realigner_Target_Creator )
        echo "$(basename $0): Creating targets files for realignment..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Realigner_Target_Creator.sh
        checkDependencies Realigner_Target_Creator_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${RTC_BAM_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkBaiIndex "${RTC_BAM_LIST}" # Check to see if our samples are indexed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If they're not indexed, exit out with error
        GATK_JAR=$(checkGATK "${GATK_JAR}") # Check to make sure GATK is installed
        if [[ "${GATK_JAR}" == 1 ]]; then exit 1; fi # If we don't have GATK, exit with error
        if checkVersion gatk 4.0 ; then   # check gatk version
	        echo "Indel realignment functionality is no longer available in GATK 4. Please use GATK 3." >&2
	    elif checkVersion gatk 3.8.0; then
	        gatkVer=3
	    else
	        echo "Please install GATK 3.8.0 or newer" >&2
	        exit 1
	    fi
        RTC_MEM=$(getMemory "${RTC_QSUB}") # Figure out memory requirements based on the Qsub settings
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            RTC_MEM=$(getMemory "${RTC_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            RTC_MEM=$(getMemory "${RTC_SBATCH}" "slurm")
        fi
        checkDict "${REF_GEN}" "${PICARD_JAR}" # Check to make sure our reference genome has a dict file
        declare -a RTC_LIST=($(grep -E ".bam" "${RTC_BAM_LIST}")) # Create an array of the BAM files
        SINGLE_ARRAY_LIMIT=$[${#RTC_LIST[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
                echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
                exit 1
        elif [[ "${USE_PBS}" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Realigner_Target_Creator.sh && Realigner_Target_Creator ${RTC_BAM_LIST} ${OUT_DIR} ${GATK_JAR} ${REF_GEN} ${RTC_MEM} ${FIX_QUALITY_SCORES} ${gatkVer}" | qsub -t "${job_array_index_range}" -q "${RTC_QUEUE}" -l "${RTC_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Realigner_Target_Creator
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            sbatch --job-name=${PROJECT}_Realigner_Target_Creator --array=${job_array_index_range} ${RTC_SBATCH} -o ${ERROR}/${PROJECT}_Realigner_Target_Creator.%A_%a.out -e ${ERROR}/${PROJECT}_Realigner_Target_Creator.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Realigner_Target_Creator.job ${CONFIG} ${SEQUENCE_HANDLING} ${RTC_MEM} ${gatkVer}
        fi
        ;;
    17 | Indel_Realigner )
        echo "$(basename $0): Realigning reads around insertions and deletions..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/Indel_Realigner.sh
        checkDependencies Indel_Realigner_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${IR_BAM_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkBaiIndex "${IR_BAM_LIST}" # Check to see if our samples are indexed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If they're not indexed, exit out with error
        GATK_JAR=$(checkGATK "${GATK_JAR}") # Check to make sure GATK is installed
        if [[ "${GATK_JAR}" == 1 ]]; then exit 1; fi # If we don't have GATK, exit with error
        if checkVersion gatk 4.0 ; then   # check gatk version
	        echo "Indel realignment functionality is no longer available in GATK 4. Please use GATK 3." >&2
	    elif checkVersion gatk 3.8.0; then
	        gatkVer=3
	    else
	        echo "Please install GATK 3.8.0 or newer" >&2
	        exit 1
	    fi
        # Figure out memory requirements based on the Qsub or Slurm settings
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ ${USE_PBS} == true ]]; then
            IR_MEM=$(getMemory "${IR_QSUB}" "pbs")
        elif [[ ${USE_SLURM} == true ]]; then
            IR_MEM=$(getMemory "${IR_SBATCH}" "slurm")
        fi
        checkDict "${REF_GEN}" "${PICARD_JAR}" # Check to make sure our reference genome has a dict file
        declare -a IR_LIST=($(grep -E ".bam" "${IR_BAM_LIST}")) # Create an array of the BAM files
        SINGLE_ARRAY_LIMIT=$[${#IR_LIST[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "$USE_PBS" == true ]] && [[ "$USE_SLURM" == true ]]; then
            echo "USE_PBS and USE_SLURM are both set to true in the Config. Only one of the two can be true, please set one of the two variables to false depending on the workload manager your cluster uses. Exiting..."
            exit 1
        elif [[ "${USE_PBS}" == true ]]; then
            echo "PBS is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/Indel_Realigner.sh && Indel_Realigner ${IR_BAM_LIST} ${OUT_DIR} ${GATK_JAR} ${REF_GEN} ${IR_MEM} ${IR_TARGETS} ${LOD_THRESHOLD} ${ENTROPY_THRESHOLD} ${FIX_QUALITY_SCORES} ${gatkVer} ${MAX_READS_IN_MEM}" | qsub -t "${job_array_index_range}" -q "${IR_QUEUE}" -l "${IR_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Indel_Realigner
        elif [[ "${USE_SLURM}" == true ]]; then
            echo "Slurm is our workload manager/job scheduler."
            if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]; then
                # If CUSTOM_JOB_ARR is NOT set (only gets set if user uses -t flag)
                if [[ -z ${CUSTOM_JOB_ARR+x} ]]
                then
                    echo "Using default array limit from 0 to ${SINGLE_ARRAY_LIMIT}"
                    job_array_index_range="0-${SINGLE_ARRAY_LIMIT}"
                else
                    # Use job arrays following -t flag
                    echo "Using custom job arrays following -t flag: ${CUSTOM_JOB_ARR}" >&2
                    job_array_index_range="${CUSTOM_JOB_ARR}"
                fi
            else # If we are creating a single gendb workspace
                echo "Using a single job"
                job_array_index_range="0"
            fi
            sbatch --job-name=${PROJECT}_Indel_Realigner --array=${job_array_index_range} ${IR_SBATCH} -o ${ERROR}/${PROJECT}_Indel_Realigner.%A_%a.out -e ${ERROR}/${PROJECT}_Indel_Realigner.%A_%a.err ${SEQUENCE_HANDLING}/SlurmJobScripts/Indel_Realigner.job ${CONFIG} ${SEQUENCE_HANDLING} ${IR_MEM} ${gatkVer}
        fi
        ;;
    18 | Freebayes_Variant_Calling )
        echo "$(basename $0): Variant calling using Freebayes..." >&2
        echo "Freebayes_Variant_Calling is not yet functional, exiting..." >&2
        exit 1
        source "${SEQUENCE_HANDLING}"/Handlers/Freebayes_VC.sh
        checkDependencies Freebayes_VC_Dependencies[@]
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        ;;
    1NP | NP_Quality_Assessment )
        echo "$(basename $0): Assessing quality..." >&2
        source "${SEQUENCE_HANDLING}"/Handlers/NP_Quality_Assessment.sh
        checkDependencies NP_Quality_Assessment_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        #if [[ ! -d "${FAST5_DIRECTORY}" ]]; then echo "Error: Failed to find FAST5 directory at "${FAST5_DIRECTORY}", exiting..." >&2; exit 1; fi # If the sample directory doesn't exist, exit
        checkSamples "${QA_SAMPLE_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        declare -a NP_QA_ARRAY=($(grep -E ".fastq|.fastq.gz" "${QA_SAMPLE_LIST}")) # Create an array of the files
        SINGLE_ARRAY_LIMIT=$[${#NP_QA_ARRAY[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
        then # If we have enough samples for a task array
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Quality_Assessment.sh && NP_Quality_Assessment ${QA_SAMPLE_LIST} ${OUT_DIR} ${PROJECT}" | qsub -t 0-"${SINGLE_ARRAY_LIMIT}" -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Quality_Assessment
        else # If we have only one sample
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Quality_Assessment.sh && NP_Quality_Assessment ${QA_SAMPLE_LIST} ${OUT_DIR} ${PROJECT}" | qsub -t 0 -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Quality_Assessment
        fi
        ;;
    2NP | NP_Adapter_Trimming )
        echo "$(basename $0): Trimming adapters..." >&2
        source "${SEQUENCE_HANDLING}/Handlers/NP_Adapter_Trimming.sh"
        checkDependencies NP_Adapter_Trimming_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${AT_SAMPLE_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        declare -a AT_ARRAY=($(grep -E ".fastq|.fastq.gz|.fasta|.fasta.gz|.fa|.fq|.fa.gz|.fq.gz" "${AT_SAMPLE_LIST}")) # Create an array of the files
        SINGLE_ARRAY_LIMIT=$[${#AT_ARRAY[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
        then # If we have enough samples for a task array
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Adapter_Trimming.sh && NP_Adapter_Trimming ${AT_ARRAY} ${OUT_DIR} ${PROJECT}" | qsub -t 0-"${SINGLE_ARRAY_LIMIT}" -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Adapter_Trimming
        else # If we have only one sample
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Adapter_Trimming.sh && NP_Adapter_Trimming ${AT_ARRAY} ${OUT_DIR} ${PROJECT}" | qsub -t 0 -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Adapter_Trimming
        fi
        ;;
    3NP | NP_Read_Mapping )
        echo "$(basename $0): Mapping reads..." >&2
        source "${SEQUENCE_HANDLING}/Handlers/NP_Read_Mapping.sh"
        checkDependencies NP_Read_Mapping_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${RM_SAMPLE_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        declare -a RM_ARRAY=($(grep -E ".fastq|.fastq.gz|.fasta|.fasta.gz|.fa|.fq|.fa.gz|.fq.gz" "${RM_SAMPLE_LIST}")) # Create an array of the files
        SINGLE_ARRAY_LIMIT=$[${#RM_ARRAY[@]} - 1] # Get the maximum number of Torque tasks (# in array - 1)
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
        then # If we have enough samples for a task array
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Read_Mapping.sh && NP_Read_Mapping ${RM_SAMPLE_LIST} ${OUT_DIR} ${OUTPUT_FORMAT} ${REF_GEN} ${REFERENCE_MINIMAP_INDEX} ${BANDWIDTH} ${MATCHING_SCORE} ${MISMATCH_PENALTY} ${GAP_OPEN_PENALTY} ${GAP_EXTENSION_PENALTY} ${Z_DROP_SCORE} ${MINIMAL_PEAK_DP_SCORE} ${THREADS}" | qsub -t 0-"${SINGLE_ARRAY_LIMIT}" -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Read_Mapping
        else # If we have only one sample
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_Read_Mapping.sh && NP_Read_Mapping ${RM_SAMPLE_LIST} ${OUT_DIR} ${OUTPUT_FORMAT} ${REF_GEN} ${REFERENCE_MINIMAP_INDEX} ${BANDWIDTH} ${MATCHING_SCORE} ${MISMATCH_PENALTY} ${GAP_OPEN_PENALTY} ${GAP_EXTENSION_PENALTY} ${Z_DROP_SCORE} ${MINIMAL_PEAK_DP_SCORE} ${THREADS}" | qsub -t 0 -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_Read_Mapping
        fi
        ;;
    4NP | NP_SAM_Processing )
        echo "$(basename $0): Converting SAM files into BAM files..." >&2
        source "${SEQUENCE_HANDLING}/Handlers/NP_SAM_Processing.sh"
        checkDependencies NP_SAM_Processing_Dependencies[@] # Check to see if the dependencies are installed
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a dependency, exit out with error
        checkSamples "${SP_SAMPLE_LIST}" # Check to see if samples and sample list exists
        if [[ "$?" -ne 0 ]]; then exit 1; fi # If we're missing a sample or our list, exit out with error
        checkVersion 'samtools' '1.3' # Check SAMtools version 1.3 or higher
        if [[ "$?" -ne 0 ]]; then echo "Please use SAMtools version 1.3 or higher" >&2; exit 1; fi
        #   Create the header for the mapping stats summary file
        mkdir -p "${OUT_DIR}/SAM_Processing/Statistics"
        echo -e "Sample name\tTotal reads\tPercent mapped" > "${OUT_DIR}/SAM_Processing/Statistics/${PROJECT}_mapping_summary.tsv"
        #   Run SAM_Processing using a task array
        declare -a SAM_LIST=($(grep -E ".sam" "${SP_SAMPLE_LIST}"))
        SINGLE_ARRAY_LIMIT=$[${#SAM_LIST[@]} - 1] # Get the maximum number of Torque tasks we're doing for SAM samples
        echo "Max array index is ${SINGLE_ARRAY_LIMIT}...">&2
        if [[ "${SINGLE_ARRAY_LIMIT}" -ne 0 ]]
        then # If we have enough samples for a task array
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_SAM_Processing.sh && NP_SAM_Processing ${SP_SAMPLE_LIST} ${OUT_DIR} ${REF_GEN} ${PROJECT}" | qsub -t 0-"${SINGLE_ARRAY_LIMIT}" -l "${SP_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_SAM_Processing
        else # If we have only one sample
            echo "source ${CONFIG} && source ${SEQUENCE_HANDLING}/Handlers/NP_SAM_Processing.sh && NP_SAM_Processing ${SP_SAMPLE_LIST} ${OUT_DIR} ${REF_GEN} ${PROJECT}" | qsub -t 0 -l "${SP_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_NP_SAM_Processing
        fi
        ;;
    * )
        Usage
        ;;
esac
