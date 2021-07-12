#!/bin/bash

#   This script performs preliminary SNP calling
#   on a single BAM sample using GATK.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_HaplotypeCaller.job

set -e
set -o pipefail

#   What are the dependencies for Haplotype_Caller?
declare -a Haplotype_Caller_Dependencies=(java samtools parallel)

function Fix_Qscores() {
    local sample_list="$1" # What is our sample list?
	local qscores="$2"
	# Check if out directories exist
	mkdir -p ${out_dir}
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array

    if [[ "${USE_PBS}" == true ]] || [[ "${USE_SLURM}" == true ]]; then
		if [[ "${USE_PBS}" == true ]]; then
			local sample=${sample_arr[${PBS_ARRAYID}]} # Which sample file are we working on
		elif [[ "${USE_SLURM}" == true ]]; then
			local sample=${sample_arr[${SLURM_ARRAY_TASK_ID}]} # Which sample file are we working on
		fi
		#local fixedBQ=false # It doesn't look like this gets used anywhere, will comment out for now
		if [[ "${qscores}" == true ]]; then
			gatk FixMisencodedBaseQualityReads -I ${sample} -O ${sample}-fixed.bam
			if [[ "$?" -ne 0 ]]; then
				echo "WARN: gatk FixMisencodedBaseQualityReads didn't run correctly. Using the original: $sample" >&2
				rm -f ${sample}-fixed.bam
				ln -s ${sample} ${sample}-fixed.bam
			fi
		fi
    else # No PBS
		# I'm not using this parallel since I have to deal with the case when it fails
		# printf '%s\n' "${sample_array[@]}" | parallel "gatk FixMisencodedBaseQualityReads -I {} -O {}-fixed.bam"
		for index in "${!sample_array[@]}"; do
			gatk FixMisencodedBaseQualityReads -I ${sample_array[${index}]} -O "${sample_array[$index]}-fixed.bam"
			if [[ "$?" -ne 0 ]]; then
				echo "WARN: gatk FixMisencodedBaseQualityReads didn't run correctly. Using the original: ${sample_array[$index]}" >&2
				rm -f ${sample}-fixed.bam
				ln -s ${sample} ${sample}-fixed.bam
			fi
		done
    fi
}

export -f Fix_Qscores

#   A function to run the SNP calling
function Haplotype_Caller_GATK4() {
	local sample_list="$1" # What is our sample list?
    local out="$2"/Haplotype_Caller # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local memory="$6" # How much memory can java use?
    local num_threads="$7" # How many threads can we use?
    local qscores="$8" # Do we fix quality scores?
    local trim_active="$9" # Do we trim the active regions?
    local force_active="${10}" # Do we force all bases to be active?
    #local gatkVer="${11}" # Either 3 or 4
    local parallelize="${11}" # Are we parallelizing across regions?
    local custom_intervals="${12}" # List of custom intervals
    local scaffolds="${13}" # List of scaffolds or sequences not covered by chromosomes
    local tmp="${14}" # temp directory
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
	# Perform the most basic checks
	mkdir -p "${out}" # Make sure the out directory exists
    # Check if tmp variable is empty, if not empty make directory
    if [ -z $test_var ]; then
        # tmp variable is empty
        echo "Not using temp directory, proceeding..."
    else
        echo "Making temp directory"
        mkdir -p "${tmp}"
    fi
    
	# Check if we need to fix quality scores
	if [[ ${qscores} == true ]]
	then
		echo "Fixing quality scores before running Haplotype Caller function..."
		Fix_Qscores ${sample_list} ${qscores}
	fi

	# Index the reference is needed
    if ! [[ -f "${reference}.fai" ]]; then samtools faidx "${reference}"; fi
    #   Determine the GATK settings based on the parameters set in the Config file
    local settings='' # Make an empty string
	# Options changes in GATK4
	# Options gone in GATK4: -nct -variant_index_type -variant_index_parameter
	# -o became -O
	# -genotyping-mode ('-' instead of '_')
	# -ERC spelled differently.
	# Need to use FixMisencodedBaseQualityReads instead of --fix_misencoded_quality_scores
	if [[ "${trim_active}" == true ]]; then settings="${settings} --dont-trim-active-regions "; fi
	if [[ "${force_active}" == true ]]; then settings="${settings} --active-probability-threshold 0.000 "; fi

	# Prepare intervals list filepath. If no custom intervals are provided this is ignored.
	if [[ "${USE_PBS}" == "true" ]]; then
		## With PBS, multiple processes might write to this file, so
		## making unique file for each job.  Is PBS_JOBID better?
		local intervals_filepath=$(echo "${out}/intervals-${PBS_ARRAYID}.list")
	elif [[ "${USE_SLURM}" == "true" ]]; then
		#intervals_filepath=$(echo "${out}/intervals-${test_arrayid}.list") # testing
		local intervals_filepath=$(echo "${out}/intervals-${SLURM_ARRAY_TASK_ID}.list")
	else
		local intervals_filepath=$(echo "${out}/intervals.list")
	fi

	# Check if we have custom intervals
	if [[ "${custom_intervals}" == false ]]; then
		echo "No custom intervals were provided, proceeding to run Haplotype Caller..."
	elif [[ "${custom_intervals}" != false ]]; then
		# We have a list of custom intervals, check the file extension of the intervals
		if [[ "${custom_intervals}" == *.bed ]]; then
			custom_intervals="${out_dir}/Haplotype_Caller/intervals.list"
			# converting bed file to GATK .intervals/.list
			sort -k1,1 -k2,2n "${custom_intervals}" | \
				perl -ne 'chomp; @a=split/\t/;print $a[0], ":", $a[1]+1,"-", $a[2], "\n"' > "${intervals_filepath}"
			# With the perl 1-liner
			# each line of *.bed becomes firstColumn:(2nd column+1)-(3rd column)
			# NOTE: https://software.broadinstitute.org/gatk/documentation/article?id=1319
			#   GATK interval uses 1-based (first position in the genome is position 1, not position 0).
			#   But bed file has the start position 0-based, and end position 1-based, so 1 is added to 2nd column
		elif [[ "${custom_intervals}" == *.intervals ]] || [[ "${custom_intervals}" == *.list ]]; then
			# check that first position does not start at 0.
			if [ $(grep ":0-" "${custom_intervals}" | wc -l) != 0 ]; then
				# For GATK-style .list or .intervals file in the form <chr>:<start>-<stop>
				# GATK 4 will throw an error if the first position starts at 0.
				# This could indicate a user error in interval conversion from .bed to .list
				echo "GATK-style .list or .intervals file uses 1-based coordinates. Please check your custom intervals file to make sure there are no positions that start at 0. If so, change so the starting position is 1. See docs: https://software.broadinstitute.org/gatk/documentation/article?id=11009. Exiting..." >&2
				exit 30
			fi
			cut -f 1 "${custom_intervals}" | sort -V | uniq > "${intervals_filepath}"
		else
			echo "ERROR: In Config, HC_CUSTOM_INTERVALS is set to ${custom_intervals}. It needs following file extensions: .intervals, .list, or .bed" >&2
			exit 30
		fi
	fi

	# Notes on intended behavior and cases to handle:
	#	1. custom_intervals=false, scaffolds=false, parallelize=false: Don't need to use -L option, provide bam file as input to -I.
	#	2. custom_intervals=false, scaffolds=path, parallelize=false: In this case, we won't parallelize across regions but the scaffolds will be included in the run.
	#	3. custom_intervals=path, scaffolds=false, parallelize=true: Use -L option for custom_intervals, each interval will run as its own job array. custom_intervals here is one interval per line and can include either just chromosome names or chromosome names plus positions.
	#	4. custom_intervals=path, scaffolds=path, parallelize=true: Same as case 2 except we will have one extra array index added for the scaffolds (IMPORTANT: these are NOT covered by the custom_intervals) and the list of scaffolds will be passed to the -L option. We don't want one array for each scaffold since often times there are thousands of these.
	#	5. custom_intervals=false, scaffolds=false, parallelize=true: Invalid combination, we only want to parallelize if at least custom_intervals list is provided.
	#	6. custom_intervals=false, scaffolds=path, parallelize=true: Invalid combination, in this case we only want to parallelize if both custom_intervals and scaffolds are provided.

	## Create arrays of input filenames, regions, output filenames.
	## Easier to deal with PBS and non-PBS parallelization with the same code.
	## Also, this method doesn't have problems if the path/filename contains '-'.
	## Finally, scaffolds list was only with parallelizing across the region. But it is fixed.
	## scaffolds list isn't used if custom_intervals aren't set.  Is this OK?
	sample_file_arr=()  # path/filenames of input bam
	intvl_arr=()        # each element is either "-L chr1" or "-L chr1:1-10000" etc, the option which can be passed to gatk
	out_filename_arr=() # output path/filenames
	
	# Handle case where we are NOT parallelizing.
	#	Case 1: custom_intervals=false, scaffolds=false, parallelize=false
	#	Case 2: custom_intervals=false, scaffolds=path, parallelize=false
	if [[ "${custom_intervals}" == false ]] && [[ "${parallelize}" == "false" ]]
	then
		echo "Not parallelizing across regions, preparing arrays..."
		# Need sample_file_arr and out_filename_arr filled, intvl_arr should be empty
		for i in ${sample_array[@]}
		do
			sample_file_arr+=( "${i}" )
			this_sample_name=$(basename ${i} .bam)
			out_filename_arr+=( "${out}/${this_sample_name}_RawGLs.g.vcf" )
		done
	# Handle cases when parallelize is set to true.
	#	Case 3: custom_intervals=path, scaffolds=false, parallelize=true
	#	Case 4: custom_intervals=path, scaffolds=path, parallelize=true
	elif [[ "${custom_intervals}" != false ]] && [[ "${parallelize}" == "true" ]]
	then
		echo "Parallelizing across regions and/or scaffolds, preparing arrays..."
		# The length of arrays become (number of samples) * (number of intervals)
		for i in ${sample_array[@]}
		do
			#local this_sample_name=$(basename ${i} .bam)
			this_sample_name=$(basename ${i} .bam)
			for x in $(cat "${intervals_filepath}")
			do
				sample_file_arr+=( ${i} )
				intvl_arr+=( "${x}" )
				# output name is sample_interval_RawGLs.g.vcf
				out_filename_arr+=( "${out}/${this_sample_name}_${x}_RawGLs.g.vcf" )
			done
		done
		# Check if we have scaffolds (Case 3) or not (Case 2)
		if [[ "${scaffolds}" == "false" ]]
		then
			echo "We don't have scaffolds to add, proceeding to next step..."
		elif [[ "${scaffolds}" != "false" ]]
		then
			# We have scaffolds or sequences not covered by the chromosomes,
			# append filename of scaffolds list to array
			for i in ${sample_array[@]}
			do
				#local this_sample_name=$(basename ${i} .bam)
				this_sample_name=$(basename ${i} .bam)
				sample_file_arr+=( ${i} )
				intvl_arr+=( "${scaffolds}" )
				out_filename_arr+=( "${out}/${this_sample_name}_additional_intervals_RawGLs.g.vcf" )
			done
		fi
	# Handle invalid cases.
	#	Case 5: custom_intervals=false, scaffolds=false, parallelize=true
	#	Case 6: custom_intervals=false, scaffolds=path, parallelize=true
	elif [[ "${custom_intervals}" == false ]] && [[ "${parallelize}" == "true" ]]
	then
		echo "Invalid combination of input variables for custom intervals, scaffolds, and parallelize. If parallelize is set to true, a custom intervals list must be provided. Please modify your config and resubmit the handler. Exiting..."
		exit 30
	fi

	if [[ "${qscores}" == true ]]; then
		# Appending "-fixed.bam" to each element of the array
		# Previous to this Haplotype_Caller() function, Fix_Qscores() created *-fixed.bam, which should be used as the input
		sample_file_arr=( "${sample_file_arr[@]/%/-fixed.bam}" )
	fi

	# Run Haplotype Caller function
	if [[ "${USE_PBS}" == true ]] || [[ "${USE_SLURM}" == true ]]
	then
		if [[ "${USE_PBS}" == true ]]
		then
			sample="${sample_file_arr[${PBS_ARRAYID}]}" # Which sample file are we working on
			out_filename="${out_filename_arr[${PBS_ARRAYID}]}"
			intvl_opt="${intvl_arr[${PBS_ARRAYID}]}"
		elif [[ "${USE_SLURM}" == true ]]
		then
			sample="${sample_file_arr[${SLURM_ARRAY_TASK_ID}]}" # Which sample file are we working on
			out_filename="${out_filename_arr[${SLURM_ARRAY_TASK_ID}]}"
			intvl_opt="${intvl_arr[${SLURM_ARRAY_TASK_ID}]}"
		fi
		# Check if we will add an intervals option
		if [[ "${intvl_opt}" != "" ]]
		then
			intvl_opt=('-L' "${intvl_opt}")
		else
			intvl_opt=()
		fi
		# Run HaplotypeCaller
		# Check if we have a tmp directory specified
		if [[ -z "${tmp}" ]]; then
			# No tmp directory specified
			set -x
			gatk --java-options "-Xmx${memory}" \
				HaplotypeCaller \
				-R "${reference}" \
				-I "${sample}" \
				-O "${out_filename}" \
				"${intvl_opt[@]}" \
				--heterozygosity "${heterozygosity}" \
				--native-pair-hmm-threads "${num_threads}" \
				--emit-ref-confidence GVCF \
				${settings}
			set +x
		else
			# tmp directory is specified
			set -x
			gatk --java-options "-Xmx${memory}" \
				HaplotypeCaller \
				-R "${reference}" \
				-I "${sample}" \
				-O "${out_filename}" \
				"${intvl_opt[@]}" \
				--tmp-dir ${tmp} \
				--heterozygosity "${heterozygosity}" \
				--native-pair-hmm-threads "${num_threads}" \
				--emit-ref-confidence GVCF \
				${settings}
			set +x
		fi
	else # Not running on cluster
		if [[ -z "${HAPLOTYPE_CALLER_THREADS}" ]]; then
			HAPLOTYPE_CALLER_THREADS='100%'   # Use all cores if not defined
		fi
		local temp_sample_filepath="${out}/temp_sample.txt"
		local temp_out_filepath="${out}/temp_out.txt"
		local temp_intvl_filepath="${out}/temp_intvl.txt"
		printf '%s\n' "${sample_file_arr[@]}" > "${temp_sample_filepath}"
		printf '%s\n' "${out_filename_arr[@]}" > "${temp_out_filepath}"
		printf '%s\n' "${intvl_arr[@]}" > "${temp_intvl_filepath}"
		
		# Naoki's Note: I wanted to make addition of -L automatic based on whether ${intvl_arr[]} is "" or not
		# But I can't make parallel behaves well... Fix here if possible.
		if [ "${custom_intervals}" == false ]; then
			set -x
			parallel --jobs ${HAPLOTYPE_CALLER_THREADS} \
				gatk --java-options -Xmx${memory} \
				HaplotypeCaller \
				-R "${reference}" \
				-I {1} \
				-O {2} \
				--heterozygosity "${heterozygosity}" \
				--emit-ref-confidence GVCF \
				"${settings}" \
				:::: "${temp_sample_filepath}" ::::+ "${temp_out_filepath}" ::::+ "${temp_intvl_filepath}"
				# ::: "${sample_file_arr[@]}" :::+ "${out_filename_arr[@]}" :::+ "${intvl_arr[@]}"
			set +x
		else
			set -x
			parallel --jobs ${HAPLOTYPE_CALLER_THREADS} \
				gatk --java-options -Xmx${memory} \
				HaplotypeCaller \
				-R "${reference}" \
				-I {1} \
				-O {2} \
				-L $( printf -- '%s' {3} )  \
				--heterozygosity "${heterozygosity}" \
				--emit-ref-confidence GVCF \
				"${settings}" \
				:::: "${temp_sample_filepath}" ::::+ "${temp_out_filepath}" ::::+ "${temp_intvl_filepath}"
				# ::: "${sample_file_arr[@]}" :::+ "${out_filename_arr[@]}" :::+ "${intvl_arr[@]}"
			set +x
		fi
		rm -f "${temp_sample_filepath}" "${temp_out_filepath}" "${temp_intvl_filepath}"
	fi
	# Cleanup
	if [ "${custom_intervals}" != false ]; then
		echo "Cleaning ${intervals_filepath}"
		rm -f "${intervals_filepath}"
	fi
}

export -f Haplotype_Caller_GATK4

#	Separate GATK 3 code since moving forwards we will only be modifying GATK 4 code
#	This ensures any changes we make in the future don't affect the original GATK 3 functionallity
function Haplotype_Caller_GATK3() {
	local sample_list="$1" # What is our sample list?
    local out="$2"/Haplotype_Caller # Where are we storing our results?
    local gatk="$3" # Where is the GATK jar?
    local reference="$4" # Where is the reference sequence?
    local heterozygosity="$5" # What is the nucleotide diversity/bp?
    local memory="$6" # How much memory can java use?
    local num_threads="$7" # How many threads can we use?
    local qscores="$8" # Do we fix quality scores?
    local trim_active="$9" # Do we trim the active regions?
    local force_active="${10}" # Do we force all bases to be active?
    #local gatkVer="${11}" # Either 3 or 4
    local parallelize="${11}" # Are we parallelizing across regions?
    local custom_intervals="${12}" # List of custom intervals
    local scaffolds="${13}" # List of scaffolds or sequences not covered by chromosomes
    local tmp="${14}" # temp directory
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
    mkdir -p "${out}" # Make sure the out directory exists
	# Index the reference is needed
    if ! [[ -f "${reference}.fai" ]]; then samtools faidx "${reference}"; fi
    #   Determine the GATK settings based on the parameters set in the Config file
    local settings='' # Make an empty string
	if [[ "${qscores}" == true ]]; then settings="--fix_misencoded_quality_scores "; fi
	if [[ "${trim_active}" == true ]]; then settings="${settings} --dontTrimActiveRegions "; fi
	if [[ "${force_active}" == true ]]; then settings="${settings} --forceActive "; fi
	
	#   Run GATK using the parameters given
	if [[ "${USE_PBS}" == true ]] || [[ "${USE_SLURM}" == true ]]; then
		# Check if we are using PBS or Slurm
		if [[ "${USE_PBS}" == true ]]; then
			local sample="${sample_array[${PBS_ARRAYID}]}" # Which sample are we working on currently?
		elif [[ "${USE_SLURM}" == true ]]; then
			local sample="${sample_array[${SLURM_ARRAY_TASK_ID}]}" # Which sample are we working on currently?
		fi
		local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
		set -x # prints extra info for troubleshooting/debugging
		java -Xmx"${memory}" -jar "${gatk}" \
			-T HaplotypeCaller \
			-R "${reference}" \
			-I "${sample}" \
			-o "${out}/${sample_name}_RawGLs.g.vcf" \
			-nct 8 \
			--genotyping_mode DISCOVERY \
			--heterozygosity "${heterozygosity}" \
			--emitRefConfidence GVCF \
			-variant_index_type LINEAR \
			-variant_index_parameter 128000 \
			${settings}
		set +x # ends debugging mode
	else
		printf '%s\n' "${sample_array[@]}" | parallel "(set -x;  java -Xmx${memory} -jar ${gatk} \
			-T HaplotypeCaller \
			-R ${reference} \
			-I {} \
			-o ${out}/{/.}_RawGLs.g.vcf \
			-nct 8 \
			--genotyping_mode DISCOVERY \
			--heterozygosity ${heterozygosity} \
			--emitRefConfidence GVCF \
			-variant_index_type LINEAR \
			-variant_index_parameter 128000 \
			${settings})"
	fi
}

export -f Haplotype_Caller_GATK3
