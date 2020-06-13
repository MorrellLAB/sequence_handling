#!/bin/bash

#   This script performs preliminary SNP calling
#   on a single BAM sample using GATK.

#   This code is modified from code written by Tom Kono at:
#   https://github.com/MorrellLAB/Deleterious_GP/blob/master/Job_Scripts/Seq_Handling/GATK_HaplotypeCaller.job

set -e
set -o pipefail

#   What are the dependencies for Haplotype_Caller?
declare -a Haplotype_Caller_Dependencies=(java samtools parallel)

#   A function to run the SNP calling
function Haplotype_Caller() {
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
    local gatkVer="${11}" # Either 3 or 4
    local parallelize="${12}" # Are we parallelizing across regions?
    local custom_intervals="${13}" # List of custom intervals
    local scaffolds="${14}" # List of scaffolds or sequences not covered by chromosomes
    local tmp="${15}" # temp directory
    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array
    mkdir -p "${out}" # Make sure the out directory exists
    # Index the reference is needed
    if ! [[ -f "${reference}.fai" ]]; then samtools faidx "${reference}"; fi
    #   Determine the GATK settings based on the parameters set in the Config file
    local settings='' # Make an empty string
    if [[ "${gatkVer}" == 3 ]]; then
		if [[ "${qscores}" == true ]]; then settings="--fix_misencoded_quality_scores "; fi
		if [[ "${trim_active}" == true ]]; then settings="${settings} --dontTrimActiveRegions "; fi
		if [[ "${force_active}" == true ]]; then settings="${settings} --forceActive "; fi
		#   Run GATK using the parameters given
		if [[ "${USE_PBS}" == true ]]; then
			local sample="${sample_array[${PBS_ARRAYID}]}" # Which sample are we working on currently?
			local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
			(set -x; java -Xmx"${memory}" -jar "${gatk}" \
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
			${settings})
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
    else
	# Options changes in GATK4
	# Options gone in GATK4: -nct -variant_index_type -variant_index_parameter
	# -o became -O
	# -genotyping-mode ('-' instead of '_')
	# -ERC spelled differently.
	# Need to use FixMisencodedBaseQualityReads instead of --fix_misencoded_quality_scores
	if [[ "${trim_active}" == true ]]; then settings="${settings} --dont-trim-active-regions "; fi
	if [[ "${force_active}" == true ]]; then settings="${settings} --active-probability-threshold 0.000 "; fi

	if [[ "${USE_PBS}" == "true" ]]; then
            ## With PBS, multiple processes might write to this file, so
            ## making unique file for each job.  Is PBS_JOBID better?
            local intervals_filepath=$(echo "${out}/intervals-${PBS_ARRAYID}.list")
	else
            local intervals_filepath=$(echo "${out}/intervals.list")
	fi
	if [ "${custom_intervals}" != false ]; then
	    if [[ "${custom_intervals}" == *.bed ]]; then
		custom_intervals="${out_dir}/Haplotype_Caller/intervals.list"
		# converting bed file to GATK .intervals/.list
		sort -k1,1 -k2,2n "${custom_intervals}" | \
		    perl -ne '@a=split/\t/;print $a[0], ":", $a[1]+1,"-", $a[2], "\n"' > "${intervals_filepath}"
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
		echo "ERROR: In Config, HC_CUSTOME_INTERVALS is set to ${custom_intervals}.  It needs following file extensions: .intervals, .list, or .bed" >&2
		exit 30
	    fi
	fi

	## Create arrays of input filenames, regions, output filenames.
	## Easier to deal with PBS and non-PBS parallelization with the same code.
	## Also, this method doesn't have problems if the path/filename contains '-'.
	## Finally, scaffolds list was only with parallelizing across the region. But it is fixed.
	## scaffolds list isn't used if custom_intervals aren't set.  Is this OK?
	sample_file_arr=()  # path/filenames of input bam
	intvl_arr=()        # each element is "-L chr1:1-10000" etc, the option which can be passed to gatk
	out_filename_arr=() # output path/filenames
        if [ "${parallelize}" == "true" ] && [ "${custom_intervals}" != false ]; then
	    # parallelizing across regions
	    # the length of arrays become (number of samples) * (number of intervals)
            for i in ${sample_array[@]}
            do
		local this_sample_name=$(basename ${i} .bam)
                for x in $(cat "${intervals_filepath}")
                do
		    sample_file_arr+=( ${i} )
		    intvl_arr+=( "${x}" )
		    # output name is sample_interval_RawGLs.g.vcf
		    out_filename_arr+=( "${out}/${this_sample_name}_${x}_RawGLs.g.vcf" )
                done
            done

	    # If we have scaffolds or sequences not covered by the chromosomes,
	    # append filename of scaffolds list to array
	    if [[ "${scaffolds}" != "false" ]]; then
		for i in ${sample_array[@]}
		do
		    local this_sample_name=$(basename ${i} .bam)
		    sample_file_arr+=( ${i} )
                    intvl_arr+=( "${scaffolds}" )
                    out_filename_arr+=( "${out}/${this_sample_name}_additional_intervals_RawGLs.g.vcf" )
		done
	    fi
	else
	    # no parallelizing across regions
	    if [[ "${custom_intervals}" != false && "${scaffolds}" != "false" ]]; then
		cat "${scaffolds}" >> "${intervals_filepath}"
	    fi
            for i in "${sample_array[@]}"
            do
		sample_file_arr+=( "${i}" )
		if [ "${parallelize}" == "false" ] && [ "${custom_intervals}" != false ]; then
		    intvl_arr+=( "${intervals_filepath}" )
		else
		    # if HC_CUSTOM_INTERVALS=false, HC_SCAFFOLDS is ignored.  Is this ok?
		    intvl_arr+=( "" )
		fi

		local this_sample_name=$(basename ${i} .bam)
		out_filename_arr+=( "${out}/${this_sample_name}_RawGLs.g.vcf" )
	    done
	fi

	if [[ "${qscores}" == true ]]; then
	    # Appending "-fixed.bam" to each element of the array
	    # Previous to this Haplotype_Caller() function, Fix_Qscores() created *-fixed.bam, which should be used as the input
	    sample_file_arr=( "${sample_file_arr[@]/%/-fixed.bam}" )
	fi

	if [[ "${USE_PBS}" == true ]]; then
	    local sample="${sample_file_arr[${PBS_ARRAYID}]}" # Which sample file are we working on
	    local out_filename="${out_filename_arr[${PBS_ARRAYID}]}"
	    local intvl_opt="${intvl_arr[${PBS_ARRAYID}]}"

	    if [[ "${intvl_opt}" != "" ]]; then
	       intvl_opt=('-L' "${intvl_opt}")
	    else
		intvl_opt=()
	    fi
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
	else # No PBS
	    if [[ -z "${HAPLOTYPE_CALLER_THREADS}" ]]; then
		HAPLOTYPE_CALLER_THREADS=0   # Use all cores if not defined
	    fi
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
			 "${settings}" ::: "${sample_file_arr[@]}" :::+ "${out_filename_arr[@]}" :::+ "${intvl_arr[@]}"
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
			 "${settings}" ::: "${sample_file_arr[@]}" :::+ "${out_filename_arr[@]}" :::+ "${intvl_arr[@]}"
		set +x
	    fi
	fi
	if [ "${custom_intervals}" != false ]; then
	    echo "Cleaning ${intervals_filepath}"
	    rm -f "${intervals_filepath}"
	fi
    fi
}

#   Export the function
export -f Haplotype_Caller

function Fix_Qscores() {
    local sample_list="$1" # What is our sample list?

    declare -a sample_array=($(grep -E ".bam" "${sample_list}")) # Turn the list into an array

    if [[ "${USE_PBS}" == true ]]; then
	local sample=${sample_arr[${PBS_ARRAYID}]} # Which sample file are we working on
	local fixedBQ=false
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
	    gatk FixMisencodedBaseQualityReads -I ${sample_array[${index}]} -O "${saample_array[$index]}-fixed.bam"
	    if [[ "$?" -ne 0 ]]; then
		echo "WARN: gatk FixMisencodedBaseQualityReads didn't run correctly. Using the original: ${sample_array[$index]}" >&2
		rm -f ${sample}-fixed.bam
		ln -s ${sample} ${sample}-fixed.bam
	    fi
	done
    fi
}

export -f Fix_Qscores
