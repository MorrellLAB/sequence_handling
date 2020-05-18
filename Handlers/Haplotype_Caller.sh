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
		if [[ "${USE_PBS}" == true ]]; then
			local sample="${sample_array[${PBS_ARRAYID}]}" # Which sample are we working on
            local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
			local fixedBQ=false
			if [[ "${qscores}" == true ]]; then
				gatk FixMisencodedBaseQualityReads -I ${sample} -O ${sample}-fixed.bam
				if [[ "$?" -eq 0 ]]; then
					sample="${sample}-fixed.bam"
					fixedBQ=true
				else
					echo "WARN: gatk FixMisencodedBaseQualityReads didn't run correctly. Using the original: $sample" >&2
				fi
			fi
            # Check if we parallelizing across regions
            if [ "${parallelize}" == "true" ] && [ "${custom_intervals}" != false ]; then
                # Create array of regions
                new_arr=()
                for i in $(cat ${sample_list})
                do
                    for x in $(cat ${custom_intervals})
                    do
                        temp=$(printf "${i}-${x}\n")
                        new_arr+=( ${temp} )
                    done
                done
				out_name_arr=()
				for i in ${new_arr[@]}
				do
					temp_intvl=$(basename $i | cut -d'-' -f 2)
					out_name_arr+=( ${temp_intvl} )
				done
				# If we have scaffolds or sequences not covered by the chromosomes,
				# append scaffolds list to array
				if [[ "${scaffolds}" != "false" ]]; then
					for i in $(cat ${sample_list})
					do
						temp=$(printf "${i}-${scaffolds}\n")
						new_arr+=("${temp}")
						out_name_arr+=($(printf "additional_intervals\n"))
					done
				fi
                local sample=$(echo ${new_arr[${PBS_ARRAYID}]} | cut -d'-' -f 1) # Which sample are we working on
                local current_intvl=$(echo ${new_arr[${PBS_ARRAYID}]} | cut -d'-' -f 2) # Current region we are processing
				local current_intvl_name="${out_name_arr[${PBS_ARRAYID}]}"
                local sample_name=$(basename ${sample} .bam) # What is the sample name without the suffix?
                #regions_arr=($(cat ${custom_intervals}))
                set -x
                #parallel --verbose --tmpdir ${tmp} gatk --java-options -Xmx${memory} HaplotypeCaller -R ${reference} -I ${sample} -O ${out}/${sample_name}_{}_RawGLs.g.vcf -L {} --heterozygosity ${heterozygosity} --native-pair-hmm-threads ${num_threads} --emit-ref-confidence GVCF ${settings} ::: ${regions_arr[@]}
                gatk --java-options "-Xmx${memory}" \
                        HaplotypeCaller \
                        -R "${reference}" \
                        -I "${sample}" \
                        -O "${out}/${sample_name}_${current_intvl_name}_RawGLs.g.vcf" \
                        -L "${current_intvl}" \
                        --heterozygosity "${heterozygosity}" \
                        --native-pair-hmm-threads "${num_threads}" \
                        --emit-ref-confidence GVCF \
                        ${settings}
                set +x
            else
                # We are not parallelizing across regions
                # Check if we are using custom intervals
                if [ "${parallelize}" == "false" ] && [ "${custom_intervals}" != false ]; then
					set -x
                    gatk --java-options "-Xmx${memory}" \
                        HaplotypeCaller \
                        -R "${reference}" \
                        -I "${sample}" \
                        -O "${out}/${sample_name}_RawGLs.g.vcf" \
                        -L "${custom_intervals}" \
                        --heterozygosity "${heterozygosity}" \
                        --native-pair-hmm-threads "${num_threads}" \
                        --emit-ref-confidence GVCF \
                        ${settings}
                    set +x
                else
                    set -x
                    gatk --java-options "-Xmx${memory}" \
                        HaplotypeCaller \
                        -R "${reference}" \
                        -I "${sample}" \
                        -O "${out}/${sample_name}_RawGLs.g.vcf" \
                        --heterozygosity "${heterozygosity}" \
                        --native-pair-hmm-threads "${num_threads}" \
                        --emit-ref-confidence GVCF \
                        ${settings}
                    set +x
                fi
            fi
			if [[ "$fixedBQ" == true ]]; then rm -f ${sample}; fi # remove the temp file
		else # No PBS
			if [[ "${qscores}" == true ]]; then
				# printf '%s\n' "${sample_array[@]}" | parallel "gatk FixMisencodedBaseQualityReads -I {} -O {}-fixed.bam"
				declare -a fixedBQ_arr=()
				for index in "${!sample_array[@]}"; do
					gatk FixMisencodedBaseQualityReads -I ${sample_array[${index}]} -O "${saample_array[$index]}-fixed.bam"
					if [[ "$?" -eq 0 ]]; then
						sample_array[$index]="${sample_array[$index]}-fixed.bam"
						fixedBQ_arr[${index}]=true
					else
						echo "WARN: gatk FixMisencodedBaseQualityReads didn't run correctly. Using the original: ${sample_array[$index]}" >&2
						fixedBQ_arr[${index}]=false
					fi
				done
			fi
			if [[ -z "${HAPLOTYPE_CALLER_THREADS}" ]]; then
				HAPLOTYPE_CALLER_THREADS=0   # Use all cores if not defined
			fi
			printf '%s\n' "${sample_array[@]}" | parallel --jobs ${HAPLOTYPE_CALLER_THREADS} "(set -x; \
					gatk --java-options -Xmx${memory} \
					HaplotypeCaller \
					-R ${reference} \
					-I {} \
					-O ${out}/{/.}_RawGLs.g.vcf \
					--genotyping-mode DISCOVERY \
					--heterozygosity ${heterozygosity} \
					--emit-ref-confidence GVCF \
					${settings})"
			# Clean temp files *-fixed.bam
			for index in "${!sample_array[@]}"; do if [[ "${fixedBQ_arr[${index}]}" == "true" ]]; then rm -f ${sample_array[${index}]}; fi; done
		fi
    fi
}

#   Export the function
export -f Haplotype_Caller
