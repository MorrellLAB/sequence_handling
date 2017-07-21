#!/bin/bash

set -e
set -o pipefail

#   Check to make sure we have our argument
if [ "$#" -ne 1 ]
then
    echo -e "\
Usage:  ./$0 sample_list.txt \n\
where:  'sample_list.txt' is a file containing a list of samples \n\
\n\
$0 is designed to check a sample list to ensure that it is properly formatted \n\
for use with sequence_handling. There are two checks run in this script: \n\
    First, we check to make sure files exist. The best way to ensure that a \n\
        sample list passes this check is to use full file paths for each sample \n\
        in the list.
    Second, we check to make sure all samples have unique names. We use the names given \n\
        to each sample to name subsequent files; if any have the same name, outfiles will be \n\
        overwritten. We check this to make sure no files are overwritten while using sequence_handling \n\
" >&2
    exit 1
fi

#   Assign our argument to a variable
SAMPLE_INFO=$1


#   Check to make sure files exist
TIME=`date +%m-%d-%y-%H.%M.%S` # Figure out what the time is so that the file with missing samples isn't one messy file
declare -a MISSING # Set up an array to hold missing samples
counter=0 # Start a counter for adding to the bash array
for sample in `cat "${SAMPLE_INFO}"`
do
    if [[ ! -f "$sample" ]] # If this sample doesn't exist
    then
        echo "Cannot find $sample" # Say which is missing
        MISSING["$counter"]="$sample" # Add the sample to our missing array
        let "counter += 1" # Increment the counter for the next index
    fi
done

if [[ ! -z "${MISSING}" ]] # If we're missing something
then
    #   Write to a missing_samples file
    echo "Cannot find:" > missing_samples_"${TIME}".txt
    for missing in "${MISSING[@]}"
    do
        echo "$missing" >> missing_samples_"${TIME}".txt
    done
    #   Say where the missing samples are
    echo "A list of samples missing can be found at `pwd`/missing_samples_${TIME}.txt" # Say where the file of missing samples is
    exit 1 # Exit
fi

#   Make sure sample names are unique
declare -a sample_names=() # Set up an array to hold these sample names
for i in `seq 0 "$(( $( wc -l < ${SAMPLE_INFO} ) - 1 ))"`
do
    sample=`basename $( head -"$(( $i + 1 ))" "${SAMPLE_INFO}" | tail -1 )`
    sample_names["$i"]="$sample"
done

oldIFS="$IFS" # Save the IFS variable
sorted_samples=($(sort <<< "${sample_names[@]}")) # Sort the sample_names array
IFS="$oldIFS" # Restore the IFS variable

declare -a unique_names=(`tr ' ' '\n' <<< "${sorted_samples[@]}" | sort -u | tr '\n' ' '`) # Create an array of unique sample names

if [[ "${#sorted_samples[@]}" -ne "${#unique_names[@]}" ]]
then
    echo "$(( ${#sorted_samples[@]} - ${#unique_names[@]} )) duplicate sample name(s) found!"
    oldIFS="$IFS" # Save the IFS variable
    IFS=$'\n\t' # Set a new IFS variable to trick 'comm' into working with arrays
    Differences=($( comm --nocheck-order -3 <(echo "${sorted_samples[@]}") <(echo "${unique_names[@]}") ) )
    IFS="$oldIFS"
    declare -p Differences
    exit 5
fi

echo "No problems found."