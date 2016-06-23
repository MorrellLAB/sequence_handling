#!/bin/bash

#   Check to make sure our samples exist
function checkSamples() {
    local sample_list="$1" # Sample ist
    if [[ -f "$sample_list" ]] # If the sample list exists
    then
        for sample in `cat "$sample_list"` # For each sample in the sample list
        do
            if ! [[ -f "$sample" ]] # If the sample doesn't exist
            then
                echo "$sample doesn't exist, exiting..." >&2 # Exit out with error
                return 1
            fi
        done
    else # If the sample list doesn't exist
        echo "$sample_list doesn't exist, exiting..." >&2 # Exit out with error
        return 1
    fi
}

#   Export the function to be used elsewhere
export -f checkSamples

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

#   Check versions of tools
function checkVersion() {
    local tool="$1"
    local version="$2"
    "${tool}" --version | grep "${version}" > /dev/null 2> /dev/null || return 1
}

#   Export the function to be used elsewhere
export -f checkVersion
