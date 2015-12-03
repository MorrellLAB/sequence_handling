#!/bin/bash

#   Check to make sure our samples exist
function checkSamples() {
    local sample_list="$1"
    if [[ -f "$sample_list" ]]
    then
        for sample in `cat "$sample_list"`
        do
            if ! [[ -f "$sample" ]]
            then
                echo "$sample doesn't exist, exiting..." >&2
                return 1
            fi
        done
    else
        echo "$sample_list doesn't exist, exiting..." >&2
        return 1
    fi
}

export -f checkSamples

#   Check to make sure our dependencies are installed
function checkDependencies() {
    local dependencies=("${!1}") # BASH array to hold dependencies
    for dep in "${dependencies[@]}"
    do
        if ! `command -v "$dep" > /dev/null 2> /dev/null`
        then
            echo "Failed to find $dep installation, exiting..." >&2
            return 1
        fi
    done
}

export -f checkDependencies
