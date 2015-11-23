#!/bin/bash

function checkSamples() {
    local sample_list="$1"
    if [[ -f "$sample_list" ]]
    then
        for sample in `cat "$sample_list"`
        do
            if ! [[ -f "$sample" ]]
            then
                echo "$sample doesn't exist, exiting..."
                return 1
            fi
        done
    else
        echo "$sample_list doesn't exist, exiting..."
        return 1
    fi
}

export -f checkSamples
