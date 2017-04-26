#!/bin/bash

#   This script creates a reference that is
#   compatible with Long Ranger.

set -o pipefail

#   What are the dependencies of Make_Reference?
declare -a Make_Reference_Dependencies=(longranger)

#   A function to run make reference
function Make_Reference() {
    local reference="$1" # What is our reference FASTA file?
    local outDir="$2" # Where are we storing our results?
    mkdir -p "${out}" # Make our output directory
    cd "${out}" # Long Ranger creates new directory in current working directory
    longranger mkref "${reference}" # Create compatible reference
}

export -f Make_Reference
