#!/bin/bash

#   A function to do gzipping while preserving the original file
function gzip_parts() {
    local sample="$1" # The vcf file to be gzipped
    local out="$2" # Where to put the gzipped file
    local name=$(basename ${sample} .vcf) # The name of the sample
    gzip -c "${sample}" > "${out}/${name}.vcf.gz" # Perform the gzipping without altering the original file
}

#   Export the function
export -f gzip_parts