#!/bin/bash

#   Uses awk to remove 33 (or 64) from each header in the seqqs quality
#   score matrix output.
#   Written by Tom Kono

MATRIX=$1
OFFSET=$2

#   Create a temporary file
#   This does NOT work on OS X
TEMPFILE=`mktemp`

#   Uncomment this line for OS X
#TEMPFILE=`mktemp -t /tmp`

awk -v offset=$OFFSET '
NR==1 {
    OFS="\t"
    gsub("Q", "")
    for(i=1; i<=NF; i++)
        $i = $i - offset
    print
    }
NR>1 {
    print
}' $MATRIX > $TEMPFILE

mv $TEMPFILE ${MATRIX}_adj
