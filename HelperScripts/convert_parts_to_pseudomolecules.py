#!/usr/bin/env python3

#   Code written by Tom Kono
"""A simple script to convert the pseudomolecule parts VCF files into full
pseudomolecule VCF files. Remember that VCF is 1-based, which makes the math
a lot nicer."""

import sys

#   Store the lengths of the pseudomolecule parts in a dictionary. We only
#   really need to store the first parts, since we will just add the part2
#   positions to it.
part_lengths = {
    'chr1H_part1': 312837513,
    'chr2H_part1': 393532674,
    'chr3H_part1': 394310633,
    'chr4H_part1': 355061206,
    'chr5H_part1': 380865482,
    'chr6H_part1': 294822070,
    'chr7H_part1': 325797516
    }

#   Then iterate through the VCF and print out the modified lines
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line.strip())
        else:
            tmp = line.strip().split('\t')
            #   Modify the chromosome to not have the part
            chrom = tmp[0].split('_')[0]
            #   And check if we have to modify the position. If the chromosome
            #   name has '_part2' in it, then we have to modify it.
            if '_part2' in tmp[0]:
                offset = part_lengths[chrom + '_part1']
                newpos = str(int(tmp[1]) + offset)
            else:
                newpos = tmp[1]
            #   Print out the VCF line with the new position and chromosome name
            toprint = [chrom, newpos] + tmp[2:]
            print('\t'.join(toprint))
