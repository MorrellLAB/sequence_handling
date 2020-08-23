#!/usr/bin/env python3
#   Code originally written by Tom Kono
"""A simple script to convert the pseudomolecule parts VCF files into full
pseudomolecule VCF files. Remember that VCF is 1-based, which makes the math
a lot nicer."""

import argparse
import gzip

#   Store the lengths of the barley pseudomolecule parts in a dictionary. We only
#   really need to store the first parts, since we will just add the part2
#   positions to it.
#   Sizes listed below are for Morex v2
PARTS_SIZES = {
    'chr1H_part1': 205502676,
    'chr2H_part1': 305853815,
    'chr3H_part1': 271947776,
    'chr4H_part1': 282386439,
    'chr5H_part1': 205989812,
    'chr6H_part1': 260041240,
    'chr7H_part1': 328847863
    }

#   Positions below are for Morex v1
# PARTS_SIZES = {
#     'chr1H_part1': 312837513,
#     'chr2H_part1': 393532674,
#     'chr3H_part1': 394310633,
#     'chr4H_part1': 355061206,
#     'chr5H_part1': 380865482,
#     'chr6H_part1': 294822070,
#     'chr7H_part1': 325797516
#     }

def parse_args():
    """Set up an argument parser, and actually parse the args."""
    parser = argparse.ArgumentParser(
        description='Convert parts coordinates to pseudomolecule coordinates.',
        add_help=True
    )
    # Add mutually exclusive group for the file format. This is so that only
    # one format can be specified.
    fmt = parser.add_mutually_exclusive_group(required=True)
    fmt.add_argument(
        '--vcf',
        '-v',
        action='store_true',
        help='Input file is VCF or VCF.gz'
    )
    fmt.add_argument(
        '--bed',
        '-b',
        action='store_true',
        help='Input file is BED'
    )
    parser.add_argument(
        'intervals',
        metavar='file',
        help='File with intervals to convert'
    )
    a = parser.parse_args()
    return a


def set_format(args):
    """Set the format that is being used for conversion."""
    arg_dict = vars(args)
    if arg_dict['vcf']:
        return 'VCF'
    elif arg_dict['bed']:
        return 'BED'
    else:
        return None


def vcf_conv(intervals, parts_sizes):
    """Convert VCF. Remember that VCF is 1-based, which makes the math a lot nicer.
    Handles both VCF and gzipped VCFs."""
    if '.gz' in intervals:
        # If lines start with '#', print them
        with gzip.open(intervals, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    print(line.strip())
                else:
                    tmp = line.strip().split()
                    # modify the chromosome to not have the part
                    chrom = tmp[0].split('_')[0]
                    # And check if we have to modify the position. If the chromosome
                    # name has '_part2' in it, then we have to modify it.
                    if '_part2' in tmp[0]:
                        offset = parts_sizes[chrom + '_part1']
                        newpos = str(int(tmp[1]) + offset)
                    else:
                        newpos = tmp[1]
                    # then print out the VCF line with the new position and chromosome
                    # name
                    toprint = [chrom, newpos] + tmp[2:]
                    print('\t'.join(toprint))
    else:
        # If lines start with '#', print them
        with open(intervals, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    print(line.strip())
                else:
                    tmp = line.strip().split()
                    # modify the chromosome to not have the part
                    chrom = tmp[0].split('_')[0]
                    # And check if we have to modify the position. If the chromosome
                    # name has '_part2' in it, then we have to modify it.
                    if '_part2' in tmp[0]:
                        offset = parts_sizes[chrom + '_part1']
                        newpos = str(int(tmp[1]) + offset)
                    else:
                        newpos = tmp[1]
                    # then print out the VCF line with the new position and chromosome
                    # name
                    toprint = [chrom, newpos] + tmp[2:]
                    print('\t'.join(toprint))


def bed_conv(intervals, parts_sizes):
    """Convert BED. Remember that BED is 0-based (i.e., the first base is 0)."""
    # If lines start with '#', print them
    with open(intervals, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
                continue
            else:
                tmp = line.strip().split()
                # modify the chromosome to not have th part
                chrom = tmp[0].split('_')[0]
                # And check if we have to modify the position. If the chromosome
                # name has '_part2' in it, then we have to modify it.
                if '_part2' in tmp[0]:
                    offset = parts_sizes[chrom + '_part1']
                    newstartpos = str(int(tmp[1]) + offset)
                    newendpos = str(int(tmp[2]) + offset)
                else:
                    newstartpos = tmp[1]
                    newendpos = tmp[2]
                # Check if standard BED file format with 3 columns
                if len(tmp) == 3:
                    # Then print out BED line with chromosome name and new start and end positions
                    toprint = [chrom, newstartpos, newendpos]
                    print('\t'.join(toprint))
                elif len(tmp) > 3:
                    # Print out BED line with chromosome name and new start and end positions
                    # along with remaining columns
                    toprint = [chrom, newstartpos, newendpos] + tmp[3:]
                    print('\t'.join(toprint))
                else:
                    print('Too few columns, please check input file is in BED format.')
                    exit(1)


def main():
    """Driver function."""
    args = parse_args()
    f = set_format(args)
    # Then, use the format and the supplied file to convert the coordinates
    if f == 'VCF':
        vcf_conv(args.intervals, PARTS_SIZES)
    elif f == 'BED':
        bed_conv(args.intervals, PARTS_SIZES)
    else:
        print('Unrecognized format, weird.')
        exit(1)


main() # Run the program
