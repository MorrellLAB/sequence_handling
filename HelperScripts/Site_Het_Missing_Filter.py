#!/usr/bin/env python
"""Filter a VCF for sites that have a proportion of heterozygous genotypes or
proportion of missing calls above given thresholds. Takes six arguments:
    1) Input VCF (.vcf or .vcf.gz)
    2) Maximum proportion of heterozygous genotypes per site
    3) Maximum proportion of missing calls per site
    4) Quality
    5) Genotype Quality (GQ)
    6) DP per sample (DP)

Note: Arguments 4-6 are for the purposes of adding a complete header line for
the Create_HC_Subset handler and automated checks within the handler. These
arguments were used in the filtering step prior to this script.
"""

# This is a modified version of code written by Tom Kono (https://github.com/MorrellLAB/Selective_Sweeps/blob/master/VCF_Filtering/Site_Het_Missing_Filter.py)
# The modifications are:
#   1) allows the script to handle both uncompressed and compressed vcf files
#   2) handles hets split by "/" and "|"
#   3) Adds header line to filtered VCF containing Create_HC_Subset handler cutoffs used

import sys
import gzip

try:
    vcf_filename = sys.argv[1]
    max_het_cutoff = float(sys.argv[2])
    max_miss_cutoff = float(sys.argv[3])
    assert max_het_cutoff >= 0 and max_het_cutoff <= 1
    assert max_miss_cutoff >= 0 and max_miss_cutoff <= 1
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    exit(1)
except ValueError:
    sys.stderr.write('Please provide positive floating point numbers for the '
                     'proportion of heterozygous and missing calls. '
                     'Proportions must be between 0 and 1.\n')
    exit(1)
except AssertionError:
    sys.stderr.write('Please provide positive foating point numbers for the '
                     'proportion of heterozygous and missing calls. '
                     'Proportions must be between 0 and 1.\n')
    exit(1)

def filter_sites(f, max_het, max_miss, qual, gq, dp):
    """Read in VCF file and filter vcf for sites that have a proportion
    of heterozygous genotypes or proportion of missing calls above a given
    threshold."""
    for line in f:
        if line.startswith('##'):
            print(line.strip())
        elif line.startswith('#CHROM'):
            # Add in line to track cutoffs used for filtering
            handler_line = "##Create_HC_Subset_filter_cutoffs=" + "Quality:" + qual + ",Max_het:" + str(max_het) + ",Max_miss:" + str(max_miss) + ",Genotype_Quality:" + gq + ",DP_per_sample:" + dp
            print(handler_line)
            # Write out #CHROM line
            print(line.strip())
        else:
            tmp = line.strip().split('\t')
            # Figure out which column has the genotype in it
            fmt = tmp[8].split(':')
            gt_col = fmt.index('GT')
            # Count up the number of heterozygous and missing genotypes
            gts = [x.split(':')[gt_col] for x in tmp[9:]]
            # Accept / or | and count how many hets we have
            n_het = gts.count('0/1') + gts.count('1/0') + gts.count('0|1') + gts.count('1|0')
            n_miss = gts.count('.') + gts.count('./.')
            # Convert them to proportions
            p_het = n_het / len(gts)
            p_miss = n_miss / len(gts)
            # Skip the site if the het or miss proportion exceeds the given
            # values.
            if p_het > max_het or p_miss > max_miss:
                continue
            else:
                print(line.strip())


if '.gz' in vcf_filename:
    with gzip.open(vcf_filename, 'rt') as vf:
        filter_sites(vf, max_het_cutoff, max_miss_cutoff, str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]))
else:
    with open(vcf_filename, 'rt') as vf:
        filter_sites(vf, max_het_cutoff, max_miss_cutoff, str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]))
