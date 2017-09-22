#!/usr/bin/env python

#   This code is based on code originally written by Tom Kono
#   Used to be named Filter_VCF.py

#   A script to apply various arbitrary filters to a VCF
#       Filters out indels and sites with more than 2 alleles
#       Filters out sites with GT:AD:DP:GQ:PGT:PID:PL
#       If the quality score is missing or the site quality score is too low, filters out the site
#       If too many samples are het, missing, low quality, or low depth, filters out the site

#   Some versions of GATK may produce incompatible output for this script
#   This script writes the filtered VCF lines to standard output

import sys

#   The QUAL cutoff
quality_cutoff = sys.argv[2]
#   The most # of samples that can be heterozygous and still keep the site
het_cutoff = sys.argv[3]
#   The most # of samples that can be missing and still keep the site
missing_cutoff = sys.argv[4]
#   The genotype quality cutoff
gt_cutoff = sys.argv[5]
#   The most # of samples that can have low GQ and still keep the site
n_gt_cutoff = sys.argv[6]
#   Our coverage cutoff
per_sample_coverage_cutoff = sys.argv[7]
#   The most # of samples that can have low DP and still keep the site
n_low_coverage_cutoff = sys.argv[8]

#   Read the file in line-by-line
with open(sys.argv[1]) as f:
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            tmp = line.strip().split('\t')
            #   If indel or more than 2 alleles, skip the site - don't write it
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            format = tmp[8]
            sample_information = tmp[9:]
            #   Start counting up # of hets, # low quality genotypes, and # low coverage sites
            nhet = 0
            low_qual_gt = 0
            low_coverage = 0
            missing_data = 0
            for s in sample_information:
                #   For the GATK HaplotypeCaller, the per-sample information follows the form GT:AD:DP:GQ:PL
                info = s.split(':')
                if len(info) != 5:
                    break
                gt = info[0]
                #   We have to check for missing data first, because if it is missing, then the other fields are not filled in
                if '.' in gt:
                    missing_data += 1
                else:
                    dp = info[2]
                    gq = info[3]
                    #   Split on / (we have unphased genotypes) and count how many hets we have
                    if len(set(gt.split('/'))) > 1:
                        nhet += 1
                    if dp == '.' or int(dp) < per_sample_coverage_cutoff:
                        low_coverage += 1
                    if gq == '.' or int(gq) < gt_cutoff:
                        low_qual_gt += 1
            #   The quality score is the sixth element
            #   If the quality score is missing or the site quality score is too low or any of the counts are above the cutoff, don't print the site
            if tmp[5] == '.' or float(tmp[5]) < quality_cutoff or nhet > het_cutoff or low_qual_gt > n_gt_cutoff or low_coverage > n_low_coverage_cutoff or missing_data > missing_cutoff:
                continue
            else:
                sys.stdout.write(line)
