#!/usr/bin/env python3

#   This code is based on code originally written by Tom Kono
#   Used to be named HeterozogotesVcfFilter.py

#   A script to apply various arbitrary filters to a VCF after Variant Recalibrator
#       - minimum # reads
#       - maximum # reads
#       - read balance in heterozygotes
#       - GQ
#       - DP
#   Some versions of GATK may produce incompatible output for this script
#   This script writes the filtered VCF lines to standard output

import sys

#   Minimum number of reads needed to support a genotype
mindp = float(sys.argv[2])
#   Maximum number of reads allowed to support a genotype (too many = gene duplication problems)
maxdp = float(sys.argv[3])
#   The maximum percent deviation from 50/50 ref/alt reads allowed 
#   For example, mindev = 0.1 allows 60/40 ref/alt and also 40/60 ref/alt but not 70/30 ref/alt reads
mindev = float(sys.argv[4])
#   Minimum genotyping quality (GQ - stored as a PHRED-scaled probability)
gt_cutoff = float(sys.argv[5])
#   Minimum number of reads
per_sample_coverage_cutoff = float(sys.argv[6])

#   Read the file in line-by-line
with open(sys.argv[1]) as f:
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            tmp = line.strip().split('\t')
            #   We don't want indels or multiple alleles
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue #  Skip those lines
            else:
                #   We want to keep the PASS SNPs
                if tmp[6] == 'PASS':
                    genotypes = tmp[9:]
                    #   enumerate is the function to iterate the index and the values for the list:
                    for geno_index, s in enumerate(genotypes):
                        #   GT:AD:DP:GQ:PL
                        gt_metadata = s.split(':')
                        gt = gt_metadata[0]
                        dp = gt_metadata[2]
                        ad = gt_metadata[1].split(',')
                        gq = gt_metadata[3]
                        if  dp == '.' or gq == '.' :
                            continue    #   If there is no depth or genotyping quality info, don't change anything
                        else:
                            #   Apply filters for DP and GQ
                            if int(dp) < per_sample_coverage_cutoff or int(gq) < gt_cutoff:
                                tmp[9+geno_index] = ':'.join(['./.'] + tmp[9+geno_index].split(':')[1:])
                            else:                           
                                if gt == '0/1':
                                    ref = float(ad[0])
                                    alt = float(ad[1]) 
                                    if ref+alt !=0:
                                        balance = ref/(ref+alt)                         
                                        if dp != '.':
                                            #   Apply filters for balance and max and min depth
                                            if int(dp) < mindp or int(dp) > maxdp or abs(0.5 - balance) > mindev:
                                                 tmp[9+geno_index] = ':'.join(['./.'] + tmp[9+geno_index].split(':')[1:])
                                    else:
                                        tmp[9+geno_index] = ':'.join(['./.'] + tmp[9+geno_index].split(':')[1:])
                sys.stdout.write('\t'.join(tmp) + '\n')