#!/usr/bin/env python3

#   Original code written by Tom Kono
#   A python program to read in a VCF file and output a Hudson table-like
#   format. Not a "true" Hudson table, since hets are not handled properly.
#   Usage:
#       VCF_To_Htable.py [VCF file] > [Htable.txt]

import sys
import gzip

# User provided input argument
vcf_fp = sys.argv[1]

#   If the "minor genotype frequency" falls below this threshhold, then we
#   omit the site.
MAFThreshhold = 0.00

#   A function to calculate the minor allele frequency
def MAF(x):
    #   get the set of genotypes in the list
    genotypes = set(x)
    #   start counting up the frequencies
    freqs = []
    for g in genotypes:
        freqs.append(x.count(g)/float(len(x)))
    return min(freqs)


def vcf_to_htable(vcf_fp, f, loci, g_matrix):
    """Convert vcf to htable format."""
    for index, line in enumerate(f):
        if index%10000 == 0:
            sys.stderr.write('Read ' + str(index) + ' sites.\n')
        if line.startswith('##'):
            continue
        #   #CHROM line is the one that contains the sample information
        elif line.startswith('#CHROM'):
            #   Split up the line on tabs
            tmp = line.strip().split('\t')
            #   Look for the 'FORMAT' field
            format_field = tmp.index('FORMAT')
            #   And get the samples out of the list
            #   Add 1 to the format_field because we don't actually want to
            #   include 'FORMAT' in the sample info
            samples = tmp[format_field + 1:]
            #   Write a little diagnostic message
            sys.stderr.write(vcf_fp + ' has ' + str(len(samples)) + ' samples.\n')
        #   Now that we have the number and names of the samples, we print the
        #   genotype data
        else:
            #   we split the line as above
            tmp = line.strip().split('\t')
            #   assign the variables for clarity
            scaffold = tmp[0]
            pos = tmp[1]
            var_id = tmp[2]
            ref_allele = tmp[3]
            #   If there are multiple alternate alleles, they are listed with a comma between them
            alt_alleles = tmp[4].split(',')
            qual = tmp[5]
            filt = tmp[6]
            info = tmp[7]
            format = tmp[8]
            #   And then the genotypes
            genotypes = tmp[9:]
            #   The locus will be the scaffold number and then the bp position
            #locus = scaffold + '_' + pos
            locus = pos
            #   then we parse the genotypes
            #   we create a list here that will be a single column of the genotype matrix
            g_column = []
            for g in genotypes:
                #   In the genotype string, the first element (separated by :) is the actual genotype call
                call = g.split(':')[0]
                #   These are diploid calls, and we are assuming they are unphased
                #   the are listed in the form allele1/allele2
                #   with 0 = ref, 1 = alt1, 2 = alt2, and so on...
                alleles = call.split('/')
                individual_call = ''
                for x in alleles:
                    if x == '.':
                        individual_call += 'N'
                    else:
                        #   cast to integer so we can use it in slicing a list
                        c = int(x)
                        #   if it's 0, we just tack on the reference state
                        if c == 0:
                            individual_call += ref_allele
                        else:
                            #   Otherwise, we use it to slice the list of alternate alleles
                            individual_call += alt_alleles[c-1]
                #   Then append the individual call to the column for the genotype matrix
                g_column.append(individual_call)
            #   Then, append that column to the genotype matrix
            #   If there is no variation in genotype calls (that is, all lines have the same genotype)
            #   then we don't care about it
            unique_calls = set(g_column)
            if len(unique_calls) <= 1:
                continue
            else:
                if MAF(g_column) > MAFThreshhold:
                    g_matrix.append(g_column)
                    loci.append(locus)
                else:
                    continue


#   Empty lists for the genotype matrix and the loci
loci = []
g_matrix = []
#   start reading through the file. We will skip any lines that start wtih '##'
if "gz" in vcf_fp:
    with gzip.open(vcf_fp, 'rt') as f:
        vcf_to_htable(vcf_fp, f, loci, g_matrix)
else:
    with open(vcf_fp, 'r') as f:
        vcf_to_htable(vcf_fp, f, loci, g_matrix)

#   Now, we have to transpose the genotype matrix
g_matrix_t = zip(*g_matrix)
#   print the number of samples and the number of loci
print(str(len(samples)) + '\t' + str(len(loci)))
#   print the loci
print('\t' + '\t'.join(loci))
#   Print the line for unknown ancestral state
print('anc\t' + '?\t'*(len(loci)-1) + '?')
#   then print the transposed genotype matrix
for index, g in enumerate(g_matrix_t):
    print(samples[index] + '\t' + '\t'.join(g))
