#!/usr/bin/env Rscript

# bedtools genomecov function histogram output produces a file with 5 columns
# see http://www.quinlanlab.org/pdf/bedtools.protocols.pdf
# pages 11.12.9 - 11.12.10
# 1) chromosome
# 2) depth
# 3) number of base pairs with depth = column2
# 4) size of the chromosome
# 5) fraction of base pairs with depth = column2

#   To take command line arguments
args <- commandArgs(TRUE)
#   This creates a vector of character strings for arguments
#   we will just take two arguments here, the stats directory
#   and the sample name
genomeinput <- args[1]
exoninput <- args[2]
geneinput <- args[3]
outdir <- args[4]
samplename <- args[5]

#   Read the three tables
cov <- read.table(genomeinput, header = FALSE)
ecov <- read.table(exoninput, header = FALSE)
gcov <- read.table(geneinput, header = FALSE)

#   Create a list of the types of graphs we're making
typelist <- list('genome', 'exons', 'genes')

#   A label for the x-axis
lab <- 'Depth'
#   Three labels for the y-axis
frac <- paste('Fraction of', typelist, 'at depth')

#   What are we calling each of our output plots?
output <- paste0(outdir, samplename, '_outplot_', typelist, '.pdf')

#   Make the first graph: genome coverage
pdf(file = output[1], width = 6, height = 6)
options(scipen=5)
plot(cov[1:51, 2], cov[1:51, 5], type='h', col='blue', lwd=5, xlab = lab, ylab = frac[1], xlim = c(0, max(cov[1:51, 2])))
dev.off()

#   Make the second graph: exon coverage
pdf(file = output[2], width = 6, height = 6)
options(scipen=5)
plot(ecov[1:101, 10], ecov[1:101, 13], type = 'h', col = 'green' , lwd = 5, xlab = lab, ylab = frac[2], xlim = c(0, max(ecov[1:101, 10])))
dev.off()

#   Make the third graph: gene coverage
pdf(file = output[3], width = 6, height = 6)
options(scipen=5)
plot(gcov[1:101, 10], gcov[1:101, 13], type = 'h', col = 'red' , lwd = 5, xlab = lab, ylab = frac[3], xlim = c(0, max(gcov[1:101, 10])))
dev.off()