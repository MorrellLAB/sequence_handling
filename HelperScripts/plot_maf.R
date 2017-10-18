#!/usr/bin/env Rscript

#   This script plots a MAF histogram from a text file
#   Usage: Rscript plot_maf.R <input_maf.txt> <vcf name> <output_maf.pdf>

#   A function to run the script
runScript <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    input <- args[1]
    name <- args[2]
    output <- args[3]
    title <- paste0(name, " Minor Allele Frequency Distribution")
    maf_table <- read.table(file=input, sep="\t", header=TRUE)
    pdf(output)
    hist(maf_table$MAF, breaks=50, main=title, xlab="Frequency", ylab="Number of Minor Alleles", xlim=c(0,0.5))
    dev.off()
}

runScript() # Run program