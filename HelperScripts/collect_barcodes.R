#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

#   A function to read in data from the keyfile
readKeyfile <- function(filename) {
    key <- read.table(file = filename, header = TRUE, as.is = TRUE, sep = '\t')
    return(key)
}

#   A function to find all possible barcode lengths
collectLengths <- function(barcodes) {
    allLengths <- sapply(X = barcodes, FUN = nchar) # Get all lengths
    uniqueLengths <- unique(x = allLengths) # Get the unique ones
    return(uniqueLengths)
}

#   A function to write the barcodes of the same length to a file
writeBarcodeFile <- function(outdir, project, barcodes, length) {
    theseBarcodes <- barcodes[which(x = nchar(barcodes) == length)]
    outfile <- paste0(outdir, '/', project, '_barcodes_length_', length, '.barcode')
    for(i in seq(length(x = theseBarcodes))) {
        cat('BC', i, '\t', theseBarcodes[i], '\n', file = outfile, sep = '', append = TRUE)
    }
}

main <- function() {
    #   Collect the arguments
    keyfile <- args[1]
    outdirectory <- args[2]
    project <- args[3]
    #   Read in the keyfile
    key <- readKeyfile(filename = keyfile)
    #   Barcodes are always column 3
    barcodes <- as.vector(x = key[[3]])
    #   Get the lengths
    lengths <- collectLengths(barcodes = barcodes)
    #   Write the files
    sapply(X = lengths, FUN = writeBarcodeFile, outdir = outdirectory, project = project, barcodes = barcodes)
}

main()
