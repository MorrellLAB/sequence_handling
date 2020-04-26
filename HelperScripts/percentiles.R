#!/usr/bin/env Rscript

#   This script writes a percentiles table of the GQ and DP per sample from a VCF file
#   Usage: Rscript percentiles.R <GQ matrix> <DP per sample matrix> <output_directory> <project_name> <filtered_status>

#   If running this on MSI and there are problems with the installation above, please see here for
#   instructions to install this package in your home directory: https://www.msi.umn.edu/sw/r
#library(data.table)
#library(bit) # Required by ff
#library(ff)
#library(ffbase)
library(bigmemory)

#   A function to run the script
runScript <- function() {
    options(warn=2)
    args <- commandArgs(trailingOnly = TRUE)
    #   Set the arguments
    # Use absolute filepath
    GQ_matrix <- args[1] # Note: bigmemory will return error if "~" is in filepath
    DP_matrix <- args[2] # Note: bigmemory will return error if "~" is in filepath
    out <- args[3]
    project <- args[4]
    status <- args[5]
    #   Make the GQ table
    GQ_path <- paste0(out, "/", project, "_", status, "_GQ.txt")
    #GQ <- read.table(file=GQ_matrix, sep="\t", na.strings=".")
    #GQ <- fread(file=GQ_matrix, sep="\t", na.strings=".") # Faster way to read in file
    #GQ <- read.table.ffdf(file=GQ_matrix, sep="\t", na.string=".")
    GQ <- read.big.matrix(GQ_matrix, sep="\t", type="integer") # Better handles large files
    #GQ_quantile <- quantile(as.numeric(as.matrix(as.vector(GQ))),probs=seq(0,1,0.01),na.rm = TRUE)
    GQ_quantile <- quantile(GQ[, 1], probs=seq(0,1,0.01))
    write.table(x=GQ_quantile, file=GQ_path, sep="\t", quote=FALSE, col.names=FALSE)
    # Remove GQ to free up memory
    rm(GQ, GQ_quantile)

    #   Make the DP table
    DP_path <- paste0(out, "/", project, "_", status, "_DP_per_sample.txt")
    #DP <- read.table(file=DP_matrix, sep="\t", na.strings="NA")
    #DP <- fread(file=DP_matrix, sep="\t", na.strings="NA") # Faster way to read in file
    #DP <- read.table.ffdf(file=DP_matrix, sep="\t", na.strings="NA")
    DP <- read.big.matrix(file=DP_matrix, sep="\t", type="integer") # Better handles large files
    #DP_quantile <- quantile(as.numeric(as.matrix(as.vector(DP))),probs=seq(0,1,0.01),na.rm = TRUE)
    DP_quantile <- quantile(DP[, 1], probs=seq(0,1,0.01))
    write.table(x=DP_quantile, file=DP_path, sep="\t", quote=FALSE, col.names=FALSE)
}

runScript() # Run program
