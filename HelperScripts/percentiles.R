#!/usr/bin/env Rscript

#   This script writes a percentiles table of the GQ and DP per sample from a VCF file
#   Usage: Rscript percentiles.R <GQ matrix> <DP per sample matrix> <output_directory> <project_name> <filtered_status>

#   If running this on MSI and there are problems with the installation above, please see here for
#   instructions to install this package in your home directory: https://www.msi.umn.edu/sw/r
.libPaths() # Check, remove after debugging
library(data.table)

#   A function to run the script
runScript <- function() {
    options(warn=2)
    args <- commandArgs(trailingOnly = TRUE)
    #   Set the arguments
    GQ_matrix <- args[1]
    DP_matrix <- args[2]
    out <- args[3]
    project <- args[4]
    status <- args[5]
    #   Make the GQ table
    GQ_path <- paste0(out, "/", project, "_", status, "_GQ.txt")
    #GQ <- read.table(file=GQ_matrix, sep="\t", na.strings=".")
    GQ <- fread(file=GQ_matrix, sep="\t", na.strings=".") # Faster way to read in file
    #GQ_quantile <- quantile(as.numeric(as.matrix(as.vector(GQ))),probs=seq(0,1,0.01),na.rm = TRUE)
    GQ_quantile <- quantile(GQ$V1, probs=seq(0,1,0.01))
    write.table(x=GQ_quantile, file=GQ_path, sep="\t", quote=FALSE, col.names=FALSE)
    # Remove GQ to free up memory
    rm(GQ, GQ_quantile)

    #   Make the DP table
    DP_path <- paste0(out, "/", project, "_", status, "_DP_per_sample.txt")
    #DP <- read.table(file=DP_matrix, sep="\t", na.strings="NA")
    DP <- fread(file=DP_matrix, sep="\t", na.strings="NA") # Faster way to read in file
    #DP_quantile <- quantile(as.numeric(as.matrix(as.vector(DP))),probs=seq(0,1,0.01),na.rm = TRUE)
    DP_quantile <- quantile(DP$V1, probs=seq(0,1,0.01))
    write.table(x=DP_quantile, file=DP_path, sep="\t", quote=FALSE, col.names=FALSE)
}

runScript() # Run program
