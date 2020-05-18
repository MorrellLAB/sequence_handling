#!/usr/bin/env Rscript

#   This script writes a percentiles table of the GQ and DP per sample from a VCF file
#   Usage: Rscript percentiles.R <GQ matrix> <DP per sample matrix> <output_directory> <project_name> <filtered_status>

#   If running this on MSI and there are problems with the installation, please see here for
#   instructions to install this package in your home directory: https://www.msi.umn.edu/sw/r
library(bigmemory)

#   A function to run the script
runGQ <- function(GQ_matrix, GQ_path) {
    #   Make the GQ table
    print("Reading in GQ matrix...")
    GQ <- read.big.matrix(GQ_matrix, sep="\t", type="integer") # Better handles large files
    print("Finished reading in GQ matrix.")
    GQ_quantile <- quantile(GQ[, 1], probs=seq(0,1,0.01))
    print("Finished calculating percentiles for GQ matrix.")
    write.table(x=GQ_quantile, file=GQ_path, sep="\t", quote=FALSE, col.names=FALSE)
    print("Finished saving GQ percentiles table.")
    # Remove GQ to free up memory
    rm(GQ, GQ_quantile)
}

runDP <- function(DP_matrix, DP_path) {
    #   Make the DP table
    print("Reading in DP matrix...")
    DP <- read.big.matrix(file=DP_matrix, sep="\t", type="integer") # Better handles large files
    print("Finished reading in DP matrix.")
    DP_quantile <- quantile(DP[, 1], probs=seq(0,1,0.01))
    print("Finished calculating percentiles for DP matrix.")
    write.table(x=DP_quantile, file=DP_path, sep="\t", quote=FALSE, col.names=FALSE)
    print("Finished saving DP percentiles table.")
}

runScript <- function() {
    # Driver function
    options(warn=2)
    args <- commandArgs(trailingOnly = TRUE)
    #   Set the arguments
    # Use absolute filepath
    GQ_matrix <- args[1] # Note: bigmemory will return error if "~" is in filepath
    DP_matrix <- args[2] # Note: bigmemory will return error if "~" is in filepath
    out <- args[3]
    project <- args[4]
    status <- args[5]
    # Prepare output filepaths
    GQ_path <- paste0(out, "/", project, "_", status, "_GQ.txt")
    DP_path <- paste0(out, "/", project, "_", status, "_DP_per_sample.txt")
    Log_file_path <- paste0(out, "/", project, "_", status, ".log")

    # Get file sizes
    # If the files are larger than 160G (equivalent to 171798691840 bytes), we will run into memory
    #   issues when calculating percentiles. Please use another method to identify the appropriate
    #   threshold to use.
    GQ_file_size <- file.info(GQ_matrix)$size
    DP_file_size <- file.info(DP_matrix)$size

    # Save all output messages to a log file
    sink(file = Log_file_path) # Start sinking (start writing output messages to file)

    # Process GQ matrix if file does not exist already
    if (file.exists(GQ_path)) {
        print("GQ percentiles table exists, proceeding to next step.")
    } else {
        # File does not exist, process GQ matrix if file is not too large
        if (GQ_file_size < 171798691840) {
            runGQ(GQ_matrix, GQ_path)
        } else {
            print("GQ file is larger than 160GB and is too large to process. Please use another approach to determine appropriate GQ cutoff. Proceeding to next step.")
        }
    }

    # Process DP matrix
    if (file.exists(DP_path)) {
        print("DP percentiles table exists, we are done generating percentiles tables.")
    } else {
        # File does not exist, process DP matrix
        if (DP_file_size < 171798691840) {
            runDP(DP_matrix, DP_path)
        } else {
            print("DP file is larger than 160GB and is too large to process. Please use another approach to determine appropriate DP cutoff. Proceeding to next step.")
        }
    }

    sink() # Stop sinking (stop writing output to file)
}

runScript() # Run program
