#!/usr/bin/env Rscript

#   This script writes a percentiles table of the GQ and DP per sample from a VCF file
#   Usage: Rscript percentiles.R <GQ matrix> <DP per sample matrix> <output_directory> <project_name> <filtered_status>

#   A function to run the script
runScript <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    #   Set the arguments
    GQ_matrix <- args[1]
    DP_matrix <- args[2]
    out <- args[3]
    project <- args[4]
    status <- args[5]
    #   Make the GQ table
    GQ_path <- paste0(out, "/", project, "_", status, "_GQ.txt")
    GQ <- read.table(file=GQ_matrix, sep="\t", na.strings=".")
    GQ_quantile <- quantile(as.numeric(as.matrix(as.vector(GQ))),probs=seq(0,1,0.01),na.rm = TRUE)
    write.table(x=GQ_quantile,file=GQ_path,sep="\t", quote=FALSE, col.names=FALSE)
    #   Make the DP table
    DP_path <- paste0(out, "/", project, "_", status, "_DP_per_sample.txt")
    DP <- read.table(file=DP_matrix, sep="\t", na.strings="NA")
    DP_quantile <- quantile(as.numeric(as.matrix(as.vector(DP))),probs=seq(0,1,0.01),na.rm = TRUE)
    write.table(x=DP_quantile,file=DP_path,sep="\t", quote=FALSE, col.names=FALSE)
}

runScript() # Run program