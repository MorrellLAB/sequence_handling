#!/usr/bin/env Rscript

# This script visualizes the amount of missingness present in your VCF file

library(ggplot2)
library(gridExtra)

# Functions
read_indv_miss <- function(filename) {
  imiss <- read.delim(filename, sep="\t")
  return(imiss)
}

read_lmiss <- function(filename) {
  lmiss <- read.delim(filename, sep="\t")
  return(lmiss)
}

plot_imiss <- function(imiss_df) {
  outplot <- ggplot(imiss_df, aes(F_MISS)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light() +
    ggtitle("Missing data per individual") +
    theme(plot.title = element_text(hjust = 0.5)
  return(outplot)
}

plot_lmiss <- function(lmiss_df) {
  outplot <- ggplot(lmiss_df, aes(F_MISS)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    theme_light() +
    ggtitle("Missing data per site") +
    theme(plot.title = element_text(hjust = 0.5)
  return(outplot)
}

# Driver function
main <- function() {
  # User provided arguments
  args <- commandArgs(trailingOnly=TRUE)
  imiss_fp <- args[1]
  lmiss_fp <- args[2]
  out_dir <- args[3]
  # Read in data
  imiss_df <- read_indv_miss(imiss_fp)
  lmiss_df <- read_lmiss(lmiss_fp)
  # Plot missingness
  imiss_plot <- plot_imiss(imiss_df)
  lmiss_plot <- plot_lmiss(lmiss_df)
  # Arrange plots into multi-panel figure
  combo_plot <- arrangeGrob(imiss_plot, lmiss_plot, nrow=2)
  # Save to file
  outfp <- paste0(out_dir, "/missingness_plots.png")
  ggsave(outfp, plot=combo_plot, device="png")
}

main() # Run the program
