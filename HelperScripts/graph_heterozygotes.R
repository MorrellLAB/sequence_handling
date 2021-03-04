#!/usr/bin/env Rscript

# This script visualizes the amount of heterozygous sites present in your VCF file
#   This script works with the GATK VariantsToTable function output.
#   See more details about what each annotation means here:
#     https://gatk.broadinstitute.org/hc/en-us/articles/360036732151-VariantsToTable
#     https://gatk.broadinstitute.org/hc/en-us/articles/360037261311-ExcessHet
#     https://gatk.broadinstitute.org/hc/en-us/articles/360036366852-InbreedingCoeff

library(ggplot2)
library(gridExtra)
library(egg)
library(tidyverse)

# Functions
read_variants_table <- function(filename) {
  variants_table <- read.delim(filename, sep="\t", header=TRUE)
  # Calculate proportion of heterozygous sites, save to new column
  variants_table$het_prop <- variants_table$HET / variants_table$NCALLED
  return(variants_table)
}

plot_het_prop <- function(vtdf) {
  # Histogram
  outplot1 <- ggplot(vtdf, aes(het_prop)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins = 100) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Het Proportion (HET / NCALLED)") +
    facet_wrap(~TYPE)
  
  # Boxplot
  outplot2 <- ggplot(vtdf, aes(x="", y=het_prop)) +
    geom_boxplot(fill = "dodgerblue1", alpha = 0.6, outlier.colour = "red", outlier.alpha = 0.1) +
    coord_flip() +
    theme_light() +
    xlab("") +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab("Het Proportion (HET / NCALLED)") + # ylab because of coord_flip
    ggtitle("Boxplot (Red = Outliers)") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    facet_wrap(~TYPE)
    
  # Arrange plots into one
  outplot <- egg::ggarrange(outplot1, outplot2, heights = 2:1, widths = 10)
  return(outplot)
}

plot_excess_het <- function(vtdf) {
  # Histogram
  outplot1 <- ggplot(vtdf, aes(ExcessHet)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins = 100) +
    theme_light() +
    xlab("Excess Heterozygosity") +
    # Add text to top right corner
    annotate("text", x=Inf, y=Inf, 
             label = "*Caveats: Only for diploid samples and assumes\nthe samples are unrelated. If the samples are\nrelated, a pedigree file should be provided.", 
             vjust=1.2, hjust=1.05, size=2)  +
    facet_wrap(~TYPE)
  
  # Boxplot
  outplot2 <- ggplot(vtdf, aes(x="", y=ExcessHet)) +
    geom_boxplot(fill = "dodgerblue1", alpha = 0.6, outlier.colour = "red", outlier.alpha = 0.1) +
    coord_flip() +
    theme_light() +
    ylab("Excess Heterozygosity") +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("Boxplot (Red = Outliers)") +
    theme(plot.title = element_text(size = 8, hjust = 0.5))  +
    facet_wrap(~TYPE)
  
  # Arrange plots into one
  outplot <- egg::ggarrange(outplot1, outplot2, heights = 2:1, widths = 10)
  return(outplot)
}

plot_inbreeding_coeff <- function(vtdf) {
  # Density plot
  outplot1 <- ggplot(vtdf, aes(InbreedingCoeff)) +
    geom_density(fill=alpha('dodgerblue1', 0.3), na.rm=TRUE) +
    theme_light() +
    xlab("Inbreeding Coefficient") +
    annotate("text", x=Inf, y=Inf,
             label = "Note: Calculated as F = 1 - (# of het genotypes)/(# of het genotypes expected under HWE)",
             vjust=1.2, hjust=1.05, size=2) +
    facet_wrap(~TYPE)
  
  # Boxplot
  outplot2 <- ggplot(vtdf, aes(x="", InbreedingCoeff)) +
    geom_boxplot(fill = "dodgerblue1", alpha = 0.6, outlier.colour = "red", outlier.alpha = 0.1, na.rm=TRUE) +
    coord_flip() +
    theme_light() +
    ylab("Inbreeding Coefficient") +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("Boxplot (Red = Outliers)") +
    theme(plot.title = element_text(size = 8, hjust = 0.5))  +
    facet_wrap(~TYPE)
  
  # Arrange plots into one
  outplot <- egg::ggarrange(outplot1, outplot2, heights = 2:1, widths = 4)
  return(outplot)
}

# Driver function
main <- function() {
  # User provided arguments
  args <- commandArgs(trailingOnly=TRUE)
  variants_table_fp <- args[1]
  out_dir <- args[2]
  # Read in data
  vtdf <- read_variants_table(variants_table_fp)
  # Plot amount of heterozygous sites
  het_prop_plot <- plot_het_prop(vtdf)
  excesshet_plot <- plot_excess_het(vtdf)
  # Check if we were able to calculate inbreeding coefficients
  if ("InbreedingCoeff" %in% colnames(vtdf)) {
    # Make plot of inbreeding coefficients
    inbCoef_plot <- plot_inbreeding_coeff(vtdf)
    # Arrange all plots into multi-panel figure
    combo_plot <- arrangeGrob(het_prop_plot, excesshet_plot, inbCoef_plot, ncol=2, nrow=2)
    # Prep outfile name
    outfp <- paste0(out_dir, "/Het-ExcessHet-InbreedingCoef_plots.png")
  } else {
    # Only Het Prop and Excess Het plots generated
    # Arrange all plots into multi-panel figure
    combo_plot <- arrangeGrob(het_prop_plot, excesshet_plot, nrow=2)
    # Prep outfile name
    outfp <- paste0(out_dir, "/Het-ExcessHet_plots.png")
  }
  # Save to file
  ggsave(outfp, plot=combo_plot, device="png")
}

main() # Run the program
