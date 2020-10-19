#!/usr/bin/env Rscript

### A script to plot distributions of Variant Annotations
### Adapted material from: https://evodify.com/gatk-in-non-model-organism/

library(ggplot2)
library(gridExtra)

options(warn=1)

##########################
######## Setup_df ########
##########################

Setup_DF <- function(Raw_table, HC_table) {
  raw_variants <- read.csv(Raw_table, 
                           header = T, na.strings=c("","NA"), sep = "\t")
  hc_variants <- read.csv(HC_table, 
                           header = T, na.strings=c("","NA"), sep = "\t")
  raw_variants$Chrom_Pos <- paste0(raw_variants$CHROM, sep = "_", raw_variants$POS, sep = "_", raw_variants$TYPE)
  hc_variants$Chrom_Pos <- paste0(hc_variants$CHROM, sep = "_", hc_variants$POS, sep = "_", hc_variants$TYPE)
  raw_variants$Cat <- ifelse(raw_variants$Chrom_Pos %in% hc_variants$Chrom_Pos, "Passed_filtering", "Failed_filtering")
  types <- split(raw_variants, raw_variants$TYPE)
  return (types)
}

##########################
# SETUP LABELS + CUTOFFS #
##########################

Cutoffs <- list("QD" = 2, "FS" = 60, 
                "SOR" = 3, "MQ" = 40, 
                "MQRankSum" = -12.5, 
                "ReadPosRankSum" = -8.0, "QUAL" = 40, 
                "DP" = 5)

X_labels <- list("QD" = c("Quality by Depth"), "FS" = c("Fisher Strand"), 
                "SOR" = c("Strand Odds Ratio"), "MQ" = c("Mapping Quality"),
                "MQRankSum" = c("Mapping Quality Rank Sum Test"), 
                "ReadPosRankSum" = c("Read Position Rank Sum Test"), "QUAL" = c("Quality Score"), 
                "DP" = c("Depth"))

X_limits <- list("QD" = c(0,40), "FS" = c(0,200), 
                 "SOR" = c(0,20), "MQ" = c(0,100),
                 "MQRankSum" = c(-12,12), 
                 "ReadPosRankSum" = c(-8,8), "QUAL" = c(0,1000), 
                 "DP" = c(0,10000))

# might want to change MQRankSum (-5,5) and ReadPosRankSum limits (-3,3) ?

#########################
#### FUNCTION TO PLOT ###
#########################

Ann_plot <- function(dataset, vector, lab, cutoff, xrange) {
  plot <- ggplot(dataset, aes(x=vector)) + geom_density(aes(fill=Cat)) + 
    scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
    geom_vline(xintercept=cutoff, size=0.7, lty=2) +
    xlab(lab) +
    xlim(xrange[1], xrange[2]) +
    theme(legend.title=element_blank())
  return (plot)
}

#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Directory <- args[1]
RawVariants <- args[2]
HCVariants <- args[3]

Mydf <- Setup_DF(RawVariants, HCVariants)

#########################
####### PLOT SNPs #######
#########################

SNP_plots <- lapply(names(Mydf$SNP[,4:11]), 
                function(x) {Ann_plot(Mydf$SNP, Mydf$SNP[[x]], X_labels[[x]], Cutoffs[[x]], X_limits[[x]])})
names(SNP_plots) <- names(Mydf$SNP[,4:11])


### Fisher's Strand on log scale
# add 1 so zero's are included
SNP_FS <- ggplot(Mydf$SNP, aes(x=FS + 1)) + geom_density(aes(fill=Cat)) +
  scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
  geom_vline(xintercept=60, size=0.7, lty=2) + 
  scale_x_log10() +
  theme(legend.title=element_blank()) +
  xlab("Fisher Strand (log-scaled)")

# Build output filepath
outfile_snp_fp <- paste0(Directory, "/Percentile_Tables/SNP_distributions.png")
SNP_plot <- grid.arrange(SNP_plots$QD, SNP_FS, SNP_plots$SOR, SNP_plots$MQ, SNP_plots$MQRankSum, SNP_plots$ReadPosRankSum, 
             SNP_plots$DP, SNP_plots$QUAL, nrow=4)
ggsave(outfile_snp_fp,
       plot = SNP_plot,
       device = "png")

print("Finished making SNP plot")


#########################
##### PLOT INDELS #######
#########################

Indel_plots <- lapply(names(Mydf$INDEL[,4:11]), 
                    function(x) {Ann_plot(Mydf$INDEL, Mydf$INDEL[[x]], X_labels[[x]], Cutoffs[[x]], X_limits[[x]])})
names(Indel_plots) <- names(Mydf$INDEL[,4:11])

print("Created Indel Plots")

### Fisher's Strand on log scale
# add 1 so zero's are included
Indel_FS <- ggplot(Mydf$INDEL, aes(x=FS + 1)) + geom_density(aes(fill=Cat)) +
  scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
  geom_vline(xintercept=60, size=0.7, lty=2) + 
  scale_x_log10() +
  theme(legend.title=element_blank()) +
  xlab("Fisher Strand (log-scaled)")

print("Created FS plot for Indels")

# Build output filepath
outfile_indel_fp <- paste0(Directory, "/Percentile_Tables/INDEL_distributions.png")
Indel_plot <- grid.arrange(Indel_plots$QD, Indel_FS, Indel_plots$SOR, Indel_plots$MQ, Indel_plots$MQRankSum, Indel_plots$ReadPosRankSum, 
                           Indel_plots$DP, Indel_plots$QUAL, nrow=4)
​
print("Finished making Indel Plot, saving...")
ggsave(outfile_indel_fp, 
       plot = Indel_plot,
       device = "png")

