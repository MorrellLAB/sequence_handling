#!/usr/bin/env Rscript

### A script to plot distributions of Variant Annotations
### Adapted material from: https://evodify.com/gatk-in-non-model-organism/

library(ggplot2)
library(gridExtra)

options(warn=1)

# added this line to prevent a device window from opening
pdf(NULL)

##########################
######## Setup_df ########
##########################

load_variants <- function(vcf_table) {
    variants <- read.csv(vcf_table, header=T, na.strings=c("","NA"), sep = "\t")
    # Prepare columns
    variants$Chrom_Pos <- paste0(variants$CHROM, sep = "_", variants$POS, sep = "_", variants$TYPE)
    return(variants)
}

# VCFs to plot
Setup_DF <- function(vcf1_table, vcf2_table, plot_type) {
    if (plot_type == "single") {
        # Load single variant table
        variants <- load_variants(vcf1_table)
        types <- split(variants, variants$TYPE)
    } else if (plot_type == "pair") {
        # Load two variant tables
        variants_set1 <- load_variants(vcf1_table)
        variants_set2 <- load_variants(vcf2_table)
        # Setup data
        variants_set1$Cat <- ifelse(variants_set1$Chrom_Pos %in% variants_set2$Chrom_Pos,
                                    "Passed_filtering", "Failed_filtering")
        types <- split(variants_set1, variants_set1$TYPE)
    }
    return(types)
}

# Prepare xlimits based on data
# Note: Check that this works for "pair" case as well
# This helps set an appropriate xlim for certain annotation value ranges
Find_Limit <- function(lower_bound, upper_bound) {
    limit_value <- pmax(abs(lower_bound), abs(upper_bound))
    return(limit_value)
}

Get_X_limits <- function(Mydf, Cutoffs) {
    # Prepare MQRankSum xlim range
    mqranksum_limit <- ceiling(Find_Limit(min(Mydf$SNP[["MQRankSum"]], na.rm=T), max(Mydf$SNP[["MQRankSum"]], na.rm=T)))
    if (mqranksum_limit <= abs(Cutoffs$MQRankSum)) {
        mqranksum_limit <- abs(Cutoffs$MQRankSum) + 1
    } else {
        mqranksum_limit <- mqranksum_limit + 2
    }
    # GATK suggests filtering out ReadPosRankSum values less than -8.0
    # So, we'll adjust the limits so we can see this cutoff
    readposranksum_limit <- ceiling(Find_Limit(min(Mydf$SNP[["ReadPosRankSum"]], na.rm=T), max(Mydf$SNP[["ReadPosRankSum"]], na.rm=T)))
    if (readposranksum_limit <= abs(Cutoffs$ReadPosRankSum)) {
        readposranksum_limit <- abs(Cutoffs$ReadPosRankSum) + 1
    } else {
        readposranksum_limit <- readposranksum_limit + 2
    }
    # Since QUAL can exceed 31000, we will set a xlim cutoff for display purposes if the max exceeds it
    qual_max <- ceiling(max(Mydf$SNP[["QUAL"]]))
    if (qual_max > 2000) {
        qual_limit <- 2000
    } else {
        qual_limit <- qual_max
    }
    # Same for DP, it can exceed 27000
    dp_max <- ceiling(max(Mydf$SNP[["DP"]]))
    if (dp_max > 6000) {
        dp_limit <- 6000
    } else {
        dp_limit <- dp_max
    }
    # Get the final xlim ranges
    X_limits <- list("QD" = c(0, ceiling(max(Mydf$SNP[["QD"]])) + 5),
                     "FS" = c(0, ceiling(max(Mydf$SNP[["FS"]])) + 5), 
                     "SOR" = c(0, ceiling(max(Mydf$SNP[["SOR"]])) + 5),
                     "MQ" = c(0, ceiling(max(Mydf$SNP[["MQ"]])) + 10),
                     "MQRankSum" = c(-mqranksum_limit, mqranksum_limit),
                     "ReadPosRankSum" = c(-readposranksum_limit, readposranksum_limit),
                     "QUAL" = c(0, qual_limit + 10),
                     "DP" = c(0, dp_limit + 10))
    return(X_limits)
}

###############################
#### FUNCTIONS FOR PLOTTING ###
###############################

# Annotation plots for SNPs and Indels
Ann_plot <- function(dataset, vector, lab, cutoff, xrange, plot_type, curr_ann) {
    if (plot_type == "single") {
        # Plot single variant table
        plot <- ggplot(dataset, aes(x=vector)) +
            geom_density(fill=alpha('#F4CCCA', 0.6), na.rm=TRUE) +
            theme(legend.position = "none")
    } else if (plot_type == "pair") {
        # Plot two variant tables
        plot <- ggplot(dataset, aes(x=vector)) + 
            geom_density(aes(fill=Cat), na.rm=TRUE) + 
            scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
            theme(legend.title=element_blank())
    }
    
    # Add cutoffs, text labels, and other shared plotting options
    plot <- plot + 
        geom_vline(xintercept=cutoff, size=0.7, lty=2, na.rm=TRUE) +
        geom_text(x=cutoff, y=0, aes(label=paste0("Cutoff = ", cutoff)), 
                  hjust=-0.1, vjust=1.5, nudge_x=-1, nudge_y=-1,
                  colour=alpha("red", 0.9), size=2, angle=90) +
        xlab(lab) +
        xlim(xrange[1], xrange[2])
    
    # If the current cutoff is either QUAL or DP, add an additional label for the max
    #   since the max value most likely won't be plotted
    if (curr_ann == "QUAL" || curr_ann == "DP") {
        print("Adding max value to plot")
        # Add max label and value to the plot
        max_value <- max(vector)
        text_label <- paste0("Max (not shown) = ", max_value)
        plot <- plot + geom_text(aes(label=text_label), x=Inf, y=Inf,
                                 hjust=1.1, vjust=1.5, nudge_x=-1, nudge_y=-1,
                                 colour="grey50", size=2)
    }
    return (plot)
}

#--------
# SNPs
### Fisher's Strand on log scale
# add 1 so zero's are included
SNP_plot_FS <- function(Mydf, Cutoffs, plot_type) {
    if (plot_type == "single") {
        # Plot single variant table
        plot <- ggplot(Mydf$SNP, aes(x=FS + 1)) +
            geom_density(fill=alpha('#F4CCCA', 0.6), na.rm=TRUE) +
            theme(legend.position = "none")
    } else if (plot_type == "pair") {
        # Plot two variant tables
        plot <- ggplot(Mydf$SNP, aes(x=FS + 1)) +
            geom_density(aes(fill=Cat), na.rm=TRUE) +
            scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
            theme(legend.title=element_blank())
    }
    # Add cutoffs, text labels, and other sahred plotting options
    plot <- plot +
        geom_vline(xintercept=Cutoffs$FS, size=0.7, lty=2) +
        geom_text(x=Cutoffs$FS, y=0, aes(label=paste0("Cutoff = ", Cutoffs$FS)), 
                  hjust=-0.1, vjust=1.5, nudge_x=-1, nudge_y=-1,
                  colour=alpha("red", 0.9), size=2, angle=90) +
        scale_x_log10() +
        xlab("Fisher Strand (log-scaled)")
    return (plot)
}

# Do the plotting and rearrange output plots
Plot_SNP_Only <- function(Mydf, X_labels, Cutoffs, X_limits, plot_type) {
    ### SNPs only
    SNP_plots <- lapply(
        names(Mydf$SNP[,4:11]),
        function(x) {
            curr_ann <- x
            Ann_plot(Mydf$SNP, Mydf$SNP[[x]], X_labels[[x]], Cutoffs[[x]], X_limits[[x]], plot_type, curr_ann)
        })
    names(SNP_plots) <- names(Mydf$SNP[,4:11])
    # Plot SNPs, Fisher's strand
    SNP_FS <- SNP_plot_FS(Mydf, Cutoffs, plot_type)
    
    # Arrange plots into multi-panel plot
    SNP_plot <- arrangeGrob(SNP_plots$QD, SNP_FS, SNP_plots$SOR, SNP_plots$MQ,
                            SNP_plots$MQRankSum, SNP_plots$ReadPosRankSum, 
                            SNP_plots$DP, SNP_plots$QUAL, nrow=4)
    return (SNP_plot)
}

#--------
# INDELS
### Fisher's Strand on log scale
# add 1 so zero's are included
Indel_plot_FS <- function(Mydf, Cutoffs, plot_type) {
    if (plot_type == "single") {
        # Single variant table
        plot <- ggplot(Mydf$INDEL, aes(x=FS + 1)) +
            geom_density(fill=alpha('#F4CCCA', 0.6), na.rm=TRUE) +
            theme(legend.position = "none")
    } else if (plot_type == "pair") {
        # Two variant tables
        plot <- ggplot(Mydf$INDEL, aes(x=FS + 1)) +
            geom_density(aes(fill=Cat), na.rm=TRUE) +
            scale_fill_manual(values=c(alpha('#F4CCCA', 0.6), alpha('#A9E2E4', 0.6))) +
            geom_vline(xintercept=Cutoffs$FS, size=0.7, lty=2) + 
            scale_x_log10() +
            theme(legend.title=element_blank())
    }
    # Add cutoffs, text labels, and other shared plotting options
    plot <- plot +
        geom_vline(xintercept=Cutoffs$FS, size=0.7, lty=2) + 
        scale_x_log10() +
        xlab("Fisher Strand (log-scaled)")
    return (plot)
}

# Do the plotting and rearrange output plots
Plot_INDEL_Only <- function(Mydf, X_labels, Cutoffs, X_limits, plot_type) {
    Indel_plots <- lapply(
        names(Mydf$INDEL[,4:11]), 
        function(x) {
            curr_ann <- x
            Ann_plot(Mydf$INDEL, Mydf$INDEL[[x]], X_labels[[x]], Cutoffs[[x]], X_limits[[x]], plot_type, curr_ann)
        })
    names(Indel_plots) <- names(Mydf$INDEL[,4:11])
    # Plot indels, Fisher's strand
    Indel_FS <- Indel_plot_FS(Mydf, Cutoffs, plot_type)
    # Arrange plots into multi-panel plot
    Indel_plot <- arrangeGrob(Indel_plots$QD, Indel_FS, Indel_plots$SOR, Indel_plots$MQ, 
                              Indel_plots$MQRankSum, Indel_plots$ReadPosRankSum, 
                              Indel_plots$DP, Indel_plots$QUAL, nrow=4)
    return (Indel_plot)
}

#----------
# Build output filepath and save plots to file
Save_SNP_Plots <- function(Directory, Prefix, Suffix, SNP_plot) {
    outfile_snp_fp <- paste0(Directory, "/Percentile_Tables/", Prefix, "_SNP_distributions", Suffix, ".png")
    print(paste0("SNP graph filename path is ", outfile_snp_fp))
    ggsave(outfile_snp_fp, plot = SNP_plot, device = "png")
}

Save_Indel_Plots <- function(Directory, Prefix, Suffix, Indel_plot) {
    outfile_indel_fp <- paste0(Directory, "/Percentile_Tables/", Prefix, "_INDEL_distributions", Suffix, ".png")
    print(paste0("Indel graph filename  path is ", outfile_indel_fp))
    ggsave(outfile_indel_fp, plot = Indel_plot, device = "png")
}

#############################
##### DRIVER FUNCTION #######
#############################

main <- function() {
    # Driver function
    # Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    # User provided command line arguments
    Directory <- args[1]
    Variants_Set1 <- args[2]
    # If we are only plotting for a single VCF, put "NA" as the input here
    Variants_Set2 <- args[3]
    Suffix <- args[4]
    # Type of plot we want to make. Valid options include:
    #   1) single
    #   2) pair
    plot_type <- args[5]
    Variants_Set1_Prefix <- args[6] # Prefix will be used in output file naming
    Variants_Set2_Prefix <- args[7] # Put "NA" if we are only plotting a single VCF file
    
    ##########################
    # SETUP LABELS + CUTOFFS #
    ##########################
    
    # Set default cutoffs to use if user did not provide cutoffs
    Cutoffs <- list("QD" = 2, "FS" = 60, 
                    "SOR" = 3, "MQ" = 40, 
                    "MQRankSum" = -12.5, 
                    "ReadPosRankSum" = -8.0, "QUAL" = 30, 
                    "DP" = 5)
    
    X_labels <- list("QD" = c("Quality by Depth"), "FS" = c("Fisher Strand"), 
                     "SOR" = c("Strand Odds Ratio"), "MQ" = c("Mapping Quality"),
                     "MQRankSum" = c("Mapping Quality Rank Sum Test"), 
                     "ReadPosRankSum" = c("Read Position Rank Sum Test"), "QUAL" = c("Quality Score"), 
                     "DP" = c("Depth"))
    
    # X_limits will be set later on
    
    ############
    # PLOTTING #
    ############
    
    if (plot_type == "single") {
        ##### SINGLE plots
        print("single plot mode")
        # Read in only Variants_Set1
        Mydf <- Setup_DF(Variants_Set1, Variants_Set2, plot_type)
        # Prepare xlims
        X_limits <- Get_X_limits(Mydf, Cutoffs)
        
        # Are we plotting SNPs only or also Indels?
        if ("SNP" %in% names(Mydf)) {
            ### Plot SNPs
            SNP_plot <- Plot_SNP_Only(Mydf, X_labels, Cutoffs, X_limits, plot_type)
            Save_SNP_Plots(Directory, Variants_Set1_Prefix, Suffix, SNP_plot)
        }
        
        if ("INDEL" %in% names(Mydf)) {
            ### Plot INDELs
            INDEL_plot <- Plot_INDEL_Only(Mydf, X_labels, Cutoffs, X_limits, plot_type)
            Save_Indel_Plots(Directory, Variants_Set1_Prefix, Suffix, Indel_plot)
        }
    } else if (plot_type == "pair") {
        ##### PAIR plots
        print("pair plot mode")
        # Read in two variants tables
        Mydf <- Setup_DF(Variants_Set1, Variants_Set2, plot_type)
        # Prepare xlims
        X_limits <- Get_X_limits(Mydf, Cutoffs)
        
        # Prepare output file prefix
        Prefix <- paste0(Variants_Set1_Prefix, "_vs_", Variants_Set2_Prefix)
        # Generate plots
        if ("SNP" %in% names(Mydf)) {
            ### Plot SNPs
            SNP_plot <- Plot_SNP_Only(Mydf, X_labels, Cutoffs, X_limits, plot_type)
            Save_SNP_Plots(Directory, Prefix, Suffix, SNP_plot)
        }
        
        if ("INDEL" %in% names(Mydf)) {
            ### Plot INDELs
            INDEL_plot <- Plot_INDELs_Only(Mydf, X_labels, Cutoffs, X_limits, plot_type)
            Save_Indel_Plots(Directory, Prefix, Suffix, Indel_plot)
        }
    }
}

main() # Run the program
