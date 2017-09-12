#!/usr/bin/env Rscript

#   This script is currently broken, but will hopefully be fixed for a future version

#   To take command line arguments
args <- commandArgs(TRUE)

#   A list of required packages
PACKAGE_LIST <- c('Hmisc')

#   Test to see if required packages are installed
pkgTest <- function(package) {
    if(package %in% rownames(installed.packages()) == FALSE) { # check to see if a packages is available
        install.packages(package) # if not, install it
    }
}

#   Load the packages
batchInstall <- function(pkgList) {
    options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
    # for(dep in pkgList) {
    #     pkgTest(package = dep) # test to see if the package is installed
    # }
    lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used
}

#   A function to read in the coverage map
read.coverage.map <- function(infile) {
    cat('Reading in', infile, '\n', sep = ' ', file = stderr())
    #   Read in our coverage map
    map <- read.table(file = infile, header = FALSE, as.is = TRUE, fill = TRUE)
    #   Return the map
    return(map)
}

#   A function to subset the coverage maps into genome, gene, and exon space
subset.coverage <- function(map, type, position, depth, coverage) {
    cat('Subsetting by', type, '\n', sep = ' ')
    #   Subset our coverage map
    #   eval(expr = as.symbol(x = var)) forces variable expression
    sub <- subset(x = map, subset = eval(expr = as.symbol(x = position)) == type, select = c(depth, coverage))
    #   Name the columns
    names(x = sub) <- c('Depth', 'Fraction')
    #   Return the subset
    return(sub)
}

#   A function to plot coverage
plot.coverage <- function(map, nrow = nrow(x = map), type = 'Coverage', col = 'peachpuff4') {
    #   Create a title, using Hmisc to capitalize the type
    main <- paste(capitalize(string = type), 'Coverage', sep = ' ')
    #   Create a y label for the graph
    ylab <- paste('Fraction of', type, 'at depth', sep = ' ')
    #   Set the limit for scientific notation
    options(scipen = 5)
    #   Make the plot
    plot(
        x = map[1:nrow, 'Depth'], # Depth up to a point
        y = map[1:nrow, 'Fraction'], # Coverage to the same point
        type = 'h', # Histogram
        main = main, # main is our title
        xlab = 'Depth', # x axis is depth
        ylab = ylab, # y axis is our fraction of coverage
        col = col, # Color the graph
        lwd = 5, # Make our bars thicker
        las = 1 # All axis tickmark labels are horizontal
    )
}

#   The driver function
main <- function() {
    #   Collect the arguments
    covfile <- args[1]
    outdir <- args[2]
    samplename <- args[3]
    sequencehandling <- args[4]
    #   Create and set directory for R packages
    library.directory <- file.path(sequencehandling, '.RLibs', fsep = '/')
    dir.create(path = library.directory, showWarnings = FALSE)
    .libPaths(new = library.directory)
    #   Check for and load required packages
    batchInstall(pkgList = PACKAGE_LIST)
    #   Read in the coverage map
    coverage <- read.coverage.map(infile = covfile)
    #   Subset our data into genome, gene, and exon maps
    genome <- subset.coverage(map = coverage, type = 'all', position = 'V1', depth = 'V2', coverage = 'V5')
    gene <- subset.coverage(map = coverage, type = 'gene', position = 'V3', depth = 'V10', coverage = 'V13')
    exon <- subset.coverage(map = coverage, type = 'exon', position = 'V3', depth = 'V10', coverage = 'V13')
    #   Plot our coverage
    outfile <- paste0(outdir, '/', samplename, '_outplot.pdf')
    cat('Creating plots at', outfile, '\n', sep = ' ', file = stderr())
    pdf(file = outfile, width = 6, height = 6)
    plot.coverage(map = genome, nrow = 201, type = 'genome', col = 'navy')
    plot.coverage(map = gene, nrow = 151, type = 'gene', col = 'cornflowerblue')
    plot.coverage(map = exon, nrow = 151, type = 'exon', col = 'violetred1')
    dev.off()
}

#   Run the program
main()
