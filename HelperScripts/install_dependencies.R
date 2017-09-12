#!/usr/bin/env Rscript

#   This is currently untested and unused, but might be added into a future version

#   Collect arguments, leave trailingOnly to 'FALSE' to collect the filename
args <- commandArgs(trailingOnly = FALSE)

#   Collect the filename from the argument list
#   Doing it this way allows for other arguments to be passed to Rscript rather
#   than depending on the filename to be at commandArgs[4]
FILEARG <- args[grep(pattern = 'file', x = args)]
FILENAME <- unlist(x = strsplit(x = FILEARG, split = '='))[2]

#   A function to check for and install packages
pkgTest <- function(package) {
    #   A function to say if something failed
    failed <- function(cond) { # Needs a 'cond' argument to take the error/warngings from R
        cat('Failed to install', package, '\n', file=stderr())
        quit(save = 'default', status = 1)
    }
    #   If the package is not installed already
    if(!(package %in% rownames(x = installed.packages()))) {
        tryCatch(
            expr = { # Try to install the package and say that we did it
                install.packages(package)
                cat('Successfully installed', package, '\n', file=stderr())
            },
            warning = failed(w), # Warnings run the failed function
            error = failed(e) # Errors run the failed function
        )
    } else { # If we have the package already, let us know
        cat(package, 'is already installed\n', file=stderr())
    }
}

#   A function to create a directory
createLibrary <- function(libbase){
    library.directory <- file.path(libbase, '.RLibs', fsep = '/') # Create a filepath
    cat('Creating R Library at', library.directory, '\n', file=stderr())
    dir.create(path = library.directory, showWarnings = FALSE) # Make the file path
    return(library.directory) # Return the file path
}

#   A usage message that quits afterwards
usage <- function(filename=FILENAME) {
    line1 <- paste('Usage:', filename, '<RLib>, <dependency, [dependency, ...]>')
    line2 <- paste('Where:\t<RLib> is the directory to place the R Library and')
    line3 <- paste('\t<dependency, [dependency, ...]> is a list of one or more packages to install')
    cat(line1, line2, line3, sep = '\n', file = stderr())
    quit(save = 'no', status = 1)
}

#   The driver function
main <- function() {
    #   Usage if no arguments
    if(!('--args' %in% args)) {
        usage()
    }
    #   Parse our arguments
    argstart <- grep(pattern = '--args', x = args)
    rlib <- args[argstart + 1] # R Library base is fisrt
    deps <- args[seq(from = argstart + 2, to = length(args))] # Dependencies come after
    #   Clean up the dependencies vector
    deps <- deps[!(is.na(x = deps))]
    deps <- deps[!(rlib == deps)]
    if(length(x = deps) < 1) {
        usage()
    }
    #   Create our library
    library.directory <- createLibrary(libbase = rlib)
    .libPaths(new = library.directory)
    #   Install our dependencies
    for(dependency in deps) {
        pkgTest(package = dependency)
    }
}

#   Run it
main()
