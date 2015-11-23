# sequence_handling
#### A series of scripts to automate DNA sequence aligning and quality control workflows via list-based batch submission and parallel processing
___
___
## Introduction

> For greater detail about everything, please see the [wiki](https://github.com/MorrellLAB/sequence_handling/wiki) for this repository

### What is `sequence_handling` for?

`sequence_handling` is a series of scripts to automate and speed up DNA sequence aligning and quality control through the use of our workflow outlined here. This repository contains two general kinds of scripts: *Shell Scripts* and *Batch Submission Scripts*, with one exception.

The former group is designed to be run directly from the command line. These serve as partial dependency installers, a way to generate a list for batch submission, QSub starters, and others that have issues with either running in parallel or using the [_Portable Batch System_](http://www.pbsworks.com/) due to memory issues. Running any of these scripts without any arguments generates a usage message for more details. Each script is named entirely in lower-case letters.

The latter group is designed to run the workflow in batch and in parallel. These scripts use a list of sequences, with full sequence paths, as their input and utilize [_GNU Parallel_](http://www.gnu.org/software/parallel/) to speed up the analysis and work they are designed for. Due to the length of time and resources needed for these scripts to run, they are designed to be submitted to a job scheduler, specifically the [_Portable Batch System_](http://www.pbsworks.com/). Each script is named using capital and lower-case letters.

Finally, there is one script that is neither designed to run directly from the shell nor submitted to a job scheduler. This script, `plot_cov.R` is designed to be called by `Plot_Coverage.sh` for creating coverage plots. This is done automatically; one does not need to change this script unless they wish to change the graphing parameters.

> **NOTE:** the latter group of scripts and `read_mapping_start.sh` are designed to use the [_Portable Batch System_](http://www.pbsworks.com/) and run on the [Minnesota Supercomputing Institute](https://www.msi.umn.edu). Heavy modifications will need to be made if not using these systems.

### Why use list-based batch submission?

Piping one sample alone through this workflow can take over 12 hours to completely run. Most sequence handling jobs are not dealing with one sample, so the amount of time to run this workflow increases drastically. List-based batch submission simplifies the amount of typing that one has to do, and enables parallel processing to decrease time spent waiting for samples to finish. An example list is shown velow

>/home/path\_to\_sample/sample\_001\_R1.fastq.gz

>/home/path\_to\_sample/sample\_001\_R2.fastq.gz

>/home/path\_to\_sample/sample\_003_R1.fastq.gz

>/home/path\_to\_sample/sample\_003\_R2.fastq.gz


### Why use parallel processing?

Parallel processing decreases the amount of time by running multiple jobs at once and keeping track of which are done, which are running, and which have yet to be run. This workflow, with the list-based batch submissions and parallel processing, both simplifies and quickens the process of sequence handling.

### Do I have to use the entire workflow as is?

No, with the one exception of `Plot_Coverage.sh` and `plot_cov.R`, no two scripts are entirely dependent on one another. While all these scripts are designed to easily use the output from one to the next, these scripts are not required to achive the end result of `sequence_handling`. If you prefer tools other than the ones used within this workflow, you can modify or replace any or all of the scripts offered in `sequence_handling`. This creates a pseudo-modularity for the entire workflow that allows for customization for each and every user.

### Dependencies

Due to the pseudo-modularity of this workflow, specific dependencies for each individual script are listed below. Some general dependencies for the workflow as a whole are listed here:

 - A quality trimmer, such as [_Seqqs_](https://github.com/morrelllab.seqqs), [_Sickle_](https://github.com/vsbuffalo/sickle), and [_Scythe_](https://github.com/vsbuffalo/scythe)
 - Tools for plotting results, such as [_R_](http://cran.r-project.org/)
 - SAM file processing utilities, such as [_SAMTools_](http://www.htslib.org/) and [_Picard_](http://broadinstitute.github.io/picard/)
 - A quality control mechanism, such as [_FastQC_](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - A read mapper, such as [_The Burrows-Wheeler Aligner_](http://bio-bwa.sourceforge.net/) (BWA)
 - [_GNU Parallel_](http://www.gnu.org/software/parallel/)

Please note that this is not a complete list of dependencies. Check below for specific dependencies for each desired script.

When running these scripts on the Minnesota Supercomputing Institute's (MSI) resources, most dependencies are included through MSI's module system. These modules are set to be automatically called by each script that calls upon them. However, some dependencies are not available through MSI; please check each script for which dependencies need to be installed separately.
___

## Shell Scripts

**NOTE: Running any of these scripts without arguments generates a usage message for greater detail about how to use them**

### installer.sh

The `installer.sh` script installs [_Seqqs_](https://github.com/morrelllab.seqqs), [_Sickle_](https://github.com/vsbuffalo/sickle), and [_Scythe_](https://github.com/vsbuffalo/scythe) for use with the `Quality_Triming.sh` script. It also has options for installing [_Bioawk_](https://github.com/lh3/bioawk), [_SAMTools_](http://www.htslib.org/) and [_R_](http://cran.r-project.org/), all dependencies for various scripts within this package.

##### dependencies

The `installer.sh` script depends on [_Git_](http://www.git-scm.com/), [_Wget_](http://www.gnu.org/software/wget/), the [_GNU Compiler Collection_](https://gcc.gnu.org/) (GCC), and [_GNU Make_](http://www.gnu.org/software/make/) to run.

### sample\_list\_generator.sh

The `sample_list_generator.sh` script creates a list of samples using a directory tree for its searching. This will find **all** samples in a given directory and its subdirectories. Only use this if you are using all samples within a directory tree. `sample_list_generator.sh` is designed to be run from the command line directly.

##### dependencies

The `sample_list_generator.sh` script has no external dependencies.

### read\_counts.sh

The `read_counts.sh` script calls [_Bioawk_](https://github.com/lh3/bioawk) to get accurate counts for read number for a list of samples. Output is written to a tab-delimited file file with sample name drawn from the file name for the list of samples.

##### dependencies

The `read_counts.sh` script depends on [_Bioawk_](https://github.com/lh3/bioawk) to run.

### read\_mapping\_start.sh

The `read_mapping_start.sh` script generates a series of QSub submissions for use with the [_Portable Batch System_](http://www.pbsworks.com/) on MSI's resources. starts a series of [BWA](http://bio-bwa.sourceforge.net/) sessions to map reads back to a reference genome.

##### dependencies

The `read_mapping_start.sh` script depends on the [_Portable Batch System_](http://www.pbsworks.com/) and [BWA](http://bio-bwa.sourceforge.net/) to run.
___

## Batch Submission Scripts

**NOTE: Each of these scripts contains usage information within the script itself. Furthermore, all values for these scripts are hard-coded into the script itself. Please open each script using your favourite text editor (ex. [_Vim_](http://www.vim.org), [_Sublime Text_](http://www.sublimetext.com), [_Visual Studio Code_](http://code.visualstudio.com), etc.) to read usage information and set values**

### Assess_Quality.sh

The `Assess_Quality.sh` script runs [_FastQC_](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the command line on a series of samples organized in a project directory for quality control. In addition, a list of all output zip files will be generated for use with the `Read_Depths.sh` script. Our recommendation is using this both before and after quality trimming and before read mapping. This script is designed to be run using the [_Portable Batch System_](http://www.pbsworks.com/).

##### dependencies

The `Assess_Quality.sh` script depends on [_FastQC_](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### Read\_Depths.sh

The `Read_Depths.sh` script utilizes the output from [_FastQC_](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to calculate the read depths for a batch of samples and outputs them into one convenient text file.

##### dependencies

The `Read_Depths.sh` script depends on the [_Portable Batch System_](http://www.pbsworks.com/) and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### Quality\_Trimming.sh

The `Quality_Trimming.sh` script runs `trim_autoplot.sh` (part of the [_Seqqs_](https://github.com/morrelllab.seqqs) repository on GitHub) on a series of samples organized in a project directory.. In addition to requiring _Seqqs_ to be installed, this also requires [GNU Parallel](http://www.gnu.org/software/parallel/) to be installed on the system.

##### dependencies

The `Quality_Trimming.sh` script depends on [_Sickle_](https://github.com/vsbuffalo/sickle), [_Scythe_](https://github.com/vsbuffalo/scythe), [_Seqqs_](https://github.com/morrelllab.seqqs), [_R_](http://cran.r-project.org/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### SAM\_Processing\_SAMTools.sh

The `SAM_Processing_SAMTools.sh` script converts the SAM files from read mapping with [BWA](http://bio-bwa.sourceforge.net/) to the BAM format using [_SAMTools_](http://www.htslib.org/). In the conversion process, it will sort and deduplicate the data for the finished BAM file, also using [_SAMTools_](http://www.htslib.org/). Alignment statistics will also be generated for both raw and finished BAM files. A list of finished BAM files will be generated at the end of this script.

##### dependencies

The `SAM_Processing_SAMTools.sh` script depends on [_SAMTools_](http://www.htslib.org/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### SAM\_Processing\_Picard.sh

The `SAM_Processing_Picard.sh` script converts the SAM files from read mapping with [BWA](http://bio-bwa.sourceforge.net/) to the BAM format using [_SAMTools_](http://www.htslib.org/). In the conversion process, it will sort and deduplicate the data for the finished BAM file, using [_Picard_](http://broadinstitute.github.io/picard/). Alignment statistics will also be generated for both raw and finished BAM files. A list of finished BAM files will be generated at the end of this script.

**NOTE: This script is extremely resource intensive, please use with caution.**

**NOTE: This script has not been tested, use with caution**

##### dependencies

The `SAM_Processing_Picard.sh` script depends on [_SAMTools_](http://www.htslib.org/), [_Picard_](http://broadinstitute.github.io/picard/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### Coverage\_Map.sh

The `Coverage_Map.sh` script generates coverage maps from BAM files using [_BEDTools_](http://bedtools.readthedocs.org/en/latest/). This map is in text format and is used for making coverage plots. In addition to generating coverage maps, this script will create a list of all the coverage maps generated for use in other scripts.

##### dependencies

The `Coverage_Map.sh` script depends on [_BEDTools_](http://bedtools.readthedocs.org/en/latest/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

### Plot\_Coverage.sh

The `Plot_Coverage.sh` script creates plots using [_R_](http://cran.r-project.org/) based off of coverage maps. It will generate three plots: one showing coverage across the genome, one showing coverage across exons, and one showing coverage across genes. This script uses `plot_cov.R` to generate the plots.

##### dependencies

The `Plot_Coverage.sh` script depends on the `plot_cov.R` script, [_R_](http://cran.r-project.org/), the [_Portable Batch System_](http://www.pbsworks.com/), and [_GNU Parallel_](http://www.gnu.org/software/parallel/) to run.

___

## Other Scripts

### plot\_cov.R

The `plot_cov.R` script is the graphical brains behind the `Plot_Coverage.sh` script. The latter will automatically call upon the former to create the coverage plots based off coverage maps. It is not necessary to open this script directly, except for making modifications to the graphical parameters.

##### dependencies

The `plot_cov.R` script has no external dependencies.

___

## TODO

 - ~~Generalize `read_counts.sh` for any project.~~ DONE!
 - ~~Add better list-out methods~~ DONE!
 - ~~Fix memory issues with `Read_Mapping.sh`~~ ~~Redesign read mapping scripts~~ DONE!
 - ~~Add coverage map script to workflow~~ ~~Finish integrating `Coverage_Map.sh` with the rest of the pipeline~~ DONE!
 - ~~Get `Plot_Coverage.sh` and `plot_cov.R` integrated into the pipeline~~ DONE!
 - ~~Add information about `plot_cov.R` to the README~~ DONE!
 - ~~Add script to easily convert SAM files from `Read_Mapping.sh` to BAM files for `Coverage_Map.sh`~~ ~~DONE!~~ ~~ish...~~ DONE!
 - ~~Add Deduplication script~~ Get ~~`Deduplication.sh`~~ `SAM_Processing_Picard.sh` working
 - ~~Add read mapping statistics via `samtools flagstat`~~ DONE! This is integrated into `SAM_Processing_SAMTools.sh`
 - Add config file and move all argument passing to that
 - Parallelize BWA
 - Create master `sequence_handling` script
 - Incorporate variant calling scripts into the pipeline
 - keep README updated
