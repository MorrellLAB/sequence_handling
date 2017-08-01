# `sequence_handling`
#### A pipeline to automate DNA sequence aligning and quality control workflows via list-based batch submission and parallel processing
___

> For greater detail, usage information, and troubleshooting please see the [`sequence_handling` wiki](https://github.com/MorrellLAB/sequence_handling/wiki).

## What is `sequence_handling` for?

`sequence_handling` is a series of scripts, called handlers, that automate and speed up DNA sequence alignment and quality control. Currently, `sequence_handling` is designed to work with Illumina paired-end whole-genome or exome capture sequences. Work is underway to expand `sequence_handling` to accept GBS sequences.

The workflow is intended to be 100% reproducible provided that you have the Config file and version of `sequence_handling` that was used. It is also intended to be easy for beginner UNIX users to configure and run independently, given that they read the [wiki](https://github.com/MorrellLAB/sequence_handling/wiki). 

> **NOTE:** This workflow is designed to use the [Portable Batch System](http://www.pbsworks.com/) and run on the [Minnesota Supercomputing Institute](https://www.msi.umn.edu). Heavy modifications will need to be made if not using these systems.

## Dependencies

Due to the pseudo-modularity of this workflow, dependencies for each individual handler are listed below. Some general dependencies for the workflow as a whole are also listed here:

 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - A quality control mechanism, such as [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - An adapter trimmer, such as [Scythe](https://github.com/vsbuffalo/scythe)
 - A read mapper, such as [The Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA)
 - SAM file processing utilities, such as [SAMTools](http://www.htslib.org/) and/or [Picard](http://broadinstitute.github.io/picard/)
 - Tools for plotting coverage, such as [R](http://cran.r-project.org/)

Please note that this is not a complete list of dependencies. Check the [dependencies wiki page](https://github.com/MorrellLab/sequence_handling/wiki/Dependencies) for dependencies for each handler.

## Basic Usage

A guide on how to install, set up, and run `sequence_handling` for the first time can be found on the [Getting Started](https://github.com/MorrellLab/sequence_handling/wiki) page of the wiki.

A brief usage message can be viewed by passing no arguments to `sequence_handling`:

```shell
./sequence_handling
```

To run `sequence_handling`, use the following command, assuming you are in the `sequence_handling` directory:

```shell
./sequence_handling <handler> Config
```

Where `<handler>` is one of the handlers listed below and `Config` is the full file path to the configuration file. 

## Recommended Workflow

![Workflow](https://raw.githubusercontent.com/MorrellLAB/sequence_handling/master/.Sequence_Handling_Workflow.png)

#### 1. [Quality\_Assessment](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Assessment)

To start, run Quality_Assessment on your raw FastQ files. The Quality_Assessment handler runs [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on a series of samples and outputs metrics used for quality control. It accepts FASTQ, SAM, and BAM files as input and outputs a summary table and individual HTML files for visualization. The Quality_Assessment handler depends on [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 2. [Adapter\_Trimming](https://github.com/MorrellLab/sequence_handling/wiki/Adapter_Trimming)

The Adapter_Trimming handler uses [Scythe](https://github.com/vsbuffalo/scythe) to trim specific adapter sequences from FastQ files. This handler differentiates between forward, reverse, and single-end FastQ files automatically. The Adapter_Trimming handler depends on [Scythe](https://github.com/vsbuffalo/scythe) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 3. [Quality\_Assessment](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Assessment)

After Adapter_Trimming, it is recommended to run Quality_Assessment again on the trimmed FastQ files to ensure that all adapter contamination was properly removed. 

#### 4. [Read\_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Read_Mapping)

The Read_Mapping handler maps sequence reads to a reference genome using [BWA-MEM](http://bio-bwa.sourceforge.net/). This handler uses Torque Task Arrays, part of the [Portable Batch System](http://www.pbsworks.com/). The Read_Mapping handler depends on the [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/).

#### 5. [SAM\_Processing with Picard](https://github.com/MorrellLAB/sequence_handling/wiki/SAM_Processing)

The SAM_Processing handler converts the SAM files from read mapping with [BWA](http://bio-bwa.sourceforge.net/) to the BAM format using [SAMTools](http://www.htslib.org/). In the conversion process, it will sort and deduplicate the data for the finished BAM file, also using [SAMTools](http://www.htslib.org/). Alignment statistics will also be generated for both raw and finished BAM files. The SAM_Processing handler depends on [SAMTools](http://www.htslib.org/) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 6. [Coverage_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Coverage_Mapping)

The Coverage_Mapping handler generates coverage histograms and summary statistics from BAM files using [BEDTools](http://bedtools.readthedocs.org/en/latest/). Plots of coverage are generated using [R](http://cran.r-project.org/) based on coverage maps. The Coverage_Mapping handler depends on [BEDTools](http://bedtools.readthedocs.org/en/latest/), [R](http://cran.r-project.org/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

## Troubleshooting

Please check the [FAQ and Troubleshooting](https://github.com/MorrellLAB/sequence_handling/wiki/FAQ) page on the wiki for help and visit the [wiki](https://github.com/MorrellLAB/sequence_handling/wiki) page for the handler(s) you're using. If you are still having difficulties, please submit an issue through the [Git issues page](https://github.com/MorrellLab/sequence_handling/issues) for this repository.

## Future Handlers

The following handlers are planned for future versions of `sequence_handling`.

#### GBS_Demultiplexing

The GBS_Demultiplexing handler will demulitplex raw GBS reads into split FastQ files using the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). Code for GBS_Demultiplexing is already present in `sequence_handling` for the purposes of collaborative alpha testing, but should not be used since it's extremely buggy.

#### SNP Calling with GATK

#### Coverage Mapping plots using R

#### Other suggestions

If you feel like there are tools or alternative processing techniques missing from `sequence_handling`, please submit a feature request through the Git issues page.

## Citation

`sequence_handling` can be cited like:

> Hoffman PJ, Wyant SR, Kono TJY, Morrell PL. (2016). sequence_handling: A pipeline to automate sequence workflows [Software]. Available at https://github.com/MorrellLAB/sequence_handling.
