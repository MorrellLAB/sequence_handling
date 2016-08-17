# sequence_handling
#### A pipeline to automate DNA sequence aligning and quality control workflows via list-based batch submission and parallel processing
___

## For greater detail, usage information, and troubleshooting please see the [sequence_handling wiki](https://github.com/MorrellLAB/sequence_handling/wiki).

## What is `sequence_handling` for?

`sequence_handling` is a series of scripts, or handlers, to automate and speed up DNA sequence alignment and quality control. Currently, `sequence_handling` is designed to work with Illumina paired-end whole-genome and exome capture data. Work is underway to expand `sequence_handling` to accept single-end and GBS data.

The workflow is designed to process samples in batches and in parallel. It is also intended to be easy for users to configure. The handlers use a list of sequences, with full sequence paths as input. `sequence_handling` uses [GNU Parallel](http://www.gnu.org/software/parallel/) to speed up  analysis.The handlers are designed to run as jobs submitted to a job scheduler, specifically the [Portable Batch System](http://www.pbsworks.com/).

> **NOTE:** This workflow is designed to use the [Portable Batch System](http://www.pbsworks.com/) and run on the [Minnesota Supercomputing Institute](https://www.msi.umn.edu). Heavy modifications will need to be made if not using these systems.

## Configuration File

The included configuration file, `Config`, provides information needed to run each of the handlers within it. No other information is needed as `sequence_handling` pulls all necessary information from `Config`. Variables that are used by more than one handler are located at the top of `Config`, followed by handler-specific variables, ending with software definitions. Please read `Config` or visit [the wiki page](https://github.com/MorrellLAB/sequence_handling/wiki/Configuration) for more usage information.

`Config` is broken up into several sections. The first section, at the top of `Config` contains variables that are used by more than one handler. Each section below is headed by a block of hash (`#`) marks and contains variables for one specific handler only. For example, the section headed by

```shell
############################################
##########    Adapter_Trimming    ##########
############################################
```

contains variables for Adapter_Trimming only. These variables are completely ignored by other handlers.

Please note, some of the variables are pre-defined in `Config`. These have been set for using the entirety of `sequence_handling`, and follows naming conventions used by all of the handlers. If you choose to not use some of the handlers in your analyses (See *Do I have to use the entire workflow as is?* below), please modify variables as needed.

## Why use list-based batch submission?

Piping one sample alone through this workflow can take over 12 hours to run to completion. Most sequence handling jobs are not dealing with one sample, so the amount of time to run this workflow increases drastically. List-based batch submission simplifies the amount of typing that one has to do, and enables parallel processing to decrease time spent waiting for samples to finish. An example list is shown below.

>/home/path\_to\_sample/sample\_001\_R1.fastq.gz

>/home/path\_to\_sample/sample\_001\_R2.fastq.gz

>/home/path\_to\_sample/sample\_003_R1.fastq.gz

>/home/path\_to\_sample/sample\_003\_R2.fastq.gz


## Why use parallel processing?

Parallel processing decreases the amount of time by running multiple jobs at once and keeping track of which are done, which are running, and which have yet to be run. This workflow, with the list-based batch submissions and parallel processing, both simplifies and quickens the process of sequence handling.

## Do I have to use the entire workflow as is?

No. No two handlers are entirely dependent on one another. While all these handlers are designed to easily use the output from one to the next, these handlers are not required to achieve the end result of `sequence_handling`. If you prefer tools other than the ones used within this workflow, you can modify or replace any or all of the handlers offered in `sequence_handling`. This creates a pseudo-modularity for the entire workflow that allows for customization for each and every user. 

## Dependencies

Due to the pseudo-modularity of this workflow, specific dependencies for each individual handler are listed below. Some general dependencies for the workflow as a whole are also listed here:

 - An adapter trimmer, such as [Scythe](https://github.com/vsbuffalo/scythe)
 - A quality trimmer, such as the [Seqqs](https://github.com/morrelllab.seqqs)/[Sickle](https://github.com/najoshi/sickle) combo
 - Tools for plotting results, such as [R](http://cran.r-project.org/)
 - SAM file processing utilities, such as [SAMTools](http://www.htslib.org/) and/or [Picard](http://broadinstitute.github.io/picard/)
 - A quality control mechanism, such as [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - A read mapper, such as [The Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

Please note that this is not a complete list of dependencies. Check the [dependencies wiki page](https://github.com/MorrellLab/sequence_handling/wiki/Dependencies) for specific dependencies for each handler.

___

## Basic Usage

To run `sequence_handling`, use the following command, assuming you are in the `sequence_handling` directory:

```shell
./sequence_handling <handler> Config
```

Where `<handler>` is one of the handlers listed below and `Config` is the full file path to the configuration file.

A brief usage message can be viewed by passing no arguments to `sequence_handling`:

```shell
./sequence_handling
```

## Handlers

> **For greater detail about each handler, visit the [sequence_handling wiki](https://github.com/MorrellLAB/sequence_handling/wiki).**

### [Quality\_Assessment](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Assessment)

The Quality_Assessment handler runs [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on a series of samples organized in a project directory for quality control. In addition, a list of all output zip files will be generated for use with the Read_Depths handler. Our recommendation is to perform quality assesment on your raw samples and before and after quality trimming. This script is designed to be run using the [Portable Batch System](http://www.pbsworks.com/). The Quality_Assessment handler depends on [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), the [Portable Batch System](http://www.pbsworks.com/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

### [Read\_Depths](https://github.com/MorrellLab/sequence_handling/wiki/Read_Depths)

The Read_Depths handler utilizes the output from [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to estimate the read depths for a batch of samples and outputs them into one convenient text file. The Read_Depths handler depends on the [Portable Batch System](http://www.pbsworks.com/) and [GNU Parallel](http://www.gnu.org/software/parallel/) to run.

### [Adapter\_Trimming](https://github.com/MorrellLab/sequence_handling/wiki/Adapter_Trimming)

The Adapter_Trimming handler uses [Scythe](https://github.com/vsbuffalo/scythe) to trim adapter sequences from FastQ files. This handler differentiates between forward, reverse, and single-end FastQ files automatically. The Adapter_Trimming handler depends on [Scythe](https://github.com/vsbuffalo/scythe), the [Portable Batch System](http://www.pbsworks.com/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

### [Quality\_Trimming](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Trimming)

The Quality_Trimming handler uses [Sickle](https://github.com/najoshi/sickle) to trim sequences based on quality. It also uses [Seqqs](https://github.com/vsbuffalo/seqqs) to generate quality statistics for raw and trimmed sequences. This handler differentiates between forward, reverse, and single-end FastQ files when given proper suffixes for each class of FastQ file. The Quality_Trimming handler depends on [Sickle](https://github.com/najoshi/sickle), [Seqqs](https://github.com/vsbuffalo/seqqs), [R](http://cran.r-project.org/), the [Portable Batch System](http://www.pbsworks.com/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

### [Read\_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Read_Mapping)

The Read_Mapping handler maps reads to a reference genome using [BWA-MEM](http://bio-bwa.sourceforge.net/). This handler uses Torque Task Arrays, part of the [Portable Batch System](http://www.pbsworks.com/). The Read_Mapping handler depends on the [Portable Batch System](http://www.pbsworks.com/) and [BWA](http://bio-bwa.sourceforge.net/).

### [SAM\_Processing](https://github.com/MorrellLAB/sequence_handling/wiki/SAM_Processing)

The SAM_Processing handler converts the SAM files from read mapping with [BWA](http://bio-bwa.sourceforge.net/) to the BAM format using [SAMTools](http://www.htslib.org/). In the conversion process, it will sort and deduplicate the data for the finished BAM file, also using [SAMTools](http://www.htslib.org/). Alignment statistics will also be generated for both raw and finished BAM files. The SAM_Processing handler depends on [SAMTools](http://www.htslib.org/), the [Portable Batch System](http://www.pbsworks.com/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

### [Coverage_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Coverage_Mapping)

The Coverage_Mapping handler generates coverage maps from BAM files using [BEDTools](http://bedtools.readthedocs.org/en/latest/). This map is in text format and is used for making coverage plots. Plots of coverage are generated using [R](http://cran.r-project.org/) based on coverage maps. Plots for coverage over genome, exon, and gene space are generated. The Coverage_Mapping handler depends on [BEDTools](http://bedtools.readthedocs.org/en/latest/), [R](http://cran.r-project.org/), the [Portable Batch System](http://www.pbsworks.com/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

### [Indel_Realignment](https://github.com/MorrellLab/sequence_handling/wiki/Indel_Realignment)

The Indel_Realignment handler realigns the BAM files using [GATK's](http://gatkforums.broadinstitute.org/gatk/discussion/38/local-realignment-around-indels) RealignerTargetCreator and IndelRealigner. Indel_Realignment depends on [Java](https://www.java.com/en/), [GATK](https://www.broadinstitute.org/gatk/), and the [Portable Batch System](http://www.pbsworks.com/). Indel_Realignment is the first step in SNP calling for many pipelines, and thus may be removed in the coming weeks.

## Future Handlers

The following handlers are not yet implemented, but will come online in the coming weeks (from August 17th, 2016).

### [GBS_Demultiplexing](https://github.com/MorrellLab/sequence_handling/wiki/GBS_Demultiplexing)

The GBS_Demultiplexing handler will demulitplex raw GBS reads into split FastQ files using the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). Code for GBS_Demultiplexing is already present in `sequence_handling` for the purposes of collaborative alpha testing, but should not be used since it's extremely buggy.

### Other suggestions

If you feel like there are tools or alternative processing techniques missing from `sequence_handling`, you may submit a feature request through the Git issues page.

## Troubleshooting

Please check the [FAQ and Troubleshooting](https://github.com/MorrellLAB/sequence_handling/wiki/FAQ) page on the wiki for help and visit the wiki page for the handler(s) you're using. If you are still having difficulties, please submit an issue through the [Git issues page](https://github.com/MorrellLab/sequence_handling/issues) for this repository.

## Citation

`sequence_handling` can be cited like:

> Hoffman PJ, Wyant SR, Kono TJY, Morrell PL. (2016). sequence_handling: A pipeline to automate sequence workflows [Software]. Available at https://github.com/MorrellLAB/sequence_handling.