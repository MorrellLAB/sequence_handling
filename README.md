# `sequence_handling`
[![DOI](https://zenodo.org/badge/36831916.svg)](https://zenodo.org/badge/latestdoi/36831916)
#### A pipeline to automate DNA sequence aligning and quality control workflows via list-based batch submission and parallel processing

For the latest updates and to chat with our team, please join our [Slack workspace: sequencehandling.slack.com](https://join.slack.com/t/sequencehandling/shared_invite/enQtNzQxNTc3NDE2NDUwLTlhYmFkOTM1NWQ1YjBjOGY5MmFlMmQ5MDVlNmM4NThhZjEzNzdkZjg1YWFlYjA1NTNmMTgwMjQ2YjY1NGE1ZDc).

___

> For greater detail, usage information, and troubleshooting please see the [`sequence_handling` wiki](https://github.com/MorrellLAB/sequence_handling/wiki).

## What is `sequence_handling` for?

`sequence_handling` is a series of scripts, called handlers, that automate and speed up DNA sequence alignment and quality control. Currently, `sequence_handling` is designed to work with Illumina paired-end whole-genome or exome capture sequences. Some parts of `sequence_handling` can accept GBS sequences.

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
 - Tools for variant calling, such as the [Genome Analysis Toolkit](https://www.broadinstitute.org/gatk/)
 - Tools for filtering and manipulating VCF files, such as [VCFtools](https://vcftools.github.io/man_latest.html) and [vcflib](https://github.com/vcflib/vcflib)

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

For any handler that utilizes PBS job arrays, there is an optional flag, `-t custom_array_indices` that you can use to re-run specific job array indices that errored out or were aborted. The `custom_array_indices` is a range of arrays and/or comma separated list of specific arrays to run WITHOUT spaces in between. Currently, the `-t` flag must be provided as the 3rd argument on the command line (see example below). If left blank (you do not use the `-t` flag), the DEFAULT runs all samples in your sample list. This is helpful if only some of your jobs arrays fail and you need to re-run only those.

Here is an example using the -t flag and SAM_Processing handler:

```bash
./sequence_handling SAM_Processing /path/to/config -t 1-5,10,12
```

## Recommended Workflow

![Workflow](https://github.com/MorrellLAB/sequence_handling/blob/master/.workflow_images/Sequence_Handling_Workflow.png)

#### 1. [Quality\_Assessment](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Assessment)

To start, run Quality_Assessment on your raw FastQ files. The Quality_Assessment handler runs [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on a series of samples and outputs metrics used for quality control. It accepts FASTQ, SAM, and BAM files as input and outputs a summary table and individual HTML files for visualization. The Quality_Assessment handler depends on [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 2. [Adapter\_Trimming](https://github.com/MorrellLab/sequence_handling/wiki/Adapter_Trimming)

The Adapter_Trimming handler uses [Scythe](https://github.com/vsbuffalo/scythe) to trim specific adapter sequences from FastQ files. This handler differentiates between forward, reverse, and single-end FastQ files automatically. The Adapter_Trimming handler depends on [Scythe](https://github.com/vsbuffalo/scythe) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### [Quality\_Assessment](https://github.com/MorrellLab/sequence_handling/wiki/Quality_Assessment)

After Adapter_Trimming, it is recommended to run Quality_Assessment again on the trimmed FastQ files to ensure that all adapter contamination was properly removed. 

#### 3. [Read\_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Read_Mapping)

The Read_Mapping handler maps sequence reads to a reference genome using [BWA-MEM](http://bio-bwa.sourceforge.net/). This handler uses Torque Task Arrays, part of the [Portable Batch System](http://www.pbsworks.com/). The Read_Mapping handler depends on the [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/).

#### 4. [SAM\_Processing with Picard](https://github.com/MorrellLAB/sequence_handling/wiki/SAM_Processing)

The SAM_Processing handler converts the SAM files from read mapping with [BWA](http://bio-bwa.sourceforge.net/) to the BAM format using [SAMTools](http://www.htslib.org/). In the conversion process, it will sort and deduplicate the data for the finished BAM file, also using [SAMTools](http://www.htslib.org/). Alignment statistics will also be generated for both raw and finished BAM files. The SAM_Processing handler depends on [SAMTools](http://www.htslib.org/) and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 5. [Coverage_Mapping](https://github.com/MorrellLab/sequence_handling/wiki/Coverage_Mapping)

The Coverage_Mapping handler generates coverage histograms and summary statistics from BAM files using [BEDTools](http://bedtools.readthedocs.org/en/latest/). Plots of coverage are generated using [R](http://cran.r-project.org/) based on coverage maps. The Coverage_Mapping handler depends on [BEDTools](http://bedtools.readthedocs.org/en/latest/), [R](http://cran.r-project.org/), and [GNU Parallel](http://www.gnu.org/software/parallel/).

#### 6. [Haplotype_Caller](https://github.com/MorrellLab/sequence_handling/wiki/Haplotype_Caller)

To begin the variant discovery process from your finished BAM files, the Haplotype_Caller handler uses [GATK](https://software.broadinstitute.org/gatk/) to generate genomic VCF files for each sample.

#### 7. [Genomics_DB_Import](https://github.com/MorrellLAB/sequence_handling/wiki/Genomics_DB_Import)

Import GVCF files output from Haplotype_Caller into a GenomicsDB workspace. Genomics_DB_Import in GATK 4 merges GVCF files from multiple samples before joint genotyping (Note: it has the same functionality as CombineGVCFs in previous versions of GATK). Because this step pools together all of your samples into one file, it is **essential that all samples are included for this step**. Automatically breaking the process into chromosome parts or smaller regions for each chromosome allows the job to be "parallelized across regions" by running as a task array and speeds up computing time.

#### 8. [Genotype_GVCFs](https://github.com/MorrellLab/sequence_handling/wiki/Genotype_GVCFs)

The Genotype_GVCFs hander converts the GVCF files for the entire dataset into VCF files broken up by chromosome or chromosome part using [GATK](https://software.broadinstitute.org/gatk/). Breaking the output into chromosome parts allows the process to be split into a task array and greatly speeds up processing time.

#### 9. [Create_HC_Subset](https://github.com/MorrellLab/sequence_handling/wiki/Create_HC_Subset)

The Create_HC_Subset handler creates a single VCF file that contains only the high-confidence sites for your samples. This filtering is performed in multiple steps using several different user-defined parameters and before-and-after percentile tables are generated. Create_HC_Subset depends on [VCFtools](https://vcftools.github.io/man_latest.html) and [vcflib](https://github.com/vcflib/vcflib) for manipulating the VCF file. 

#### 10. [Variant_Recalibrator](https://github.com/MorrellLAB/sequence_handling/wiki/Variant_Recalibrator)

The Variant_Recalibrator handler uses the [GATK](https://software.broadinstitute.org/gatk/) and user-provided prior sets of "truth" variants to create a model that attempts to separate true variants from false positives. An unfiltered VCF file the the FILTER field annotated is generated.

#### 11. [Variant_Filtering](https://github.com/MorrellLab/sequence_handling/wiki/Variant_Filtering)

The Variant_Filtering handler creates a single variant call format (VCF) file that contains only high-quality sites and genotypes for your samples. This filtering is performed in multiple steps using several different user-defined parameters and before-and-after percentile tables are generated. Variant_Filtering depends on [VCFtools](https://vcftools.github.io/man_latest.html) and [vcflib](https://github.com/vcflib/vcflib) for manipulating the VCF file. 

#### 12. [Variant_Analysis](https://github.com/MorrellLab/sequence_handling/wiki/Variant_Analysis)

The Variant_Analysis handler uses a variety of dependencies to produce statistics about the input VCF file. Information generated by the handler includes heterozygosity summaries, missing-ness summaries, a minor allele frequency histogram, the Ts/Tv ratio, and the raw count of SNPs. Additional information is output for barley samples. Variant_Analysis depends on [VCFtools](https://vcftools.github.io/man_latest.html), [vcflib](https://github.com/vcflib/vcflib), [molpopgen](https://github.com/molpopgen/analysis), [Python3](https://www.python.org/), [GNU Parallel](http://www.gnu.org/software/parallel/), [BCFtools](https://samtools.github.io/bcftools/bcftools.html), [R](https://www.r-project.org/), [TeX Live](https://www.tug.org/texlive/), and the [Enthought Python Distribution](https://www.enthought.com/product/enthought-python-distribution/). 

## Troubleshooting

Please check the [Common Problems and Errors](https://github.com/MorrellLAB/sequence_handling/wiki/FAQ) page on the wiki for help and visit the [wiki](https://github.com/MorrellLAB/sequence_handling/wiki) page for the handler(s) you're using. If you are still having difficulties, please submit an issue through the [Git issues page](https://github.com/MorrellLab/sequence_handling/issues) for this repository.

## Future Handlers

The following handlers are planned for future versions of `sequence_handling`.

#### GBS_Demultiplexing

The GBS_Demultiplexing handler will demulitplex raw GBS reads into split FastQ files using the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). Code for GBS_Demultiplexing is already present in `sequence_handling` for the purposes of collaborative alpha testing, but should not be used since it's extremely buggy.

#### Coverage Mapping histograms using R

#### Longranger support for 10x reads

#### Tools for QC and mapping long nanopore reads with minimap2

#### Other suggestions

If you feel like there are tools or alternative processing techniques missing from `sequence_handling`, please submit a feature request through the Git issues page.

## Citation

`sequence_handling` can be cited like:

**Version 2:**

> Paul J. Hoffman, Skylar R. Wyant, Thomas J.Y. Kono, & Peter L. Morrell. (2018, June 1). MorrellLAB/sequence_handling: Release v2.0: SNP calling with GATK 3.8 (Version v2.0). Zenodo. http://doi.org/10.5281/zenodo.1257692

**Version 1:**

> Hoffman PJ, Wyant SR, Kono TJY, Morrell PL. (2018). sequence_handling: A pipeline to automate sequence workflows. Zenodo. https://doi.org/10.5281/zenodo.1257692
