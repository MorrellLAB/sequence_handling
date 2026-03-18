# pixy all-sites VCF generation

This guide describes how to generate an all-sites VCF (variant + invariant sites) compatible with pixy v2.0 using:

- `HelperScripts/generate_pixy_allsites_vcf.sh`

The script supports both the bcftools and GATK workflows described in the pixy documentation.

## Requirements

- bash
- bcftools (required for both output indexing and bcftools mode)
- GATK (required for gatk mode)
- Reference FASTA file

## Quick start

Run from the repository root:

```bash
chmod +x HelperScripts/generate_pixy_allsites_vcf.sh
```

### Method 1: bcftools mode (from BAMs)

Inputs:

- BAM list file (one BAM path per line)
- Reference FASTA

Command:

```bash
HelperScripts/generate_pixy_allsites_vcf.sh bcftools \
  --bam-list /path/to/bams.txt \
  --reference /path/to/reference.fa \
  --output /path/to/output/all_sites.bcftools.vcf.gz \
  --threads 8 \
  --min-mapq 20 \
  --min-baseq 20 \
  --min-qual 20
```

Output:

- bgzipped VCF: `all_sites.bcftools.vcf.gz`
- tabix index: `all_sites.bcftools.vcf.gz.tbi`

### Method 2A: GATK mode from GenomicsDB workspace

Inputs:

- GenomicsDB workspace path
- Reference FASTA

Command:

```bash
HelperScripts/generate_pixy_allsites_vcf.sh gatk \
  --reference /path/to/reference.fa \
  --gendb-workspace /path/to/combinedDB/gendb_wksp \
  --output /path/to/output/all_sites.gatk.vcf.gz \
  --gatk /path/to/gatk
```

### Method 2B: GATK mode from gVCF list

Inputs:

- gVCF list file (one gVCF path per line)
- Reference FASTA

Command:

```bash
HelperScripts/generate_pixy_allsites_vcf.sh gatk \
  --reference /path/to/reference.fa \
  --gvcf-list /path/to/gvcfs.txt \
  --output /path/to/output/all_sites.gatk.vcf.gz \
  --gatk /path/to/gatk
```

## Notes

- For GATK mode, the script uses `GenotypeGVCFs --include-non-variant-sites true`.
- For bcftools mode, the script uses `mpileup` and `call -A` to retain non-variant sites.
- If `--gatk` is omitted, the script uses `GATK_JAR` from the environment when set, otherwise it falls back to `gatk` on PATH.
- The script exits early with clear errors when required files or commands are missing.

## Help

```bash
HelperScripts/generate_pixy_allsites_vcf.sh --help
```
