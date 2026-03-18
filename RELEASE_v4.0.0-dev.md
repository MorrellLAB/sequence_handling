# sequence_handling v4.0.0-dev Release Notes

Status: draft (dev branch)
Date: 2026-03-18

This document summarizes major updates currently staged on the `dev` branch for the next `sequence_handling` release.

## Highlights

- Migrated active fastp configuration files from GATK 4.1 to GATK 4.6.
- Added a new accessory helper for generating all-sites VCFs for pixy v2.0 workflows.
- Updated non-config scripts to remove hardcoded GATK 4.1.8 executable paths where applicable.
- Verified that gzipped FASTQ concatenation in the active long-read trimming path uses `zcat`.

## Detailed Changes

### 1) GATK update in active fastp configs

The following config files were updated to use GATK 4.6:

- `Config_fastp`
  - `GATK_JAR=/panfs/jay/groups/9/morrellp/public/Software/gatk-4.6.0.0/gatk`
  - `module load gatk/4.6.0`
- `Config_fastp_test`
  - `GATK_JAR=/panfs/jay/groups/9/morrellp/public/Software/gatk-4.6.0.0/gatk`
  - `module load gatk/4.6.0`

Note: Older legacy config files were intentionally left unchanged.

### 2) Removed hardcoded legacy GATK path in handler logic

`Handlers/Genomics_DB_Import.sh` was updated to use `${GATK_JAR}` rather than a hardcoded `/panfs/.../gatk-4.1.8.0/gatk` path in all execution branches.

This makes GenomicsDBImport behavior consistent with config-driven reproducibility and version pinning.

### 3) Updated helper module pin

`HelperScripts/combine_and_sort_split_vcf.sh` now loads:

- `module load gatk/4.6.0`

instead of `gatk/4.1.2`.

`HelperScripts/Intervals_at_Ns.sh` was also updated to a GATK 4.6 default module/path example.

### 4) Added pixy all-sites accessory script

New file:

- `HelperScripts/generate_pixy_allsites_vcf.sh`

Supported methods:

- `bcftools` mode: runs `mpileup -> call -A -> filter -> norm` and outputs indexed `.vcf.gz`
- `gatk` mode: runs `GenotypeGVCFs --include-non-variant-sites true` from either:
  - a GenomicsDB workspace (`gendb://...`), or
  - a list of gVCFs

The script validates inputs, prints tool metadata, and creates a tabix index for output.

### 5) FASTQ concatenation behavior

In the active long-read fastp workflow (`Handlers/Fastplong.sh`), gzipped FASTQ concatenation is done with `zcat` into a FIFO stream.

Additionally, a repository search found no remaining `cat ...fastq.gz` pattern in `Handlers/`.

### 6) GATK 4.6 indel/SV assessment

- Current germline SNP/indel workflow remains valid (`Haplotype_Caller -> Genomics_DB_Import -> Genotype_GVCFs`).
- No new mandatory germline indel caller integration is required for this release.
- GATK 4.6 includes HaplotypeCaller bug fixes (including long-deletion edge cases), so no handler redesign is needed, but regression checks are recommended for representative long-indel regions.
- GATK 4.6 SV updates are primarily in supporting SV tooling (annotation/concordance), not a drop-in replacement for a full SV calling branch in this pipeline.
- Conclusion: no release-blocking workflow changes required for indel/SV based on GATK 4.6 updates; optional SV branch design can be planned separately.

## Compatibility Notes

- Ensure `gatk/4.6.0` and the `GATK_JAR` path are available on your cluster module stack.
- The pixy helper requires:
  - `bcftools` for indexing in both modes
  - `gatk` for GATK mode

## Migration Notes

- Existing runs pinned to older versions remain reproducible with archived historical configs.
- New dev-branch runs should archive:
  - the exact config file,
  - this release note,
  - and the git commit/tag used for execution.
