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

### 1) New long-read variant calling with Clair3 (Handler 13)

New handler: `Handlers/Clair3_Variant_Calling.sh`

- Enables SNP/indel calling on long-read BAM files (ONT, PacBio HiFi, ONT Q20)
- Automatically selects Clair3 model based on `SEQ_PLATFORM` and optional `CLAIR3_BASECALLER` variables
- Outputs indexed VCF suitable for downstream filtering and joint calling with GLnexus
- Supports basecaller tracking for ONT reads (guppy, dorado)
- Includes per-sample logging and error handling
- SLURM job script: `SlurmJobScripts/Clair3_Variant_Calling.job`

Config parameters:
- `CLAIR3_BASECALLER` (optional, for ONT reads)
- `CLAIR3_GENERATE_ALL_SITES` (boolean, optional, for pixy AllSites VCF generation)

### 2) SNP-only variant calling with bcftools mpileup (Handler 14)

New handler: `Handlers/Bcftools_Mpileup_Variant_Calling.sh`

- Lightweight alternative for SNP-only variant calling from short-read BAMs
- Runs after BAM files are sorted and indexed
- Uses `bcftools mpileup` followed by filtering to SNPs only
- Outputs indexed VCF ready for joint calling or downstream filtering
- Suitable for workflows where indel calling is not required
- SLURM job script: `SlurmJobScripts/Bcftools_Mpileup_Variant_Calling.job`

Config parameters:
- `FINISHED_BAM_LIST` (path to file listing BAM files to process)

### 3) Joint variant calling with GLnexus (Handler 15)

New handler: `Handlers/GLnexus_Joint_Calling.sh`

- Consolidates single-sample VCFs from Clair3, bcftools mpileup, or other callers into joint calls
- Converts VCF inputs to BCF, runs GLnexus, outputs joint VCF
- Recommended for:
  - Clair3 long-read variant consolidation (instead of GATK Genomics_DB_Import)
  - bcftools mpileup multi-sample consolidation
  - Ultima Genomics UG100 joint calling via sequence_accessories
- Supports configurable GLnexus config (default: DeepVariant)
- SLURM job script: `SlurmJobScripts/GLnexus_Joint_Calling.job`

Config parameters:
- `VCF_LIST` (path to file listing single-sample VCFs)
- `GLNEXUS_CONFIG` (GLnexus config name, optional)

### 4) Updated sequence_handling_fastp menu and handler numbering

Cleaned up and reorganized the workflow:

- **Removed deprecated handlers:**
  - GBS_Demultiplex (was 13)
  - All Nanopore Workflow options (1NP, 2NP, 3NP, 4NP) — now integrated into main workflow
  
- **Added new handler tiers:**
  - Alternative variant calling branch (Clair3, bcftools mpileup, GLnexus)
  
- **Updated handler numbering:**
  - Handlers 1-12: Standard short-read GATK workflow (unchanged)
  - Handlers 13-15: Alternative variant calling (new)
  - Handlers 16-18: Other utilities (Quality_Trimming, Realigner_Target_Creator, Indel_Realigner)

### 5) GATK update in active fastp configs

The following config files were updated to use GATK 4.6:

- `Config_fastp`
  - `GATK_JAR=/panfs/jay/groups/9/morrellp/public/Software/gatk-4.6.0.0/gatk`
  - `module load gatk/4.6.0`
- `Config_fastp_test`
  - `GATK_JAR=/panfs/jay/groups/9/morrellp/public/Software/gatk-4.6.0.0/gatk`
  - `module load gatk/4.6.0`

Note: Older legacy config files were intentionally left unchanged.

### 6) Removed hardcoded legacy GATK path in handler logic

`Handlers/Genomics_DB_Import.sh` was updated to use `${GATK_JAR}` rather than a hardcoded `/panfs/.../gatk-4.1.8.0/gatk` path in all execution branches.

This makes GenomicsDBImport behavior consistent with config-driven reproducibility and version pinning.

### 7) Updated helper module pin

`HelperScripts/combine_and_sort_split_vcf.sh` now loads:

- `module load gatk/4.6.0`

instead of `gatk/4.1.2`.

`HelperScripts/Intervals_at_Ns.sh` was also updated to a GATK 4.6 default module/path example.

### 8) Added pixy all-sites accessory script

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
