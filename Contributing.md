* * *

# `sequence_handling` Design Principles

This document records the guiding principles for developing and extending the [`sequence_handling`](https://www.google.com/search?q=%5Bhttps://github.com/MorrellLAB/sequence_handling%5D\(https://github.com/MorrellLAB/sequence_handling\)) pipeline. It is intended as an enduring reference for contributors, particularly during the ongoing modernization of the pipeline to support current GATK versions, long-read sequencing technologies (ONT and PacBio HiFi), and updated tooling.

These principles apply to all new handlers, to revisions of existing handlers, and to supporting infrastructure such as config files and Slurm job scripts.

* * *

## 1. Reproducibility Through Configuration and Versioning

The combination of a config file and a tagged version of `sequence_handling` should be sufficient to reproduce any analysis exactly. This means:

- All file paths, tool parameters, reference genome locations, and sample lists must come from the config file. Hard-coded paths are not permitted.
- Every release should be tagged (e.g., `v3.0`) and archived (e.g., via Zenodo) so that a specific version can always be retrieved.
- The config file should be treated as a scientific record and archived alongside results, just as raw data are.
- Tool version numbers used during a run should be captured in logs (see Principle 4).

Reproducibility is the baseline expectation of cumulative science: if a result cannot be reproduced from the config and version alone, something is missing.

* * *

## 2. UNIX Philosophy: One Handler, One Job

Each handler should do one thing and do it well. This is the UNIX philosophy applied to genomic workflows.

- A handler encapsulates a single logical step in the pipeline (e.g., read mapping, variant calling). It should not silently perform steps that belong to another handler.
- Avoid handler creep: if new functionality is needed, evaluate whether it belongs in an existing handler or warrants a new one.
- Handlers should accept well-defined inputs (paths from the config) and produce well-defined outputs (files in a predictable location and format).
- Composability matters: handlers should chain cleanly, with the outputs of one serving as the inputs of the next without manual intervention.

This principle also constrains scope: if another handler already covers a tool, do not duplicate it. The pipeline's modularity is a feature, not a limitation.

* * *

## 3. Architecture and Terminology

Consistent terminology prevents ambiguity across documentation, code, and issues.

- **Handlers** (in `Handlers/`) are the scripts that run computational jobs. Each handler corresponds to one pipeline step and is submitted to the cluster.
- **SlurmJobScripts** (in `SlurmJobScripts/`) are the scripts that construct Slurm job headers (resource requests, partition, memory, time limits, etc.) and wrap handler submission. They do not contain analysis logic or values for memory and time limits; these are in the `Config` files.
- **Config files** (e.g., `Config`) define all parameters for a run. A handler should source a single config and rely entirely on the variables defined there.
- **HelperScripts** (in `HelperScripts`) are scripts for handling multiple samples.
- **Sequence_Accessories** (e.g., `PanDepthCoverage.sh`) are tools that may be used occasionally or that supplement `sequence_handling`.
- This separation of concerns -- job logic in Handlers, scheduling logic in SlurmJobScripts, parameters in Config -- should be preserved as the pipeline grows.

When adding support for new sequencing platforms or tools, follow this multi-layer structure. Do not mix scheduling directives into handler logic.

* * *

## 4. Shell Quality and Safety

All shell code must be safe, readable, and linter-compliant.

- Use strict shell options at the top of every handler: `set -euo pipefail`. This ensures the script exits on errors, treats unset variables as errors, and catches failures in pipelines.
- All handlers should pass [`shellcheck`](https://www.google.com/search?q=%5Bhttps://www.shellcheck.net/%5D\(https://www.shellcheck.net/\)) without warnings. `shellcheck` compliance is the baseline standard.
- Where applicable, code should also satisfy [DeepSource](https://deepsource.com/) static analysis alerts.
- Variables that come from the config should be validated before use (e.g., check that a file path is non-empty and the file exists before passing it to a tool).
- Avoid `eval`, unquoted variable expansions in paths, and other patterns that are fragile or unsafe in HPC environments.

**Variable Scoping**

In Bash, variables are global by default. To prevent variables from leaking and causing unintended side effects, scope them locally within functions:

Bash
    
    
    # Define variables inside functions as local
    function parse_bam() {
        local input_bam="$1"
        local sample_name="$2"
        # ...
    }
    

**Exception Handling and Cleanup**

Bash lacks standard `try/catch` blocks, so errors must be handled explicitly. Catch potential failures, provide informative error messages to `stderr`, and use `trap` to ensure temporary files are cleaned up even if the script fails prematurely.

Bash
    
    
    # Explicit error handling for tool execution
    samtools index "${input_bam}" || {
        echo "[ERROR] Failed to index ${input_bam}. Exiting..." >&2
        exit 1
    }
    
    # Ensure temporary directories are removed on exit or failure
    local temp_dir
    temp_dir=$(mktemp -d)
    trap 'rm -rf "${temp_dir}"' EXIT
    

Code that runs silently and produces incorrect results is worse than code that fails loudly. Prefer explicit failure over silent success.

* * *

## 5. Logging, Metadata, and Observability

Handlers should produce output that is useful for both humans and AI-assisted troubleshooting.

- Log the tool version, key parameters, input files, and output files at the start and end of each handler run. This information should appear in both stdout and the Slurm log.
- Capture metrics that are meaningful for QC: read counts, mapping rates, coverage depth, duplication rates, variant counts, and so forth -- whichever are appropriate for the handler's task.
- Use structured, grep-friendly log lines where possible (e.g., `[Fastplong] Sample: ${SAMPLE} | Reads before: ${N_BEFORE} | Reads after: ${N_AFTER}`).
- Do not suppress stderr from tools unless you have explicitly handled the error conditions. Suppressed errors are a common source of silent failures in pipelines.
- Metadata captured at runtime (tool versions, parameters, timestamps) should be written to a per-sample or per-run log file that persists after the job completes.

The goal is that any run -- successful or failed -- leaves enough of a record to diagnose what happened without re-running it.
* * *

## 6. Dependency and Compatibility Checking

Known incompatibilities between tools, modules, or versions should be caught early — before a job completes and produces unusable output.
- Where known version constraints exist between tools (e.g., GATK version and Java version, samtools and htslib, model compatibility in Clair3), add an explicit check at the start of the handler and exit with an informative error if the constraint is violated.
- Dependency checks should test for both tool availability (`command -v tool || { echo "Error: tool not found"; exit 1; }`) and, where feasible, version compatibility.
- Document known incompatibilities as comments in the handler code and in the relevant config stanza.
- As GATK and other tools release breaking changes (e.g., GATK 4.x argument syntax changes), handlers should be updated and the incompatibility documented in the changelog.

Failing fast with a clear error message saves far more time than diagnosing a corrupted VCF hours into a run.

* * *

## 7. Development Workflow and Releases

All active development occurs on the `dev` branch.

- New handlers, refactored handlers, and bug fixes should be developed on `dev` (or on feature branches that merge into `dev`). The `main` branch reflects the last stable release.
- `dev` is the integration target: it should remain functional and pass basic testing at all times, even if individual features are incomplete.
- When a meaningful set of changes has accumulated -- for example, the addition of long-read support, or a GATK version bump -- prepare a new versioned release. Update the changelog, tag the release, and archive via Zenodo so the version is citable.
- Pull requests into `dev` should include: updated config stanzas (if new parameters are introduced), handler-level comments explaining non-obvious logic, and a brief description in the PR of what changed and why.
- Breaking changes to the config format or handler interface should be clearly flagged in the changelog and in a migration note for existing users.

The `dev` branch is where the next version of `sequence_handling` is being built. Treat it accordingly.

* * *

## 8. Documentation

- Documentation of the workflow is provided in the `README.md` file. More specifics are available in a [wiki](https://github.com/MorrellLAB/sequence_handling/wiki).
- [`Dependencies`](https://www.google.com/search?q=%5Bhttps://github.com/MorrellLAB/sequence_handling/wiki/Dependencies%5D\(https://github.com/MorrellLAB/sequence_handling/wiki/Dependencies\)) are listed on a dedicated page, which should be updated so users can easily identify the tools needed for successful execution.
- Where necessary, clarifying comments should be included in code, particularly for more complex operations.
* * *

## Relationship to the Original Goals

These principles extend and operationalize the original design goals of `sequence_handling` as stated in the project README:

> _The workflow is intended to be 100% reproducible, provided that you have the Config file and the version of `sequence_handling` that was used. It is also intended to be easy for beginner UNIX users to configure and run independently._

The principles above are the mechanisms by which those goals are achieved and maintained as the pipeline evolves to support new sequencing platforms (ONT, PacBio HiFi), updated variant callers (GATK 4.x, Clair3, pbsv, Sniffles2), and the cluster environments in which it runs.

* * *

_Last updated: March 2026. To be revised as the pipeline evolves._

* * *

## Next steps

- Update from GATK v4.1 to GATK v4.6.

- Add an accessory script to generate an AllSites VCF for [pixy v2.0](https://github.com/ksamuk/pixy). The specific examples for using bcftools mpileup and GATK are [here](https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html).

- Document all changes to the new version in a Release file similar to that for [v3.0.0](https://github.com/MorrellLAB/sequence_handling/releases/tag/v3.0.0). Recent updates to the dev branch include replacing' fastp' and' fastplong' with `fastp` and `fastplong` for quality assessment and adapter trimming. This replaces the full front end of the workflow. We have also added `minimap2` for long-read mapping and updated the `config` to include read-mapping presets for PacBio (HiFi reads) and ONT (Q20 reads).

- Make sure that the concatenation of gzipped fastq files uses `zcat` rather than `cat`. They don't produce the same results.


* * *
