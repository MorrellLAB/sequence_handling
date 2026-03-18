#!/usr/bin/env bash

# Generate an all-sites VCF (variant + invariant sites) suitable for pixy v2.0.
# Supports two methods:
#   1) bcftools mpileup/call (from BAMs)
#   2) GATK GenotypeGVCFs with --include-non-variant-sites (from gVCFs or GenomicsDB)

set -euo pipefail

function usage() {
    cat <<'EOF'
Usage:
  generate_pixy_allsites_vcf.sh bcftools \
    --bam-list BAM_LIST \
    --reference REF_FASTA \
    --output OUT_VCF_GZ \
    [--threads 4] [--min-mapq 20] [--min-baseq 20] [--min-qual 20]

  generate_pixy_allsites_vcf.sh gatk \
    --reference REF_FASTA \
    --output OUT_VCF_GZ \
    [--gatk /path/to/gatk] \
    (--gendb-workspace GENOMICSDB_PATH | --gvcf-list GVCF_LIST)

Notes:
  - Output is bgzipped VCF (.vcf.gz) and indexed with bcftools index -t.
  - For bcftools mode, BAM list must contain one BAM path per line.
  - For gatk mode, gVCF list must contain one gVCF path per line.
EOF
}

function require_cmd() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    }
}

function ensure_readable_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || {
        echo "[ERROR] ${label} not found or empty: ${path}" >&2
        exit 1
    }
}

function ensure_readable_lines() {
    local list_file="$1"
    local label="$2"
    ensure_readable_file "${list_file}" "${label} list"
    while IFS= read -r line || [[ -n "${line}" ]]; do
        [[ -z "${line}" ]] && continue
        [[ -s "${line}" ]] || {
            echo "[ERROR] ${label} path from list is not readable: ${line}" >&2
            exit 1
        }
    done < "${list_file}"
}

[[ $# -ge 1 ]] || { usage; exit 1; }

method="$1"
shift

threads=4
min_mapq=20
min_baseq=20
min_qual=20

bam_list=""
reference=""
output=""
gendb_workspace=""
gvcf_list=""
gatk_bin="${GATK_JAR:-gatk}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-list)
            bam_list="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --min-mapq)
            min_mapq="$2"
            shift 2
            ;;
        --min-baseq)
            min_baseq="$2"
            shift 2
            ;;
        --min-qual)
            min_qual="$2"
            shift 2
            ;;
        --gendb-workspace)
            gendb_workspace="$2"
            shift 2
            ;;
        --gvcf-list)
            gvcf_list="$2"
            shift 2
            ;;
        --gatk)
            gatk_bin="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

[[ -n "${reference}" ]] || { echo "[ERROR] --reference is required" >&2; usage; exit 1; }
[[ -n "${output}" ]] || { echo "[ERROR] --output is required" >&2; usage; exit 1; }
ensure_readable_file "${reference}" "Reference FASTA"

mkdir -p "$(dirname "${output}")"

case "${method}" in
    bcftools)
        [[ -n "${bam_list}" ]] || { echo "[ERROR] --bam-list is required for bcftools mode" >&2; usage; exit 1; }
        ensure_readable_lines "${bam_list}" "BAM"
        require_cmd bcftools

        echo "[pixy-allsites] Method: bcftools"
        echo "[pixy-allsites] Reference: ${reference}"
        echo "[pixy-allsites] BAM list: ${bam_list}"
        echo "[pixy-allsites] Output: ${output}"
        echo "[pixy-allsites] bcftools version: $(bcftools --version | head -n 1)"

        # Emit all sites (-A), then retain the full-site representation for pixy input.
        bcftools mpileup \
            --threads "${threads}" \
            -f "${reference}" \
            -b "${bam_list}" \
            -a AD,DP \
            -q "${min_mapq}" \
            -Q "${min_baseq}" \
            -Ou | \
        bcftools call \
            --threads "${threads}" \
            -m \
            -A \
            -Ou | \
        bcftools filter \
            --threads "${threads}" \
            -i "QUAL>=${min_qual}" \
            -Ou | \
        bcftools norm \
            --threads "${threads}" \
            -f "${reference}" \
            -Oz \
            -o "${output}"

        bcftools index --threads "${threads}" -t "${output}"
        ;;
    gatk)
        require_cmd bcftools

        [[ -n "${gendb_workspace}" || -n "${gvcf_list}" ]] || {
            echo "[ERROR] For gatk mode, provide --gendb-workspace or --gvcf-list" >&2
            usage
            exit 1
        }

        if [[ -n "${gendb_workspace}" && -n "${gvcf_list}" ]]; then
            echo "[ERROR] Use either --gendb-workspace or --gvcf-list, not both" >&2
            exit 1
        fi

        if [[ -n "${gvcf_list}" ]]; then
            ensure_readable_lines "${gvcf_list}" "gVCF"
        fi

        if [[ -x "${gatk_bin}" ]]; then
            :
        else
            require_cmd "${gatk_bin}"
        fi

        echo "[pixy-allsites] Method: gatk"
        echo "[pixy-allsites] Reference: ${reference}"
        echo "[pixy-allsites] Output: ${output}"
        echo "[pixy-allsites] GATK binary: ${gatk_bin}"

        if [[ -n "${gendb_workspace}" ]]; then
            [[ -d "${gendb_workspace}" ]] || {
                echo "[ERROR] GenomicsDB workspace not found: ${gendb_workspace}" >&2
                exit 1
            }

            "${gatk_bin}" GenotypeGVCFs \
                -R "${reference}" \
                -V "gendb://${gendb_workspace}" \
                --include-non-variant-sites true \
                -O "${output}"
        else
            declare -a gvcf_args=()
            while IFS= read -r gvcf || [[ -n "${gvcf}" ]]; do
                [[ -z "${gvcf}" ]] && continue
                gvcf_args+=( -V "${gvcf}" )
            done < "${gvcf_list}"

            "${gatk_bin}" GenotypeGVCFs \
                -R "${reference}" \
                "${gvcf_args[@]}" \
                --include-non-variant-sites true \
                -O "${output}"
        fi

        bcftools index --threads "${threads}" -t "${output}"
        ;;
    *)
        echo "[ERROR] Unknown method: ${method}" >&2
        usage
        exit 1
        ;;
esac

echo "[pixy-allsites] Completed: ${output}"
