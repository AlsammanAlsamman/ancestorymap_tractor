#!/bin/bash
set -euo pipefail

INPUT_VCF=""
SAMPLE_LIST=""
OUTPUT_VCF=""
BCFTOOLS_MODULE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-vcf)
            INPUT_VCF="$2"
            shift 2
            ;;
        --sample-list)
            SAMPLE_LIST="$2"
            shift 2
            ;;
        --output-vcf)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        --bcftools-module)
            BCFTOOLS_MODULE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "$INPUT_VCF" || -z "$SAMPLE_LIST" || -z "$OUTPUT_VCF" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "Input VCF not found: $INPUT_VCF" >&2
    exit 1
fi

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "Sample list not found: $SAMPLE_LIST" >&2
    exit 1
fi

if [[ ! -f "${INPUT_VCF}.tbi" && ! -f "${INPUT_VCF}.csi" ]]; then
    echo "Input index not found (.tbi or .csi): ${INPUT_VCF}" >&2
    exit 1
fi

have_working_bcftools=false
if command -v bcftools >/dev/null 2>&1; then
    if bcftools --version >/dev/null 2>&1; then
        have_working_bcftools=true
    fi
fi

if [[ "$have_working_bcftools" == false ]]; then
    if type module >/dev/null 2>&1 && [[ -n "$BCFTOOLS_MODULE" ]]; then
        module load "$BCFTOOLS_MODULE"
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        echo "bcftools not found in PATH (module attempted: ${BCFTOOLS_MODULE:-none})" >&2
        exit 1
    fi

    if ! bcftools --version >/dev/null 2>&1; then
        echo "bcftools exists but is not runnable (likely shared library mismatch)." >&2
        echo "Module attempted: ${BCFTOOLS_MODULE:-none}" >&2
        echo "Try updating configs/software.yml bcftools.module to a compatible version (e.g., bcftools/1.15)." >&2
        exit 1
    fi
fi

mkdir -p "$(dirname "$OUTPUT_VCF")"

bcftools view \
    -S "$SAMPLE_LIST" \
    -Oz \
    -o "$OUTPUT_VCF" \
    "$INPUT_VCF"

bcftools index -t "$OUTPUT_VCF"

echo "Created: $OUTPUT_VCF"
