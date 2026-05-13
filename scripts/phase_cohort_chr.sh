#!/bin/bash
set -euo pipefail

INPUT_VCF=""
REFERENCE_VCF=""
GENETIC_MAP=""
CHR=""
OUTPUT_VCF=""
THREADS="2"
BCFTOOLS_MODULE=""
SHAPEIT5_MODULE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-vcf)
            INPUT_VCF="$2"
            shift 2
            ;;
        --reference-vcf)
            REFERENCE_VCF="$2"
            shift 2
            ;;
        --genetic-map)
            GENETIC_MAP="$2"
            shift 2
            ;;
        --chr)
            CHR="$2"
            shift 2
            ;;
        --output-vcf)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --bcftools-module)
            BCFTOOLS_MODULE="$2"
            shift 2
            ;;
        --shapeit5-module)
            SHAPEIT5_MODULE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "$INPUT_VCF" || -z "$REFERENCE_VCF" || -z "$GENETIC_MAP" || -z "$CHR" || -z "$OUTPUT_VCF" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

for f in "$INPUT_VCF" "$REFERENCE_VCF" "$GENETIC_MAP"; do
    if [[ ! -f "$f" ]]; then
        echo "Required file not found: $f" >&2
        exit 1
    fi
done

if [[ ! -f "${INPUT_VCF}.tbi" && ! -f "${INPUT_VCF}.csi" ]]; then
    echo "Input cohort VCF index missing for: ${INPUT_VCF}" >&2
    exit 1
fi

if [[ ! -f "${REFERENCE_VCF}.tbi" && ! -f "${REFERENCE_VCF}.csi" ]]; then
    echo "Reference VCF index missing for: ${REFERENCE_VCF}" >&2
    exit 1
fi

have_working_bcftools=false
if command -v bcftools >/dev/null 2>&1; then
    if bcftools --version >/dev/null 2>&1; then
        have_working_bcftools=true
    fi
fi

have_working_shapeit5=false
SHAPEIT5_BIN=""
if command -v phase_common_static >/dev/null 2>&1; then
    if phase_common_static --help >/dev/null 2>&1; then
        have_working_shapeit5=true
        SHAPEIT5_BIN="phase_common_static"
    fi
fi

if [[ "$have_working_shapeit5" == false ]] && command -v phase_common >/dev/null 2>&1; then
    if phase_common --help >/dev/null 2>&1; then
        have_working_shapeit5=true
        SHAPEIT5_BIN="phase_common"
    fi
fi

if [[ "$have_working_bcftools" == false || "$have_working_shapeit5" == false ]]; then
    if type module >/dev/null 2>&1; then
        if [[ "$have_working_bcftools" == false && -n "$BCFTOOLS_MODULE" ]]; then
            module load "$BCFTOOLS_MODULE"
        fi
        if [[ "$have_working_shapeit5" == false && -n "$SHAPEIT5_MODULE" ]]; then
            module load "$SHAPEIT5_MODULE"
        fi
    fi
fi

if ! command -v bcftools >/dev/null 2>&1 || ! bcftools --version >/dev/null 2>&1; then
    echo "bcftools not runnable (PATH/module issue). Module requested: ${BCFTOOLS_MODULE:-none}" >&2
    exit 1
fi

if command -v phase_common_static >/dev/null 2>&1 && phase_common_static --help >/dev/null 2>&1; then
    SHAPEIT5_BIN="phase_common_static"
elif command -v phase_common >/dev/null 2>&1 && phase_common --help >/dev/null 2>&1; then
    SHAPEIT5_BIN="phase_common"
else
    echo "SHAPEIT5 phasing executable not runnable (expected phase_common_static or phase_common). Module requested: ${SHAPEIT5_MODULE:-none}" >&2
    exit 1
fi

OUT_DIR="$(dirname "$OUTPUT_VCF")"
mkdir -p "$OUT_DIR"

WORK_DIR="$(mktemp -d -p "$OUT_DIR" ".tmp_phase_chr${CHR}.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

TMP_TAGGED="$WORK_DIR/chr${CHR}.prephase.vcf.gz"
TMP_BCF="$WORK_DIR/chr${CHR}.phased.bcf"
PHASE_INPUT="$INPUT_VCF"

if bcftools +fill-tags "$INPUT_VCF" -Oz -o "$TMP_TAGGED" -- -t "AN,AC"; then
    bcftools index -t "$TMP_TAGGED"
    PHASE_INPUT="$TMP_TAGGED"
else
    echo "Warning: bcftools +fill-tags failed; proceeding with original input VCF." >&2
fi

"$SHAPEIT5_BIN" \
    --input "$PHASE_INPUT" \
    --reference "$REFERENCE_VCF" \
    --region "$CHR" \
    --map "$GENETIC_MAP" \
    --output "$TMP_BCF" \
    --thread "$THREADS"

bcftools convert -Oz -o "$OUTPUT_VCF" "$TMP_BCF"
bcftools index -t "$OUTPUT_VCF"

echo "Created phased VCF: $OUTPUT_VCF"
