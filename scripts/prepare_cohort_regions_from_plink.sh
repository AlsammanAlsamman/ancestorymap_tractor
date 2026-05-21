#!/bin/bash
set -euo pipefail

PLINK_PREFIX=""
LOCI_FILE=""
CHR=""
OUTPUT_VCF=""
PLINK_MODULE=""
BCFTOOLS_MODULE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --plink-prefix)
            PLINK_PREFIX="$2"
            shift 2
            ;;
        --loci-file)
            LOCI_FILE="$2"
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
        --plink-module)
            PLINK_MODULE="$2"
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

if [[ -z "$PLINK_PREFIX" || -z "$LOCI_FILE" || -z "$CHR" || -z "$OUTPUT_VCF" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

for ext in bed bim fam; do
    if [[ ! -f "${PLINK_PREFIX}.${ext}" ]]; then
        echo "Missing PLINK file: ${PLINK_PREFIX}.${ext}" >&2
        exit 1
    fi
done

if [[ ! -f "$LOCI_FILE" ]]; then
    echo "Loci file not found: $LOCI_FILE" >&2
    exit 1
fi

if type module >/dev/null 2>&1; then
    if [[ -n "$PLINK_MODULE" ]]; then
        module load "$PLINK_MODULE"
    fi
    if [[ -n "$BCFTOOLS_MODULE" ]]; then
        module load "$BCFTOOLS_MODULE"
    fi
fi

if ! command -v plink2 >/dev/null 2>&1; then
    echo "plink2 not found in PATH" >&2
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "bcftools not found in PATH" >&2
    exit 1
fi

OUT_DIR="$(dirname "$OUTPUT_VCF")"
mkdir -p "$OUT_DIR"

WORK_DIR="$(mktemp -d -p "$OUT_DIR" ".tmp_region_chr${CHR}.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

RANGE_FILE="$WORK_DIR/chr${CHR}.ranges.txt"
EXTRACT_IDS="$WORK_DIR/chr${CHR}.extract_ids.txt"
TMP_PREFIX="$WORK_DIR/chr${CHR}.subset"
TMP_OUT_PREFIX="$WORK_DIR/chr${CHR}.region"
TMP_EXPORTED_VCF="${TMP_OUT_PREFIX}.vcf.gz"

awk -v chr="$CHR" 'BEGIN{FS="[ \t]+"; OFS="\t"}
NR==1 { next }
{
    locus=$1;
    chrom=$2;
    gsub(/^chr/, "", chrom);
    start=$3;
    end=$4;
    if (chrom != chr) next;
    if (start == "" || end == "") next;
    if (start > end) {
        tmp = start;
        start = end;
        end = tmp;
    }
    if (start < 1) start = 1;
    print chrom, start, end, locus;
}' "$LOCI_FILE" > "$RANGE_FILE"

if [[ ! -s "$RANGE_FILE" ]]; then
    echo "No loci boundaries found for chr${CHR} in $LOCI_FILE" >&2
    exit 1
fi

# Build variant ID list from BIM for any position falling within one of the loci boundaries.
awk -v chr="$CHR" 'BEGIN{FS="[ \t]+"; OFS="\t"}
FNR==NR {
    n++;
    start[n]=$2;
    end[n]=$3;
    next
}
{
    chrom=$1;
    gsub(/^chr/, "", chrom);
    if (chrom != chr) next;
    pos=$4;
    for (i=1; i<=n; i++) {
        if (pos>=start[i] && pos<=end[i]) {
            print $2;
            break;
        }
    }
}' "$RANGE_FILE" "${PLINK_PREFIX}.bim" | sort -u > "$EXTRACT_IDS"

if [[ ! -s "$EXTRACT_IDS" ]]; then
    echo "No variants found in PLINK BIM for chr${CHR} within requested loci boundaries" >&2
    exit 1
fi

plink2 \
    --bfile "$PLINK_PREFIX" \
    --chr "$CHR" \
    --extract "$EXTRACT_IDS" \
    --make-bed \
    --out "$TMP_PREFIX"

plink2 \
    --bfile "$TMP_PREFIX" \
    --export vcf bgz \
    --out "$TMP_OUT_PREFIX"

bcftools sort -Oz -o "$OUTPUT_VCF" "$TMP_EXPORTED_VCF"
bcftools index -t "$OUTPUT_VCF"

echo "Created loci-based region VCF: $OUTPUT_VCF"
