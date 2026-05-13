#!/bin/bash
set -euo pipefail

# Build a dedicated chr6 no-HLA VCF from the harmonized cohort subset.
INPUT_VCF="results_hispanic/cohort_region_subset_harmonized/chr6.loci.harmonized.vcf.gz"
OUTPUT_VCF="results_hispanic/cohort_region_subset_harmonized/chr6.loci.harmonized.noHLA.vcf.gz"
HLA_REGION="6:25000000-34000000"

mkdir -p logs "$(dirname "$OUTPUT_VCF")"

if type module >/dev/null 2>&1; then
    module load bcftools || true
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools is not available in PATH." >&2
    exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF" >&2
    exit 1
fi

echo "Extracting chr6 variants excluding HLA region $HLA_REGION from $INPUT_VCF"
bcftools view -t "^$HLA_REGION" -Oz -o "$OUTPUT_VCF" "$INPUT_VCF"
bcftools index -t "$OUTPUT_VCF"

N_VARIANTS=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
HLA_LEFT=$(bcftools view -r "$HLA_REGION" -H "$OUTPUT_VCF" | wc -l)
echo "Done. Wrote $OUTPUT_VCF with $N_VARIANTS variants outside HLA."
echo "Validation: variants still in HLA interval = $HLA_LEFT"
