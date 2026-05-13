#!/usr/bin/env bash
set -euo pipefail

cd /s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor

if type ml >/dev/null 2>&1; then
  ml bcftools || true
fi

RESULTS_DIR="results_hispanic_mod"
CHR="11"
POS="128605169"

echo "RESULTS_DIR=${RESULTS_DIR}"

echo "-- cohort_region_subset chr${CHR}.loci --"
if [[ -f "${RESULTS_DIR}/cohort_region_subset/chr${CHR}.loci.vcf.gz" ]]; then
  bcftools view -r "${CHR}:${POS}-${POS}" "${RESULTS_DIR}/cohort_region_subset/chr${CHR}.loci.vcf.gz" | grep -v '^#' || true
else
  echo "missing ${RESULTS_DIR}/cohort_region_subset/chr${CHR}.loci.vcf.gz"
fi

echo "-- reference_panel_subset chr${CHR}.* --"
shopt -s nullglob
for f in "${RESULTS_DIR}/reference_panel_subset/chr${CHR}."*.vcf.gz; do
  echo "file=${f}"
  bcftools view -r "${CHR}:${POS}-${POS}" "$f" | grep -v '^#' || true
done
shopt -u nullglob

echo "-- cohort_phasing chr${CHR}.phased --"
if [[ -f "${RESULTS_DIR}/cohort_phasing/chr${CHR}.phased.vcf.gz" ]]; then
  bcftools view -r "${CHR}:${POS}-${POS}" "${RESULTS_DIR}/cohort_phasing/chr${CHR}.phased.vcf.gz" | grep -v '^#' || true
else
  echo "missing ${RESULTS_DIR}/cohort_phasing/chr${CHR}.phased.vcf.gz"
fi

echo "-- per-source MAF rows --"
for t in \
  "${RESULTS_DIR}/maf_summary/cohort.maf.tsv" \
  "${RESULTS_DIR}/maf_summary/AFR.maf.tsv" \
  "${RESULTS_DIR}/maf_summary/AMR.maf.tsv" \
  "${RESULTS_DIR}/maf_summary/EUR.maf.tsv" \
  "${RESULTS_DIR}/maf_summary/EAS.maf.tsv"; do
  [[ -f "$t" ]] || continue
  echo "table=$t"
  rg -n "^${CHR}\t${POS}\t" "$t" || true
done

echo "-- merged.annotated.maf.tsv row --"
if [[ -f "${RESULTS_DIR}/maf_summary/merged.annotated.maf.tsv" ]]; then
  rg -n "^${CHR}\t${POS}\trs1236176\b" "${RESULTS_DIR}/maf_summary/merged.annotated.maf.tsv" || true
fi

echo "-- GWAS pair presence --"
shopt -s nullglob
found=0
for g in "${RESULTS_DIR}/tractor_gwas/merged.model_"*.gwas.tsv; do
  if rg -n "^${CHR}\t${POS}\trs1236176\b" "$g" >/tmp/snp_hit.txt 2>/dev/null; then
    echo "present in $(basename "$g")"
    cat /tmp/snp_hit.txt
    found=1
  fi
done
[[ $found -eq 0 ]] && echo "not present in any merged.model_*.gwas.tsv"
shopt -u nullglob

echo "-- merged.maf_gwas_summary.tsv row --"
if [[ -f "${RESULTS_DIR}/maf_gwas_summary/merged.maf_gwas_summary.tsv" ]]; then
  rg -n "^${CHR}\t${POS}\trs1236176\b" "${RESULTS_DIR}/maf_gwas_summary/merged.maf_gwas_summary.tsv" || true
fi
