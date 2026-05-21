#!/bin/bash
set -euo pipefail

ROOT="/s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor"
OUT="${ROOT}/results_hispanic/snp_count_diagnostic/chr10_15_step_counts.tsv"
mkdir -p "$(dirname "$OUT")"

cd "$ROOT"

# Load modules requested by user
if type ml >/dev/null 2>&1; then
  ml slurm >/dev/null 2>&1 || true
  ml bcftools >/dev/null 2>&1 || true
fi
if type module >/dev/null 2>&1; then
  module load slurm >/dev/null 2>&1 || true
  module load bcftools >/dev/null 2>&1 || true
fi

echo -e "chr\toriginal_bim\tharmonized_vcf\tphased_vcf\trfmix_msp_rows\tdosage_txt_rows_anc0\tpvar_rows_AFR\tpvar_rows_AMR\tpvar_rows_EAS\tpvar_rows_EUR" > "$OUT"

# Precompute original BIM SNP counts for chr10-15 in one pass (faster than rescanning each loop).
declare -A ORIG_BY_CHR
while IFS=$'\t' read -r c n; do
  ORIG_BY_CHR["$c"]="$n"
done < <(awk '$1>=10 && $1<=15{n[$1]++} END{for(i=10;i<=15;i++) printf "%d\t%d\n", i, n[i]+0}' inputs/Hispanic.bim)

count_vcf_records() {
  local f="$1"
  # If CSI/TBI exists, this is near-instant. Fallback to full scan otherwise.
  if [[ -f "${f}.csi" || -f "${f}.tbi" ]]; then
    bcftools index -n "$f"
  else
    bcftools view -H "$f" | wc -l
  fi
}

echo "Starting per-chromosome diagnostics..."

for CHR in 10 11 12 13 14 15; do
  ORIG="${ORIG_BY_CHR[$CHR]:-0}"
  echo "[chr${CHR}] counting records..."

  HARM="results_hispanic/cohort_region_subset_harmonized/chr${CHR}.loci.harmonized.vcf.gz"
  PHASED="results_hispanic/cohort_phasing/chr${CHR}.phased.vcf.gz"
  MSP="results_hispanic/rfmix/chr${CHR}.deconvoluted.msp.tsv"
  DOSAGE="results_hispanic/tractor/chr${CHR}.phased.anc0.dosage.txt"

  HARM_N=$(count_vcf_records "$HARM")
  PHASED_N=$(count_vcf_records "$PHASED")
  MSP_N=$(grep -v '^#' "$MSP" | awk 'NF>0{n++} END{print n+0}')
  DOSAGE_N=$(grep -v '^#' "$DOSAGE" | awk 'NR>1 && NF>0{n++} END{print n+0}')

  AFR_N=$(awk 'END{print NR-1}' "results_hispanic/ancestry_plink_dosage/AFR/chr${CHR}.dosage.pvar")
  AMR_N=$(awk 'END{print NR-1}' "results_hispanic/ancestry_plink_dosage/AMR/chr${CHR}.dosage.pvar")
  EAS_N=$(awk 'END{print NR-1}' "results_hispanic/ancestry_plink_dosage/EAS/chr${CHR}.dosage.pvar")
  EUR_N=$(awk 'END{print NR-1}' "results_hispanic/ancestry_plink_dosage/EUR/chr${CHR}.dosage.pvar")

  echo -e "${CHR}\t${ORIG}\t${HARM_N}\t${PHASED_N}\t${MSP_N}\t${DOSAGE_N}\t${AFR_N}\t${AMR_N}\t${EAS_N}\t${EUR_N}" >> "$OUT"
done

echo "Wrote $OUT"
