#!/usr/bin/env bash
set -euo pipefail

cd /s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor
if type ml >/dev/null 2>&1; then
  ml bcftools || true
fi

echo "cohort"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -r 11:128605169-128605169 results_hispanic_mod/cohort_region_subset/chr11.loci.vcf.gz || true

echo "reference"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -r 11:128605169-128605169 results_hispanic_mod/reference_panel_subset/chr11.YRI_PEL_CEU_CHB.vcf.gz || true

echo "phased"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -r 11:128605169-128605169 results_hispanic_mod/cohort_phasing/chr11.phased.vcf.gz || true
