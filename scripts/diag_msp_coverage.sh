#!/bin/bash
#SBATCH --job-name=diag_msp
#SBATCH --output=logs/diag_msp_%j.out
#SBATCH --error=logs/diag_msp_%j.err
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00

set -uo pipefail
cd /s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor
ml bcftools 2>/dev/null || true

echo "=== MSP window coverage per chromosome ==="
printf "%-4s  %-10s  %-14s  %-14s  %-14s\n" "chr" "msp_windows" "bp_covered" "phased_snps" "snps_in_windows"

for CHR in 10 11 12 13 14 15; do
    MSP="results_hispanic/rfmix/chr${CHR}.deconvoluted.msp.tsv"
    PHASED_VCF="results_hispanic/cohort_phasing/chr${CHR}.phased.vcf.gz"
    TRACTOR_D="results_hispanic/tractor/chr${CHR}.phased.anc0.dosage.txt"

    # Count MSP windows and bp covered (skip header lines starting with # or 'chm')
    if [[ -f "$MSP" ]]; then
        read WIN_COUNT BP_COVERED < <(awk 'NR>1 && $1!~/^#/ && $1!="chm" {w++; b+=($3-$2)} END{print w+0, b+0}' "$MSP")
    else
        WIN_COUNT="MISSING"; BP_COVERED="MISSING"
    fi

    # Count phased VCF SNPs (bcftools index -s gives chr\tlen\tnrecords)
    if [[ -f "${PHASED_VCF}.tbi" ]]; then
        PHASED_SNPS=$(bcftools index -s "$PHASED_VCF" 2>/dev/null | awk '{s+=$3} END{print s+0}')
    elif [[ -f "$PHASED_VCF" ]]; then
        PHASED_SNPS=$(bcftools view -H "$PHASED_VCF" 2>/dev/null | wc -l)
    else
        PHASED_SNPS="MISSING"
    fi

    # Count tractor dosage rows (= SNPs that actually made it through)
    if [[ -f "$TRACTOR_D" ]]; then
        TRACTOR_SNPS=$(awk 'NR>1{c++} END{print c+0}' "$TRACTOR_D")
    else
        TRACTOR_SNPS="MISSING"
    fi

    printf "%-4s  %-10s  %-14s  %-14s  %-14s\n" \
        "$CHR" "$WIN_COUNT" "$BP_COVERED" "$PHASED_SNPS" "$TRACTOR_SNPS"
done

echo ""
echo "=== MSP first 5 data lines for chr10 (shows coordinate range) ==="
MSP10="results_hispanic/rfmix/chr10.deconvoluted.msp.tsv"
if [[ -f "$MSP10" ]]; then
    awk 'NR<=6' "$MSP10" | cut -f1-5
    echo "..."
    awk 'END{print "last line:"; print}' "$MSP10" | cut -f1-5
fi

echo ""
echo "=== Phased VCF chr10: first 3 and last 3 variant positions ==="
VCF10="results_hispanic/cohort_phasing/chr10.phased.vcf.gz"
if [[ -f "$VCF10" ]]; then
    bcftools view -H "$VCF10" 2>/dev/null | awk 'NR<=3{print $1,$2} END{print "...last:", $1,$2}'
fi
