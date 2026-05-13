#!/bin/bash
#SBATCH --job-name=diag_snp_loss
#SBATCH --output=logs/diag_snp_loss_%j.out
#SBATCH --error=logs/diag_snp_loss_%j.err
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

set -euo pipefail

WDIR=/s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor
cd "$WDIR"

# Load tools
ml bcftools 2>/dev/null || true

CHROMS="10 11 12 13 14 15"

LOCI_DIR=results_hispanic/cohort_region_subset
HARM_DIR=results_hispanic/cohort_region_subset_harmonized
PHASE_DIR=results_hispanic/cohort_phasing
TRACTOR_DIR=results_hispanic/tractor
PLINK_DIR=results_hispanic/ancestry_plink_dosage/AFR   # all ancestries identical

OUT=results_hispanic/snp_count_diagnostic/diag_snp_loss_per_stage.tsv

mkdir -p results_hispanic/snp_count_diagnostic logs

printf "%-4s\t%-14s\t%-14s\t%-14s\t%-18s\t%-14s\t%s\n" \
    "chr" "loci_vcf" "harmonized_vcf" "phased_vcf" "tractor_dosage_anc0" "plink_pvar_AFR" "note" \
    > "$OUT"

for CHR in $CHROMS; do

    # --- Stage 1: loci VCF (bcftools index -s → col3 = nrecords) ---
    F1="${LOCI_DIR}/chr${CHR}.loci.vcf.gz"
    if [[ -f "${F1}.tbi" ]]; then
        N1=$(bcftools index -s "$F1" 2>/dev/null | awk '{s+=$3} END{print s+0}')
    elif [[ -f "$F1" ]]; then
        N1=$(bcftools view -H "$F1" 2>/dev/null | wc -l)
    else
        N1="MISSING"
    fi

    # --- Stage 2: harmonized VCF ---
    F2="${HARM_DIR}/chr${CHR}.loci.harmonized.vcf.gz"
    if [[ -f "${F2}.tbi" ]]; then
        N2=$(bcftools index -s "$F2" 2>/dev/null | awk '{s+=$3} END{print s+0}')
    elif [[ -f "$F2" ]]; then
        N2=$(bcftools view -H "$F2" 2>/dev/null | wc -l)
    else
        N2="MISSING"
    fi

    # --- Stage 3: phased VCF ---
    F3="${PHASE_DIR}/chr${CHR}.phased.vcf.gz"
    if [[ -f "${F3}.tbi" ]]; then
        N3=$(bcftools index -s "$F3" 2>/dev/null | awk '{s+=$3} END{print s+0}')
    elif [[ -f "$F3" ]]; then
        N3=$(bcftools view -H "$F3" 2>/dev/null | wc -l)
    else
        N3="MISSING"
    fi

    # --- Stage 4: Tractor dosage anc0 (line count minus header) ---
    F4="${TRACTOR_DIR}/chr${CHR}.phased.anc0.dosage.txt"
    if [[ -f "$F4" ]]; then
        N4=$(awk 'NR>1{c++} END{print c+0}' "$F4")
    else
        N4="MISSING"
    fi

    # --- Stage 5: PLINK pvar AFR (non-header lines) ---
    F5="${PLINK_DIR}/chr${CHR}.dosage.pvar"
    if [[ -f "$F5" ]]; then
        N5=$(grep -c -v '^#' "$F5" || true)
    else
        N5="MISSING"
    fi

    # Detect which step the big drop occurs
    NOTE=""
    if [[ "$N1" =~ ^[0-9]+$ ]] && [[ "$N2" =~ ^[0-9]+$ ]]; then
        DIFF12=$(( N1 - N2 ))
        if (( DIFF12 > 10000 )); then NOTE="${NOTE}DROP_AT_HARMONIZE(${DIFF12}) "; fi
    fi
    if [[ "$N2" =~ ^[0-9]+$ ]] && [[ "$N3" =~ ^[0-9]+$ ]]; then
        DIFF23=$(( N2 - N3 ))
        if (( DIFF23 > 10000 )); then NOTE="${NOTE}DROP_AT_PHASE(${DIFF23}) "; fi
    fi
    if [[ "$N3" =~ ^[0-9]+$ ]] && [[ "$N4" =~ ^[0-9]+$ ]]; then
        DIFF34=$(( N3 - N4 ))
        if (( DIFF34 > 10000 )); then NOTE="${NOTE}DROP_AT_TRACTOR(${DIFF34}) "; fi
    fi
    if [[ "$N4" =~ ^[0-9]+$ ]] && [[ "$N5" =~ ^[0-9]+$ ]]; then
        DIFF45=$(( N4 - N5 ))
        if (( DIFF45 > 10000 )); then NOTE="${NOTE}DROP_AT_PLINK_EXPORT(${DIFF45}) "; fi
    fi
    [[ -z "$NOTE" ]] && NOTE="no_large_drop_detected"

    printf "%-4s\t%-14s\t%-14s\t%-14s\t%-18s\t%-14s\t%s\n" \
        "$CHR" "$N1" "$N2" "$N3" "$N4" "$N5" "$NOTE" \
        >> "$OUT"

    echo "chr${CHR}: loci=$N1  harm=$N2  phased=$N3  tractor_anc0=$N4  plink_AFR=$N5  NOTE=$NOTE"
done

echo ""
echo "=== FULL TABLE ==="
column -t "$OUT"
echo ""
echo "Output saved to: $OUT"
