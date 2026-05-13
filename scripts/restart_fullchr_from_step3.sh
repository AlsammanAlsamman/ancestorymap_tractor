#!/bin/bash
#SBATCH --job-name=restart_fullchr_s3
#SBATCH --output=logs/restart_fullchr_s3_%j.out
#SBATCH --error=logs/restart_fullchr_s3_%j.err
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00

set -euo pipefail

cd /s/nath-lab/alsamman/____MyCodes____/ancestorymap_tractor

# Required by user workflow
ml slurm 2>/dev/null || true

# Rebuild from Step 3 onward using full chromosome loci ranges
./submit.sh --snakefile rules/prepare_cohort_regions_from_plink.smk prepare_cohort_regions_from_plink_done --jobs 22
./submit.sh --snakefile rules/harmonize_cohort_alleles.smk harmonize_cohort_alleles_done --jobs 22
./submit.sh --snakefile rules/calculate_maf.smk calculate_maf_done --jobs 22
./submit.sh --snakefile rules/phase_cohort_for_rfmix.smk phase_cohort_for_rfmix_done --jobs 22
./submit.sh --snakefile rules/run_rfmix.smk run_rfmix_done --jobs 22
./submit.sh --snakefile rules/merge_rfmix_outputs.smk merge_rfmix_outputs_done --jobs 1
./submit.sh --snakefile rules/extract_tracts.smk extract_tracts_done --jobs 22
./submit.sh --snakefile rules/export_ancestry_dosage_plink.smk ancestry_dosage_plink_done --subjob-time 72:00:00 --jobs 22
./submit.sh --snakefile rules/snp_count_diagnostic.smk --force snp_count_diagnostic_done --jobs 1
