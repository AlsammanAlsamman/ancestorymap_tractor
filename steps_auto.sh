#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${RESULTS_DIR:-results_hispanic_mod}"
SLEEP_SECONDS="${SLEEP_SECONDS:-900}"
MAX_WAIT_CYCLES="${MAX_WAIT_CYCLES:-0}"

has_submission_error() {
    local submit_log="$1"

    grep -Eq \
        'Full Traceback|Traceback \(most recent call last\)|Error in rule|Exception in line|FileNotFoundError in line|Snakemake workflow failed for target\(s\)' \
        "$submit_log"
}

wait_for_done() {
    local done_file="$1"
    local step_name="$2"
    local cycles=0

    until [[ -f "$done_file" ]]; do
        printf '[%s] Waiting for %s to finish. Missing done file: %s\n' "$(date '+%F %T')" "$step_name" "$done_file"
        cycles=$((cycles + 1))
        if [[ "$MAX_WAIT_CYCLES" -gt 0 && "$cycles" -ge "$MAX_WAIT_CYCLES" ]]; then
            printf '[%s] ERROR: Timed out waiting for %s after %s checks.\n' "$(date '+%F %T')" "$step_name" "$cycles" >&2
            printf '[%s] Tip: increase MAX_WAIT_CYCLES or inspect logs for the failing step.\n' "$(date '+%F %T')" >&2
            exit 1
        fi
        sleep "$SLEEP_SECONDS"
    done

    printf '[%s] Confirmed done file for %s: %s\n' "$(date '+%F %T')" "$step_name" "$done_file"
}

run_step() {
    local step_name="$1"
    local snakefile="$2"
    local target="$3"
    local jobs="$4"
    local done_file="$5"
    local previous_done_file="${6:-}"
    local previous_step_name="${7:-previous step}"

    if [[ -n "$previous_done_file" ]]; then
        wait_for_done "$previous_done_file" "$previous_step_name"
    fi

    if [[ -f "$done_file" ]]; then
        printf '[%s] Skipping %s because done file already exists: %s\n' "$(date '+%F %T')" "$step_name" "$done_file"
        return
    fi

    printf '[%s] Submitting %s\n' "$(date '+%F %T')" "$step_name"

    mkdir -p "$PROJECT_DIR/logs/steps_auto"
    local safe_name
    safe_name="$(printf '%s' "$step_name" | tr '[:upper:]' '[:lower:]' | tr ' :' '__')"
    local submit_log="$PROJECT_DIR/logs/steps_auto/${safe_name}_$(date '+%Y%m%d_%H%M%S').log"

    if [[ "$jobs" == "" ]]; then
        ./submit.sh --snakefile "$snakefile" "$target" 2>&1 | tee "$submit_log"
    else
        ./submit.sh --snakefile "$snakefile" "$target" --jobs "$jobs" 2>&1 | tee "$submit_log"
    fi

    if has_submission_error "$submit_log"; then
        printf '[%s] ERROR: %s submission reported a Snakemake error.\n' "$(date '+%F %T')" "$step_name" >&2
        printf '[%s] See submission log: %s\n' "$(date '+%F %T')" "$submit_log" >&2
        exit 1
    fi

    wait_for_done "$done_file" "$step_name"
}

cd "$PROJECT_DIR"

STEP1_DONE="$RESULTS_DIR/reference_panel_prep/build_reference_sample_lists.done"
STEP2_DONE="$RESULTS_DIR/reference_panel_subset/subset_reference_panel.done"
STEP3_DONE="$RESULTS_DIR/cohort_region_subset/prepare_cohort_regions_from_plink.done"
STEP3C_DONE="$RESULTS_DIR/cohort_region_subset_harmonized/harmonize_cohort_alleles.done"
STEP3B_DONE="$RESULTS_DIR/maf_summary/calculate_maf.done"
STEP4_DONE="$RESULTS_DIR/cohort_phasing/phase_cohort_for_rfmix.done"
STEP5_DONE="$RESULTS_DIR/rfmix/run_rfmix.done"
STEP7_DONE="$RESULTS_DIR/tractor/extract_tracts.done"
STEP8_DONE="$RESULTS_DIR/tractor_gwas/tractor_gwas_pairwise.done"
STEP9_DONE="$RESULTS_DIR/maf_gwas_summary/merge_maf_gwas_summary.done"
STEP10_DONE="$RESULTS_DIR/locus_ancestry_report/locus_ancestry_report.done"
STEP11_DONE="$RESULTS_DIR/admixture_validation/validate_admixture_prediction.done"
STEP12_DONE="$RESULTS_DIR/admixture_importance/report_admixture_importance.done"

run_step \
    "Step 1: build reference sample lists" \
    "rules/example_complete.smk" \
    "build_reference_sample_lists" \
    "" \
    "$STEP1_DONE"

run_step \
    "Step 2: subset reference panel" \
    "rules/subset_reference_panel.smk" \
    "subset_reference_panel_done" \
    "22" \
    "$STEP2_DONE" \
    "$STEP1_DONE" \
    "Step 1"

run_step \
    "Step 3: prepare cohort regions from plink" \
    "rules/prepare_cohort_regions_from_plink.smk" \
    "prepare_cohort_regions_from_plink_done" \
    "22" \
    "$STEP3_DONE" \
    "$STEP2_DONE" \
    "Step 2"

run_step \
    "Step 3c: harmonize cohort alleles" \
    "rules/harmonize_cohort_alleles.smk" \
    "harmonize_cohort_alleles_done" \
    "22" \
    "$STEP3C_DONE" \
    "$STEP3_DONE" \
    "Step 3"

run_step \
    "Step 3b: calculate maf" \
    "rules/calculate_maf.smk" \
    "calculate_maf_done" \
    "22" \
    "$STEP3B_DONE" \
    "$STEP3C_DONE" \
    "Step 3c"

run_step \
    "Step 4: phase cohort for rfmix" \
    "rules/phase_cohort_for_rfmix.smk" \
    "phase_cohort_for_rfmix_done" \
    "22" \
    "$STEP4_DONE" \
    "$STEP3C_DONE" \
    "Step 3c"

run_step \
    "Step 5: run rfmix" \
    "rules/run_rfmix.smk" \
    "run_rfmix_done" \
    "22" \
    "$STEP5_DONE" \
    "$STEP4_DONE" \
    "Step 4"

run_step \
    "Step 7: extract tracts" \
    "rules/extract_tracts.smk" \
    "extract_tracts_done" \
    "22" \
    "$STEP7_DONE" \
    "$STEP5_DONE" \
    "Step 5"

run_step \
    "Step 8: tractor gwas pairwise" \
    "rules/tractor_gwas_pairwise.smk" \
    "tractor_gwas_pairwise_done" \
    "40" \
    "$STEP8_DONE" \
    "$STEP7_DONE" \
    "Step 7"

run_step \
    "Step 9: merge maf and gwas summary" \
    "rules/merge_maf_gwas_summary.smk" \
    "merge_maf_gwas_summary_done" \
    "1" \
    "$STEP9_DONE" \
    "$STEP8_DONE" \
    "Step 8"

run_step \
    "Step 10: build locus ancestry report" \
    "rules/report_locus_ancestry.smk" \
    "build_locus_ancestry_report" \
    "1" \
    "$STEP10_DONE" \
    "$STEP9_DONE" \
    "Step 9"

run_step \
    "Step 11: validate admixture prediction" \
    "rules/validate_admixture_prediction.smk" \
    "validate_admixture_prediction_done" \
    "1" \
    "$STEP11_DONE" \
    "$STEP10_DONE" \
    "Step 10"

run_step \
    "Step 12: report admixture importance" \
    "rules/report_admixture_importance.smk" \
    "report_admixture_importance_done" \
    "1" \
    "$STEP12_DONE" \
    "$STEP11_DONE" \
    "Step 11"

printf '[%s] All configured pipeline steps completed successfully.\n' "$(date '+%F %T')"