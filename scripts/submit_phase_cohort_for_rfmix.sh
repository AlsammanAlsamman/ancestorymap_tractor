#!/bin/bash
set -euo pipefail

TARGET=${1:-results/cohort_phasing/phase_cohort_for_rfmix.done}
shift || true

./submit.sh --snakefile rules/phase_cohort_for_rfmix.smk "$TARGET" "$@"
