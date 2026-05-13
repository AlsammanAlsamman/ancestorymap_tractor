#!/bin/bash
set -euo pipefail

TARGET=${1:-results/cohort_region_subset/prepare_cohort_regions_from_plink.done}
shift || true

./submit.sh --snakefile rules/prepare_cohort_regions_from_plink.smk "$TARGET" "$@"
