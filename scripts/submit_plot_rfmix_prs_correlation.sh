#!/bin/bash
set -euo pipefail

TARGET=${1:-results/rfmix_prs_correlation/plot_rfmix_prs_correlation.done}
shift || true

./submit.sh --snakefile rules/plot_rfmix_prs_correlation.smk "$TARGET" "$@"
