#!/bin/bash
set -euo pipefail

TARGET=${1:-results/rfmix_bin_case_control/plot_rfmix_bins_case_control.done}
shift || true

./submit.sh --snakefile rules/plot_rfmix_bins_case_control.smk "$TARGET" "$@"
