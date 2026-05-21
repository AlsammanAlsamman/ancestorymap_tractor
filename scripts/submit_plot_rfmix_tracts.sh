#!/bin/bash
set -euo pipefail

TARGET=${1:-results/rfmix_plots/plot_rfmix_tracts.done}
shift || true

./submit.sh --snakefile rules/plot_rfmix_tracts.smk "$TARGET" "$@"
