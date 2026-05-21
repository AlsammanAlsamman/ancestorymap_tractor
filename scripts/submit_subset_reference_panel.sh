#!/bin/bash
set -euo pipefail

TARGET=${1:-results/reference_panel_subset/subset_reference_panel.done}
shift || true

./submit.sh --snakefile rules/subset_reference_panel.smk "$TARGET" "$@"
