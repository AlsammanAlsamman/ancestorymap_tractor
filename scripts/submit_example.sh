#!/bin/bash
set -euo pipefail

TARGET=${1:-results/reference_panel_prep/build_reference_sample_lists.done}
shift || true

./submit.sh --snakefile rules/example_complete.smk "$TARGET" "$@"
