#!/bin/bash
set -euo pipefail

TARGET=${1:-results/maf_summary/calculate_maf.done}
shift || true

./submit.sh --snakefile rules/calculate_maf.smk "$TARGET" "$@"
