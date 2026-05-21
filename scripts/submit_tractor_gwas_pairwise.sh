#!/bin/bash
set -euo pipefail

TARGET=${1:-results/tractor_gwas/tractor_gwas_pairwise.done}
shift || true

./submit.sh --snakefile rules/tractor_gwas_pairwise.smk "$TARGET" "$@"
