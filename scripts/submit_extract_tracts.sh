#!/bin/bash
set -euo pipefail

TARGET=${1:-results/tractor/extract_tracts.done}
shift || true

./submit.sh --snakefile rules/extract_tracts.smk "$TARGET" "$@"
