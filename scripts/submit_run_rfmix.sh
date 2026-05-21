#!/bin/bash
set -euo pipefail

TARGET=${1:-results/rfmix/run_rfmix.done}
shift || true

./submit.sh --snakefile rules/run_rfmix.smk "$TARGET" "$@"
