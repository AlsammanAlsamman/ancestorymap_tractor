#!/bin/bash
set -euo pipefail

PORT="${1:-8765}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cat <<EOF
Serving the project root at:
  http://127.0.0.1:${PORT}/scripts/maf_gwas_summary_viewer.html

Use the "Load default results file" button in the page,
or upload:
  results/maf_gwas_summary/merged.maf_gwas_summary.tsv
EOF

cd "$PROJECT_ROOT"
python3 -m http.server "$PORT"
