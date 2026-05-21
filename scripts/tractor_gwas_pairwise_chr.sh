#!/bin/bash
set -euo pipefail

CHR=""
ANC_HAPCOUNT=""
ANC_I_DOSAGE=""
ANC_J_DOSAGE=""
LABEL_I=""
LABEL_J=""
PHENOTYPE=""
PHENOTYPE_FORMAT=""
FAM_CASE_CODE="2"
FAM_CONTROL_CODE="1"
TABLE_IID_COLUMN="IID"
TABLE_PHENO_COLUMN="PHENO"
PCA_ENABLED="false"
PCA_FILE=""
PCA_FORMAT="plink_eigenvec"
PCA_HAS_HEADER="false"
PCA_FID_COLUMN="FID"
PCA_IID_COLUMN="IID"
PCA_N_PCS="10"
MIN_PARTITIONS="32"
OUTPUT_TSV=""
SCRIPT_PATH=""
PYTHON_MODULE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chr)
            CHR="$2"
            shift 2
            ;;
        --anc-hapcount)
            ANC_HAPCOUNT="$2"
            shift 2
            ;;
        --anc-i-dosage)
            ANC_I_DOSAGE="$2"
            shift 2
            ;;
        --anc-j-dosage)
            ANC_J_DOSAGE="$2"
            shift 2
            ;;
        --label-i)
            LABEL_I="$2"
            shift 2
            ;;
        --label-j)
            LABEL_J="$2"
            shift 2
            ;;
        --phenotype)
            PHENOTYPE="$2"
            shift 2
            ;;
        --phenotype-format)
            PHENOTYPE_FORMAT="$2"
            shift 2
            ;;
        --fam-case-code)
            FAM_CASE_CODE="$2"
            shift 2
            ;;
        --fam-control-code)
            FAM_CONTROL_CODE="$2"
            shift 2
            ;;
        --table-iid-column)
            TABLE_IID_COLUMN="$2"
            shift 2
            ;;
        --table-phenotype-column)
            TABLE_PHENO_COLUMN="$2"
            shift 2
            ;;
        --pca-enabled)
            PCA_ENABLED="$2"
            shift 2
            ;;
        --pca-file)
            PCA_FILE="$2"
            shift 2
            ;;
        --pca-format)
            PCA_FORMAT="$2"
            shift 2
            ;;
        --pca-has-header)
            PCA_HAS_HEADER="$2"
            shift 2
            ;;
        --pca-fid-column)
            PCA_FID_COLUMN="$2"
            shift 2
            ;;
        --pca-iid-column)
            PCA_IID_COLUMN="$2"
            shift 2
            ;;
        --pca-n-pcs)
            PCA_N_PCS="$2"
            shift 2
            ;;
        --min-partitions)
            MIN_PARTITIONS="$2"
            shift 2
            ;;
        --output-tsv)
            OUTPUT_TSV="$2"
            shift 2
            ;;
        --script-path)
            SCRIPT_PATH="$2"
            shift 2
            ;;
        --python-module)
            PYTHON_MODULE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "$CHR" || -z "$ANC_HAPCOUNT" || -z "$ANC_I_DOSAGE" || -z "$ANC_J_DOSAGE" || -z "$LABEL_I" || -z "$LABEL_J" || -z "$PHENOTYPE" || -z "$PHENOTYPE_FORMAT" || -z "$OUTPUT_TSV" || -z "$SCRIPT_PATH" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

for f in "$ANC_HAPCOUNT" "$ANC_I_DOSAGE" "$ANC_J_DOSAGE" "$PHENOTYPE" "$SCRIPT_PATH"; do
    if [[ ! -f "$f" ]]; then
        echo "Required file not found: $f" >&2
        exit 1
    fi
done

if [[ "$PCA_ENABLED" == "true" || "$PCA_ENABLED" == "1" || "$PCA_ENABLED" == "yes" ]]; then
    if [[ -z "$PCA_FILE" || ! -f "$PCA_FILE" ]]; then
        echo "PCA covariates are enabled but the PCA file was not found: $PCA_FILE" >&2
        exit 1
    fi
fi

mkdir -p "$(dirname "$OUTPUT_TSV")"

if ! type module >/dev/null 2>&1; then
    if [[ -f /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    elif [[ -f /usr/share/Modules/init/bash ]]; then
        source /usr/share/Modules/init/bash
    fi
fi

if [[ -n "$PYTHON_MODULE" ]] && type module >/dev/null 2>&1; then
    if ! module load "$PYTHON_MODULE"; then
        echo "Warning: failed to load requested python module '$PYTHON_MODULE'; falling back to python in PATH" >&2
    fi
fi

PYTHON_BIN=""
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
fi

if [[ -z "$PYTHON_BIN" ]]; then
    echo "python3/python is required but was not found in PATH" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEPS_ROOT="${HAIL_USERBASE:-$HOME/.local}"
LOCK_DIR="$PROJECT_ROOT/.deps/.hail_install_lock"
HAIL_INSTALLER="$PROJECT_ROOT/scripts/install_hail_local.sh"

mkdir -p "$DEPS_ROOT" "$PROJECT_ROOT/.deps"

export PYTHONUSERBASE="$DEPS_ROOT"
export PATH="$PYTHONUSERBASE/bin:$PATH"

USER_SITE="$($PYTHON_BIN - <<'PY'
import site
print(site.getusersitepackages())
PY
)"

if [[ -n "$USER_SITE" ]]; then
    if [[ -n "${PYTHONPATH:-}" ]]; then
        export PYTHONPATH="$USER_SITE:$PYTHONPATH"
    else
        export PYTHONPATH="$USER_SITE"
    fi
fi

PYVER="$($PYTHON_BIN - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
HAIL_MARKER="$PROJECT_ROOT/.deps/hail_installed_py${PYVER}.ok"

if ! "$PYTHON_BIN" -c "import hail" >/dev/null 2>&1; then
    while ! mkdir "$LOCK_DIR" 2>/dev/null; do
        sleep 2
    done
    trap 'rmdir "$LOCK_DIR" >/dev/null 2>&1 || true' EXIT

    if [[ ! -f "$HAIL_MARKER" ]] || ! "$PYTHON_BIN" -c "import hail" >/dev/null 2>&1; then
        echo "hail not found in $PYTHON_BIN; installing once into $PYTHONUSERBASE" >&2
        bash "$HAIL_INSTALLER" --python-module "$PYTHON_MODULE" --user-base "$PYTHONUSERBASE" --skip-init-check
        "$PYTHON_BIN" -c "import hail" >/dev/null 2>&1
        date +"%Y-%m-%dT%H:%M:%S" > "$HAIL_MARKER"
    fi

    rmdir "$LOCK_DIR" >/dev/null 2>&1 || true
    trap - EXIT
fi

if ! command -v java >/dev/null 2>&1; then
    if type module >/dev/null 2>&1; then
        module load java >/dev/null 2>&1 || module load jdk >/dev/null 2>&1 || module load openjdk >/dev/null 2>&1 || true
    fi
fi

if command -v java >/dev/null 2>&1; then
    JAVA_BIN_REAL="$(readlink -f "$(command -v java)")"
    if [[ "$JAVA_BIN_REAL" == */bin/java ]]; then
        export JAVA_HOME="${JAVA_BIN_REAL%/bin/java}"
    fi
fi

if [[ -z "${JAVA_HOME:-}" ]] || [[ ! -x "${JAVA_HOME}/bin/java" ]]; then
    echo "JAVA_HOME is not set and java was not auto-detected. Hail/Spark requires Java on compute nodes." >&2
    exit 1
fi

"$PYTHON_BIN" "$SCRIPT_PATH" \
    --chr "$CHR" \
    --anc-hapcount "$ANC_HAPCOUNT" \
    --anc-i-dosage "$ANC_I_DOSAGE" \
    --anc-j-dosage "$ANC_J_DOSAGE" \
    --label-i "$LABEL_I" \
    --label-j "$LABEL_J" \
    --phenotype "$PHENOTYPE" \
    --phenotype-format "$PHENOTYPE_FORMAT" \
    --fam-case-code "$FAM_CASE_CODE" \
    --fam-control-code "$FAM_CONTROL_CODE" \
    --table-iid-column "$TABLE_IID_COLUMN" \
    --table-phenotype-column "$TABLE_PHENO_COLUMN" \
    --pca-enabled "$PCA_ENABLED" \
    --pca-file "$PCA_FILE" \
    --pca-format "$PCA_FORMAT" \
    --pca-has-header "$PCA_HAS_HEADER" \
    --pca-fid-column "$PCA_FID_COLUMN" \
    --pca-iid-column "$PCA_IID_COLUMN" \
    --pca-n-pcs "$PCA_N_PCS" \
    --min-partitions "$MIN_PARTITIONS" \
    --output-tsv "$OUTPUT_TSV"
