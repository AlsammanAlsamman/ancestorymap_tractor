#!/bin/bash
set -euo pipefail

PYTHON_MODULE="python/3.10.2"
USER_BASE="${HOME}/.local"
HAIL_SPEC="hail==0.2.137"
SKIP_INIT_CHECK=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python-module)
            PYTHON_MODULE="$2"
            shift 2
            ;;
        --user-base)
            USER_BASE="$2"
            shift 2
            ;;
        --hail-spec)
            HAIL_SPEC="$2"
            shift 2
            ;;
        --skip-init-check)
            SKIP_INIT_CHECK=true
            shift
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

# Make the environment modules command available in non-login shells.
if ! type module >/dev/null 2>&1; then
    if [[ -f /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    elif [[ -f /usr/share/Modules/init/bash ]]; then
        source /usr/share/Modules/init/bash
    fi
fi

if type module >/dev/null 2>&1 && [[ -n "$PYTHON_MODULE" ]]; then
    module load "$PYTHON_MODULE"
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

mkdir -p "$USER_BASE"
export PYTHONUSERBASE="$USER_BASE"
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

if ! "$PYTHON_BIN" -m pip --version >/dev/null 2>&1; then
    "$PYTHON_BIN" -m ensurepip --user >/dev/null 2>&1 || true
fi

echo "Using Python: $($PYTHON_BIN --version 2>&1)"
echo "Installing into user base: $PYTHONUSERBASE"
echo "Installing package spec: $HAIL_SPEC"

"$PYTHON_BIN" -m pip install --user --upgrade pip setuptools wheel
"$PYTHON_BIN" -m pip install --user --upgrade-strategy only-if-needed "$HAIL_SPEC"

"$PYTHON_BIN" - <<'PY'
import hail
print(f"Hail import OK: version {hail.__version__}")
PY

if [[ "$SKIP_INIT_CHECK" == false ]]; then
    if ! command -v java >/dev/null 2>&1; then
        if type module >/dev/null 2>&1; then
            module load java >/dev/null 2>&1 || module load jdk >/dev/null 2>&1 || module load openjdk >/dev/null 2>&1 || true
        fi
    fi

    if command -v java >/dev/null 2>&1; then
        "$PYTHON_BIN" - <<'PY'
import hail as hl
hl.init(quiet=True, log='/tmp/hail_init_test.log')
print('Hail init OK')
hl.stop()
PY
    else
        echo "Warning: java not found; skipped 'hl.init()' validation." >&2
    fi
fi

echo "Local Hail installation is ready."
