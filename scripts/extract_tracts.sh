#!/bin/bash
set -euo pipefail

VCF=""
MSP=""
NUM_ANCS=""
OUTPUT_DIR=""
SCRIPT_PATH=""
PYTHON_MODULE=""
OUTPUT_VCF="False"
COMPRESS_OUTPUT="False"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)
            VCF="$2"
            shift 2
            ;;
        --msp)
            MSP="$2"
            shift 2
            ;;
        --num-ancs)
            NUM_ANCS="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
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
        --output-vcf)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        --compress-output)
            COMPRESS_OUTPUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "$VCF" || -z "$MSP" || -z "$NUM_ANCS" || -z "$OUTPUT_DIR" || -z "$SCRIPT_PATH" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

for f in "$VCF" "$MSP" "$SCRIPT_PATH"; do
    if [[ ! -f "$f" ]]; then
        echo "Required file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

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

CMD=("$PYTHON_BIN" "$SCRIPT_PATH" --vcf "$VCF" --msp "$MSP" --num-ancs "$NUM_ANCS" --output-dir "$OUTPUT_DIR")

if [[ "$OUTPUT_VCF" == "True" || "$OUTPUT_VCF" == "true" ]]; then
    CMD+=(--output-vcf)
fi

if [[ "$COMPRESS_OUTPUT" == "True" || "$COMPRESS_OUTPUT" == "true" ]]; then
    CMD+=(--compress-output)
fi

"${CMD[@]}"
