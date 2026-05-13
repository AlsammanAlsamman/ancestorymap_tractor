#!/bin/bash
set -euo pipefail

COHORT_PHASED_VCF=""
REFERENCE_VCF=""
SAMPLE_MAP=""
GENETIC_MAP=""
CHR=""
OUTPUT_PREFIX=""
RFMIX_MODULE=""
BCFTOOLS_MODULE=""
THREADS="20"
MIN_COMMON_SNPS="2"
LABEL_CODE_MAP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort-phased-vcf)
            COHORT_PHASED_VCF="$2"
            shift 2
            ;;
        --reference-vcf)
            REFERENCE_VCF="$2"
            shift 2
            ;;
        --sample-map)
            SAMPLE_MAP="$2"
            shift 2
            ;;
        --genetic-map)
            GENETIC_MAP="$2"
            shift 2
            ;;
        --chr)
            CHR="$2"
            shift 2
            ;;
        --output-prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --rfmix-module)
            RFMIX_MODULE="$2"
            shift 2
            ;;
        --bcftools-module)
            BCFTOOLS_MODULE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --min-common-snps)
            MIN_COMMON_SNPS="$2"
            shift 2
            ;;
        --label-code-map)
            LABEL_CODE_MAP="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "$COHORT_PHASED_VCF" || -z "$REFERENCE_VCF" || -z "$SAMPLE_MAP" || -z "$GENETIC_MAP" || -z "$CHR" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Missing required arguments" >&2
    exit 2
fi

if [[ -z "$RFMIX_MODULE" ]]; then
    echo "Missing required --rfmix-module argument from workflow config" >&2
    exit 2
fi

if [[ -z "$LABEL_CODE_MAP" ]]; then
    echo "Missing required --label-code-map argument from analysis config" >&2
    exit 2
fi

if ! [[ "$MIN_COMMON_SNPS" =~ ^[0-9]+$ ]]; then
    echo "--min-common-snps must be a non-negative integer: $MIN_COMMON_SNPS" >&2
    exit 2
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
    echo "--threads must be a positive integer: $THREADS" >&2
    exit 2
fi

for f in "$COHORT_PHASED_VCF" "$REFERENCE_VCF" "$SAMPLE_MAP" "$GENETIC_MAP"; do
    if [[ ! -f "$f" ]]; then
        echo "Required file not found: $f" >&2
        exit 1
    fi
done

if [[ ! -f "${COHORT_PHASED_VCF}.tbi" && ! -f "${COHORT_PHASED_VCF}.csi" ]]; then
    echo "Missing index for cohort phased VCF: ${COHORT_PHASED_VCF}" >&2
    exit 1
fi

if [[ ! -f "${REFERENCE_VCF}.tbi" && ! -f "${REFERENCE_VCF}.csi" ]]; then
    echo "Missing index for reference VCF: ${REFERENCE_VCF}" >&2
    exit 1
fi

# Validate sample map format (sample and label fields exist)
if ! awk 'NF>=2 {good=1} END{exit !good?1:0}' "$SAMPLE_MAP"; then
    echo "Sample map has invalid format (expected: sample_id label): $SAMPLE_MAP" >&2
    exit 1
fi

# Ensure module command is available in non-login batch shells.
if ! type module >/dev/null 2>&1; then
    if [[ -f /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    elif [[ -f /usr/share/Modules/init/bash ]]; then
        source /usr/share/Modules/init/bash
    fi
fi

if ! type module >/dev/null 2>&1; then
    echo "Environment module system is unavailable; cannot load required module: $RFMIX_MODULE" >&2
    exit 1
fi

if ! module load "$RFMIX_MODULE"; then
    echo "Failed to load required RFMix module from software.yml: $RFMIX_MODULE" >&2
    exit 1
fi

have_working_bcftools=false
if command -v bcftools >/dev/null 2>&1 && bcftools --version >/dev/null 2>&1; then
    have_working_bcftools=true
fi

if [[ "$have_working_bcftools" == false ]]; then
    if [[ -n "$BCFTOOLS_MODULE" ]]; then
        module load "$BCFTOOLS_MODULE" || true
    fi
fi

if ! command -v bcftools >/dev/null 2>&1 || ! bcftools --version >/dev/null 2>&1; then
    echo "bcftools is required to pre-check cohort/reference overlap before running RFMix" >&2
    echo "Configured bcftools module: ${BCFTOOLS_MODULE:-<unset>}" >&2
    exit 1
fi

module list 2>&1 || true

RFMIX_BIN=""
RFMIX_BIN_MODE="native"
for candidate in rfmix rfmix2 RunRFMix RFMix RunRFMix.py; do
    if command -v "$candidate" >/dev/null 2>&1; then
        RFMIX_BIN="$candidate"
        if [[ "$candidate" == "RunRFMix.py" ]]; then
            RFMIX_BIN_MODE="python"
        fi
        break
    fi
done

if [[ -z "$RFMIX_BIN" ]]; then
    for root in "${RFMIX_HOME:-}" "${EBROOTRFMIX:-}"; do
        if [[ -n "$root" ]]; then
            for path_candidate in \
                "$root/rfmix" \
                "$root/rfmix2" \
                "$root/RunRFMix" \
                "$root/RFMix" \
                "$root/RunRFMix.py" \
                "$root/bin/rfmix" \
                "$root/bin/rfmix2" \
                "$root/bin/RunRFMix" \
                "$root/bin/RFMix" \
                "$root/bin/RunRFMix.py"; do
                if [[ -x "$path_candidate" ]]; then
                    RFMIX_BIN="$path_candidate"
                    case "$path_candidate" in
                        *.py)
                            RFMIX_BIN_MODE="python"
                            ;;
                        *)
                            RFMIX_BIN_MODE="native"
                            ;;
                    esac
                    break
                fi
            done
        fi
        if [[ -n "$RFMIX_BIN" ]]; then
            break
        fi
    done
fi

if [[ -z "$RFMIX_BIN" ]]; then
    echo "PATH after module load: $PATH" >&2
    echo "RFMIX_HOME=${RFMIX_HOME:-<unset>}" >&2
    echo "EBROOTRFMIX=${EBROOTRFMIX:-<unset>}" >&2
    command -v rfmix >/dev/null 2>&1 && command -v rfmix >&2 || true
    command -v rfmix2 >/dev/null 2>&1 && command -v rfmix2 >&2 || true
    command -v RunRFMix >/dev/null 2>&1 && command -v RunRFMix >&2 || true
    command -v RFMix >/dev/null 2>&1 && command -v RFMix >&2 || true
    command -v RunRFMix.py >/dev/null 2>&1 && command -v RunRFMix.py >&2 || true
    echo "RFMix executable not found in PATH after loading module: $RFMIX_MODULE (tried: rfmix, rfmix2, RunRFMix, RFMix, RunRFMix.py)" >&2
    exit 1
fi

OUT_DIR="$(dirname "$OUTPUT_PREFIX")"
mkdir -p "$OUT_DIR"

WORK_DIR="$(mktemp -d -p "$OUT_DIR" ".tmp_rfmix_chr${CHR}.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

NORMALIZED_MAP="$WORK_DIR/chr${CHR}.rfmix.map"
NORMALIZED_SAMPLE_MAP="$WORK_DIR/chr${CHR}.rfmix.sample_map.txt"

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required to normalize the genetic map format for RFMix" >&2
    exit 1
fi

python3 - "$GENETIC_MAP" "$NORMALIZED_MAP" "$CHR" <<'PY'
import gzip
import sys

input_map, output_map, chr_target = sys.argv[1], sys.argv[2], sys.argv[3]


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def clean_chr(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        return value[3:]
    return value


open_fn = gzip.open if input_map.endswith('.gz') else open
written = 0

with open_fn(input_map, 'rt', encoding='utf-8') as src, open(output_map, 'w', encoding='utf-8') as dst:
    for line in src:
        row = line.strip()
        if not row:
            continue
        parts = row.split()
        if len(parts) < 3:
            continue

        p1, p2, p3 = parts[0], parts[1], parts[2]

        # Skip known header patterns.
        if p1.lower() in {"pos", "position", "bp"}:
            continue

        if is_number(p1) and not is_number(p2):
            # pos chr cM
            pos, chr_val, cm = p1, clean_chr(p2), p3
        elif not is_number(p1) and is_number(p2):
            # chr pos cM
            chr_val, pos, cm = clean_chr(p1), p2, p3
        elif is_number(p1) and is_number(p2):
            # likely pos chr cM with numeric chr column
            pos, chr_val, cm = p1, clean_chr(p2), p3
        else:
            continue

        if chr_val == chr_target:
            dst.write(f"{chr_val}\t{pos}\t{cm}\n")
            written += 1

if written == 0:
    print(
        f"No usable genetic map positions found for chromosome {chr_target} in {input_map}",
        file=sys.stderr,
    )
    sys.exit(3)
PY

python3 - "$SAMPLE_MAP" "$NORMALIZED_SAMPLE_MAP" "$LABEL_CODE_MAP" <<'PY'
import sys

sample_map_in, sample_map_out, label_code_map_str = sys.argv[1], sys.argv[2], sys.argv[3]

label_to_code = {}
for pair in label_code_map_str.split(','):
    item = pair.strip()
    if not item:
        continue
    if ':' not in item:
        print(f"Invalid --label-code-map entry (expected LABEL:CODE): {item}", file=sys.stderr)
        sys.exit(4)
    label, code = item.split(':', 1)
    label = label.strip()
    code = code.strip()
    if not label or not code:
        print(f"Invalid --label-code-map entry (empty label or code): {item}", file=sys.stderr)
        sys.exit(4)
    if label in label_to_code:
        print(f"Duplicate label in --label-code-map: {label}", file=sys.stderr)
        sys.exit(4)
    if code in label_to_code.values():
        print(f"Duplicate numeric code in --label-code-map: {code}", file=sys.stderr)
        sys.exit(4)
    label_to_code[label] = code

if not label_to_code:
    print("--label-code-map produced no ancestry mappings", file=sys.stderr)
    sys.exit(4)

# First pass: extract all unique labels from sample map
labels_found = set()
with open(sample_map_in, 'r', encoding='utf-8') as src:
    for line in src:
        row = line.strip()
        if not row:
            continue
        parts = row.split()
        if len(parts) >= 2:
            labels_found.add(parts[1])

expected_labels = set(label_to_code.keys())
if labels_found != expected_labels:
    print(
        f"Sample map labels {sorted(labels_found)} do not match expected labels {sorted(expected_labels)}",
        file=sys.stderr,
    )
    sys.exit(4)

# Second pass: write normalized sample map
written = 0
with open(sample_map_in, 'r', encoding='utf-8') as src, open(sample_map_out, 'w', encoding='utf-8') as dst:
    for line in src:
        row = line.strip()
        if not row:
            continue
        parts = row.split()
        if len(parts) < 2:
            continue
        sample, label = parts[0], parts[1]
        if label not in label_to_code:
            print(f"Unsupported ancestry label in sample map: {label}", file=sys.stderr)
            sys.exit(4)
        dst.write(f"{sample}\t{label_to_code[label]}\n")
        written += 1

if written == 0:
    print("Sample map normalization produced no rows", file=sys.stderr)
    sys.exit(5)
PY

count_common_snps() {
    local cohort_sites="$WORK_DIR/chr${CHR}.cohort.sites.tsv"
    local reference_sites="$WORK_DIR/chr${CHR}.reference.sites.tsv"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$COHORT_PHASED_VCF" > "$cohort_sites"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$REFERENCE_VCF" > "$reference_sites"

    python3 - "$cohort_sites" "$reference_sites" <<'PY'
import sys

cohort_path, reference_path = sys.argv[1], sys.argv[2]

with open(cohort_path, 'r', encoding='utf-8') as handle:
    cohort_sites = {line.rstrip('\n') for line in handle if line.strip()}

count = 0
seen = set()
with open(reference_path, 'r', encoding='utf-8') as handle:
    for raw_line in handle:
        line = raw_line.rstrip('\n')
        if not line or line in seen:
            continue
        seen.add(line)
        if line in cohort_sites:
            count += 1

print(count)
PY
}

COMMON_SNPS=$(count_common_snps)
echo "Common SNPs between cohort and reference for chr${CHR}: ${COMMON_SNPS}"

if (( COMMON_SNPS < MIN_COMMON_SNPS )); then
    reason="#SKIPPED chr${CHR} insufficient_common_snps=${COMMON_SNPS} threshold=${MIN_COMMON_SNPS}"
    for suffix in msp.tsv fb.tsv rfmix.Q sis.tsv; do
        printf '%s\n' "$reason" > "${OUTPUT_PREFIX}.${suffix}"
    done
    echo "$reason"
    exit 0
fi

if [[ "$RFMIX_BIN_MODE" == "python" ]]; then
    if ! command -v python3 >/dev/null 2>&1; then
        echo "python3 not found, required for RunRFMix.py launcher" >&2
        exit 1
    fi
    RFMIX_CMD=(python3 "$RFMIX_BIN")
else
    RFMIX_CMD=("$RFMIX_BIN")
fi

SUPPORTS_N_THREADS=0
if [[ "$RFMIX_BIN_MODE" == "python" ]]; then
    if python3 "$RFMIX_BIN" --help 2>&1 | grep -q -- "--n-threads"; then
        SUPPORTS_N_THREADS=1
    fi
else
    if "$RFMIX_BIN" --help 2>&1 | grep -q -- "--n-threads"; then
        SUPPORTS_N_THREADS=1
    fi
fi

RFMIX_ARGS=(
    -f "$COHORT_PHASED_VCF"
    -r "$REFERENCE_VCF"
    -m "$NORMALIZED_SAMPLE_MAP"
    -g "$NORMALIZED_MAP"
    -o "$OUTPUT_PREFIX"
    --chromosome="$CHR"
)

if (( SUPPORTS_N_THREADS == 1 )); then
    RFMIX_ARGS+=(--n-threads="$THREADS")
else
    echo "Warning: RFMix binary does not advertise --n-threads; running without explicit thread flag" >&2
fi

"${RFMIX_CMD[@]}" "${RFMIX_ARGS[@]}"

echo "Created RFMix outputs for chr${CHR} with prefix: $OUTPUT_PREFIX (threads=${THREADS})"
