import os
import sys
from itertools import combinations

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("maf_gwas_summary.output_dir")
DONE_DIR = OUTPUT_DIR

MAF_DONE = os.path.join(RESULTS_DIR, "maf_summary", "calculate_maf.done")
GWAS_DONE = os.path.join(RESULTS_DIR, "tractor_gwas", "tractor_gwas_pairwise.done")
MAF_TABLE = os.path.join(RESULTS_DIR, "maf_summary", "merged.annotated.maf.tsv")
OUTPUT_TABLE = os.path.join(OUTPUT_DIR, "merged.maf_gwas_summary.tsv")

# Build dynamic pair names from ancestry labels (same logic as tractor_gwas_pairwise.smk)
_ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
_INDICES = sorted(_ANCESTRY_LABELS.keys(), key=lambda x: int(x))
_FOUND_CODES = [int(index) for index in _INDICES]
_EXPECTED_CODES = list(range(len(_INDICES)))
if _FOUND_CODES != _EXPECTED_CODES:
    raise ValueError(
        "merge_maf_gwas_summary requires contiguous ancestry codes 0..N-1. "
        f"Found: {_FOUND_CODES}; Expected: {_EXPECTED_CODES}"
    )

_ORDERED_ANCESTRIES = [str(_ANCESTRY_LABELS[index]) for index in _INDICES]
if len(set(_ORDERED_ANCESTRIES)) != len(_ORDERED_ANCESTRIES):
    raise ValueError("merge_maf_gwas_summary requires unique ancestry labels")

GWAS_PAIR_NAMES = []
for _i, _j in combinations(_INDICES, 2):
    _la = str(_ANCESTRY_LABELS[_i])
    _lb = str(_ANCESTRY_LABELS[_j])
    GWAS_PAIR_NAMES.append(f"{_la}_{_lb}")

GWAS_FILES = {
    pair: os.path.join(RESULTS_DIR, "tractor_gwas", f"merged.model_{pair}.gwas.tsv")
    for pair in GWAS_PAIR_NAMES
}

PYTHON_MODULE = get_software_module("python")


rule merge_maf_gwas_summary:
    input:
        maf_done=MAF_DONE,
        gwas_done=GWAS_DONE,
        maf=MAF_TABLE,
        gwas=list(GWAS_FILES.values()),
    output:
        tsv=OUTPUT_TABLE,
    log:
        os.path.join(LOG_DIR, "merge_maf_gwas_summary.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/merge_maf_gwas_summary.py",
        python_module=PYTHON_MODULE,
        pair_names=" ".join(GWAS_PAIR_NAMES),
        gwas_dir=os.path.join(RESULTS_DIR, "tractor_gwas"),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" || true
        fi
        GWAS_ARGS=""
        for pair in {params.pair_names}; do
            GWAS_ARGS="$GWAS_ARGS --gwas $pair={params.gwas_dir}/merged.model_${{pair}}.gwas.tsv"
        done
        python3 "{params.script}" \
            --maf "{input.maf}" \
            $GWAS_ARGS \
            --output "{output.tsv}" \
            > "{log}" 2>&1
        """


rule merge_maf_gwas_summary_done:
    input:
        tsv=OUTPUT_TABLE,
    output:
        done=os.path.join(DONE_DIR, "merge_maf_gwas_summary.done"),
    log:
        os.path.join(LOG_DIR, "merge_maf_gwas_summary_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:10:00"),
        cores=get_default_resource("cores", 1),
    params:
        checker="scripts/check_merge_maf_gwas_summary_integrity.py",
        ancestry_labels=" ".join(_ORDERED_ANCESTRIES),
        pair_names=" ".join(GWAS_PAIR_NAMES),
    shell:
        """
        mkdir -p {DONE_DIR}
        python3 "{params.checker}" \
            --summary "{input.tsv}" \
            --ancestries {params.ancestry_labels} \
            --pairs {params.pair_names} \
            > "{log}" 2>&1
        touch "{output.done}"
        echo "MAF + GWAS summary complete (ancestry integrity verified)" >> "{log}"
        """
