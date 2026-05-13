import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
STEP5_DIR = get_analysis_value("step5.output_dir")
STEP5_DONE = os.path.join(STEP5_DIR, "run_rfmix.done")
CHROMOSOMES = get_step5_successful_chromosomes()
LABEL_CODES = get_analysis_value("step5.sample_map_label_codes")

if not isinstance(LABEL_CODES, dict) or not LABEL_CODES:
    raise ValueError("step5.sample_map_label_codes must be a non-empty mapping for merge_rfmix_outputs")

LABEL_CODES = {str(label): str(code) for label, code in LABEL_CODES.items()}
if len(set(LABEL_CODES.values())) != len(LABEL_CODES):
    raise ValueError("step5.sample_map_label_codes contains duplicate numeric codes")

SORTED_LABELS = sorted(LABEL_CODES, key=lambda label: int(LABEL_CODES[label]))
LABEL_CODE_MAP_STRING = ",".join(f"{label}:{LABEL_CODES[label]}" for label in SORTED_LABELS)

MERGED_Q = os.path.join(STEP5_DIR, "merged.deconvoluted.rfmix.Q.tsv")
MERGED_MSP = os.path.join(STEP5_DIR, "merged.deconvoluted.msp.tsv")
MERGE_DONE = os.path.join(STEP5_DIR, "merge_rfmix_outputs.done")


rule merge_rfmix_outputs_done:
    input:
        step5_done=STEP5_DONE,
        q_files=expand(os.path.join(STEP5_DIR, "chr{chr}.deconvoluted.rfmix.Q"), chr=CHROMOSOMES),
        msp_files=expand(os.path.join(STEP5_DIR, "chr{chr}.deconvoluted.msp.tsv"), chr=CHROMOSOMES),
    output:
        merged_q=MERGED_Q,
        merged_msp=MERGED_MSP,
        done=MERGE_DONE,
    log:
        os.path.join(LOG_DIR, "merge_rfmix_outputs.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "04:00:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/merge_rfmix_outputs.py",
        label_code_map=LABEL_CODE_MAP_STRING,
    shell:
        """
        mkdir -p {STEP5_DIR} {LOG_DIR}
        python3 {params.script} \
            --q-files {input.q_files} \
            --msp-files {input.msp_files} \
            --label-code-map '{params.label_code_map}' \
            --output-q {output.merged_q} \
            --output-msp {output.merged_msp} \
            > {log} 2>&1
        touch {output.done}
        echo "Merged RFMix Q and MSP outputs complete" >> {log}
        """