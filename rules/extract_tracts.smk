import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "tractor")

STEP4_DONE = os.path.join(RESULTS_DIR, "cohort_phasing", "phase_cohort_for_rfmix.done")
STEP5_DONE = os.path.join(RESULTS_DIR, "rfmix", "run_rfmix.done")

CHROMOSOMES = get_step5_successful_chromosomes()
VCF_PATTERN = get_analysis_value("tractor.vcf_pattern")
MSP_PATTERN = get_analysis_value("tractor.msp_pattern")
NUM_ANCS = int(get_analysis_value("tractor.num_ancs"))
OUTPUT_DIR = get_analysis_value("tractor.output_dir")
OUTPUT_VCF = bool(get_analysis_value("tractor.output_vcf"))
COMPRESS_OUTPUT = bool(get_analysis_value("tractor.compress_output"))
LABEL_CODES = get_analysis_value("step5.sample_map_label_codes")

if not isinstance(LABEL_CODES, dict) or not LABEL_CODES:
    raise ValueError("step5.sample_map_label_codes must be a non-empty mapping for extract_tracts")

LABEL_CODES = {str(k): str(v) for k, v in LABEL_CODES.items()}
EXPECTED_CODES = sorted(int(code) for code in LABEL_CODES.values())
EXPECTED_CODE_SEQUENCE = list(range(NUM_ANCS))
if EXPECTED_CODES != EXPECTED_CODE_SEQUENCE:
    raise ValueError(
        "Tractor extract_tracts requires ancestry codes to be contiguous 0..N-1. "
        f"Found codes: {EXPECTED_CODES}; Expected: {EXPECTED_CODE_SEQUENCE}"
    )

ANCESTRY_INDICES = [str(code) for code in EXPECTED_CODES]
LABEL_CODE_MAP_STRING = ",".join(
    f"{label}:{LABEL_CODES[label]}" for label in sorted(LABEL_CODES, key=lambda label: int(LABEL_CODES[label]))
)

EXPECTED_VCF_PATTERN = os.path.join(get_analysis_value("step4.output_dir"), "chr{chr}.phased.vcf.gz")
EXPECTED_MSP_PATTERN = os.path.join(get_analysis_value("step5.output_dir"), "chr{chr}.deconvoluted.msp.tsv")

if VCF_PATTERN != EXPECTED_VCF_PATTERN:
    raise ValueError(
        "tractor.vcf_pattern is not compatible with step4 outputs. "
        f"Expected: {EXPECTED_VCF_PATTERN}; Found: {VCF_PATTERN}"
    )

if MSP_PATTERN != EXPECTED_MSP_PATTERN:
    raise ValueError(
        "tractor.msp_pattern is not compatible with step5 outputs. "
        f"Expected: {EXPECTED_MSP_PATTERN}; Found: {MSP_PATTERN}"
    )

PYTHON_MODULE = get_software_module("python")


def resolve_pattern(pattern, chr_name):
    return pattern.replace("{chr}", str(chr_name))


rule extract_tracts_chr:
    input:
        step4_done=STEP4_DONE,
        step5_done=STEP5_DONE,
        vcf=lambda wildcards: resolve_pattern(VCF_PATTERN, wildcards.chr),
        msp=lambda wildcards: resolve_pattern(MSP_PATTERN, wildcards.chr),
    output:
        dosage=expand(os.path.join(OUTPUT_DIR, "chr{{chr}}.phased.anc{anc}.dosage.txt"), anc=ANCESTRY_INDICES),
        hapcount=expand(os.path.join(OUTPUT_DIR, "chr{{chr}}.phased.anc{anc}.hapcount.txt"), anc=ANCESTRY_INDICES),
    log:
        os.path.join(LOG_DIR, "extract_tracts_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        helper="scripts/extract_tracts.sh",
        script="scripts/extract_tracts.py",
        num_ancs=NUM_ANCS,
        python_module=PYTHON_MODULE,
        output_vcf=OUTPUT_VCF,
        compress_output=COMPRESS_OUTPUT,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --vcf {input.vcf} \
            --msp {input.msp} \
            --num-ancs {params.num_ancs} \
            --output-dir {OUTPUT_DIR} \
            --script-path {params.script} \
            --python-module {params.python_module} \
            --output-vcf {params.output_vcf} \
            --compress-output {params.compress_output} \
            > {log} 2>&1
        """


rule extract_tracts_done:
    input:
        step4_done=STEP4_DONE,
        step5_done=STEP5_DONE,
        dosage=expand(os.path.join(OUTPUT_DIR, "chr{chr}.phased.anc{anc}.dosage.txt"), chr=CHROMOSOMES, anc=ANCESTRY_INDICES),
        hapcount=expand(os.path.join(OUTPUT_DIR, "chr{chr}.phased.anc{anc}.hapcount.txt"), chr=CHROMOSOMES, anc=ANCESTRY_INDICES),
    output:
        done=os.path.join(DONE_DIR, "extract_tracts.done"),
    log:
        os.path.join(LOG_DIR, "extract_tracts_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        checker="scripts/check_extract_tracts_integrity.py",
        output_dir=OUTPUT_DIR,
        num_ancs=NUM_ANCS,
        label_code_map=LABEL_CODE_MAP_STRING,
    shell:
        """
        mkdir -p {DONE_DIR}
        python3 {params.checker} \
            --output-dir {params.output_dir} \
            --num-ancs {params.num_ancs} \
            --label-code-map '{params.label_code_map}' \
            --chromosomes {CHROMOSOMES} \
            > {log} 2>&1
        touch {output.done}
        echo "Tractor extract_tracts step complete (ancestry integrity verified)" >> {log}
        """
