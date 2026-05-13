import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_loci_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "rfmix")

STEP4_DONE = os.path.join(RESULTS_DIR, "cohort_phasing", "phase_cohort_for_rfmix.done")

CHROMOSOMES = get_loci_chromosomes()
COHORT_PHASED_PATTERN = get_analysis_value("step5.cohort_phased_pattern")
REFERENCE_PATTERN = get_analysis_value("step5.reference_vcf_pattern")
GENETIC_MAP_PATTERN = get_analysis_value("step5.genetic_map_pattern")
SAMPLE_MAP_FILE = get_analysis_value("step5.sample_map_file")
OUTPUT_DIR = get_analysis_value("step5.output_dir")
LABEL_CODES = get_analysis_value("step5.sample_map_label_codes")

if not isinstance(LABEL_CODES, dict) or not LABEL_CODES:
    raise ValueError("step5.sample_map_label_codes must be a non-empty mapping")

LABEL_CODES = {str(k): str(v) for k, v in LABEL_CODES.items()}
if len(set(LABEL_CODES.values())) != len(LABEL_CODES):
    raise ValueError("step5.sample_map_label_codes contains duplicate numeric codes")

_sorted_labels = sorted(LABEL_CODES.keys(), key=lambda lbl: int(LABEL_CODES[lbl]))
LABEL_CODE_MAP_STRING = ",".join([f"{lbl}:{LABEL_CODES[lbl]}" for lbl in _sorted_labels])
EXPECTED_CODES = sorted({int(v) for v in LABEL_CODES.values()})

EXPECTED_COHORT_PATTERN = os.path.join(get_analysis_value("step4.output_dir"), "chr{chr}.phased.vcf.gz")
if COHORT_PHASED_PATTERN != EXPECTED_COHORT_PATTERN:
    raise ValueError(
        "step5.cohort_phased_pattern is not compatible with step4 outputs. "
        f"Expected: {EXPECTED_COHORT_PATTERN}; Found: {COHORT_PHASED_PATTERN}"
    )

RFMIX_MODULE = get_software_module("rfmix")
BCFTOOLS_MODULE = get_software_module("bcftools")
MIN_COMMON_SNPS = int(get_analysis_value("step5.min_common_snps"))


def resolve_pattern(pattern, chr_name):
    return pattern.replace("{chr}", str(chr_name))


rule run_rfmix_chr:
    input:
        step4_done=STEP4_DONE,
        cohort_phased=lambda wildcards: resolve_pattern(COHORT_PHASED_PATTERN, wildcards.chr),
        reference_vcf=lambda wildcards: resolve_pattern(REFERENCE_PATTERN, wildcards.chr),
        genetic_map=lambda wildcards: resolve_pattern(GENETIC_MAP_PATTERN, wildcards.chr),
        sample_map=SAMPLE_MAP_FILE,
    output:
        msp=os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.msp.tsv"),
        fb=os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.fb.tsv"),
        q=os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.rfmix.Q"),
        sis=os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.sis.tsv"),
    log:
        os.path.join(LOG_DIR, "run_rfmix_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 128000),
        time=get_default_resource("time", "72:00:00"),
        cores=get_default_resource("cores", 20),
    params:
        checker="scripts/check_rfmix_ancestry_integrity.py",
        sample_map=SAMPLE_MAP_FILE,
        label_code_map=LABEL_CODE_MAP_STRING,
        helper="scripts/run_rfmix_chr.sh",
        rfmix_module=RFMIX_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
        min_common_snps=MIN_COMMON_SNPS,
        threads=lambda wildcards, resources: resources.cores,
        output_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, f"chr{wildcards.chr}.deconvoluted"),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --cohort-phased-vcf {input.cohort_phased} \
            --reference-vcf {input.reference_vcf} \
            --sample-map {input.sample_map} \
            --genetic-map {input.genetic_map} \
            --chr {wildcards.chr} \
            --output-prefix {params.output_prefix} \
            --rfmix-module {params.rfmix_module} \
            --bcftools-module {params.bcftools_module} \
            --threads {params.threads} \
            --min-common-snps {params.min_common_snps} \
            --label-code-map '{params.label_code_map}' \
            > {log} 2>&1
        """


rule run_rfmix_done:
    input:
        step4_done=STEP4_DONE,
        msp=expand(os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.msp.tsv"), chr=CHROMOSOMES),
        fb=expand(os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.fb.tsv"), chr=CHROMOSOMES),
        q=expand(os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.rfmix.Q"), chr=CHROMOSOMES),
        sis=expand(os.path.join(OUTPUT_DIR, "chr{chr}.deconvoluted.sis.tsv"), chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "run_rfmix.done"),
    log:
        os.path.join(LOG_DIR, "run_rfmix_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "72:00:00"),
        cores=get_default_resource("cores", 2),
    params:
        checker="scripts/check_rfmix_ancestry_integrity.py",
        sample_map=SAMPLE_MAP_FILE,
        label_code_map=LABEL_CODE_MAP_STRING,
    shell:
        """
        mkdir -p {DONE_DIR}
        python3 - <<'PY'
import os

output_dir = {OUTPUT_DIR!r}
chromosomes = {CHROMOSOMES!r}
successful = []
skipped = []

for chrom in chromosomes:
    msp_path = os.path.join(output_dir, f"chr{{chrom}}.deconvoluted.msp.tsv")
    if not os.path.exists(msp_path):
        raise FileNotFoundError(f"Missing expected RFMix MSP output: {{msp_path}}")
    with open(msp_path, "r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    if first_line.startswith("#SKIPPED"):
        skipped.append((chrom, first_line))
    else:
        successful.append(chrom)

success_path = os.path.join(output_dir, "successful_chromosomes.txt")
skip_path = os.path.join(output_dir, "skipped_chromosomes.txt")
with open(success_path, "w", encoding="utf-8") as handle:
    for chrom in successful:
        print(f"{{chrom}}", file=handle)
with open(skip_path, "w", encoding="utf-8") as handle:
    for chrom, reason in skipped:
        print(f"{{chrom}}\t{{reason}}", file=handle)

if not successful:
    raise RuntimeError("All chromosomes were skipped in Step 5; no usable RFMix outputs were produced")
PY
        SUCCESS_FILES=$(python3 - <<'PY'
import os

output_dir = {OUTPUT_DIR!r}
success_path = os.path.join(output_dir, "successful_chromosomes.txt")
paths = []
with open(success_path, "r", encoding="utf-8") as handle:
    for raw_line in handle:
        chrom = raw_line.strip()
        if chrom:
            paths.append(os.path.join(output_dir, f"chr{{chrom}}.deconvoluted.msp.tsv"))
print(" ".join(paths))
PY
)
        python3 {params.checker} \
            --sample-map {params.sample_map} \
            --label-code-map '{params.label_code_map}' \
            --msp-files $SUCCESS_FILES
        touch {output.done}
        echo "RFMix step complete (ancestry integrity verified on successful chromosomes)" > {log}
        """
