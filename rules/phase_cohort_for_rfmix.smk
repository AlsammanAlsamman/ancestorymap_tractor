import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "cohort_phasing")

STEP3_DONE = os.path.join(RESULTS_DIR, "cohort_region_subset", "prepare_cohort_regions_from_plink.done")
STEP3C_DONE = os.path.join(get_analysis_value("step3c.output_dir"), "harmonize_cohort_alleles.done")
STEP2_DONE = os.path.join(RESULTS_DIR, "reference_panel_subset", "subset_reference_panel.done")

LOCI_FILE = get_analysis_value("step3.loci_file")
INPUT_PATTERN = get_analysis_value("step4.cohort_input_pattern")
REFERENCE_PATTERN = get_analysis_value("step4.reference_vcf_pattern")
MAP_PATTERN = get_analysis_value("step4.genetic_map_pattern")
OUTPUT_DIR = get_analysis_value("step4.output_dir")


def _chr_sort_key(value):
    value = str(value).replace("chr", "")
    return (0, int(value)) if value.isdigit() else (1, value)


def load_chromosomes_from_loci(loci_file):
    chromosomes = set()
    with open(loci_file, "r", encoding="utf-8") as handle:
        next(handle, None)
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chromosomes.add(str(parts[1]).replace("chr", ""))
    if not chromosomes:
        raise ValueError(f"No chromosomes found in loci file: {loci_file}")
    return sorted(chromosomes, key=_chr_sort_key)


CHROMOSOMES = load_chromosomes_from_loci(LOCI_FILE)

EXPECTED_INPUT_PATTERN = os.path.join(
    get_analysis_value("step3c.output_dir"),
    "chr{chr}." + get_analysis_value("step3c.output_suffix") + ".vcf.gz",
)

if INPUT_PATTERN != EXPECTED_INPUT_PATTERN:
    raise ValueError(
        "step4.cohort_input_pattern is not compatible with step3c outputs. "
        f"Expected: {EXPECTED_INPUT_PATTERN}; Found: {INPUT_PATTERN}"
    )

BCFTOOLS_MODULE = get_software_module("bcftools")
SHAPEIT5_MODULE = get_software_module("shapeit5")


def resolve_pattern(pattern, chr_name):
    return pattern.replace("{chr}", str(chr_name))


def phase_mem_mb(chr_name):
    return get_default_resource("mem_mb", 256000 if str(chr_name) == "6" else 128000)


def phase_cores(chr_name):
    return get_default_resource("cores", 8 if str(chr_name) == "6" else 20)


rule phase_cohort_chr:
    input:
        step2_done=STEP2_DONE,
        step3_done=STEP3_DONE,
        step3c_done=STEP3C_DONE,
        cohort_vcf=lambda wildcards: resolve_pattern(INPUT_PATTERN, wildcards.chr),
        reference_vcf=lambda wildcards: resolve_pattern(REFERENCE_PATTERN, wildcards.chr),
        genetic_map=lambda wildcards: resolve_pattern(MAP_PATTERN, wildcards.chr),
    output:
        vcf=os.path.join(OUTPUT_DIR, "chr{chr}.phased.vcf.gz"),
        index=os.path.join(OUTPUT_DIR, "chr{chr}.phased.vcf.gz.tbi"),
    log:
        os.path.join(LOG_DIR, "phase_cohort_chr{chr}.log")
    resources:
        mem_mb=lambda wildcards: phase_mem_mb(wildcards.chr),
        time=get_default_resource("time", "72:00:00"),
        cores=lambda wildcards: phase_cores(wildcards.chr)
    params:
        helper="scripts/phase_cohort_chr.sh",
        bcftools_module=BCFTOOLS_MODULE,
        shapeit5_module=SHAPEIT5_MODULE,
        threads=lambda wildcards, resources: resources.cores,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --input-vcf {input.cohort_vcf} \
            --reference-vcf {input.reference_vcf} \
            --genetic-map {input.genetic_map} \
            --chr {wildcards.chr} \
            --output-vcf {output.vcf} \
            --threads {params.threads} \
            --bcftools-module {params.bcftools_module} \
            --shapeit5-module {params.shapeit5_module} \
            > {log} 2>&1
        """


rule phase_cohort_for_rfmix_done:
    input:
        step2_done=STEP2_DONE,
        step3_done=STEP3_DONE,
        step3c_done=STEP3C_DONE,
        vcfs=expand(os.path.join(OUTPUT_DIR, "chr{chr}.phased.vcf.gz"), chr=CHROMOSOMES),
        indexes=expand(os.path.join(OUTPUT_DIR, "chr{chr}.phased.vcf.gz.tbi"), chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "phase_cohort_for_rfmix.done")
    log:
        os.path.join(LOG_DIR, "phase_cohort_for_rfmix_done.log")
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "72:00:00"),
        cores=get_default_resource("cores", 2)
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "Cohort phasing complete" > {log}
        """
