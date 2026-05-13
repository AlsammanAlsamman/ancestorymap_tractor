import os
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_loci_chromosomes, get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("step3c.output_dir")
OUTPUT_SUFFIX = get_analysis_value("step3c.output_suffix")
DONE_DIR = OUTPUT_DIR

STEP2_DONE = os.path.join(RESULTS_DIR, "reference_panel_subset", "subset_reference_panel.done")
STEP3_DONE = os.path.join(RESULTS_DIR, "cohort_region_subset", "prepare_cohort_regions_from_plink.done")

LOCI_FILE = get_analysis_value("step3.loci_file")
CHROMOSOMES = get_loci_chromosomes(LOCI_FILE)

COHORT_INPUT_PATTERN = os.path.join(get_analysis_value("step3.output_dir"), f"chr{{chr}}.{get_analysis_value('step3.output_suffix')}.vcf.gz")
REFERENCE_PATTERN = os.path.join(get_analysis_value("step2.output_dir"), f"chr{{chr}}.{get_analysis_value('step2.output_suffix')}.vcf.gz")

HARMONIZED_PATTERN = os.path.join(OUTPUT_DIR, f"chr{{chr}}.{OUTPUT_SUFFIX}.vcf.gz")
REPORT_PATTERN = os.path.join(OUTPUT_DIR, "reports", "chr{chr}.harmonize_report.txt")

PYTHON_MODULE = get_software_module("python")
BCFTOOLS_MODULE = get_software_module("bcftools")


rule harmonize_cohort_alleles_chr:
    input:
        step2_done=STEP2_DONE,
        step3_done=STEP3_DONE,
        cohort_vcf=COHORT_INPUT_PATTERN,
        reference_vcf=REFERENCE_PATTERN,
    output:
        vcf=HARMONIZED_PATTERN,
        index=HARMONIZED_PATTERN + ".tbi",
        report=REPORT_PATTERN,
    log:
        os.path.join(LOG_DIR, "harmonize_cohort_alleles_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/harmonize_vcf_to_reference.py",
        python_module=PYTHON_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {OUTPUT_DIR}/reports {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" "{params.bcftools_module}" || true
        fi

        tmp_vcf="{OUTPUT_DIR}/.tmp.chr{wildcards.chr}.harmonized.vcf"
        python3 "{params.script}" \
            --cohort-vcf "{input.cohort_vcf}" \
            --reference-vcf "{input.reference_vcf}" \
            --output-vcf "$tmp_vcf" \
            --report "{output.report}" \
            > "{log}" 2>&1

        bcftools view -Oz -o "{output.vcf}" "$tmp_vcf"
        bcftools index -t "{output.vcf}"
        rm -f "$tmp_vcf"
        """


rule harmonize_cohort_alleles_done:
    input:
        vcfs=expand(HARMONIZED_PATTERN, chr=CHROMOSOMES),
        indexes=expand(HARMONIZED_PATTERN + ".tbi", chr=CHROMOSOMES),
        reports=expand(REPORT_PATTERN, chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "harmonize_cohort_alleles.done"),
    log:
        os.path.join(LOG_DIR, "harmonize_cohort_alleles_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:10:00"),
        cores=get_default_resource("cores", 1),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch "{output.done}"
        echo "Cohort/reference allele harmonization complete" > "{log}"
        """
