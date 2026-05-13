import os
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_loci_chromosomes, get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("maf_summary.output_dir")
PER_CHR_DIR = os.path.join(OUTPUT_DIR, "per_chr")
DONE_DIR = OUTPUT_DIR

STEP2_DONE = os.path.join(RESULTS_DIR, "reference_panel_subset", "subset_reference_panel.done")
STEP3C_DONE = os.path.join(get_analysis_value("step3c.output_dir"), "harmonize_cohort_alleles.done")

STEP2_OUTPUT_DIR = get_analysis_value("step2.output_dir")
STEP2_OUTPUT_SUFFIX = get_analysis_value("step2.output_suffix")
STEP3C_OUTPUT_DIR = get_analysis_value("step3c.output_dir")
STEP3C_OUTPUT_SUFFIX = get_analysis_value("step3c.output_suffix")
REF_SAMPLE_DIR = get_analysis_value("first_step.output_dir")
ANNOTATION_FILE = get_analysis_value("annotation_meta.annotation")
LOCI_FILE = get_analysis_value("step3.loci_file")
CHROMOSOMES = get_loci_chromosomes()

REF_VCF_PATTERN = os.path.join(STEP2_OUTPUT_DIR, f"chr{{chr}}.{STEP2_OUTPUT_SUFFIX}.vcf.gz")
COHORT_VCF_PATTERN = os.path.join(STEP3C_OUTPUT_DIR, f"chr{{chr}}.{STEP3C_OUTPUT_SUFFIX}.vcf.gz")

POPULATIONS = [str(item) for item in get_analysis_value("first_step.populations")]
POP_TO_LABEL = {str(pop): str(label) for pop, label in get_analysis_value("first_step.ancestry_labels").items()}
LABEL_TO_POP = {label: pop for pop, label in POP_TO_LABEL.items()}
SOURCE_LABELS = ["cohort"] + [POP_TO_LABEL[pop] for pop in POPULATIONS]
_SOURCE_LABEL_PATTERN = "|".join(SOURCE_LABELS)

REFERENCE_LABELS = [POP_TO_LABEL[pop] for pop in POPULATIONS]

PER_CHR_COHORT_PATTERN = os.path.join(PER_CHR_DIR, "cohort.chr{chr}.maf.tsv")
PER_CHR_REF_PATTERN = os.path.join(PER_CHR_DIR, "{label}.chr{chr}.maf.tsv")
SOURCE_TABLE_PATTERN = os.path.join(OUTPUT_DIR, "{source}.maf.tsv")
MERGED_TABLE = os.path.join(OUTPUT_DIR, "merged.annotated.maf.tsv")

PYTHON_MODULE = get_software_module("python")
BCFTOOLS_MODULE = get_software_module("bcftools")


def source_chr_inputs(source):
    if source == "cohort":
        return expand(PER_CHR_COHORT_PATTERN, chr=CHROMOSOMES)
    return expand(PER_CHR_REF_PATTERN, chr=CHROMOSOMES, label=[source])


rule maf_cohort_chr:
    input:
        step3c_done=STEP3C_DONE,
        vcf=COHORT_VCF_PATTERN,
    output:
        tsv=PER_CHR_COHORT_PATTERN,
    log:
        os.path.join(LOG_DIR, "maf_cohort_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/calculate_maf_from_vcf.py",
        python_module=PYTHON_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
        label="cohort",
        loci_file=LOCI_FILE,
        format_version="v2_harmonized_alt_af",
    shell:
        """
        mkdir -p {PER_CHR_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" "{params.bcftools_module}" || true
        fi
        python3 "{params.script}" \
            --vcf "{input.vcf}" \
            --regions-file "{params.loci_file}" \
            --output "{output.tsv}" \
            --label "{params.label}" \
            > "{log}" 2>&1
        echo "format={params.format_version}" >> "{log}"
        """


rule maf_reference_chr:
    input:
        step2_done=STEP2_DONE,
        vcf=REF_VCF_PATTERN,
        sample_file=lambda wildcards: os.path.join(REF_SAMPLE_DIR, f"{LABEL_TO_POP[wildcards.label]}.samples.txt"),
    output:
        tsv=PER_CHR_REF_PATTERN,
    log:
        os.path.join(LOG_DIR, "maf_{label}_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/calculate_maf_from_vcf.py",
        python_module=PYTHON_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
        label=lambda wildcards: wildcards.label,
        loci_file=LOCI_FILE,
        format_version="v2_harmonized_alt_af",
    shell:
        """
        mkdir -p {PER_CHR_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" "{params.bcftools_module}" || true
        fi
        python3 "{params.script}" \
            --vcf "{input.vcf}" \
            --regions-file "{params.loci_file}" \
            --sample-file "{input.sample_file}" \
            --output "{output.tsv}" \
            --label "{params.label}" \
            > "{log}" 2>&1
        echo "format={params.format_version}" >> "{log}"
        """


rule maf_source_table:
    wildcard_constraints:
        source=_SOURCE_LABEL_PATTERN
    input:
        lambda wildcards: source_chr_inputs(wildcards.source),
    output:
        tsv=SOURCE_TABLE_PATTERN,
    log:
        os.path.join(LOG_DIR, "maf_merge_{source}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:15:00"),
        cores=get_default_resource("cores", 1),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        set -- {input}
        first_file="$1"
        head -n 1 "$first_file" > "{output.tsv}"
        for f in "$@"; do
            tail -n +2 "$f"
        done | sort -k1,1V -k2,2n >> "{output.tsv}"
        echo "Merged chromosome tables into {output.tsv}" > "{log}"
        """


rule maf_merged_annotated:
    input:
        step2_done=STEP2_DONE,
        step3c_done=STEP3C_DONE,
        per_source=expand(SOURCE_TABLE_PATTERN, source=SOURCE_LABELS),
    output:
        merged=MERGED_TABLE,
    log:
        os.path.join(LOG_DIR, "maf_merged_annotated.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 64000),
        time=get_default_resource("time", "01:00:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/merge_annotate_maf.py",
        python_module=PYTHON_MODULE,
        annotation=ANNOTATION_FILE,
        loci_file=LOCI_FILE,
        source_args=" ".join(
            [f'--source cohort="{os.path.join(OUTPUT_DIR, "cohort.maf.tsv")}"']
            + [
                f'--source {label}="{os.path.join(OUTPUT_DIR, f"{label}.maf.tsv")}"'
                for label in REFERENCE_LABELS
            ]
        ),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" || true
        fi
        python3 "{params.script}" \
            {params.source_args} \
            --annotation "{params.annotation}" \
            --loci-file "{params.loci_file}" \
            --output "{output.merged}" \
            > "{log}" 2>&1
        """


rule calculate_maf_done:
    input:
        per_source=expand(SOURCE_TABLE_PATTERN, source=SOURCE_LABELS),
        merged=MERGED_TABLE,
    output:
        done=os.path.join(DONE_DIR, "calculate_maf.done"),
    log:
        os.path.join(LOG_DIR, "calculate_maf_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:10:00"),
        cores=get_default_resource("cores", 1),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch "{output.done}"
        echo "MAF summary complete" > "{log}"
        """
