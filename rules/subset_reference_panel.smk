import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_loci_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "reference_panel_subset")

REF_RELEASE_DIR = get_analysis_value("step2.reference_release_dir")
SAMPLE_LIST = get_analysis_value("step2.sample_list")
OUTPUT_DIR = get_analysis_value("step2.output_dir")
OUTPUT_SUFFIX = get_analysis_value("step2.output_suffix")
LOCI_FILE = get_analysis_value("step3.loci_file")
CHROMOSOMES = get_loci_chromosomes(LOCI_FILE)
STEP1_DONE = os.path.join(RESULTS_DIR, "reference_panel_prep", "build_reference_sample_lists.done")

BCFTOOLS_MODULE = get_software_module("bcftools")


def ref_input_vcf(chr_name):
    return os.path.join(
        REF_RELEASE_DIR,
        f"ALL.chr{chr_name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    )


rule subset_reference_one_chr:
    input:
        step1_done=STEP1_DONE,
        sample_list=SAMPLE_LIST,
        ref_vcf=lambda wildcards: ref_input_vcf(wildcards.chr)
    output:
        vcf=os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz"),
        index=os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz.tbi")
    log:
        os.path.join(LOG_DIR, "subset_reference_chr{chr}.log")
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2)
    params:
        helper="scripts/subset_reference_panel.sh",
        bcftools_module=BCFTOOLS_MODULE
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --input-vcf {input.ref_vcf} \
            --sample-list {input.sample_list} \
            --output-vcf {output.vcf} \
            --bcftools-module {params.bcftools_module} \
            > {log} 2>&1
        """


rule subset_reference_panel_done:
    input:
        step1_done=STEP1_DONE,
        vcfs=expand(
            os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz"),
            chr=CHROMOSOMES,
        ),
        indexes=expand(
            os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz.tbi"),
            chr=CHROMOSOMES,
        )
    output:
        done=os.path.join(DONE_DIR, "subset_reference_panel.done")
    log:
        os.path.join(LOG_DIR, "subset_reference_panel_done.log")
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2)
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "Subset reference panel complete" > {log}
        """
