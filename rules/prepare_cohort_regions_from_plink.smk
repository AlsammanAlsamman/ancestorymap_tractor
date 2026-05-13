import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "cohort_region_subset")

STEP2_DONE = os.path.join(RESULTS_DIR, "reference_panel_subset", "subset_reference_panel.done")

PLINK_PREFIX = get_analysis_value("step3.plink_prefix")
LOCI_FILE = get_analysis_value("step3.loci_file")
OUTPUT_DIR = get_analysis_value("step3.output_dir")
OUTPUT_SUFFIX = get_analysis_value("step3.output_suffix")

PLINK_MODULE = get_software_module("plink2")
BCFTOOLS_MODULE = get_software_module("bcftools")


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


rule prepare_cohort_regions_one_chr:
    input:
        step2_done=STEP2_DONE,
        bed=PLINK_PREFIX + ".bed",
        bim=PLINK_PREFIX + ".bim",
        fam=PLINK_PREFIX + ".fam",
        loci=LOCI_FILE,
    output:
        vcf=os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz"),
        index=os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz.tbi"),
    log:
        os.path.join(LOG_DIR, "prepare_cohort_regions_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        helper="scripts/prepare_cohort_regions_from_plink.sh",
        plink_prefix=PLINK_PREFIX,
        plink_module=PLINK_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --plink-prefix {params.plink_prefix} \
            --loci-file {input.loci} \
            --chr {wildcards.chr} \
            --output-vcf {output.vcf} \
            --plink-module {params.plink_module} \
            --bcftools-module {params.bcftools_module} \
            > {log} 2>&1
        """


rule prepare_cohort_regions_from_plink_done:
    input:
        step2_done=STEP2_DONE,
        vcfs=expand(os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz"), chr=CHROMOSOMES),
        indexes=expand(os.path.join(OUTPUT_DIR, "chr{chr}." + OUTPUT_SUFFIX + ".vcf.gz.tbi"), chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "prepare_cohort_regions_from_plink.done"),
    log:
        os.path.join(LOG_DIR, "prepare_cohort_regions_from_plink_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "Cohort loci-based region extraction from PLINK complete" > {log}
        """
