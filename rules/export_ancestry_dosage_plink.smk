import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "ancestry_plink_dosage")

TRACTOR_DONE = os.path.join(get_analysis_value("tractor.output_dir"), "extract_tracts.done")
TRACTOR_DOSAGE_PATTERN = get_analysis_value("admixture_validation.tractor_dosage_pattern")

CHROMOSOMES = [str(chr_name) for chr_name in get_step5_successful_chromosomes()]
ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
ANCESTRY_INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
ANCESTRY_BY_LABEL = {str(ANCESTRY_LABELS[idx]): str(idx) for idx in ANCESTRY_INDICES}
ANCESTRY_NAMES = sorted(ANCESTRY_BY_LABEL.keys())

PLINK_MODULE = get_software_module("plink2")
PYTHON_MODULE = get_software_module("python")


def dosage_path(chr_name, ancestry_label):
    anc_idx = ANCESTRY_BY_LABEL[ancestry_label]
    return TRACTOR_DOSAGE_PATTERN.replace("{chr}", str(chr_name)).replace("{anc}", str(anc_idx))


rule ancestry_dosage_to_plink_chr:
    input:
        tractor_done=TRACTOR_DONE,
        dosage=lambda wildcards: dosage_path(wildcards.chr, wildcards.ancestry),
    output:
        pgen=os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.pgen"),
        pvar=os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.pvar"),
        psam=os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.psam"),
    log:
        os.path.join(LOG_DIR, "ancestry_dosage_to_plink_chr{chr}_{ancestry}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        converter="scripts/convert_dosage_txt_to_vcf.py",
        plink_module=PLINK_MODULE,
        python_module=PYTHON_MODULE,
        out_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"chr{wildcards.chr}.dosage"),
        tmp_vcf=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"chr{wildcards.chr}.dosage.vcf.gz"),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/{wildcards.ancestry} {LOG_DIR}
        bash --login -c '
            ml "{params.python_module}" >/dev/null 2>&1 || true
            python3 "{params.converter}" \
                --dosage "{input.dosage}" \
                --output-vcf-gz "{params.tmp_vcf}" \
                --source "{wildcards.ancestry}_dosage"

            ml "{params.plink_module}" >/dev/null 2>&1 || true
            plink2 \
                --vcf "{params.tmp_vcf}" \
                --double-id \
                --allow-extra-chr \
                --make-pgen \
                --out "{params.out_prefix}"

            rm -f "{params.tmp_vcf}"
        ' > "{log}" 2>&1
        """


rule ancestry_dosage_plink_done:
    input:
        tractor_done=TRACTOR_DONE,
        pgen=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.pgen"), ancestry=ANCESTRY_NAMES, chr=CHROMOSOMES),
        pvar=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.pvar"), ancestry=ANCESTRY_NAMES, chr=CHROMOSOMES),
        psam=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "chr{chr}.dosage.psam"), ancestry=ANCESTRY_NAMES, chr=CHROMOSOMES),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_dosage_plink.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_dosage_plink_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Per-ancestry PLINK dosage export complete" > {log}
        """
