import os
import sys

sys.path.append("utils")
from bioconfigme import (
    get_results_dir,
    get_software_module,
    get_analysis_value,
    get_default_resource,
    get_step5_successful_chromosomes,
)

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
PLINK_DOSAGE_DIR = os.path.join(RESULTS_DIR, "ancestry_plink_dosage")
DIAG_DIR = os.path.join(RESULTS_DIR, "snp_count_diagnostic")

PLINK_DOSAGE_DONE = os.path.join(PLINK_DOSAGE_DIR, "ancestry_dosage_plink.done")

CHROMOSOMES = [str(c) for c in get_step5_successful_chromosomes()]

ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
ANCESTRY_INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
ANCESTRY_NAMES = sorted([str(ANCESTRY_LABELS[idx]) for idx in ANCESTRY_INDICES])

COHORT_PREFIX = get_analysis_value("cohort.plink_prefix")
ORIGINAL_BIM = os.path.join("inputs", "{}.bim".format(COHORT_PREFIX))
ORIGINAL_FAM = os.path.join("inputs", "{}.fam".format(COHORT_PREFIX))

R_MODULE = get_software_module("r")
PYTHON_MODULE = get_software_module("python")


rule snp_count_diagnostic_done:
    input:
        plink_done=PLINK_DOSAGE_DONE,
        pvar=expand(
            os.path.join(PLINK_DOSAGE_DIR, "{ancestry}", "chr{chr}.dosage.pvar"),
            ancestry=ANCESTRY_NAMES,
            chr=CHROMOSOMES,
        ),
        psam=expand(
            os.path.join(PLINK_DOSAGE_DIR, "{ancestry}", "chr{chr}.dosage.psam"),
            ancestry=ANCESTRY_NAMES,
            chr=CHROMOSOMES,
        ),
        orig_bim=ORIGINAL_BIM,
        orig_fam=ORIGINAL_FAM,
    output:
        table=os.path.join(DIAG_DIR, "snp_count_table.tsv"),
        stacked_pdf=os.path.join(DIAG_DIR, "snp_count_stacked_plot.pdf"),
        facet_pdf=os.path.join(DIAG_DIR, "snp_count_facet_plot.pdf"),
        retention_pdf=os.path.join(DIAG_DIR, "snp_retention_plot.pdf"),
        sample_pdf=os.path.join(DIAG_DIR, "sample_count_plot.pdf"),
        done=os.path.join(DIAG_DIR, "snp_count_diagnostic.done"),
    log:
        os.path.join(LOG_DIR, "snp_count_diagnostic.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 128000),
        time=get_default_resource("time", "02:00:00"),
        cores=get_default_resource("cores", 4),
    params:
        count_script="scripts/count_snps_diagnostic.py",
        plot_script="scripts/plot_snp_count_diagnostic.R",
        ancestry_plink_dir=PLINK_DOSAGE_DIR,
        ancestries=",".join(ANCESTRY_NAMES),
        chromosomes=",".join(CHROMOSOMES),
        diag_dir=DIAG_DIR,
        python_module=PYTHON_MODULE,
        r_module=R_MODULE,
    shell:
        """
        mkdir -p {params.diag_dir} {LOG_DIR}
        bash --login -c '
            # --- Python: count SNPs/samples ---
            if type module >/dev/null 2>&1; then
                module load {params.python_module} >/dev/null 2>&1 || true
            fi
            if type ml >/dev/null 2>&1; then
                ml {params.python_module} >/dev/null 2>&1 || true
            fi
            python3 {params.count_script} \
                --ancestry-plink-dir {params.ancestry_plink_dir} \
                --ancestries {params.ancestries} \
                --chromosomes {params.chromosomes} \
                --original-bim {input.orig_bim} \
                --original-fam {input.orig_fam} \
                --output-table {output.table}

            # --- R: produce plots ---
            if type module >/dev/null 2>&1; then
                module load {params.r_module} >/dev/null 2>&1 || true
            fi
            if type ml >/dev/null 2>&1; then
                ml {params.r_module} >/dev/null 2>&1 || true
            fi
            command -v Rscript >/dev/null 2>&1 || {{ echo "Rscript not found after module load"; exit 127; }}
            Rscript {params.plot_script} \
                --input {output.table} \
                --output-dir {params.diag_dir}
        ' > {log} 2>&1
        touch {output.done}
        echo "SNP count diagnostic complete" >> {log}
        """
