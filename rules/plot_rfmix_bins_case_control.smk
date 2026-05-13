import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_loci_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "rfmix_bin_case_control")

STEP5_DONE = os.path.join(RESULTS_DIR, "rfmix", "run_rfmix.done")

CHROMOSOMES = get_loci_chromosomes()
MSP_PATTERN = get_analysis_value("step7.rfmix_msp_pattern")
PLOT_OUTPUT_DIR = get_analysis_value("step7.output_dir")
PLOT_PREFIX_PATTERN = get_analysis_value("step7.output_prefix_pattern")

BIN_SIZE_BP = str(get_analysis_value("step7.bin_size_bp"))
GENE_ANNOTATION_FILE = get_analysis_value("step7.gene_annotation_file")
MAX_GENES_PER_BIN_LABEL = str(get_analysis_value("step7.max_genes_per_bin_label"))

PHENOTYPE_FILE = get_analysis_value("step7.phenotype_file")
PHENOTYPE_FORMAT = str(get_analysis_value("step7.phenotype_format"))
FAM_CASE_CODE = str(get_analysis_value("step7.fam_case_code"))
FAM_CONTROL_CODE = str(get_analysis_value("step7.fam_control_code"))
TABLE_SAMPLE_COLUMN = str(get_analysis_value("step7.table_sample_column"))
TABLE_STATUS_COLUMN = str(get_analysis_value("step7.table_status_column"))
TABLE_CASE_VALUE = str(get_analysis_value("step7.table_case_value"))
TABLE_CONTROL_VALUE = str(get_analysis_value("step7.table_control_value"))

ANCESTRY_LABELS = get_analysis_value("step7.ancestry_labels")
ANCESTRY_COLORS = get_analysis_value("step7.ancestry_colors")

LABEL0 = str(ANCESTRY_LABELS["0"])
LABEL1 = str(ANCESTRY_LABELS["1"])
LABEL2 = str(ANCESTRY_LABELS["2"])
COLOR0 = str(ANCESTRY_COLORS[LABEL0])
COLOR1 = str(ANCESTRY_COLORS[LABEL1])
COLOR2 = str(ANCESTRY_COLORS[LABEL2])

R_MODULE = get_software_module("r")


def resolve_pattern(pattern, chr_name):
    return pattern.replace("{chr}", str(chr_name))


PNG_PATTERN = os.path.join(PLOT_OUTPUT_DIR, PLOT_PREFIX_PATTERN + ".png")
PDF_PATTERN = os.path.join(PLOT_OUTPUT_DIR, PLOT_PREFIX_PATTERN + ".pdf")
TSV_PATTERN = os.path.join(PLOT_OUTPUT_DIR, PLOT_PREFIX_PATTERN + ".tsv")


rule plot_rfmix_bins_case_control_chr:
    input:
        step5_done=STEP5_DONE,
        msp=lambda wildcards: resolve_pattern(MSP_PATTERN, wildcards.chr),
        phenotype=PHENOTYPE_FILE,
        genes=GENE_ANNOTATION_FILE,
    output:
        png=PNG_PATTERN,
        pdf=PDF_PATTERN,
        tsv=TSV_PATTERN,
    log:
        os.path.join(LOG_DIR, "plot_rfmix_bins_case_control_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/plot_rfmix_bins_case_control_chr.R",
        r_module=R_MODULE,
        bin_size_bp=BIN_SIZE_BP,
        max_genes_per_bin_label=MAX_GENES_PER_BIN_LABEL,
        phenotype_format=PHENOTYPE_FORMAT,
        fam_case_code=FAM_CASE_CODE,
        fam_control_code=FAM_CONTROL_CODE,
        table_sample_column=TABLE_SAMPLE_COLUMN,
        table_status_column=TABLE_STATUS_COLUMN,
        table_case_value=TABLE_CASE_VALUE,
        table_control_value=TABLE_CONTROL_VALUE,
        label0=LABEL0,
        label1=LABEL1,
        label2=LABEL2,
        color0=COLOR0,
        color1=COLOR1,
        color2=COLOR2,
    shell:
        """
        mkdir -p {PLOT_OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.r_module}"
        fi
        Rscript "{params.script}" \
            --msp "{input.msp}" \
            --chr "{wildcards.chr}" \
            --phenotype "{input.phenotype}" \
            --phenotype-format "{params.phenotype_format}" \
            --fam-case-code "{params.fam_case_code}" \
            --fam-control-code "{params.fam_control_code}" \
            --table-sample-column "{params.table_sample_column}" \
            --table-status-column "{params.table_status_column}" \
            --table-case-value "{params.table_case_value}" \
            --table-control-value "{params.table_control_value}" \
            --genes "{input.genes}" \
            --bin-size-bp "{params.bin_size_bp}" \
            --max-genes-per-bin-label "{params.max_genes_per_bin_label}" \
            --png "{output.png}" \
            --pdf "{output.pdf}" \
            --tsv "{output.tsv}" \
            --label0 "{params.label0}" \
            --label1 "{params.label1}" \
            --label2 "{params.label2}" \
            --color0 "{params.color0}" \
            --color1 "{params.color1}" \
            --color2 "{params.color2}" \
            > "{log}" 2>&1
        """


rule plot_rfmix_bins_case_control_done:
    input:
        step5_done=STEP5_DONE,
        pngs=expand(PNG_PATTERN, chr=CHROMOSOMES),
        pdfs=expand(PDF_PATTERN, chr=CHROMOSOMES),
        tsvs=expand(TSV_PATTERN, chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "plot_rfmix_bins_case_control.done"),
    log:
        os.path.join(LOG_DIR, "plot_rfmix_bins_case_control_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "RFMix bin-level case/control plot step complete" > {log}
        """
