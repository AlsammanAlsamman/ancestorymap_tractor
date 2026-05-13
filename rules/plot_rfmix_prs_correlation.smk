import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_loci_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "rfmix_prs_correlation")

STEP5_DONE = os.path.join(RESULTS_DIR, "rfmix", "run_rfmix.done")

CHROMOSOMES = get_loci_chromosomes()
Q_PATTERN = get_analysis_value("step8.rfmix_q_pattern")
PRS_FILE = get_analysis_value("step8.prs_file")
ANCESTRY_BINS = str(get_analysis_value("step8.ancestry_contribution_bins"))
PHENOTYPE_FILE = get_analysis_value("step8.phenotype_file")
PHENOTYPE_FORMAT = str(get_analysis_value("step8.phenotype_format"))
FAM_CASE_CODE = str(get_analysis_value("step8.fam_case_code"))
FAM_CONTROL_CODE = str(get_analysis_value("step8.fam_control_code"))

OUTPUT_DIR = get_analysis_value("step8.output_dir")
PER_CHR_PREFIX_PATTERN = get_analysis_value("step8.per_chr_prefix_pattern")
MATRIX_PREFIX = get_analysis_value("step8.matrix_prefix")

PRS_IID_COLUMN = str(get_analysis_value("step8.prs_iid_column"))
PRS_SCORE_COLUMN = str(get_analysis_value("step8.prs_score_column"))
PRS_IN_REG_COLUMN = str(get_analysis_value("step8.prs_in_regression_column"))
PRS_IN_REG_KEEP = str(get_analysis_value("step8.prs_in_regression_keep_value"))

ANCESTRY_LABELS = get_analysis_value("step8.ancestry_labels")
ANCESTRY_COLORS = get_analysis_value("step8.ancestry_colors")

LABEL0 = str(ANCESTRY_LABELS["0"])
LABEL1 = str(ANCESTRY_LABELS["1"])
LABEL2 = str(ANCESTRY_LABELS["2"])
COLOR0 = str(ANCESTRY_COLORS[LABEL0])
COLOR1 = str(ANCESTRY_COLORS[LABEL1])
COLOR2 = str(ANCESTRY_COLORS[LABEL2])

R_MODULE = get_software_module("r")


def resolve_pattern(pattern, chr_name):
    return pattern.replace("{chr}", str(chr_name))


PER_CHR_PNG_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".png")
PER_CHR_PDF_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".pdf")
PER_CHR_TSV_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".tsv")
PER_CHR_CC_PNG_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".case_control.png")
PER_CHR_CC_PDF_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".case_control.pdf")
PER_CHR_CC_TSV_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".case_control.tsv")

MATRIX_PNG = os.path.join(OUTPUT_DIR, MATRIX_PREFIX + ".png")
MATRIX_PDF = os.path.join(OUTPUT_DIR, MATRIX_PREFIX + ".pdf")
MATRIX_CSV = os.path.join(OUTPUT_DIR, MATRIX_PREFIX + ".csv")


rule rfmix_prs_corr_per_chr:
    input:
        step5_done=STEP5_DONE,
        q=lambda wildcards: resolve_pattern(Q_PATTERN, wildcards.chr),
        phenotype=PHENOTYPE_FILE,
        prs=PRS_FILE,
    output:
        png=PER_CHR_PNG_PATTERN,
        pdf=PER_CHR_PDF_PATTERN,
        tsv=PER_CHR_TSV_PATTERN,
        cc_png=PER_CHR_CC_PNG_PATTERN,
        cc_pdf=PER_CHR_CC_PDF_PATTERN,
        cc_tsv=PER_CHR_CC_TSV_PATTERN,
    log:
        os.path.join(LOG_DIR, "rfmix_prs_corr_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/plot_rfmix_prs_corr_chr.R",
        r_module=R_MODULE,
        prs_iid_column=PRS_IID_COLUMN,
        prs_score_column=PRS_SCORE_COLUMN,
        prs_in_reg_column=PRS_IN_REG_COLUMN,
        prs_in_reg_keep=PRS_IN_REG_KEEP,
        ancestry_bins=ANCESTRY_BINS,
        phenotype_format=PHENOTYPE_FORMAT,
        fam_case_code=FAM_CASE_CODE,
        fam_control_code=FAM_CONTROL_CODE,
        label0=LABEL0,
        label1=LABEL1,
        label2=LABEL2,
        color0=COLOR0,
        color1=COLOR1,
        color2=COLOR2,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.r_module}"
        fi
        Rscript "{params.script}" \
            --q "{input.q}" \
            --chr "{wildcards.chr}" \
            --phenotype "{input.phenotype}" \
            --phenotype-format "{params.phenotype_format}" \
            --fam-case-code "{params.fam_case_code}" \
            --fam-control-code "{params.fam_control_code}" \
            --prs "{input.prs}" \
            --prs-iid-column "{params.prs_iid_column}" \
            --prs-score-column "{params.prs_score_column}" \
            --prs-in-regression-column "{params.prs_in_reg_column}" \
            --prs-in-regression-keep-value "{params.prs_in_reg_keep}" \
            --ancestry-bins "{params.ancestry_bins}" \
            --png "{output.png}" \
            --pdf "{output.pdf}" \
            --tsv "{output.tsv}" \
            --cc-png "{output.cc_png}" \
            --cc-pdf "{output.cc_pdf}" \
            --cc-tsv "{output.cc_tsv}" \
            --label0 "{params.label0}" \
            --label1 "{params.label1}" \
            --label2 "{params.label2}" \
            --color0 "{params.color0}" \
            --color1 "{params.color1}" \
            --color2 "{params.color2}" \
            > "{log}" 2>&1
        """


rule rfmix_prs_corr_matrix:
    input:
        step5_done=STEP5_DONE,
        per_chr_tsv=expand(PER_CHR_TSV_PATTERN, chr=CHROMOSOMES),
    output:
        png=MATRIX_PNG,
        pdf=MATRIX_PDF,
        csv=MATRIX_CSV,
    log:
        os.path.join(LOG_DIR, "rfmix_prs_corr_matrix.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/plot_rfmix_prs_corr_matrix.R",
        r_module=R_MODULE,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.r_module}"
        fi
        Rscript "{params.script}" \
            --input-tsv {input.per_chr_tsv} \
            --png "{output.png}" \
            --pdf "{output.pdf}" \
            --csv "{output.csv}" \
            > "{log}" 2>&1
        """


rule plot_rfmix_prs_correlation_done:
    input:
        step5_done=STEP5_DONE,
        chr_png=expand(PER_CHR_PNG_PATTERN, chr=CHROMOSOMES),
        chr_pdf=expand(PER_CHR_PDF_PATTERN, chr=CHROMOSOMES),
        chr_tsv=expand(PER_CHR_TSV_PATTERN, chr=CHROMOSOMES),
        chr_cc_png=expand(PER_CHR_CC_PNG_PATTERN, chr=CHROMOSOMES),
        chr_cc_pdf=expand(PER_CHR_CC_PDF_PATTERN, chr=CHROMOSOMES),
        chr_cc_tsv=expand(PER_CHR_CC_TSV_PATTERN, chr=CHROMOSOMES),
        matrix_png=MATRIX_PNG,
        matrix_pdf=MATRIX_PDF,
        matrix_csv=MATRIX_CSV,
    output:
        done=os.path.join(DONE_DIR, "plot_rfmix_prs_correlation.done"),
    log:
        os.path.join(LOG_DIR, "plot_rfmix_prs_correlation_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "RFMix-PRS correlation step complete" > {log}
        """
