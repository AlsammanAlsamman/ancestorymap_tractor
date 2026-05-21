import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_loci_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "rfmix_plots")

STEP5_DONE = os.path.join(RESULTS_DIR, "rfmix", "run_rfmix.done")

CHROMOSOMES = get_loci_chromosomes()
MSP_PATTERN = get_analysis_value("step6.rfmix_msp_pattern")
PLOT_OUTPUT_DIR = get_analysis_value("step6.output_dir")
PLOT_PREFIX_PATTERN = get_analysis_value("step6.output_prefix_pattern")

ANCESTRY_LABELS = get_analysis_value("step6.ancestry_labels")
ANCESTRY_COLORS = get_analysis_value("step6.ancestry_colors")

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


rule plot_rfmix_chr:
    input:
        step5_done=STEP5_DONE,
        msp=lambda wildcards: resolve_pattern(MSP_PATTERN, wildcards.chr),
    output:
        png=PNG_PATTERN,
        pdf=PDF_PATTERN,
    log:
        os.path.join(LOG_DIR, "plot_rfmix_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/plot_rfmix_chr.R",
        r_module=R_MODULE,
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
            --png "{output.png}" \
            --pdf "{output.pdf}" \
            --label0 "{params.label0}" \
            --label1 "{params.label1}" \
            --label2 "{params.label2}" \
            --color0 "{params.color0}" \
            --color1 "{params.color1}" \
            --color2 "{params.color2}" \
            > "{log}" 2>&1
        """


rule plot_rfmix_tracts_done:
    input:
        step5_done=STEP5_DONE,
        pngs=expand(PNG_PATTERN, chr=CHROMOSOMES),
        pdfs=expand(PDF_PATTERN, chr=CHROMOSOMES),
    output:
        done=os.path.join(DONE_DIR, "plot_rfmix_tracts.done"),
    log:
        os.path.join(LOG_DIR, "plot_rfmix_tracts_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {DONE_DIR}
        touch {output.done}
        echo "RFMix plot step complete" > {log}
        """
