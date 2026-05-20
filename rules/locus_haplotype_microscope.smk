import os
import sys

import pandas as pd

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir, get_software_module


RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "locus_haplotype_microscope")

REGIONS_FILE = get_analysis_value("validation_regions")
_regions_df = pd.read_csv(REGIONS_FILE, sep="\t", dtype=str)
if "locus" not in _regions_df.columns:
    raise ValueError(f"validation_regions file is missing required 'locus' column: {REGIONS_FILE}")
LOCUS_IDS = [str(item) for item in _regions_df["locus"].dropna().tolist()]

GWAS_DONE = os.path.join(RESULTS_DIR, "ancestry_plink_gwas", "ancestry_plink_gwas_with_freq.done")
TRACTOR_DONE = os.path.join(get_analysis_value("tractor.output_dir"), "extract_tracts.done")
ANCESTRY_LABELS = {str(code): str(label) for code, label in get_analysis_value("tractor_gwas.ancestry_labels").items()}
ANCESTRY_NAMES = sorted(ANCESTRY_LABELS.values())
PYTHON_MODULE = get_software_module("python")


rule locus_microscope_per_locus:
    input:
        gwas_done=GWAS_DONE,
        tractor_done=TRACTOR_DONE,
        regions_file=REGIONS_FILE,
    output:
        tract_plot=os.path.join(OUTPUT_DIR, "{locus}", "{locus}_tract_plot.png"),
        ancestry_composition=os.path.join(OUTPUT_DIR, "{locus}", "{locus}_ancestry_composition.png"),
        allele_freq_plot=os.path.join(OUTPUT_DIR, "{locus}", "{locus}_allele_freq_by_ancestry.png"),
        allele_counts=os.path.join(OUTPUT_DIR, "{locus}", "{locus}_allele_counts.tsv"),
        summary_card=os.path.join(OUTPUT_DIR, "{locus}", "{locus}_summary_card.xlsx"),
    log:
        os.path.join(LOG_DIR, "locus_microscope_{locus}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "01:00:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/locus_haplotype_microscope.py",
        config_file="analysis.yml",
        python_module=PYTHON_MODULE,
        output_dir=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.locus),
        n_samples=100,
    shell:
        """
        mkdir -p "{params.output_dir}" "{LOG_DIR}"
        bash --login -c '
            if type module >/dev/null 2>&1; then
                module load "{params.python_module}" || true
            fi

            python3 "{params.script}" \
                --config "{params.config_file}" \
                --regions "{input.regions_file}" \
                --locus "{wildcards.locus}" \
                --outdir "{params.output_dir}" \
                --n-samples "{params.n_samples}"
        ' > "{log}" 2>&1
        """


rule locus_microscope_merge:
    input:
        counts=expand(os.path.join(OUTPUT_DIR, "{locus}", "{locus}_allele_counts.tsv"), locus=LOCUS_IDS),
    output:
        merged_tsv=os.path.join(OUTPUT_DIR, "all_loci_allele_counts.tsv"),
        merged_excel=os.path.join(OUTPUT_DIR, "all_loci_summary_cards.xlsx"),
    log:
        os.path.join(LOG_DIR, "locus_microscope_merge.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/locus_haplotype_microscope.py",
        config_file="analysis.yml",
        python_module=PYTHON_MODULE,
        output_dir=OUTPUT_DIR,
    shell:
        """
        mkdir -p "{OUTPUT_DIR}" "{LOG_DIR}"
        bash --login -c '
            if type module >/dev/null 2>&1; then
                module load "{params.python_module}" || true
            fi

            python3 "{params.script}" \
                --config "{params.config_file}" \
                --regions "{REGIONS_FILE}" \
                --merge-only \
                --counts-dir "{params.output_dir}" \
                --outdir "{params.output_dir}"
        ' > "{log}" 2>&1
        """


rule locus_microscope_done:
    input:
        per_locus_tract=expand(os.path.join(OUTPUT_DIR, "{locus}", "{locus}_tract_plot.png"), locus=LOCUS_IDS),
        per_locus_comp=expand(os.path.join(OUTPUT_DIR, "{locus}", "{locus}_ancestry_composition.png"), locus=LOCUS_IDS),
        per_locus_freq=expand(os.path.join(OUTPUT_DIR, "{locus}", "{locus}_allele_freq_by_ancestry.png"), locus=LOCUS_IDS),
        merged_tsv=os.path.join(OUTPUT_DIR, "all_loci_allele_counts.tsv"),
        merged_excel=os.path.join(OUTPUT_DIR, "all_loci_summary_cards.xlsx"),
    output:
        done=os.path.join(OUTPUT_DIR, "locus_haplotype_microscope.done"),
    log:
        os.path.join(LOG_DIR, "locus_haplotype_microscope_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        locus_count=len(LOCUS_IDS),
    shell:
        """
        mkdir -p "{OUTPUT_DIR}" "{LOG_DIR}"
        touch "{output.done}"
        echo "Locus haplotype microscope complete for {params.locus_count} loci across {len(ANCESTRY_NAMES)} ancestries" > "{log}"
        """