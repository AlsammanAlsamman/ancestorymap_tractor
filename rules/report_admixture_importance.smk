import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("admixture_importance.output_dir")

VALIDATION_DONE = os.path.join(get_analysis_value("admixture_validation.output_dir"), "validate_admixture_prediction.done")
SAMPLE_FEATURES_TSV = get_analysis_value("admixture_importance.sample_features_tsv")
SELECTED_SNPS_TSV = get_analysis_value("admixture_importance.selected_snps_tsv")
TRACTOR_DOSAGE_PATTERN = get_analysis_value("admixture_importance.tractor_dosage_pattern")
TOP_N_LOCI = int(get_analysis_value("admixture_importance.top_n_loci_plot"))
GLOBAL_REF_ANCESTRY = str(get_analysis_value("admixture_importance.global_ancestry_reference"))

CHROMOSOMES = [str(chr_name) for chr_name in get_step5_successful_chromosomes()]
ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
ANCESTRY_INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
ANCESTRY_MAP = ",".join(f"{idx}:{ANCESTRY_LABELS[idx]}" for idx in ANCESTRY_INDICES)

R_MODULE = get_software_module("r")


rule report_admixture_importance:
    input:
        validation_done=VALIDATION_DONE,
        sample_features=SAMPLE_FEATURES_TSV,
        selected_snps=SELECTED_SNPS_TSV,
        dosage_files=expand(TRACTOR_DOSAGE_PATTERN, chr=CHROMOSOMES, anc=ANCESTRY_INDICES),
    output:
        locus_scores=os.path.join(OUTPUT_DIR, "sample_locus_dosage_scores.tsv"),
        locus_importance=os.path.join(OUTPUT_DIR, "locus_ancestry_importance.tsv"),
        locus_overall=os.path.join(OUTPUT_DIR, "locus_overall_importance.tsv"),
        ancestry_importance=os.path.join(OUTPUT_DIR, "ancestry_dosage_importance.tsv"),
        metrics=os.path.join(OUTPUT_DIR, "importance_model_metrics.tsv"),
        top_loci_png=os.path.join(OUTPUT_DIR, "top_loci_importance.png"),
        ancestry_png=os.path.join(OUTPUT_DIR, "overall_ancestry_importance.png"),
        heatmap_png=os.path.join(OUTPUT_DIR, "locus_ancestry_importance_heatmap.png"),
    log:
        os.path.join(LOG_DIR, "report_admixture_importance.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/report_admixture_importance.R",
        r_module=R_MODULE,
        chromosomes=",".join(CHROMOSOMES),
        ancestry_map=ANCESTRY_MAP,
        global_ref_ancestry=GLOBAL_REF_ANCESTRY,
        top_n_loci=TOP_N_LOCI,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            ml "{params.r_module}" >/dev/null 2>&1 || true
            Rscript "{params.script}" \
                --sample-features-tsv "{input.sample_features}" \
                --selected-snps-tsv "{input.selected_snps}" \
                --tractor-dosage-pattern "{TRACTOR_DOSAGE_PATTERN}" \
                --chromosomes "{params.chromosomes}" \
                --ancestry-map "{params.ancestry_map}" \
                --global-reference-ancestry "{params.global_ref_ancestry}" \
                --top-n-loci "{params.top_n_loci}" \
                --locus-scores-out "{output.locus_scores}" \
                --locus-importance-out "{output.locus_importance}" \
                --locus-overall-out "{output.locus_overall}" \
                --ancestry-importance-out "{output.ancestry_importance}" \
                --metrics-out "{output.metrics}" \
                --top-loci-png "{output.top_loci_png}" \
                --ancestry-png "{output.ancestry_png}" \
                --heatmap-png "{output.heatmap_png}"
        ' > "{log}" 2>&1
        """


rule report_admixture_importance_done:
    input:
        validation_done=VALIDATION_DONE,
        locus_scores=os.path.join(OUTPUT_DIR, "sample_locus_dosage_scores.tsv"),
        locus_importance=os.path.join(OUTPUT_DIR, "locus_ancestry_importance.tsv"),
        locus_overall=os.path.join(OUTPUT_DIR, "locus_overall_importance.tsv"),
        ancestry_importance=os.path.join(OUTPUT_DIR, "ancestry_dosage_importance.tsv"),
        metrics=os.path.join(OUTPUT_DIR, "importance_model_metrics.tsv"),
        top_loci_png=os.path.join(OUTPUT_DIR, "top_loci_importance.png"),
        ancestry_png=os.path.join(OUTPUT_DIR, "overall_ancestry_importance.png"),
        heatmap_png=os.path.join(OUTPUT_DIR, "locus_ancestry_importance_heatmap.png"),
    output:
        done=os.path.join(OUTPUT_DIR, "report_admixture_importance.done"),
    log:
        os.path.join(LOG_DIR, "report_admixture_importance_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Admixture importance report complete" > {log}
        """
