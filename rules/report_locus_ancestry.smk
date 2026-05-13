import os
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("locus_ancestry_report.output_dir")
DONE_DIR = OUTPUT_DIR

INPUT_DONE = os.path.join(RESULTS_DIR, "maf_gwas_summary", "merge_maf_gwas_summary.done")
INPUT_TABLE = get_analysis_value("locus_ancestry_report.input_tsv")
OUTPUT_TABLE = os.path.join(OUTPUT_DIR, "locus_ancestry_classification.tsv")
FILTERED_SNP_TABLE = os.path.join(OUTPUT_DIR, "filtered_snps_for_locus_report.tsv")
OUTPUT_EXCEL = os.path.join(OUTPUT_DIR, "locus_ancestry_classification.xlsx")
OUTPUT_DENDROGRAM = os.path.join(OUTPUT_DIR, "locus_ancestry_dendrogram.png")
OUTPUT_PANEL_TREE = os.path.join(OUTPUT_DIR, "locus_ancestry_manhattan_tree.png")
LOCI_FILE = get_analysis_value("step3.loci_file")

PYTHON_MODULE = get_software_module("python")
R_MODULE = get_software_module("r")


rule build_locus_ancestry_report:
    input:
        done=INPUT_DONE,
        tsv=INPUT_TABLE,
        loci=LOCI_FILE,
    output:
        locus=OUTPUT_TABLE,
        snps=FILTERED_SNP_TABLE,
        excel=OUTPUT_EXCEL,
        dendrogram=OUTPUT_DENDROGRAM,
        panel_tree=OUTPUT_PANEL_TREE,
        done=os.path.join(DONE_DIR, "locus_ancestry_report.done"),
    log:
        os.path.join(LOG_DIR, "build_locus_ancestry_report.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/report_locus_ancestry.py",
        r_script="scripts/report_locus_ancestry_excel_plot.R",
        python_module=PYTHON_MODULE,
        r_module=R_MODULE,
        cohort_maf_min=get_analysis_value("locus_ancestry_report.cohort_maf_min"),
        ancestry_maf_min=get_analysis_value("locus_ancestry_report.ancestry_maf_min"),
        local_ancestry_min_logp=get_analysis_value("locus_ancestry_report.local_ancestry_min_logp"),
        dosage_gwas_max_p=get_analysis_value("locus_ancestry_report.dosage_gwas_max_p"),
        min_supporting_snps=get_analysis_value("locus_ancestry_report.min_supporting_snps"),
        clustering_method=get_analysis_value("locus_ancestry_report.clustering_method"),
        clustering_metric=get_analysis_value("locus_ancestry_report.clustering_metric"),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        if type module >/dev/null 2>&1; then
            module load "{params.python_module}" || true
            module load "{params.r_module}" || true
        fi
        python3 "{params.script}" \
            --input "{input.tsv}" \
            --output-locus "{output.locus}" \
            --output-snp "{output.snps}" \
            --output-excel "{output.excel}" \
            --output-dendrogram "{output.dendrogram}" \
            --cohort-maf-min "{params.cohort_maf_min}" \
            --ancestry-maf-min "{params.ancestry_maf_min}" \
            --local-ancestry-min-logp "{params.local_ancestry_min_logp}" \
            --dosage-gwas-max-p "{params.dosage_gwas_max_p}" \
            --min-supporting-snps "{params.min_supporting_snps}" \
            --clustering-method "{params.clustering_method}" \
            --clustering-metric "{params.clustering_metric}" \
            > "{log}" 2>&1
        Rscript "{params.r_script}" \
            --input-locus "{output.locus}" \
            --input-snp "{output.snps}" \
            --output-excel "{output.excel}" \
            --output-dendrogram "{output.dendrogram}" \
            --output-panel-tree "{output.panel_tree}" \
            --loci-file "{input.loci}" \
            --cohort-maf-min "{params.cohort_maf_min}" \
            --ancestry-maf-min "{params.ancestry_maf_min}" \
            --local-ancestry-min-logp "{params.local_ancestry_min_logp}" \
            --dosage-gwas-max-p "{params.dosage_gwas_max_p}" \
            --min-supporting-snps "{params.min_supporting_snps}" \
            --clustering-method "{params.clustering_method}" \
            --clustering-metric "{params.clustering_metric}" \
            >> "{log}" 2>&1
        touch "{output.done}"
        """