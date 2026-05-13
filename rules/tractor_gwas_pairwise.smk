import os
import sys
from itertools import combinations

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "tractor_gwas")

STEP9_DONE = os.path.join(RESULTS_DIR, "tractor", "extract_tracts.done")

CHROMOSOMES = get_step5_successful_chromosomes()
EXTRACT_DIR = get_analysis_value("tractor_gwas.extract_output_dir")
PHENOTYPE_FILE = get_analysis_value("tractor_gwas.phenotype_file")
PHENOTYPE_FORMAT = str(get_analysis_value("tractor_gwas.phenotype_format"))
FAM_CASE_CODE = str(get_analysis_value("tractor_gwas.fam_case_code"))
FAM_CONTROL_CODE = str(get_analysis_value("tractor_gwas.fam_control_code"))
TABLE_IID_COLUMN = str(get_analysis_value("tractor_gwas.table_iid_column"))
TABLE_PHENO_COLUMN = str(get_analysis_value("tractor_gwas.table_phenotype_column"))
OUTPUT_DIR = get_analysis_value("tractor_gwas.output_dir")
PER_CHR_PREFIX_PATTERN = get_analysis_value("tractor_gwas.per_chr_prefix_pattern")
MERGED_PREFIX_PATTERN = get_analysis_value("tractor_gwas.merged_prefix_pattern")
MIN_PARTITIONS = int(get_analysis_value("tractor_gwas.min_partitions"))
ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
PCA_ENABLED = str(get_analysis_value("tractor_gwas.covariates.enabled")).lower() in {"1", "true", "yes", "y"}
PCA_FILE = str(get_analysis_value("tractor_gwas.covariates.pca_file"))
PCA_FORMAT = str(get_analysis_value("tractor_gwas.covariates.format"))
PCA_HAS_HEADER = str(get_analysis_value("tractor_gwas.covariates.has_header")).lower()
PCA_FID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.fid_column"))
PCA_IID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.iid_column"))
PCA_N_PCS = int(get_analysis_value("tractor_gwas.covariates.n_pcs"))

if not isinstance(ANCESTRY_LABELS, dict) or not ANCESTRY_LABELS:
    raise ValueError("tractor_gwas.ancestry_labels must be a non-empty mapping")

INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
EXPECTED_CODES = list(range(len(INDICES)))
FOUND_CODES = [int(index) for index in INDICES]
if FOUND_CODES != EXPECTED_CODES:
    raise ValueError(
        "tractor_gwas.ancestry_labels must use contiguous ancestry codes 0..N-1. "
        f"Found: {FOUND_CODES}; Expected: {EXPECTED_CODES}"
    )

LABELS_IN_ORDER = [str(ANCESTRY_LABELS[index]) for index in INDICES]
if len(set(LABELS_IN_ORDER)) != len(LABELS_IN_ORDER):
    raise ValueError("tractor_gwas.ancestry_labels must map each code to a unique ancestry label")

PAIR_INFO = []
for i, j in combinations(INDICES, 2):
    label_i = str(ANCESTRY_LABELS[i])
    label_j = str(ANCESTRY_LABELS[j])
    pair_name = f"{label_i}_{label_j}"
    PAIR_INFO.append((pair_name, i, j, label_i, label_j))

PAIR_NAMES = [item[0] for item in PAIR_INFO]
PAIR_TO_ANCESTRY = {item[0]: item for item in PAIR_INFO}

EXPECTED_EXTRACT_DIR = get_analysis_value("tractor.output_dir")
if EXTRACT_DIR != EXPECTED_EXTRACT_DIR:
    raise ValueError(
        "tractor_gwas.extract_output_dir is not compatible with tractor outputs. "
        f"Expected: {EXPECTED_EXTRACT_DIR}; Found: {EXTRACT_DIR}"
    )

PYTHON_MODULE = get_software_module("python")
R_MODULE = get_software_module("r")


def pair_indices(pair_name):
    _, i, j, _, _ = PAIR_TO_ANCESTRY[pair_name]
    return str(i), str(j)


def pair_labels(pair_name):
    _, _, _, label_i, label_j = PAIR_TO_ANCESTRY[pair_name]
    return label_i, label_j


def hapcount_path(chr_name, ancestry_index):
    return os.path.join(EXTRACT_DIR, f"chr{chr_name}.phased.anc{ancestry_index}.hapcount.txt")


def dosage_path(chr_name, ancestry_index):
    return os.path.join(EXTRACT_DIR, f"chr{chr_name}.phased.anc{ancestry_index}.dosage.txt")


def per_chr_pair_inputs(pair_name):
    return [
        os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN.format(chr=chr_name, pair=pair_name) + ".gwas.tsv")
        for chr_name in CHROMOSOMES
    ]


PER_CHR_PATTERN = os.path.join(OUTPUT_DIR, PER_CHR_PREFIX_PATTERN + ".gwas.tsv")
MERGED_PATTERN = os.path.join(OUTPUT_DIR, MERGED_PREFIX_PATTERN + ".gwas.tsv")
PLOT_PNG_PATTERN = os.path.join(OUTPUT_DIR, MERGED_PREFIX_PATTERN + ".manhattan.png")
PLOT_PDF_PATTERN = os.path.join(OUTPUT_DIR, MERGED_PREFIX_PATTERN + ".manhattan.pdf")


rule tractor_gwas_pairwise_chr:
    input:
        step9_done=STEP9_DONE,
        phenotype=PHENOTYPE_FILE,
        pca=PCA_FILE,
        anc_hapcount=lambda wildcards: hapcount_path(wildcards.chr, pair_indices(wildcards.pair)[0]),
        anc_i_dosage=lambda wildcards: dosage_path(wildcards.chr, pair_indices(wildcards.pair)[0]),
        anc_j_dosage=lambda wildcards: dosage_path(wildcards.chr, pair_indices(wildcards.pair)[1]),
    output:
        tsv=PER_CHR_PATTERN,
    log:
        os.path.join(LOG_DIR, "tractor_gwas_pairwise_chr{chr}_{pair}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        helper="scripts/tractor_gwas_pairwise_chr.sh",
        script="scripts/tractor_gwas_pairwise_chr.py",
        python_module=PYTHON_MODULE,
        phenotype_format=PHENOTYPE_FORMAT,
        fam_case_code=FAM_CASE_CODE,
        fam_control_code=FAM_CONTROL_CODE,
        table_iid_column=TABLE_IID_COLUMN,
        table_pheno_column=TABLE_PHENO_COLUMN,
        pca_enabled=str(PCA_ENABLED).lower(),
        pca_format=PCA_FORMAT,
        pca_has_header=PCA_HAS_HEADER,
        pca_fid_column=PCA_FID_COLUMN,
        pca_iid_column=PCA_IID_COLUMN,
        pca_n_pcs=PCA_N_PCS,
        min_partitions=MIN_PARTITIONS,
        anc_i_idx=lambda wildcards: pair_indices(wildcards.pair)[0],
        anc_j_idx=lambda wildcards: pair_indices(wildcards.pair)[1],
        label_i=lambda wildcards: pair_labels(wildcards.pair)[0],
        label_j=lambda wildcards: pair_labels(wildcards.pair)[1],
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --chr {wildcards.chr} \
            --anc-hapcount {input.anc_hapcount} \
            --anc-i-dosage {input.anc_i_dosage} \
            --anc-j-dosage {input.anc_j_dosage} \
            --label-i {params.label_i} \
            --label-j {params.label_j} \
            --phenotype {input.phenotype} \
            --phenotype-format {params.phenotype_format} \
            --fam-case-code {params.fam_case_code} \
            --fam-control-code {params.fam_control_code} \
            --table-iid-column {params.table_iid_column} \
            --table-phenotype-column {params.table_pheno_column} \
            --pca-enabled {params.pca_enabled} \
            --pca-file {input.pca} \
            --pca-format {params.pca_format} \
            --pca-has-header {params.pca_has_header} \
            --pca-fid-column {params.pca_fid_column} \
            --pca-iid-column {params.pca_iid_column} \
            --pca-n-pcs {params.pca_n_pcs} \
            --min-partitions {params.min_partitions} \
            --output-tsv {output.tsv} \
            --script-path {params.script} \
            --python-module {params.python_module} \
            > {log} 2>&1
        """


rule tractor_gwas_pairwise_merge_plot:
    input:
        step9_done=STEP9_DONE,
        per_chr=lambda wildcards: per_chr_pair_inputs(wildcards.pair),
    output:
        merged=MERGED_PATTERN,
        png=PLOT_PNG_PATTERN,
        pdf=PLOT_PDF_PATTERN,
    log:
        os.path.join(LOG_DIR, "tractor_gwas_pairwise_merge_{pair}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/tractor_merge_plot_pairwise.R",
        r_module=R_MODULE,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            if type module >/dev/null 2>&1; then
                module load {params.r_module} >/dev/null 2>&1 || true
            fi
            if type ml >/dev/null 2>&1; then
                ml {params.r_module} >/dev/null 2>&1 || true
            fi
            command -v Rscript >/dev/null 2>&1 || {{ echo "Rscript not found after module load"; exit 127; }}
            Rscript "{params.script}" \
                --pair "{wildcards.pair}" \
                --input-tsv "$(printf '%s,' {input.per_chr} | sed '"'"'s/,$//'"'"')" \
                --output-merged "{output.merged}" \
                --output-png "{output.png}" \
                --output-pdf "{output.pdf}"
        ' > "{log}" 2>&1
        """


rule tractor_gwas_pairwise_done:
    input:
        step9_done=STEP9_DONE,
        per_chr=expand(PER_CHR_PATTERN, chr=CHROMOSOMES, pair=PAIR_NAMES),
        merged=expand(MERGED_PATTERN, pair=PAIR_NAMES),
        png=expand(PLOT_PNG_PATTERN, pair=PAIR_NAMES),
        pdf=expand(PLOT_PDF_PATTERN, pair=PAIR_NAMES),
    output:
        done=os.path.join(DONE_DIR, "tractor_gwas_pairwise.done"),
    log:
        os.path.join(LOG_DIR, "tractor_gwas_pairwise_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        checker="scripts/check_tractor_gwas_pairwise_integrity.py",
        output_dir=OUTPUT_DIR,
        pair_names=" ".join(PAIR_NAMES),
        chromosomes=" ".join(str(chrom) for chrom in CHROMOSOMES),
    shell:
        """
        mkdir -p {DONE_DIR}
        python3 {params.checker} \
            --output-dir {params.output_dir} \
            --pairs {params.pair_names} \
            --chromosomes {params.chromosomes} \
            > {log} 2>&1
        touch {output.done}
        echo "Tractor pairwise GWAS complete (ancestry integrity verified)" >> {log}
        """
