import os
import re
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir, get_software_module, get_step5_successful_chromosomes


def _cfg(path, default=None):
    try:
        return get_analysis_value(path)
    except KeyError:
        return default


def _safe_name(value):
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        raise ValueError("Collapsed group name resolved to an empty token")
    return value


def _parse_group(group_cfg, group_key):
    if not isinstance(group_cfg, dict):
        raise ValueError(f"collapsed_tractor.{group_key} must be a mapping")

    raw_name = str(group_cfg.get("name", "")).strip()
    ancestries = group_cfg.get("ancestries", [])

    if not raw_name:
        raise ValueError(f"collapsed_tractor.{group_key}.name must be non-empty")
    if not isinstance(ancestries, list) or not ancestries:
        raise ValueError(f"collapsed_tractor.{group_key}.ancestries must be a non-empty list")

    ancestry_labels = [str(item).strip() for item in ancestries if str(item).strip()]
    if len(ancestry_labels) != len(set(ancestry_labels)):
        raise ValueError(f"collapsed_tractor.{group_key}.ancestries contains duplicates: {ancestry_labels}")

    return {
        "name_raw": raw_name,
        "name_safe": _safe_name(raw_name),
        "ancestries": ancestry_labels,
    }


def _render_pattern(template, chr_token="{chr}"):
    text = str(template)
    text = text.replace("{group1}", GROUP1["name_safe"])
    text = text.replace("{group2}", GROUP2["name_safe"])
    text = text.replace("{model}", MODEL_NAME)
    text = text.replace("{chr}", chr_token)
    return text


RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
CHROMOSOMES = [str(chrom) for chrom in get_step5_successful_chromosomes()]

GROUP1 = _parse_group(get_analysis_value("collapsed_tractor.group1"), "group1")
GROUP2 = _parse_group(get_analysis_value("collapsed_tractor.group2"), "group2")
if GROUP1["name_safe"] == GROUP2["name_safe"]:
    raise ValueError(
        "collapsed_tractor.group1.name and collapsed_tractor.group2.name must be different after filename normalization"
    )

MODEL_NAME = f"{GROUP1['name_safe']}_{GROUP2['name_safe']}"

COLLAPSED_DIR = str(_cfg("collapsed_tractor.output_dir", os.path.join(RESULTS_DIR, "tractor_collapsed")))
COLLAPSED_DONE = os.path.join(COLLAPSED_DIR, "collapse_ancestry_dosage_counts.done")
COLLAPSED_PREFIX_PATTERN = str(_cfg("collapsed_tractor.output_prefix_pattern", "chr{chr}.phased"))

PHENOTYPE_FILE = get_analysis_value("tractor_gwas.phenotype_file")
PHENOTYPE_FORMAT = str(get_analysis_value("tractor_gwas.phenotype_format"))
FAM_CASE_CODE = str(get_analysis_value("tractor_gwas.fam_case_code"))
FAM_CONTROL_CODE = str(get_analysis_value("tractor_gwas.fam_control_code"))
TABLE_IID_COLUMN = str(get_analysis_value("tractor_gwas.table_iid_column"))
TABLE_PHENO_COLUMN = str(get_analysis_value("tractor_gwas.table_phenotype_column"))

PCA_ENABLED = str(get_analysis_value("tractor_gwas.covariates.enabled")).lower() in {"1", "true", "yes", "y"}
PCA_FILE = str(get_analysis_value("tractor_gwas.covariates.pca_file"))
PCA_FORMAT = str(get_analysis_value("tractor_gwas.covariates.format"))
PCA_HAS_HEADER = str(get_analysis_value("tractor_gwas.covariates.has_header")).lower()
PCA_FID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.fid_column"))
PCA_IID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.iid_column"))
PCA_N_PCS = int(get_analysis_value("tractor_gwas.covariates.n_pcs"))

MIN_PARTITIONS = int(_cfg("collapsed_tractor.min_partitions", get_analysis_value("tractor_gwas.min_partitions")))
OUTPUT_DIR = str(_cfg("collapsed_tractor.gwas_output_dir", os.path.join(RESULTS_DIR, "tractor_gwas_collapsed")))
PER_CHR_PREFIX_PATTERN = str(_cfg("collapsed_tractor.per_chr_prefix_pattern", "chr{chr}.model_{model}"))
MERGED_PREFIX_PATTERN = str(_cfg("collapsed_tractor.merged_prefix_pattern", "merged.model_{model}"))

PYTHON_MODULE = get_software_module("python")
R_MODULE = get_software_module("r")


def collapsed_path(chr_name, group_name, kind):
    prefix = COLLAPSED_PREFIX_PATTERN.replace("{chr}", str(chr_name))
    return os.path.join(COLLAPSED_DIR, f"{prefix}.{group_name}.{kind}.txt")


def per_chr_outputs():
    output_pattern = os.path.join(OUTPUT_DIR, _render_pattern(PER_CHR_PREFIX_PATTERN, chr_token="{chr}") + ".gwas.tsv")
    return [output_pattern.replace("{chr}", str(chr_name)) for chr_name in CHROMOSOMES]


PER_CHR_TSV_PATTERN = os.path.join(OUTPUT_DIR, _render_pattern(PER_CHR_PREFIX_PATTERN, chr_token="{chr}") + ".gwas.tsv")
MERGED_PREFIX = os.path.join(OUTPUT_DIR, _render_pattern(MERGED_PREFIX_PATTERN, chr_token="{chr}"))
MERGED_TSV = MERGED_PREFIX + ".gwas.tsv"
MERGED_DOSAGE_I_MANHATTAN_PNG = MERGED_PREFIX + ".dosage_i.manhattan.png"
MERGED_DOSAGE_I_MANHATTAN_PDF = MERGED_PREFIX + ".dosage_i.manhattan.pdf"
MERGED_DOSAGE_I_QQ_PNG = MERGED_PREFIX + ".dosage_i.qq.png"
MERGED_DOSAGE_I_QQ_PDF = MERGED_PREFIX + ".dosage_i.qq.pdf"
MERGED_DOSAGE_J_MANHATTAN_PNG = MERGED_PREFIX + ".dosage_j.manhattan.png"
MERGED_DOSAGE_J_MANHATTAN_PDF = MERGED_PREFIX + ".dosage_j.manhattan.pdf"
MERGED_DOSAGE_J_QQ_PNG = MERGED_PREFIX + ".dosage_j.qq.png"
MERGED_DOSAGE_J_QQ_PDF = MERGED_PREFIX + ".dosage_j.qq.pdf"
MERGED_LOCAL_ANC_MANHATTAN_PNG = MERGED_PREFIX + ".local_ancestry.manhattan.png"
MERGED_LOCAL_ANC_MANHATTAN_PDF = MERGED_PREFIX + ".local_ancestry.manhattan.pdf"
MERGED_LOCAL_ANC_QQ_PNG = MERGED_PREFIX + ".local_ancestry.qq.png"
MERGED_LOCAL_ANC_QQ_PDF = MERGED_PREFIX + ".local_ancestry.qq.pdf"


rule tractor_gwas_collapsed_chr:
    input:
        collapsed_done=COLLAPSED_DONE,
        phenotype=PHENOTYPE_FILE,
        anc_hapcount=lambda wildcards: collapsed_path(wildcards.chr, GROUP1["name_safe"], "hapcount"),
        anc_i_dosage=lambda wildcards: collapsed_path(wildcards.chr, GROUP1["name_safe"], "dosage"),
        anc_j_dosage=lambda wildcards: collapsed_path(wildcards.chr, GROUP2["name_safe"], "dosage"),
    output:
        tsv=PER_CHR_TSV_PATTERN,
    log:
        os.path.join(LOG_DIR, "tractor_gwas_collapsed_chr{chr}.log"),
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
        pca_file=PCA_FILE,
        min_partitions=MIN_PARTITIONS,
        label_i=GROUP1["name_raw"],
        label_j=GROUP2["name_raw"],
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash {params.helper} \
            --chr "{wildcards.chr}" \
            --anc-hapcount "{input.anc_hapcount}" \
            --anc-i-dosage "{input.anc_i_dosage}" \
            --anc-j-dosage "{input.anc_j_dosage}" \
            --label-i "{params.label_i}" \
            --label-j "{params.label_j}" \
            --phenotype "{input.phenotype}" \
            --phenotype-format "{params.phenotype_format}" \
            --fam-case-code "{params.fam_case_code}" \
            --fam-control-code "{params.fam_control_code}" \
            --table-iid-column "{params.table_iid_column}" \
            --table-phenotype-column "{params.table_pheno_column}" \
            --pca-enabled "{params.pca_enabled}" \
            --pca-file "{params.pca_file}" \
            --pca-format "{params.pca_format}" \
            --pca-has-header "{params.pca_has_header}" \
            --pca-fid-column "{params.pca_fid_column}" \
            --pca-iid-column "{params.pca_iid_column}" \
            --pca-n-pcs "{params.pca_n_pcs}" \
            --min-partitions "{params.min_partitions}" \
            --output-tsv "{output.tsv}" \
            --script-path "{params.script}" \
            --python-module "{params.python_module}" \
            > {log} 2>&1
        """


rule tractor_gwas_collapsed_merge_plot:
    input:
        collapsed_done=COLLAPSED_DONE,
        per_chr=lambda wildcards: per_chr_outputs(),
    output:
        merged=MERGED_TSV,
        dosage_i_manhattan_png=MERGED_DOSAGE_I_MANHATTAN_PNG,
        dosage_i_manhattan_pdf=MERGED_DOSAGE_I_MANHATTAN_PDF,
        dosage_i_qq_png=MERGED_DOSAGE_I_QQ_PNG,
        dosage_i_qq_pdf=MERGED_DOSAGE_I_QQ_PDF,
        dosage_j_manhattan_png=MERGED_DOSAGE_J_MANHATTAN_PNG,
        dosage_j_manhattan_pdf=MERGED_DOSAGE_J_MANHATTAN_PDF,
        dosage_j_qq_png=MERGED_DOSAGE_J_QQ_PNG,
        dosage_j_qq_pdf=MERGED_DOSAGE_J_QQ_PDF,
        local_anc_manhattan_png=MERGED_LOCAL_ANC_MANHATTAN_PNG,
        local_anc_manhattan_pdf=MERGED_LOCAL_ANC_MANHATTAN_PDF,
        local_anc_qq_png=MERGED_LOCAL_ANC_QQ_PNG,
        local_anc_qq_pdf=MERGED_LOCAL_ANC_QQ_PDF,
    log:
        os.path.join(LOG_DIR, "tractor_gwas_collapsed_merge.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/tractor_merge_plot_pairwise.R",
        r_module=R_MODULE,
        pair_name=MODEL_NAME,
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
                --pair "{params.pair_name}" \
                --input-tsv "$(printf '%s,' {input.per_chr} | sed '"'"'s/,$//'"'"')" \
                --output-merged "{output.merged}" \
                --output-dosage-i-manhattan-png "{output.dosage_i_manhattan_png}" \
                --output-dosage-i-manhattan-pdf "{output.dosage_i_manhattan_pdf}" \
                --output-dosage-i-qq-png "{output.dosage_i_qq_png}" \
                --output-dosage-i-qq-pdf "{output.dosage_i_qq_pdf}" \
                --output-dosage-j-manhattan-png "{output.dosage_j_manhattan_png}" \
                --output-dosage-j-manhattan-pdf "{output.dosage_j_manhattan_pdf}" \
                --output-dosage-j-qq-png "{output.dosage_j_qq_png}" \
                --output-dosage-j-qq-pdf "{output.dosage_j_qq_pdf}" \
                --output-local-ancestry-manhattan-png "{output.local_anc_manhattan_png}" \
                --output-local-ancestry-manhattan-pdf "{output.local_anc_manhattan_pdf}" \
                --output-local-ancestry-qq-png "{output.local_anc_qq_png}" \
                --output-local-ancestry-qq-pdf "{output.local_anc_qq_pdf}"
        ' > "{log}" 2>&1
        """


rule tractor_gwas_collapsed_done:
    input:
        collapsed_done=COLLAPSED_DONE,
        per_chr=lambda wildcards: per_chr_outputs(),
        merged=MERGED_TSV,
        dosage_i_manhattan_png=MERGED_DOSAGE_I_MANHATTAN_PNG,
        dosage_i_manhattan_pdf=MERGED_DOSAGE_I_MANHATTAN_PDF,
        dosage_i_qq_png=MERGED_DOSAGE_I_QQ_PNG,
        dosage_i_qq_pdf=MERGED_DOSAGE_I_QQ_PDF,
        dosage_j_manhattan_png=MERGED_DOSAGE_J_MANHATTAN_PNG,
        dosage_j_manhattan_pdf=MERGED_DOSAGE_J_MANHATTAN_PDF,
        dosage_j_qq_png=MERGED_DOSAGE_J_QQ_PNG,
        dosage_j_qq_pdf=MERGED_DOSAGE_J_QQ_PDF,
        local_anc_manhattan_png=MERGED_LOCAL_ANC_MANHATTAN_PNG,
        local_anc_manhattan_pdf=MERGED_LOCAL_ANC_MANHATTAN_PDF,
        local_anc_qq_png=MERGED_LOCAL_ANC_QQ_PNG,
        local_anc_qq_pdf=MERGED_LOCAL_ANC_QQ_PDF,
    output:
        done=os.path.join(OUTPUT_DIR, "tractor_gwas_collapsed.done"),
    log:
        os.path.join(LOG_DIR, "tractor_gwas_collapsed_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        g1_name=GROUP1["name_safe"],
        g2_name=GROUP2["name_safe"],
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Collapsed Tractor GWAS complete: {params.g1_name} vs {params.g2_name}" > {log}
        """
