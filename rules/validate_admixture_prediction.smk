import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value, get_default_resource, get_step5_successful_chromosomes

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = get_analysis_value("admixture_validation.output_dir")
REPEATED_8020_OUTPUT_DIR = os.path.join(f"{OUTPUT_DIR}_repeated_8020")

RUN_RFMIX_DONE = os.path.join(get_analysis_value("step5.output_dir"), "run_rfmix.done")
EXTRACT_TRACTS_DONE = os.path.join(get_analysis_value("tractor.output_dir"), "extract_tracts.done")
TRACTOR_GWAS_DONE = os.path.join(get_analysis_value("tractor_gwas.output_dir"), "tractor_gwas_pairwise.done")

CHROMOSOMES = [str(chr_name) for chr_name in get_step5_successful_chromosomes()]
ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
ANCESTRY_INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
ANCESTRY_MAP = ",".join(f"{idx}:{ANCESTRY_LABELS[idx]}" for idx in ANCESTRY_INDICES)

SUMMARY_TSV = get_analysis_value("admixture_validation.summary_tsv")
RFMIX_Q_PATTERN = get_analysis_value("admixture_validation.rfmix_q_pattern")
TRACTOR_DOSAGE_PATTERN = get_analysis_value("admixture_validation.tractor_dosage_pattern")
PHENOTYPE_FILE = get_analysis_value("tractor_gwas.phenotype_file")
PHENOTYPE_FORMAT = str(get_analysis_value("tractor_gwas.phenotype_format"))
FAM_CASE_CODE = str(get_analysis_value("tractor_gwas.fam_case_code"))
FAM_CONTROL_CODE = str(get_analysis_value("tractor_gwas.fam_control_code"))
TABLE_IID_COLUMN = str(get_analysis_value("tractor_gwas.table_iid_column"))
TABLE_PHENO_COLUMN = str(get_analysis_value("tractor_gwas.table_phenotype_column"))
COVARIATES_FILE = str(get_analysis_value("tractor_gwas.covariates.pca_file"))
COVARIATES_FORMAT = str(get_analysis_value("tractor_gwas.covariates.format"))
COVARIATES_HAS_HEADER = str(get_analysis_value("tractor_gwas.covariates.has_header")).lower()
COVARIATES_FID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.fid_column"))
COVARIATES_IID_COLUMN = str(get_analysis_value("tractor_gwas.covariates.iid_column"))
N_PCS = int(get_analysis_value("tractor_gwas.covariates.n_pcs"))
INCLUDE_SEX = str(get_analysis_value("admixture_validation.include_sex")).lower()
INCLUDE_AGE = str(get_analysis_value("admixture_validation.include_age")).lower()
AGE_FILE = str(get_analysis_value("admixture_validation.age_file"))
AGE_IID_COLUMN = str(get_analysis_value("admixture_validation.age_iid_column"))
AGE_COLUMN = str(get_analysis_value("admixture_validation.age_column"))
GLOBAL_REF_ANCESTRY = str(get_analysis_value("admixture_validation.global_ancestry_reference"))
SNP_P_THRESHOLD = float(get_analysis_value("admixture_validation.snp_pvalue_threshold"))
MAX_SNPS_PER_ANCESTRY = int(get_analysis_value("admixture_validation.max_snps_per_ancestry"))
USE_TOP_SNP_PER_LOCUS = str(get_analysis_value("admixture_validation.use_top_snp_per_locus")).lower()
CV_FOLDS = int(get_analysis_value("admixture_validation.cv_folds"))
RANDOM_STATE = int(get_analysis_value("admixture_validation.random_state"))
REPEATED_8020_REPEATS = 10
REPEATED_8020_TEST_SIZE = 0.2

R_MODULE = get_software_module("r")


def resolve_pattern(pattern, chr_name=None, anc=None):
    resolved = pattern
    if chr_name is not None:
        resolved = resolved.replace("{chr}", str(chr_name))
    if anc is not None:
        resolved = resolved.replace("{anc}", str(anc))
    return resolved


rule validate_admixture_prediction:
    input:
        run_rfmix_done=RUN_RFMIX_DONE,
        extract_tracts_done=EXTRACT_TRACTS_DONE,
        tractor_gwas_done=TRACTOR_GWAS_DONE,
        summary=SUMMARY_TSV,
        phenotype=PHENOTYPE_FILE,
        covariates=COVARIATES_FILE,
        q_files=expand(RFMIX_Q_PATTERN, chr=CHROMOSOMES),
        dosage_files=expand(TRACTOR_DOSAGE_PATTERN, chr=CHROMOSOMES, anc=ANCESTRY_INDICES),
    output:
        sample_features=os.path.join(OUTPUT_DIR, "sample_prediction_features.tsv"),
        selected_snps=os.path.join(OUTPUT_DIR, "selected_admixture_snps.tsv"),
        cv_predictions=os.path.join(OUTPUT_DIR, "cv_predictions.tsv"),
        metrics=os.path.join(OUTPUT_DIR, "prediction_metrics.tsv"),
        coefficients=os.path.join(OUTPUT_DIR, "model_coefficients.tsv"),
        roc_png=os.path.join(OUTPUT_DIR, "roc_curve.png"),
        calibration_png=os.path.join(OUTPUT_DIR, "calibration_curve.png"),
    log:
        os.path.join(LOG_DIR, "validate_admixture_prediction.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/validate_admixture_prediction.R",
        r_module=R_MODULE,
        chromosomes=",".join(CHROMOSOMES),
        ancestry_map=ANCESTRY_MAP,
        phenotype_format=PHENOTYPE_FORMAT,
        fam_case_code=FAM_CASE_CODE,
        fam_control_code=FAM_CONTROL_CODE,
        table_iid_column=TABLE_IID_COLUMN,
        table_pheno_column=TABLE_PHENO_COLUMN,
        covariates_format=COVARIATES_FORMAT,
        covariates_has_header=COVARIATES_HAS_HEADER,
        covariates_fid_column=COVARIATES_FID_COLUMN,
        covariates_iid_column=COVARIATES_IID_COLUMN,
        n_pcs=N_PCS,
        include_sex=INCLUDE_SEX,
        include_age=INCLUDE_AGE,
        age_file=AGE_FILE,
        age_iid_column=AGE_IID_COLUMN,
        age_column=AGE_COLUMN,
        global_ref_ancestry=GLOBAL_REF_ANCESTRY,
        snp_p_threshold=SNP_P_THRESHOLD,
        max_snps_per_ancestry=MAX_SNPS_PER_ANCESTRY,
        use_top_snp_per_locus=USE_TOP_SNP_PER_LOCUS,
        cv_folds=CV_FOLDS,
        random_state=RANDOM_STATE,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            ml "{params.r_module}" >/dev/null 2>&1 || true
            Rscript "{params.script}" \
                --summary-tsv "{input.summary}" \
                --rfmix-q-pattern "{RFMIX_Q_PATTERN}" \
                --tractor-dosage-pattern "{TRACTOR_DOSAGE_PATTERN}" \
                --chromosomes "{params.chromosomes}" \
                --ancestry-map "{params.ancestry_map}" \
                --phenotype "{input.phenotype}" \
                --phenotype-format "{params.phenotype_format}" \
                --fam-case-code "{params.fam_case_code}" \
                --fam-control-code "{params.fam_control_code}" \
                --table-iid-column "{params.table_iid_column}" \
                --table-phenotype-column "{params.table_pheno_column}" \
                --covariates-file "{input.covariates}" \
                --covariates-format "{params.covariates_format}" \
                --covariates-has-header "{params.covariates_has_header}" \
                --covariates-fid-column "{params.covariates_fid_column}" \
                --covariates-iid-column "{params.covariates_iid_column}" \
                --n-pcs "{params.n_pcs}" \
                --include-sex "{params.include_sex}" \
                --include-age "{params.include_age}" \
                --age-file "{params.age_file}" \
                --age-iid-column "{params.age_iid_column}" \
                --age-column "{params.age_column}" \
                --global-reference-ancestry "{params.global_ref_ancestry}" \
                --snp-pvalue-threshold "{params.snp_p_threshold}" \
                --max-snps-per-ancestry "{params.max_snps_per_ancestry}" \
                --use-top-snp-per-locus "{params.use_top_snp_per_locus}" \
                --cv-folds "{params.cv_folds}" \
                --random-state "{params.random_state}" \
                --sample-features-out "{output.sample_features}" \
                --selected-snps-out "{output.selected_snps}" \
                --cv-predictions-out "{output.cv_predictions}" \
                --metrics-out "{output.metrics}" \
                --coefficients-out "{output.coefficients}" \
                --roc-png "{output.roc_png}" \
                --calibration-png "{output.calibration_png}"
        ' > "{log}" 2>&1
        """


rule validate_admixture_prediction_done:
    input:
        run_rfmix_done=RUN_RFMIX_DONE,
        extract_tracts_done=EXTRACT_TRACTS_DONE,
        tractor_gwas_done=TRACTOR_GWAS_DONE,
        sample_features=os.path.join(OUTPUT_DIR, "sample_prediction_features.tsv"),
        selected_snps=os.path.join(OUTPUT_DIR, "selected_admixture_snps.tsv"),
        cv_predictions=os.path.join(OUTPUT_DIR, "cv_predictions.tsv"),
        metrics=os.path.join(OUTPUT_DIR, "prediction_metrics.tsv"),
        coefficients=os.path.join(OUTPUT_DIR, "model_coefficients.tsv"),
        roc_png=os.path.join(OUTPUT_DIR, "roc_curve.png"),
        calibration_png=os.path.join(OUTPUT_DIR, "calibration_curve.png"),
    output:
        done=os.path.join(OUTPUT_DIR, "validate_admixture_prediction.done"),
    log:
        os.path.join(LOG_DIR, "validate_admixture_prediction_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Admixture validation step complete" > {log}
        """


rule validate_admixture_prediction_repeated_8020:
    input:
        run_rfmix_done=RUN_RFMIX_DONE,
        extract_tracts_done=EXTRACT_TRACTS_DONE,
        tractor_gwas_done=TRACTOR_GWAS_DONE,
        summary=SUMMARY_TSV,
        phenotype=PHENOTYPE_FILE,
        covariates=COVARIATES_FILE,
        q_files=expand(RFMIX_Q_PATTERN, chr=CHROMOSOMES),
        dosage_files=expand(TRACTOR_DOSAGE_PATTERN, chr=CHROMOSOMES, anc=ANCESTRY_INDICES),
    output:
        sample_features=os.path.join(REPEATED_8020_OUTPUT_DIR, "sample_prediction_features.tsv"),
        selected_snps=os.path.join(REPEATED_8020_OUTPUT_DIR, "selected_admixture_snps.tsv"),
        cv_predictions=os.path.join(REPEATED_8020_OUTPUT_DIR, "cv_predictions.tsv"),
        metrics=os.path.join(REPEATED_8020_OUTPUT_DIR, "prediction_metrics.tsv"),
        coefficients=os.path.join(REPEATED_8020_OUTPUT_DIR, "model_coefficients.tsv"),
        roc_png=os.path.join(REPEATED_8020_OUTPUT_DIR, "roc_curve.png"),
        calibration_png=os.path.join(REPEATED_8020_OUTPUT_DIR, "calibration_curve.png"),
    log:
        os.path.join(LOG_DIR, "validate_admixture_prediction_repeated_8020.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/validate_admixture_prediction.R",
        r_module=R_MODULE,
        chromosomes=",".join(CHROMOSOMES),
        ancestry_map=ANCESTRY_MAP,
        phenotype_format=PHENOTYPE_FORMAT,
        fam_case_code=FAM_CASE_CODE,
        fam_control_code=FAM_CONTROL_CODE,
        table_iid_column=TABLE_IID_COLUMN,
        table_pheno_column=TABLE_PHENO_COLUMN,
        covariates_format=COVARIATES_FORMAT,
        covariates_has_header=COVARIATES_HAS_HEADER,
        covariates_fid_column=COVARIATES_FID_COLUMN,
        covariates_iid_column=COVARIATES_IID_COLUMN,
        n_pcs=N_PCS,
        include_sex=INCLUDE_SEX,
        include_age=INCLUDE_AGE,
        age_file=AGE_FILE,
        age_iid_column=AGE_IID_COLUMN,
        age_column=AGE_COLUMN,
        global_ref_ancestry=GLOBAL_REF_ANCESTRY,
        snp_p_threshold=SNP_P_THRESHOLD,
        max_snps_per_ancestry=MAX_SNPS_PER_ANCESTRY,
        use_top_snp_per_locus=USE_TOP_SNP_PER_LOCUS,
        repeats=REPEATED_8020_REPEATS,
        holdout_test_size=REPEATED_8020_TEST_SIZE,
        random_state=RANDOM_STATE,
    shell:
        """
        mkdir -p {REPEATED_8020_OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            ml "{params.r_module}" >/dev/null 2>&1 || true
            Rscript "{params.script}" \
                --summary-tsv "{input.summary}" \
                --rfmix-q-pattern "{RFMIX_Q_PATTERN}" \
                --tractor-dosage-pattern "{TRACTOR_DOSAGE_PATTERN}" \
                --chromosomes "{params.chromosomes}" \
                --ancestry-map "{params.ancestry_map}" \
                --phenotype "{input.phenotype}" \
                --phenotype-format "{params.phenotype_format}" \
                --fam-case-code "{params.fam_case_code}" \
                --fam-control-code "{params.fam_control_code}" \
                --table-iid-column "{params.table_iid_column}" \
                --table-phenotype-column "{params.table_pheno_column}" \
                --covariates-file "{input.covariates}" \
                --covariates-format "{params.covariates_format}" \
                --covariates-has-header "{params.covariates_has_header}" \
                --covariates-fid-column "{params.covariates_fid_column}" \
                --covariates-iid-column "{params.covariates_iid_column}" \
                --n-pcs "{params.n_pcs}" \
                --include-sex "{params.include_sex}" \
                --include-age "{params.include_age}" \
                --age-file "{params.age_file}" \
                --age-iid-column "{params.age_iid_column}" \
                --age-column "{params.age_column}" \
                --global-reference-ancestry "{params.global_ref_ancestry}" \
                --snp-pvalue-threshold "{params.snp_p_threshold}" \
                --max-snps-per-ancestry "{params.max_snps_per_ancestry}" \
                --use-top-snp-per-locus "{params.use_top_snp_per_locus}" \
                --validation-design "repeated_stratified_holdout_80_20" \
                --n-repeats "{params.repeats}" \
                --holdout-test-size "{params.holdout_test_size}" \
                --random-state "{params.random_state}" \
                --sample-features-out "{output.sample_features}" \
                --selected-snps-out "{output.selected_snps}" \
                --cv-predictions-out "{output.cv_predictions}" \
                --metrics-out "{output.metrics}" \
                --coefficients-out "{output.coefficients}" \
                --roc-png "{output.roc_png}" \
                --calibration-png "{output.calibration_png}"
        ' > "{log}" 2>&1
        """


rule validate_admixture_prediction_repeated_8020_done:
    input:
        run_rfmix_done=RUN_RFMIX_DONE,
        extract_tracts_done=EXTRACT_TRACTS_DONE,
        tractor_gwas_done=TRACTOR_GWAS_DONE,
        sample_features=os.path.join(REPEATED_8020_OUTPUT_DIR, "sample_prediction_features.tsv"),
        selected_snps=os.path.join(REPEATED_8020_OUTPUT_DIR, "selected_admixture_snps.tsv"),
        cv_predictions=os.path.join(REPEATED_8020_OUTPUT_DIR, "cv_predictions.tsv"),
        metrics=os.path.join(REPEATED_8020_OUTPUT_DIR, "prediction_metrics.tsv"),
        coefficients=os.path.join(REPEATED_8020_OUTPUT_DIR, "model_coefficients.tsv"),
        roc_png=os.path.join(REPEATED_8020_OUTPUT_DIR, "roc_curve.png"),
        calibration_png=os.path.join(REPEATED_8020_OUTPUT_DIR, "calibration_curve.png"),
    output:
        done=os.path.join(REPEATED_8020_OUTPUT_DIR, "validate_admixture_prediction_repeated_8020.done"),
    log:
        os.path.join(LOG_DIR, "validate_admixture_prediction_repeated_8020_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {REPEATED_8020_OUTPUT_DIR}
        touch {output.done}
        echo "Admixture validation repeated stratified 80/20 step complete" > {log}
        """
