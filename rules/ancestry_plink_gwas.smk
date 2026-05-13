import os
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module, get_software_param, get_default_resource, get_step5_successful_chromosomes, get_gwas_validation_files

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")

DOSAGE_PLINK_DONE = os.path.join(RESULTS_DIR, "ancestry_plink_dosage", "ancestry_dosage_plink.done")
DOSAGE_PLINK_DIR = os.path.join(RESULTS_DIR, "ancestry_plink_dosage")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "ancestry_plink_gwas")

CHROMOSOMES = [str(chr_name) for chr_name in get_step5_successful_chromosomes()]
ANCESTRY_LABELS = get_analysis_value("tractor_gwas.ancestry_labels")
ANCESTRY_INDICES = sorted(ANCESTRY_LABELS.keys(), key=lambda x: int(x))
ANCESTRY_NAMES = sorted([str(ANCESTRY_LABELS[idx]) for idx in ANCESTRY_INDICES])

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
SUMMARY_TSV = get_analysis_value("admixture_validation.summary_tsv")
REFERENCE_ANCESTRY = str(get_analysis_value("admixture_validation.global_ancestry_reference"))
GENOME_BUILD = str(get_analysis_value("admixture_validation.genome_build"))

PLINK_MODULE = get_software_module("plink2")
R_MODULE = get_software_module("r")
R_LIBS_USER = str(get_software_param("r", "r_libs_user", "") or "")

GWAS_VALIDATION = get_gwas_validation_files()          # {ancestry: file_path}
VAL_ANCESTRIES  = sorted(GWAS_VALIDATION.keys())       # e.g. ["AMR", "EAS", "EUR"]
VAL_FILES       = [GWAS_VALIDATION[a] for a in VAL_ANCESTRIES]


rule ancestry_plink_gwas_per_ancestry:
    input:
        dosage_done=DOSAGE_PLINK_DONE,
        pgen=expand(os.path.join(DOSAGE_PLINK_DIR, "{{ancestry}}", "chr{chr}.dosage.pgen"), chr=CHROMOSOMES),
        pvar=expand(os.path.join(DOSAGE_PLINK_DIR, "{{ancestry}}", "chr{chr}.dosage.pvar"), chr=CHROMOSOMES),
        psam=expand(os.path.join(DOSAGE_PLINK_DIR, "{{ancestry}}", "chr{chr}.dosage.psam"), chr=CHROMOSOMES),
        phenotype=PHENOTYPE_FILE,
        covariates=COVARIATES_FILE,
    output:
        merged_pgen=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pgen"),
        merged_pvar=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pvar"),
        merged_psam=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.psam"),
        gwas_tsv=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.tsv"),
        manhattan_png=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.manhattan.png"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_{ancestry}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/run_ancestry_plink_gwas.R",
        r_module=R_MODULE,
        plink_module=PLINK_MODULE,
        chromosomes=",".join(CHROMOSOMES),
        input_dir=lambda wildcards: os.path.join(DOSAGE_PLINK_DIR, wildcards.ancestry),
        output_dir=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry),
        merged_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"{wildcards.ancestry}.merged"),
        gwas_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"{wildcards.ancestry}.gwas"),
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
    shell:
        """
        mkdir -p {params.output_dir} {LOG_DIR}
        bash --login -c '
            if ! type module >/dev/null 2>&1; then
                if [[ -f /etc/profile.d/modules.sh ]]; then
                    source /etc/profile.d/modules.sh
                elif [[ -f /usr/share/Modules/init/bash ]]; then
                    source /usr/share/Modules/init/bash
                fi
            fi

            if ! type module >/dev/null 2>&1; then
                echo "Environment module system unavailable; cannot load R/plink2 modules" >&2
                exit 1
            fi

            module load {params.r_module} {params.plink_module}

            if ! command -v Rscript >/dev/null 2>&1; then
                echo "Rscript not found after loading module: {params.r_module}" >&2
                exit 1
            fi

            if ! command -v plink2 >/dev/null 2>&1; then
                echo "plink2 not found after loading module: {params.plink_module}" >&2
                exit 1
            fi

            Rscript "{params.script}" \
                --ancestry "{wildcards.ancestry}" \
                --chromosomes "{params.chromosomes}" \
                --input-dir "{params.input_dir}" \
                --output-dir "{params.output_dir}" \
                --merged-prefix "{params.merged_prefix}" \
                --gwas-prefix "{params.gwas_prefix}" \
                --gwas-tsv "{output.gwas_tsv}" \
                --manhattan-png "{output.manhattan_png}" \
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
                --n-pcs "{params.n_pcs}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_done:
    input:
        dosage_done=DOSAGE_PLINK_DONE,
        merged_pgen=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pgen"), ancestry=ANCESTRY_NAMES),
        merged_pvar=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pvar"), ancestry=ANCESTRY_NAMES),
        merged_psam=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.psam"), ancestry=ANCESTRY_NAMES),
        gwas_tsv=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.tsv"), ancestry=ANCESTRY_NAMES),
        manhattan_png=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.manhattan.png"), ancestry=ANCESTRY_NAMES),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Per-ancestry merged PLINK GWAS complete" > {log}
        """


rule ancestry_plink_gwas_add_case_control_freq:
    input:
        merged_pgen=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pgen"),
        merged_pvar=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pvar"),
        merged_psam=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.psam"),
        gwas_tsv=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.tsv"),
        pheno=os.path.join(OUTPUT_DIR, "{ancestry}", "phenotype_for_plink.tsv"),
    output:
        case_afreq=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.case.afreq"),
        ctrl_afreq=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.control.afreq"),
        gwas_with_freq=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.with_case_control_freq.tsv"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_freq_{ancestry}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        plink_module=PLINK_MODULE,
        r_module=R_MODULE,
        merge_script="scripts/merge_gwas_case_control_freq.R",
        case_code=FAM_CASE_CODE,
        control_code=FAM_CONTROL_CODE,
        merged_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"{wildcards.ancestry}.merged"),
        case_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"{wildcards.ancestry}.case"),
        control_prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.ancestry, f"{wildcards.ancestry}.control"),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/{wildcards.ancestry} {LOG_DIR}
        bash --login -c '
            if ! type module >/dev/null 2>&1; then
                if [[ -f /etc/profile.d/modules.sh ]]; then
                    source /etc/profile.d/modules.sh
                elif [[ -f /usr/share/Modules/init/bash ]]; then
                    source /usr/share/Modules/init/bash
                fi
            fi

            if ! type module >/dev/null 2>&1; then
                echo "Environment module system unavailable; cannot load R/plink2 modules" >&2
                exit 1
            fi

            module load {params.plink_module} {params.r_module}

            if ! command -v plink2 >/dev/null 2>&1; then
                echo "plink2 not found after loading module(s): {params.plink_module}" >&2
                exit 1
            fi

            if ! command -v Rscript >/dev/null 2>&1; then
                echo "Rscript not found after loading module(s): {params.r_module}" >&2
                exit 1
            fi

            plink2 \
                --pfile "{params.merged_prefix}" \
                --allow-extra-chr \
                --pheno "{input.pheno}" \
                --pheno-name PHENO \
                --keep-if PHENO == {params.case_code} \
                --freq \
                --out "{params.case_prefix}"

            plink2 \
                --pfile "{params.merged_prefix}" \
                --allow-extra-chr \
                --pheno "{input.pheno}" \
                --pheno-name PHENO \
                --keep-if PHENO == {params.control_code} \
                --freq \
                --out "{params.control_prefix}"

            cp "{params.case_prefix}.afreq" "{output.case_afreq}"
            cp "{params.control_prefix}.afreq" "{output.ctrl_afreq}"

            Rscript "{params.merge_script}" \
                --gwas-tsv "{input.gwas_tsv}" \
                --case-afreq "{output.case_afreq}" \
                --ctrl-afreq "{output.ctrl_afreq}" \
                --output-tsv "{output.gwas_with_freq}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_with_freq_done:
    input:
        gwas_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas.done"),
        enriched=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.with_case_control_freq.tsv"), ancestry=ANCESTRY_NAMES),
        case_afreq=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.case.afreq"), ancestry=ANCESTRY_NAMES),
        ctrl_afreq=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.control.afreq"), ancestry=ANCESTRY_NAMES),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_with_freq.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_with_freq_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Per-ancestry GWAS outputs with case/control allele frequencies complete" > {log}
        """


rule ancestry_plink_gwas_effect_allele_tables_per_ancestry:
    input:
        with_freq=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.with_case_control_freq.tsv"),
        pvar=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pvar"),
    output:
        full=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.ea_standardized.tsv"),
        p5e3=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.ea_standardized.p_lt_5e3.tsv"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_ea_standardized_{ancestry}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        script="scripts/build_effect_allele_gwas_tables.py",
    shell:
        """
        mkdir -p {OUTPUT_DIR}/{wildcards.ancestry} {LOG_DIR}
        python3 "{params.script}" \
            --input-gwas "{input.with_freq}" \
            --input-pvar "{input.pvar}" \
            --output-full "{output.full}" \
            --output-significant "{output.p5e3}" \
            --p-threshold 5e-3 \
            > "{log}" 2>&1
        """


rule ancestry_plink_gwas_effect_allele_tables_done:
    input:
        with_freq_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_with_freq.done"),
        full=expand(
            os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.ea_standardized.tsv"),
            ancestry=ANCESTRY_NAMES,
        ),
        p5e3=expand(
            os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.ea_standardized.p_lt_5e3.tsv"),
            ancestry=ANCESTRY_NAMES,
        ),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_effect_allele_tables.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_effect_allele_tables_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Per-ancestry effect-allele standardized GWAS tables complete" > {log}
        """


rule ancestry_plink_gwas_harmonized_merge:
    input:
        with_freq_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_with_freq.done"),
        summary=SUMMARY_TSV,
        with_freq=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.with_case_control_freq.tsv"), ancestry=ANCESTRY_NAMES),
        pvar=expand(os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.merged.pvar"), ancestry=ANCESTRY_NAMES),
    output:
        merged_tsv=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.tsv"),
        merged_xlsx=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.significant.xlsx"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_harmonized_merge.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        r_module=R_MODULE,
        script="scripts/merge_ancestry_gwas_harmonized.R",
        ancestries=",".join(ANCESTRY_NAMES),
        gwas_map=",".join([
            f"{anc}={os.path.join(OUTPUT_DIR, anc, anc + '.gwas.with_case_control_freq.tsv')}" for anc in ANCESTRY_NAMES
        ]),
        pvar_map=",".join([
            f"{anc}={os.path.join(OUTPUT_DIR, anc, anc + '.merged.pvar')}" for anc in ANCESTRY_NAMES
        ]),
        ref_ancestry=REFERENCE_ANCESTRY,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            if ! type module >/dev/null 2>&1; then
                if [[ -f /etc/profile.d/modules.sh ]]; then
                    source /etc/profile.d/modules.sh
                elif [[ -f /usr/share/Modules/init/bash ]]; then
                    source /usr/share/Modules/init/bash
                fi
            fi

            if ! type module >/dev/null 2>&1; then
                echo "Environment module system unavailable; cannot load R module" >&2
                exit 1
            fi

            module load {params.r_module}

            if ! command -v Rscript >/dev/null 2>&1; then
                echo "Rscript not found after loading module(s): {params.r_module}" >&2
                exit 1
            fi
            Rscript "{params.script}" \
                --ancestries "{params.ancestries}" \
                --gwas-with-freq-map "{params.gwas_map}" \
                --pvar-map "{params.pvar_map}" \
                --summary-tsv "{input.summary}" \
                --reference-ancestry "{params.ref_ancestry}" \
                --output-tsv "{output.merged_tsv}" \
                --output-xlsx "{output.merged_xlsx}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_harmonized_done:
    input:
        with_freq_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_with_freq.done"),
        merged_tsv=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.tsv"),
        merged_xlsx=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.significant.xlsx"),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_harmonized.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_harmonized_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Merged harmonized all-ancestry GWAS report complete" > {log}
        """


rule ancestry_plink_gwas_sig5e8_patchwork_report:
    input:
        harmonized_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_harmonized.done"),
        harmonized_tsv=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.tsv"),
    output:
        loci_plot_dir=directory(os.path.join(OUTPUT_DIR, "loci_plots")),
        loci_excel=os.path.join(OUTPUT_DIR, "all_ancestries.loci_sig5e8.report.xlsx"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_sig5e8_patchwork_report.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        r_module=R_MODULE,
        r_libs_user=R_LIBS_USER,
        script="scripts/plot_significant_loci_patchwork.R",
        ancestries=",".join(ANCESTRY_NAMES),
        genome_build=GENOME_BUILD,
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
                        ml {params.r_module} >/dev/null 2>&1 || true
                        if [[ -n "{params.r_libs_user}" ]]; then
                            export R_LIBS_USER="{params.r_libs_user}"
                        fi
            Rscript "{params.script}" \
                --harmonized-tsv "{input.harmonized_tsv}" \
                --ancestries "{params.ancestries}" \
                --genome-build "{params.genome_build}" \
                --output-plot-dir "{output.loci_plot_dir}" \
                --output-excel "{output.loci_excel}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_sig5e8_patchwork_report_done:
    input:
        harmonized_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_harmonized.done"),
        loci_plot_dir=os.path.join(OUTPUT_DIR, "loci_plots"),
        loci_excel=os.path.join(OUTPUT_DIR, "all_ancestries.loci_sig5e8.report.xlsx"),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_sig5e8_patchwork_report.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_sig5e8_patchwork_report_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Per-locus plots in {output.loci_plot_dir} complete" > {log}
        
        """

rule ancestry_plink_gwas_validation:
    input:
        harmonized_done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_harmonized.done"),
        harmonized_tsv=os.path.join(OUTPUT_DIR, "all_ancestries.harmonized.tsv"),
        val_files=VAL_FILES,
    output:
        val_excel=os.path.join(OUTPUT_DIR, "all_ancestries.gwas_validation.xlsx"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_validation.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 16000),
        time=get_default_resource("time", "00:20:00"),
        cores=get_default_resource("cores", 1),
    params:
        r_module=R_MODULE,
        r_libs_user=R_LIBS_USER,
        script="scripts/gwas_validation.R",
        ancestries=",".join(ANCESTRY_NAMES),
        val_ancestries=",".join(VAL_ANCESTRIES),
        val_files=",".join(VAL_FILES),
    shell:
        """
        mkdir -p {OUTPUT_DIR} {LOG_DIR}
        bash --login -c '
            ml {params.r_module} >/dev/null 2>&1 || true
            if [[ -n "{params.r_libs_user}" ]]; then
                export R_LIBS_USER="{params.r_libs_user}"
            fi
            Rscript "{params.script}" \
                --harmonized-tsv "{input.harmonized_tsv}" \
                --ancestries    "{params.ancestries}" \
                --val-ancestries "{params.val_ancestries}" \
                --val-files     "{params.val_files}" \
                --output-excel  "{output.val_excel}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_validation_done:
    input:
        val_excel=os.path.join(OUTPUT_DIR, "all_ancestries.gwas_validation.xlsx"),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_validation.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_validation_done.log"),
    shell:
        """
        touch {output.done}
        echo "GWAS validation complete: {input.val_excel}" > {log}
        """


rule ancestry_plink_gwas_filter_pvalue_per_ancestry:
    input:
        gwas_tsv=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.tsv"),
    output:
        filtered_tsv=os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.p_lt_5e3.tsv"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_filter_pvalue_{ancestry}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        r_module=R_MODULE,
        script="scripts/filter_gwas_pvalue.R",
        pvalue_threshold="0.005",
    shell:
        """
        mkdir -p {OUTPUT_DIR}/{wildcards.ancestry} {LOG_DIR}
        bash --login -c '
            if ! type module >/dev/null 2>&1; then
                if [[ -f /etc/profile.d/modules.sh ]]; then
                    source /etc/profile.d/modules.sh
                elif [[ -f /usr/share/Modules/init/bash ]]; then
                    source /usr/share/Modules/init/bash
                fi
            fi

            if ! type module >/dev/null 2>&1; then
                echo "Environment module system unavailable; cannot load R module" >&2
                exit 1
            fi

            module load {params.r_module}

            if ! command -v Rscript >/dev/null 2>&1; then
                echo "Rscript not found after loading module(s): {params.r_module}" >&2
                exit 1
            fi

            Rscript "{params.script}" \
                --input "{input.gwas_tsv}" \
                --output "{output.filtered_tsv}" \
                --pvalue-threshold "{params.pvalue_threshold}"
        ' > "{log}" 2>&1
        """


rule ancestry_plink_gwas_filter_pvalue_done:
    input:
        filtered_tsv=expand(
            os.path.join(OUTPUT_DIR, "{ancestry}", "{ancestry}.gwas.p_lt_5e3.tsv"),
            ancestry=ANCESTRY_NAMES
        ),
    output:
        done=os.path.join(OUTPUT_DIR, "ancestry_plink_gwas_filter_pvalue.done"),
    log:
        os.path.join(LOG_DIR, "ancestry_plink_gwas_filter_pvalue_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "GWAS filtering (p < 5e-3) complete for all ancestries" > {log}
        """
