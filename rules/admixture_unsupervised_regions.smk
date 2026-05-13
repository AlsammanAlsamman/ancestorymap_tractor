import os
import sys

sys.path.append("utils")
from bioconfigme import (
    get_analysis_value,
    get_default_resource,
    get_loci_chromosomes,
    get_results_dir,
    get_software_module,
)

RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "admixture_unsupervised")

STEP1_DONE = os.path.join(RESULTS_DIR, "reference_panel_prep", "build_reference_sample_lists.done")
STEP2_DONE = os.path.join(RESULTS_DIR, "reference_panel_subset", "subset_reference_panel.done")
STEP3C_DONE = os.path.join(get_analysis_value("step3c.output_dir"), "harmonize_cohort_alleles.done")

LOCI_FILE = get_analysis_value("step3.loci_file")
CHROMOSOMES = [str(c) for c in get_loci_chromosomes(LOCI_FILE)]

COHORT_VCF_PATTERN = os.path.join(get_analysis_value("step3c.output_dir"), "chr{chr}." + get_analysis_value("step3c.output_suffix") + ".vcf.gz")
REF_VCF_PATTERN = os.path.join(get_analysis_value("step2.output_dir"), "chr{chr}." + get_analysis_value("step2.output_suffix") + ".vcf.gz")

REF_POPS = [str(p) for p in get_analysis_value("reference.populations")]
REF_POP_TAG = "_".join(REF_POPS)
REFERENCE_SAMPLE_MAP = os.path.join(get_analysis_value("first_step.output_dir"), f"{REF_POP_TAG}.sample_map.txt")

MIN_K = int(get_analysis_value("admixture_unsupervised.min_k"))
MAX_K = int(get_analysis_value("admixture_unsupervised.max_k"))
CV_FOLDS = int(get_analysis_value("admixture_unsupervised.cv_folds"))
COHORT_LABEL = str(get_analysis_value("admixture_unsupervised.cohort_label"))
COHORT_COLOR = str(get_analysis_value("admixture_unsupervised.cohort_color"))

if MIN_K < 1 or MAX_K < MIN_K:
    raise ValueError("admixture_unsupervised min_k/max_k are invalid")

K_VALUES = [str(k) for k in range(MIN_K, MAX_K + 1)]

LABEL_TO_COLOR = get_analysis_value("ancestry.label_to_color")
GROUP_COLORS = ",".join([f"{k}={v}" for k, v in LABEL_TO_COLOR.items()])

PLINK_MODULE = get_software_module("plink2")
BCFTOOLS_MODULE = get_software_module("bcftools")
ADMIXTURE_MODULE = get_software_module("admixture")
R_MODULE = get_software_module("r")
PYTHON_MODULE = get_software_module("python")

rule admixture_regions_merge_chr:
    input:
        step1_done=STEP1_DONE,
        step2_done=STEP2_DONE,
        step3c_done=STEP3C_DONE,
        cohort_vcf=COHORT_VCF_PATTERN,
        ref_vcf=REF_VCF_PATTERN,
    output:
        bed=os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.bed"),
        bim=os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.bim"),
        fam=os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.fam"),
        done=os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_merge_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 128000),
        time=get_default_resource("time", "24:00:00"),
        cores=get_default_resource("cores", 20),
    params:
        plink_module=PLINK_MODULE,
        bcftools_module=BCFTOOLS_MODULE,
        out_prefix=lambda wc: os.path.join(OUTPUT_DIR, "per_chr", f"chr{wc.chr}.merged"),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/per_chr {LOG_DIR}
        bash --login -c '
            set -euo pipefail
            module load {params.bcftools_module} {params.plink_module}
            tmp_vcf="{OUTPUT_DIR}/per_chr/.tmp.chr{wildcards.chr}.cohort_ref.merged.vcf.gz"
            tmp_merged_bcf="{OUTPUT_DIR}/per_chr/.tmp.chr{wildcards.chr}.cohort_ref.merged.bcf"
            tmp_norm_bcf="{OUTPUT_DIR}/per_chr/.tmp.chr{wildcards.chr}.cohort_ref.merged.norm.bcf"
            rm -f "$tmp_vcf" "$tmp_vcf.tbi" "$tmp_merged_bcf" "$tmp_norm_bcf"

            # Use explicit intermediate files instead of a long stream pipeline to avoid
            # silent truncation/missing outputs on cluster nodes.
            bcftools merge --force-samples -Ob -o "$tmp_merged_bcf" "{input.cohort_vcf}" "{input.ref_vcf}"
            test -s "$tmp_merged_bcf"

            bcftools norm -m -any -Ob -o "$tmp_norm_bcf" "$tmp_merged_bcf"
            test -s "$tmp_norm_bcf"

            bcftools view -m2 -M2 -v snps -Oz -o "$tmp_vcf" "$tmp_norm_bcf"
            test -s "$tmp_vcf"

            # Validate compressed output before PLINK conversion.
            bcftools view -h "$tmp_vcf" >/dev/null
            bcftools index -f -t "$tmp_vcf"

            plink2 --vcf "$tmp_vcf" \
                   --allow-extra-chr \
                   --double-id \
                   --make-bed \
                   --threads {resources.cores} \
                   --memory {resources.mem_mb} \
                   --out "{params.out_prefix}"

            rm -f "$tmp_vcf" "$tmp_vcf.tbi" "$tmp_merged_bcf" "$tmp_norm_bcf"
            touch "{output.done}"
        ' > "{log}" 2>&1
        """


rule admixture_regions_merge_all_chr:
    input:
        beds=expand(os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.bed"), chr=CHROMOSOMES),
        bims=expand(os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.bim"), chr=CHROMOSOMES),
        fams=expand(os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.fam"), chr=CHROMOSOMES),
    output:
        bed=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bed"),
        bim=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bim"),
        fam=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.fam"),
        done=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_merge_all_chr.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 128000),
        time=get_default_resource("time", "24:00:00"),
        cores=get_default_resource("cores", 20),
    params:
        plink_module=PLINK_MODULE,
        first_prefix=os.path.join(OUTPUT_DIR, "per_chr", f"chr{CHROMOSOMES[0]}.merged"),
        merge_list=os.path.join(OUTPUT_DIR, "merged", "merge_list.txt"),
        out_prefix=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g"),
        other_prefixes="\n".join([os.path.join(OUTPUT_DIR, "per_chr", f"chr{c}.merged") for c in CHROMOSOMES[1:]]),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/merged {LOG_DIR}
        cat > "{params.merge_list}" << 'EOF'
{params.other_prefixes}
EOF

        bash --login -c '
            module load {params.plink_module}
            if [[ -s "{params.merge_list}" ]]; then
                plink2 --bfile "{params.first_prefix}" \
                       --merge-list "{params.merge_list}" \
                       --allow-extra-chr \
                       --make-bed \
                       --threads {resources.cores} \
                       --memory {resources.mem_mb} \
                       --out "{params.out_prefix}"
            else
                plink2 --bfile "{params.first_prefix}" \
                       --allow-extra-chr \
                       --make-bed \
                       --threads {resources.cores} \
                       --memory {resources.mem_mb} \
                       --out "{params.out_prefix}"
            fi

            touch "{output.done}"
        ' > "{log}" 2>&1
        """


rule admixture_regions_sample_groups:
    input:
        fam=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.fam"),
        sample_map=REFERENCE_SAMPLE_MAP,
    output:
        groups_tsv=os.path.join(OUTPUT_DIR, "merged", "sample_groups.tsv"),
        done=os.path.join(OUTPUT_DIR, "merged", "sample_groups.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_sample_groups.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        py_module=PYTHON_MODULE,
        script="scripts/build_admixture_sample_groups.py",
        cohort_label=COHORT_LABEL,
    shell:
        """
        mkdir -p {OUTPUT_DIR}/merged {LOG_DIR}
        bash --login -c '
            module load {params.py_module}
            python3 "{params.script}" \
                --fam "{input.fam}" \
                --reference-sample-map "{input.sample_map}" \
                --cohort-label "{params.cohort_label}" \
                --output "{output.groups_tsv}"

            touch "{output.done}"
        ' > "{log}" 2>&1
        """


rule admixture_regions_run_k:
    input:
        bed=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bed"),
        bim=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bim"),
        fam=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.fam"),
    output:
        q=os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.Q"),
        p=os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.P"),
        log=os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.log"),
        done=os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.done"),
    resources:
        mem_mb=get_default_resource("mem_mb", 128000),
        time=get_default_resource("time", "24:00:00"),
        cores=get_default_resource("cores", 20),
    params:
        admixture_module=ADMIXTURE_MODULE,
        run_dir=os.path.join(OUTPUT_DIR, "admixture_k"),
        cv_folds=CV_FOLDS,
    shell:
        """
        mkdir -p {params.run_dir}
        bash --login -c '
            module load {params.admixture_module}
            cd "{params.run_dir}"
            ln -sf "{input.bed}" data.bed
            ln -sf "{input.bim}" data.bim
            ln -sf "{input.fam}" data.fam
            admixture --cv={params.cv_folds} -j{resources.cores} data.bed {wildcards.k} > "k{wildcards.k}.log" 2>&1
            cp -f data.{wildcards.k}.Q "k{wildcards.k}.Q"
            cp -f data.{wildcards.k}.P "k{wildcards.k}.P"
            touch "k{wildcards.k}.done"
        '
        """


rule admixture_regions_select_best_k:
    input:
        logs=expand(os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.log"), k=K_VALUES),
    output:
        cv_summary=os.path.join(OUTPUT_DIR, "admixture_k", "cv_summary.tsv"),
        best_k=os.path.join(OUTPUT_DIR, "admixture_k", "best_k.txt"),
        done=os.path.join(OUTPUT_DIR, "admixture_k", "select_best_k.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_select_best_k.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        py_module=PYTHON_MODULE,
        script="scripts/select_best_admixture_k.py",
        k_logs=",".join([f"{k}={os.path.join(OUTPUT_DIR, 'admixture_k', f'k{k}.log')}" for k in K_VALUES]),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/admixture_k {LOG_DIR}
        bash --login -c '
            module load {params.py_module}
            python3 "{params.script}" \
                --k-logs "{params.k_logs}" \
                --output-summary "{output.cv_summary}" \
                --output-best-k "{output.best_k}"

            touch "{output.done}"
        ' > "{log}" 2>&1
        """


rule admixture_regions_plot_best_k:
    input:
        best_k=os.path.join(OUTPUT_DIR, "admixture_k", "best_k.txt"),
        fam=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.fam"),
        groups_tsv=os.path.join(OUTPUT_DIR, "merged", "sample_groups.tsv"),
        cv_summary=os.path.join(OUTPUT_DIR, "admixture_k", "cv_summary.tsv"),
        q_files=expand(os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.Q"), k=K_VALUES),
    output:
        plot_png=os.path.join(OUTPUT_DIR, "plots", "admixture_best_k.png"),
        done=os.path.join(OUTPUT_DIR, "plots", "admixture_best_k.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_plot_best_k.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        r_module=R_MODULE,
        script="scripts/plot_admixture_best_k.R",
        q_dir=os.path.join(OUTPUT_DIR, "admixture_k"),
        cohort_label=COHORT_LABEL,
        group_colors=GROUP_COLORS,
        cohort_color=COHORT_COLOR,
    shell:
        """
        mkdir -p {OUTPUT_DIR}/plots {LOG_DIR}
        bash --login -c '
            module load {params.r_module}
            Rscript "{params.script}" \
                --best-k-file "{input.best_k}" \
                --fam-file "{input.fam}" \
                --sample-groups "{input.groups_tsv}" \
                --q-dir "{params.q_dir}" \
                --cohort-label "{params.cohort_label}" \
                --group-colors "{params.group_colors}" \
                --cohort-color "{params.cohort_color}" \
                --output-png "{output.plot_png}"

            touch "{output.done}"
        ' > "{log}" 2>&1
        """


rule admixture_regions_unsupervised_done:
    input:
        step1_done=STEP1_DONE,
        step2_done=STEP2_DONE,
        step3c_done=STEP3C_DONE,
        per_chr_done=expand(os.path.join(OUTPUT_DIR, "per_chr", "chr{chr}.merged.done"), chr=CHROMOSOMES),
        merged_done=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.done"),
        sample_groups_done=os.path.join(OUTPUT_DIR, "merged", "sample_groups.done"),
        k_done=expand(os.path.join(OUTPUT_DIR, "admixture_k", "k{k}.done"), k=K_VALUES),
        select_done=os.path.join(OUTPUT_DIR, "admixture_k", "select_best_k.done"),
        plot_done=os.path.join(OUTPUT_DIR, "plots", "admixture_best_k.done"),
        merged_bed=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bed"),
        merged_bim=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.bim"),
        merged_fam=os.path.join(OUTPUT_DIR, "merged", "all_regions_cohort_1000g.fam"),
        sample_groups=os.path.join(OUTPUT_DIR, "merged", "sample_groups.tsv"),
        cv_summary=os.path.join(OUTPUT_DIR, "admixture_k", "cv_summary.tsv"),
        best_k=os.path.join(OUTPUT_DIR, "admixture_k", "best_k.txt"),
        plot_png=os.path.join(OUTPUT_DIR, "plots", "admixture_best_k.png"),
    output:
        done=os.path.join(OUTPUT_DIR, "admixture_regions_unsupervised.done"),
    log:
        os.path.join(LOG_DIR, "admixture_regions_unsupervised_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch "{output.done}"
        echo "Unsupervised ADMIXTURE (regions, cohort+1000G) complete" > "{log}"
        """
