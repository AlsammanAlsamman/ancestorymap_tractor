"""
Common utilities for file paths, naming patterns, and ancestry mappings.
All file naming conventions and patterns are hardcoded here.
Users only provide essentials in analysis.yml
"""

def get_cohort_name(config):
    """Get cohort name from config"""
    return config['cohort']['name']

def get_result_dir(cohort):
    """Auto-generate base results directory"""
    return f"results_{cohort}"

def get_output_path(cohort, step_name):
    """Auto-generate output directory for any pipeline step"""
    return f"{get_result_dir(cohort)}/{step_name}"

##############################################################################
# STEP-SPECIFIC OUTPUT DIRECTORIES
##############################################################################

def dir_reference_panel_prep(cohort):
    return get_output_path(cohort, "reference_panel_prep")

def dir_reference_panel_subset(cohort):
    return get_output_path(cohort, "reference_panel_subset")

def dir_cohort_region_subset(cohort):
    return get_output_path(cohort, "cohort_region_subset")

def dir_cohort_phasing(cohort):
    return get_output_path(cohort, "cohort_phasing")

def dir_rfmix(cohort):
    return get_output_path(cohort, "rfmix")

def dir_rfmix_plots(cohort):
    return get_output_path(cohort, "rfmix_plots")

def dir_rfmix_bin_case_control(cohort):
    return get_output_path(cohort, "rfmix_bin_case_control")

def dir_rfmix_prs_correlation(cohort):
    return get_output_path(cohort, "rfmix_prs_correlation")

def dir_tractor(cohort):
    return get_output_path(cohort, "tractor")

def dir_tractor_gwas(cohort):
    return get_output_path(cohort, "tractor_gwas")

def dir_maf_summary(cohort):
    return get_output_path(cohort, "maf_summary")

def dir_maf_gwas_summary(cohort):
    return get_output_path(cohort, "maf_gwas_summary")

def dir_locus_ancestry_report(cohort):
    return get_output_path(cohort, "locus_ancestry_report")

def dir_admixture_validation(cohort):
    return get_output_path(cohort, "admixture_validation")

def dir_admixture_importance(cohort):
    return get_output_path(cohort, "admixture_importance")

##############################################################################
# FILE PATTERNS (with chromosome/ancestry placeholders)
##############################################################################

def ref_subset_vcf(cohort, chromosome, population_tag):
    """Reference panel VCF per chromosome"""
    return f"{dir_reference_panel_subset(cohort)}/chr{chromosome}.{population_tag}.vcf.gz"

def cohort_subset_vcf(cohort, chromosome):
    """Cohort region subset VCF"""
    return f"{dir_cohort_region_subset(cohort)}/chr{chromosome}.loci.vcf.gz"

def phased_vcf(cohort, chromosome):
    """Phased cohort VCF"""
    return f"{dir_cohort_phasing(cohort)}/chr{chromosome}.phased.vcf.gz"

def rfmix_msp(cohort, chromosome):
    """RFMix MSP file (local ancestry segments)"""
    return f"{dir_rfmix(cohort)}/chr{chromosome}.deconvoluted.msp.tsv"

def rfmix_q_file(cohort, chromosome):
    """RFMix Q file (ancestry dosages)"""
    return f"{dir_rfmix(cohort)}/chr{chromosome}.deconvoluted.rfmix.Q"

def tractor_dosage(cohort, chromosome, ancestry):
    """Tractor ancestry-specific dosage"""
    return f"{dir_tractor(cohort)}/chr{chromosome}.phased.anc{ancestry}.dosage.txt"

def rfmix_prs_per_chr(cohort, chromosome):
    """PRS correlation per chromosome"""
    return f"{dir_rfmix_prs_correlation(cohort)}/chr{chromosome}.rfmix_prs_correlation"

def rfmix_ancestry_tract_plot(cohort, chromosome):
    """RFMix ancestry tract visualization"""
    return f"{dir_rfmix_plots(cohort)}/chr{chromosome}.rfmix_ancestry_tracts"

def rfmix_case_control_bin(cohort, chromosome):
    """Case/control binned ancestry analysis"""
    return f"{dir_rfmix_bin_case_control(cohort)}/chr{chromosome}.rfmix_case_control_bins"

##############################################################################
# ANCESTRY MAPPINGS
##############################################################################

def get_ancestry_code(config, label):
    """Convert ancestry label (AFR) to code (0)"""
    mapping = {v: k for k, v in config['ancestry']['population_to_label'].items()}
    rev_mapping = {v: str(list(config['ancestry']['population_to_label'].values()).index(v)) 
                   for k, v in config['ancestry']['population_to_label'].items()}
    # Build proper code mapping: AFR→0, AMR→1, EUR→2, EAS→3
    codes = {}
    for i, (pop, lbl) in enumerate(config['ancestry']['population_to_label'].items()):
        codes[lbl] = str(i)
    return codes.get(label, label)

def get_ancestry_label(config, code):
    """Convert ancestry code (0) to label (AFR)"""
    populations = list(config['ancestry']['population_to_label'].values())
    return populations[int(code)] if int(code) < len(populations) else code

def get_ancestry_color(config, label):
    """Get color for ancestry label"""
    return config['ancestry']['label_to_color'].get(label, '#000000')

def get_population_tag(config):
    """Get population tag from reference populations (e.g., YRI_PEL_CEU_CHB)"""
    return "_".join(config['reference']['populations'])

##############################################################################
# INPUT DATA PATHS
##############################################################################

def get_plink_prefix(config):
    """Full PLINK prefix path"""
    cohort = config['cohort']['name']
    plink_file = config['cohort']['plink_prefix']
    return f"inputs/{plink_file}"

def get_loci_file(config):
    """Loci file path"""
    return f"input/{config['data']['loci_file']}"

def get_prs_file(config):
    """PRS file path"""
    return f"input/{config['data']['prs_file']}"

def get_covariates_file(config):
    """Covariates file path"""
    return f"inputs/{config['data']['covariates_file']}"

def get_phenotype_file(config):
    """Phenotype file path"""
    plink_prefix = get_plink_prefix(config)
    return f"{plink_prefix}.fam"

##############################################################################
# REFERENCE DATA PATHS (fixed locations)
##############################################################################

def ref_panel_file():
    """1000G reference panel file"""
    return "resources/1000g_phase3/panel/integrated_call_samples_v3.20130502.ALL.panel"

def ref_vcf_dir():
    """1000G VCF directory"""
    return "resources/1000g_phase3/release_20130502"

def genetic_map_file(chromosome):
    """Genetic map for chromosome"""
    return f"resources/geneticmaps37/chr{chromosome}.b37.gmap.gz"

def gene_annotation_file():
    """Gene annotation reference"""
    return "resources/NCBI37.3.gene.loc"

##############################################################################
# SAMPLE MAP FILE
##############################################################################

def rfmix_sample_map(cohort, population_tag):
    """RFMix sample map file"""
    return f"{dir_reference_panel_prep(cohort)}/{population_tag}.sample_map.txt"

def reference_sample_list(cohort, population_tag):
    """Reference panel sample list"""
    return f"{dir_reference_panel_prep(cohort)}/{population_tag}.samples.txt"

##############################################################################
# SUMMARY/OUTPUT FILES
##############################################################################

def maf_gwas_summary_merged(cohort):
    """Merged MAF/GWAS summary file"""
    return f"{dir_maf_gwas_summary(cohort)}/merged.maf_gwas_summary.tsv"

def admixture_features_tsv(cohort):
    """Admixture sample features"""
    return f"{dir_admixture_validation(cohort)}/sample_prediction_features.tsv"

def admixture_selected_snps(cohort):
    """Selected admixture SNPs"""
    return f"{dir_admixture_validation(cohort)}/selected_admixture_snps.tsv"

##############################################################################
# HELPER: Get all chromosomes
##############################################################################

def get_all_chromosomes():
    """Return all chromosome numbers"""
    return [str(i) for i in range(1, 23)]

def get_all_ancestries(config):
    """Return all ancestry codes"""
    n_anc = len(config['reference']['populations'])
    return list(range(n_anc))
