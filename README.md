# Local Ancestry and Admixture Modeling in This Repository

## What this project does
This pipeline builds an ancestry-aware analysis workflow for an admixed cohort:

1. Infer local ancestry along the genome using phased cohort + reference haplotypes.
2. Convert local ancestry calls into ancestry-specific dosage/haplotype matrices.
3. Run ancestry-pair GWAS (Tractor-style inputs and outputs).
4. Build disease-prediction models that combine PCs, global ancestry, and local-ancestry-derived burden features.

The current configuration uses four ancestry labels mapped from 1000G populations:
- AFR (YRI)
- AMR (PEL)
- EUR (CEU)
- EAS (CHB)

Reference ancestry for collinearity handling in prediction is AMR.

## End-to-end local ancestry workflow

### Step 1. Reference sample map creation
- Rule: rules/example_complete.smk
- Purpose: produce per-population sample lists and a combined RFMix sample map.
- Key output: sample map with ancestry label codes used by RFMix.

### Step 2. Reference panel subsetting
- Rule: rules/subset_reference_panel.smk
- Purpose: keep only selected reference samples in chromosome VCFs.

### Step 3 and 3c. Cohort region extraction and allele harmonization
- Rules:
  - rules/prepare_cohort_regions_from_plink.smk
  - rules/harmonize_cohort_alleles.smk
- Purpose:
  - export loci-focused cohort VCFs from PLINK
  - harmonize cohort alleles to reference before phasing

### Step 4. Cohort phasing
- Rule: rules/phase_cohort_for_rfmix.smk
- Purpose: create phased cohort VCFs compatible with RFMix.

### Step 5. Local ancestry inference (RFMix)
- Rule: rules/run_rfmix.smk
- Helper: scripts/run_rfmix_chr.sh
- Inputs per chromosome:
  - phased cohort VCF
  - subset reference VCF
  - genetic map
  - sample map (label -> numeric code)
- Outputs per chromosome:
  - chr*.deconvoluted.msp.tsv
  - chr*.deconvoluted.fb.tsv
  - chr*.deconvoluted.rfmix.Q
  - chr*.deconvoluted.sis.tsv
- Integrity behavior:
  - records successful and skipped chromosomes
  - verifies ancestry label-code consistency
  - requires at least one successful chromosome

### Step 7. Convert local ancestry to ancestry-specific dosage/hapcount
- Rule: rules/extract_tracts.smk
- Script: scripts/extract_tracts.py
- Purpose:
  - combine phased genotypes with MSP ancestry calls
  - emit ancestry-specific matrices per chromosome:
    - chr{chr}.phased.anc{anc}.dosage.txt
    - chr{chr}.phased.anc{anc}.hapcount.txt
- Notes:
  - ancestry codes must be contiguous 0..N-1
  - supports optional ancestry-specific VCF outputs

### Step 8. Pairwise ancestry GWAS
- Rule: rules/tractor_gwas_pairwise.smk
- Purpose:
  - run pairwise ancestry GWAS per chromosome
  - merge per-chromosome results by ancestry pair
  - produce Manhattan plots

### Step 9 and 10. Integrated summary and locus ancestry report
- Rules:
  - rules/merge_maf_gwas_summary.smk
  - rules/report_locus_ancestry.smk
- Purpose:
  - merge MAF + GWAS metrics
  - classify ancestry-specific/shared signal and produce report tables

## How the prediction model is created

### Rule and script
- Rule: rules/validate_admixture_prediction.smk
- Script executed by the rule: scripts/validate_admixture_prediction.R

### Inputs to prediction
- Merged summary table from Step 9.
- RFMix global ancestry proportions from .Q files.
- Tractor ancestry-specific dosage matrices.
- Phenotype from FAM (default case=2, control=1).
- PC covariates (default 5 PCs).
- Optional sex and optional age.

### Feature engineering
Three feature sets are evaluated:

1. pcs_only
- PCs (+ optional covariates such as sex/age)

2. pcs_plus_global
- pcs_only + global ancestry proportions (excluding configured reference ancestry)

3. full_admixture
- pcs_plus_global + local-ancestry burden features:
  - hap_dosage_AFR
  - hap_dosage_AMR
  - hap_dosage_EUR
  - hap_dosage_EAS

SNP selection for burden features:
- ancestry-specific p-value threshold (default 0.01)
- capped SNP count per ancestry (default 50)
- optional top SNP per locus
- per-SNP weight = log(odds_ratio)
- ancestry burden feature = weighted sum of dosage values across selected SNPs

### Logistic model form
For each model variant:

- Linear score:
  eta = b0 + b1*x1 + ... + bk*xk
- Probability:
  p(y=1|x) = 1 / (1 + exp(-eta))
- Predicted label:
  y_hat = 1 if p >= 0.5 else 0

### Standardization and fitting
Within each train/test split:
- median imputation per feature (train median)
- z-score scaling using train mean/sd
- binomial logistic regression

## Validation designs in this repository

### A) Stratified k-fold cross-validation (existing default)
- Controlled by:
  - cv_folds (default 5)
  - random_state (default 42)
- Each sample appears in one held-out fold.

### B) Exactly 10 stratified random 80/20 repeats (alternative design)
- Rule target:
  - validate_admixture_prediction_repeated_8020_done
- Implemented as:
  - validation_design = repeated_stratified_holdout_80_20
  - n_repeats = 10
  - holdout_test_size = 0.2
- Meaning:
  - run exactly 10 independent random stratified splits
  - each repeat trains on 80% and tests on 20%
  - this is not 10-fold CV

## Main prediction outputs
Output directory (default):
- results_hispanic_mod/admixture_validation

Primary files:
- sample_prediction_features.tsv
- selected_admixture_snps.tsv
- cv_predictions.tsv
- prediction_metrics.tsv
- model_coefficients.tsv
- roc_curve.png
- calibration_curve.png

Alternative 10x 80/20 outputs are written to:
- results_hispanic_mod/admixture_validation_repeated_8020

## How to run

From steps.sh:

- Default validation:
  ./submit.sh --snakefile rules/validate_admixture_prediction.smk validate_admixture_prediction_done --force --jobs 1

- Alternative validation (exactly 10 stratified random 80/20 repeats):
  ./submit.sh --snakefile rules/validate_admixture_prediction.smk validate_admixture_prediction_repeated_8020_done --force --jobs 1

## Practical interpretation
- If full_admixture outperforms pcs_only and pcs_plus_global, local ancestry burden features add predictive signal beyond global ancestry and PCs.
- model_coefficients.tsv gives direction and magnitude of feature effects.
- prediction_metrics.tsv summarizes discrimination and calibration-quality proxies (AUC, PR-AUC, Brier, sensitivity, specificity, balanced accuracy).

## Configuration knobs to review
In configs/analysis.yml:
- reference populations and ancestry labels
- reference ancestry used for collinearity handling
- SNP thresholding and max SNPs per ancestry
- cv_folds and random_state
- include_sex/include_age

## Notes
- Successful chromosomes are derived from Step 5 output checks.
- Downstream extraction and GWAS use only successful chromosomes.
- Local ancestry dosage files are the bridge between RFMix calls and prediction features.
