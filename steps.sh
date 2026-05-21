#!/bin/bash

# Step 1: Build configured reference sample lists and sample map
./submit.sh --snakefile rules/example_complete.smk build_reference_sample_lists

# Step 2: Subset 1000G reference VCFs to configured populations
./submit.sh --snakefile rules/subset_reference_panel.smk subset_reference_panel_done --jobs 22

# Step 3: Extract loci-defined regions from cohort PLINK and export chromosome-based VCFs
./submit.sh --snakefile rules/prepare_cohort_regions_from_plink.smk prepare_cohort_regions_from_plink_done --jobs 22

# Step 3c: Harmonize cohort VCF alleles to the subset reference panel before phasing
./submit.sh --snakefile rules/harmonize_cohort_alleles.smk harmonize_cohort_alleles_done --jobs 22

# Step 3d: Unsupervised ADMIXTURE on extracted regions (cohort + 1000G), K=min_k..max_k from analysis.yml
./submit.sh --snakefile rules/admixture_unsupervised_regions.smk admixture_regions_unsupervised_done --jobs 22

# Step 3b: Calculate cohort/reference-panel MAF tables and build the merged annotated SNP table
./submit.sh --snakefile rules/calculate_maf.smk calculate_maf_done --jobs 40

# Step 4: Phase cohort VCFs per chromosome for RFMix input
./submit.sh --snakefile rules/phase_cohort_for_rfmix.smk phase_cohort_for_rfmix_done --jobs 22

# Step 5: Run RFMix local ancestry per chromosome
./submit.sh --snakefile rules/run_rfmix.smk run_rfmix_done --jobs 22

# Step 5b: Merge chromosome-level RFMix Q/MSP outputs, drop comment lines, add chr, and relabel ancestry codes
./submit.sh --snakefile rules/merge_rfmix_outputs.smk merge_rfmix_outputs_done --jobs 1

# # Step 6: Plot chromosome-level ancestry tracts from RFMix MSP outputs
# ./submit.sh --snakefile rules/plot_rfmix_tracts.smk results/rfmix_plots/plot_rfmix_tracts.done --jobs 9

# Step 7: Extract ancestry-specific tracts and dosage from phased VCF + RFMix MSP
./submit.sh --snakefile rules/extract_tracts.smk extract_tracts_done --jobs 22

# Step 7a-collapsed: Collapse 4-way ancestry dosage/hapcount into two configured groups (collapsed_tractor.group1/group2)
./submit.sh --snakefile rules/collapse_ancestry_dosage_counts.smk collapse_ancestry_dosage_counts_done --jobs 22

# Step 8-collapsed: Run Tractor GWAS on collapsed two-group design and merge genome-wide outputs
./submit.sh --snakefile rules/tractor_gwas_collapsed.smk tractor_gwas_collapsed_done --jobs 22











# Step 7b: Export per-ancestry dosage matrices to separate PLINK pgen folders (AMR/EUR/EAS/AFR)
./submit.sh --snakefile rules/export_ancestry_dosage_plink.smk ancestry_dosage_plink_done --subjob-time 72:00:00 --jobs 22

# Step 7b-diag: Count SNPs and samples per ancestry per chromosome vs original Hispanic PLINK; produce diagnostic plots
./submit.sh --snakefile rules/snp_count_diagnostic.smk snp_count_diagnostic_done --jobs 1

# Step 7c: Merge ancestry PLINKs across chromosomes, run PLINK GWAS with phenotype+PCs, and plot Manhattan per ancestry
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_done --jobs 22

# Step 7c-filter: Filter GWAS results to retain only SNPs with p-value < 5e-3, sorted by chromosome and position
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_filter_pvalue_done --jobs 22

# Step 7d: Add case/control ref+alt allele frequencies to each ancestry GWAS output table
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_with_freq_done --jobs 22

# Step 7d-micro: Build per-locus tract/allele microscope plots and summary cards for validation_regions
./submit.sh --snakefile rules/locus_haplotype_microscope.smk locus_microscope_done --subjob-time 08:00:00 --jobs 22

# Step 7d2: Build per-ancestry EA/NEA standardized GWAS tables (full + P<5e-3), with EA case/control/total frequencies
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_effect_allele_tables_done  --jobs 22



# Step 8: Pairwise Tractor GWAS per chromosome, merge each pair, and plot Manhattan
./submit.sh --snakefile rules/tractor_gwas_pairwise.smk tractor_gwas_pairwise_done --jobs 22

# Step 9: Merge harmonized MAF with merged Tractor GWAS tables into one summary table
./submit.sh --snakefile rules/merge_maf_gwas_summary.smk merge_maf_gwas_summary_done --jobs 1

# Step 7e: Harmonize and merge all ancestry GWAS outputs to common EA/NEA, annotate loci/genes, and build Excel summaries
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_harmonized_done --jobs 1

# Step 7f: Plot only loci with P<5e-8 using patchwork (all ancestries together) and build formatted significant SNP Excel report
./submit.sh --snakefile rules/ancestry_plink_gwas.smk ancestry_plink_gwas_sig5e8_patchwork_report_done --jobs 1

# Step 10: Classify loci as ancestry-specific or shared and build the report outputs
./submit.sh --snakefile rules/report_locus_ancestry.smk build_locus_ancestry_report --jobs 1

# Step 11: Validate admixture-informed disease prediction using Tractor dosage, RFMix ancestry, and PCs
./submit.sh --snakefile rules/validate_admixture_prediction.smk validate_admixture_prediction_done --jobs 1

# Step 11b: Alternative validation design: exactly 10 stratified random 80/20 repeats (not k-fold CV)
./submit.sh --snakefile rules/validate_admixture_prediction.smk validate_admixture_prediction_repeated_8020_done --jobs 1

# Step 12: Quantify locus-level and ancestry-level dosage importance in the prediction model
./submit.sh --snakefile rules/report_admixture_importance.smk report_admixture_importance_done --jobs 1


