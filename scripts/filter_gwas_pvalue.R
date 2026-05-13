#!/usr/bin/env Rscript

# Filter GWAS results by p-value threshold and sort by chromosome and position
# Usage: Rscript filter_gwas_pvalue.R --input <input.tsv> --output <output.tsv> --pvalue-threshold <pval>

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Parse arguments
option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Path to input GWAS TSV file"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Path to output filtered TSV file"),
  make_option(c("--pvalue-threshold"), dest = "pvalue_threshold", type = "double", default = 0.005,
              help = "P-value threshold for filtering [default: 0.005]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both --input and --output are required")
}

# Set default pvalue threshold if not provided
if (is.null(opt$`pvalue_threshold`) || opt$`pvalue_threshold` == 0) {
  opt$`pvalue_threshold` <- 0.005
}

# Read GWAS results
cat("Reading GWAS results from:", opt$input, "\n")
gwas <- read_tsv(opt$input, show_col_types = FALSE)

# Check for required columns (common GWAS output formats)
pval_col <- NULL
chr_col <- NULL
pos_col <- NULL

# Try common p-value column names
pval_names <- c("P", "pval", "p_val", "p_value", "P.value", "p.value", "PVAL", "P_VAL")
for (name in pval_names) {
  if (name %in% names(gwas)) {
    pval_col <- name
    break
  }
}

# Try common chromosome column names
chr_names <- c("CHR", "chr", "#CHR", "CHROM", "Chromosome")
for (name in chr_names) {
  if (name %in% names(gwas)) {
    chr_col <- name
    break
  }
}

# Try common position column names
pos_names <- c("POS", "pos", "BP", "bp", "Position")
for (name in pos_names) {
  if (name %in% names(gwas)) {
    pos_col <- name
    break
  }
}

if (is.null(pval_col)) {
  stop("Could not find p-value column in input file. Tried: ", paste(pval_names, collapse = ", "))
}

if (is.null(chr_col)) {
  warning("Could not find chromosome column. Will not sort by chromosome.")
}

if (is.null(pos_col)) {
  warning("Could not find position column. Will not sort by position.")
}

# Ensure pvalue_threshold is numeric
pval_threshold <- as.numeric(opt$`pvalue_threshold`)
cat("P-value threshold:", pval_threshold, "\n")

# Filter by p-value
cat("Filtering for p-value <", pval_threshold, "\n")
cat("Total SNPs before filtering:", nrow(gwas), "\n")

gwas_filtered <- gwas %>%
  filter(!!sym(pval_col) < pval_threshold)

cat("Total SNPs after filtering:", nrow(gwas_filtered), "\n")

# Sort by chromosome and position if available
if (!is.null(chr_col) && !is.null(pos_col)) {
  cat("Sorting by chromosome and position\n")
  gwas_filtered <- gwas_filtered %>%
    mutate(!!sym(chr_col) := as.numeric(gsub("chr", "", !!sym(chr_col)))) %>%
    arrange(!!sym(chr_col), !!sym(pos_col))
} else if (!is.null(chr_col)) {
  cat("Sorting by chromosome only\n")
  gwas_filtered <- gwas_filtered %>%
    mutate(!!sym(chr_col) := as.numeric(gsub("chr", "", !!sym(chr_col)))) %>%
    arrange(!!sym(chr_col))
}

# Write output
cat("Writing filtered results to:", opt$output, "\n")
write_tsv(gwas_filtered, opt$output)

cat("Done!\n")
