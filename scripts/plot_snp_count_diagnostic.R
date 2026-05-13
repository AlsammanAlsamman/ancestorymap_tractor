#!/usr/bin/env Rscript
# plot_snp_count_diagnostic.R
# Produces stacked-bar + faceted SNP-count plots and a sample-count bar plot,
# comparing per-ancestry PLINK output (ancestry_plink) to the original Hispanic data.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
output_dir  <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir") {
    output_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(input_file) || is.null(output_dir)) {
  stop("Usage: Rscript plot_snp_count_diagnostic.R --input <tsv> --output-dir <dir>")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------- #
# Load data and pivot wide → long
# --------------------------------------------------------------------------- #
wide <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Order chromosomes numerically
chr_levels <- as.character(sort(unique(as.integer(wide$chr))))
wide$chr <- factor(wide$chr, levels = chr_levels)

# Detect ancestry names from column names (everything before "_snps" except "original")
snp_cols    <- grep("_snps$",    names(wide), value = TRUE)
sample_cols <- grep("_samples$", names(wide), value = TRUE)
sources <- sub("_snps$", "", snp_cols)   # e.g. "original", "AFR", "AMR", "EUR", "EAS"
anc_names <- setdiff(sources, "original")

# Pivot SNP counts to long
snp_long <- wide %>%
  select(chr, all_of(snp_cols)) %>%
  pivot_longer(cols = all_of(snp_cols), names_to = "source_col", values_to = "snp_count") %>%
  mutate(ancestry = sub("_snps$", "", source_col),
         source   = ifelse(ancestry == "original", "original", "ancestry_plink")) %>%
  select(source, ancestry, chr, snp_count)

# Pivot sample counts to long
sample_long <- wide %>%
  select(chr, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "source_col", values_to = "sample_count") %>%
  mutate(ancestry = sub("_samples$", "", source_col)) %>%
  select(ancestry, chr, sample_count)

df <- left_join(snp_long, sample_long, by = c("ancestry", "chr"))

orig_df <- df[df$source == "original",  ]
anc_df  <- df[df$source == "ancestry_plink", ]

# Ancestry colour palette (plus Original in grey)
pal_base   <- c(AFR = "#1f77b4", AMR = "#ff7f0e", EUR = "#2ca02c", EAS = "#d62728")
palette    <- c(pal_base[names(pal_base) %in% anc_names], Original = "#888888")

# --------------------------------------------------------------------------- #
# Plot 1: Stacked bar (all ancestries per chromosome) + original as dashed line
# --------------------------------------------------------------------------- #
p1 <- ggplot() +
  geom_col(
    data  = anc_df,
    aes(x = chr, y = snp_count, fill = ancestry)
  ) +
  geom_line(
    data  = orig_df,
    aes(x = as.integer(chr), y = snp_count, colour = "Original"),
    linetype  = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    data   = orig_df,
    aes(x = as.integer(chr), y = snp_count, colour = "Original"),
    shape  = 18,
    size   = 2.5
  ) +
  scale_fill_manual(values = palette, name = "Ancestry") +
  scale_colour_manual(values = c(Original = "black"), name = NULL) +
  labs(
    title    = "SNP Counts per Chromosome by Ancestry vs Original",
    subtitle = "Stacked bars = ancestry-specific PLINK dosage; dashed line = original Hispanic",
    x        = "Chromosome",
    y        = "SNP Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )

ggsave(file.path(output_dir, "snp_count_stacked_plot.pdf"), p1, width = 16, height = 6)
ggsave(file.path(output_dir, "snp_count_stacked_plot.png"), p1, width = 16, height = 6, dpi = 150)
cat("Saved snp_count_stacked_plot\n")

# --------------------------------------------------------------------------- #
# Plot 2: Faceted by ancestry – bar per chromosome vs original dashed line
# --------------------------------------------------------------------------- #
p2 <- ggplot() +
  geom_col(
    data  = anc_df,
    aes(x = chr, y = snp_count, fill = ancestry)
  ) +
  geom_line(
    data      = orig_df,
    aes(x = as.integer(chr), y = snp_count, colour = "Original"),
    linetype  = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    data   = orig_df,
    aes(x = as.integer(chr), y = snp_count, colour = "Original"),
    shape  = 18,
    size   = 1.8
  ) +
  scale_fill_manual(values  = palette, name = "Ancestry") +
  scale_colour_manual(values = c(Original = "black"), name = NULL) +
  facet_wrap(~ ancestry, ncol = 2, scales = "free_y") +
  labs(
    title    = "Per-Ancestry SNP Counts vs Original (per chromosome)",
    subtitle = "Bars = ancestry-specific PLINK dosage; dashed line = original Hispanic",
    x        = "Chromosome",
    y        = "SNP Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "snp_count_facet_plot.pdf"), p2, width = 14, height = 10)
ggsave(file.path(output_dir, "snp_count_facet_plot.png"), p2, width = 14, height = 10, dpi = 150)
cat("Saved snp_count_facet_plot\n")

# --------------------------------------------------------------------------- #
# Plot 3: Sample counts per ancestry vs original
# (sample count is constant across chromosomes; take the value from chr 1 or mean)
# --------------------------------------------------------------------------- #
sample_anc <- anc_df %>%
  group_by(ancestry) %>%
  summarise(sample_count = round(mean(sample_count, na.rm = TRUE)), .groups = "drop") %>%
  mutate(source = "Ancestry PLINK")

sample_orig <- data.frame(
  ancestry     = "Original",
  sample_count = orig_df$sample_count[1],
  source       = "Original",
  stringsAsFactors = FALSE
)

sample_combined <- rbind(
  as.data.frame(sample_anc),
  sample_orig
)

# Order bars: ancestries alphabetically, then Original last
anc_order <- c(sort(setdiff(sample_combined$ancestry, "Original")), "Original")
sample_combined$ancestry <- factor(sample_combined$ancestry, levels = anc_order)

pal_samples <- c(pal_base[names(pal_base) %in% anc_names], Original = "#888888")

p3 <- ggplot(sample_combined, aes(x = ancestry, y = sample_count, fill = ancestry)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sample_count), vjust = -0.4, size = 4) +
  scale_fill_manual(values = pal_samples, guide = "none") +
  labs(
    title    = "Sample Counts by Ancestry vs Original Hispanic Data",
    subtitle = "Per-ancestry values are the mean across chromosomes",
    x        = "Ancestry / Dataset",
    y        = "Number of Samples"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(output_dir, "sample_count_plot.pdf"), p3, width = 8, height = 5)
ggsave(file.path(output_dir, "sample_count_plot.png"), p3, width = 8, height = 5, dpi = 150)
cat("Saved sample_count_plot\n")

# --------------------------------------------------------------------------- #
# Plot 4: Fraction of original SNPs retained per ancestry per chromosome
# --------------------------------------------------------------------------- #
orig_lookup <- setNames(orig_df$snp_count, as.character(orig_df$chr))
anc_df$orig_snps   <- orig_lookup[as.character(anc_df$chr)]
anc_df$pct_retained <- ifelse(
  anc_df$orig_snps > 0,
  100 * anc_df$snp_count / anc_df$orig_snps,
  NA_real_
)

p4 <- ggplot(anc_df, aes(x = chr, y = pct_retained, colour = ancestry, group = ancestry)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = palette, name = "Ancestry") +
  scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Fraction of Original SNPs Retained After Ancestry-specific PLINK Export",
    subtitle = "100% = all original SNPs present; lower values indicate SNP loss",
    x        = "Chromosome",
    y        = "% SNPs Retained"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "snp_retention_plot.pdf"), p4, width = 16, height = 5)
ggsave(file.path(output_dir, "snp_retention_plot.png"), p4, width = 16, height = 5, dpi = 150)
cat("Saved snp_retention_plot\n")

cat("All plots written to", output_dir, "\n")
