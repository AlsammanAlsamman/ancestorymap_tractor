#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    stop(paste("Missing required argument", flag), call. = FALSE)
  }
  args[idx + 1]
}

pair_name <- get_arg("--pair")
input_tsv <- get_arg("--input-tsv")
output_merged <- get_arg("--output-merged")
output_png <- get_arg("--output-png")
output_pdf <- get_arg("--output-pdf")

input_files <- strsplit(input_tsv, ",", fixed = TRUE)[[1]]
input_files <- input_files[nzchar(input_files)]

if (length(input_files) == 0) {
  stop("No input TSV files provided", call. = FALSE)
}

frames <- lapply(input_files, function(path) {
  read.delim(path, stringsAsFactors = FALSE)
})

df <- do.call(rbind, frames)

if (!all(c("CHROM", "POS", "p_dosage_i", "p_dosage_j") %in% colnames(df))) {
  stop("Input table missing required columns: CHROM, POS, p_dosage_i, p_dosage_j", call. = FALSE)
}

df$p_gwas <- pmin(df$p_dosage_i, df$p_dosage_j, na.rm = TRUE)
df$p_gwas[!is.finite(df$p_gwas) | df$p_gwas <= 0] <- NA

chr_levels <- unique(df$CHROM)
chr_num <- suppressWarnings(as.numeric(chr_levels))
if (all(!is.na(chr_num))) {
  chr_levels <- chr_levels[order(chr_num)]
} else {
  chr_levels <- sort(chr_levels)
}

df$CHROM <- factor(df$CHROM, levels = chr_levels)
df <- df[order(df$CHROM, df$POS), ]

chr_lengths <- tapply(df$POS, df$CHROM, max, na.rm = TRUE)
chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
names(chr_offsets) <- names(chr_lengths)

df$BP_cum <- df$POS + chr_offsets[as.character(df$CHROM)]
axis_df <- aggregate(BP_cum ~ CHROM, data = df, FUN = function(x) mean(range(x, na.rm = TRUE)))

dir.create(dirname(output_merged), recursive = TRUE, showWarnings = FALSE)
write.table(df, file = output_merged, sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(df, aes(x = BP_cum, y = -log10(p_gwas), color = CHROM)) +
  geom_point(size = 0.8, alpha = 0.8, na.rm = TRUE) +
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$BP_cum) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(
    title = paste0("Tractor Pairwise GWAS: ", pair_name),
    x = "Chromosome",
    y = "-log10(p)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
ggsave(filename = output_png, plot = p, width = 12, height = 5, dpi = 300)
ggsave(filename = output_pdf, plot = p, width = 12, height = 5)
