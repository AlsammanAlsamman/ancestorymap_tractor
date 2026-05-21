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

get_opt_arg <- function(flag, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default)
  }
  if (idx == length(args)) {
    stop(paste("Missing value for optional argument", flag), call. = FALSE)
  }
  args[idx + 1]
}

pair_name <- get_arg("--pair")
input_tsv <- get_arg("--input-tsv")
output_merged <- get_arg("--output-merged")
output_png <- get_opt_arg("--output-png")
output_pdf <- get_opt_arg("--output-pdf")
output_qq_png <- get_opt_arg("--output-qq-png")
output_qq_pdf <- get_opt_arg("--output-qq-pdf")
output_dosage_i_manhattan_png <- get_opt_arg("--output-dosage-i-manhattan-png")
output_dosage_i_manhattan_pdf <- get_opt_arg("--output-dosage-i-manhattan-pdf")
output_dosage_i_qq_png <- get_opt_arg("--output-dosage-i-qq-png")
output_dosage_i_qq_pdf <- get_opt_arg("--output-dosage-i-qq-pdf")
output_dosage_j_manhattan_png <- get_opt_arg("--output-dosage-j-manhattan-png")
output_dosage_j_manhattan_pdf <- get_opt_arg("--output-dosage-j-manhattan-pdf")
output_dosage_j_qq_png <- get_opt_arg("--output-dosage-j-qq-png")
output_dosage_j_qq_pdf <- get_opt_arg("--output-dosage-j-qq-pdf")
output_local_manhattan_png <- get_opt_arg("--output-local-ancestry-manhattan-png")
output_local_manhattan_pdf <- get_opt_arg("--output-local-ancestry-manhattan-pdf")
output_local_qq_png <- get_opt_arg("--output-local-ancestry-qq-png")
output_local_qq_pdf <- get_opt_arg("--output-local-ancestry-qq-pdf")

input_files <- strsplit(input_tsv, ",", fixed = TRUE)[[1]]
input_files <- input_files[nzchar(input_files)]

if (length(input_files) == 0) {
  stop("No input TSV files provided", call. = FALSE)
}

frames <- lapply(input_files, function(path) {
  read.delim(path, stringsAsFactors = FALSE)
})

df <- do.call(rbind, frames)

if (!all(c("CHROM", "POS", "p_dosage_i", "p_dosage_j", "p_local_ancestry") %in% colnames(df))) {
  stop("Input table missing required columns: CHROM, POS, p_dosage_i, p_dosage_j, p_local_ancestry", call. = FALSE)
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

SIGNIF_P_MAIN <- 5e-8
SIGNIF_P_SUGGESTIVE <- 5e-5
PLOT_MAX_POINTS <- 1500000
QQ_MAX_POINTS <- 2000000

downsample_for_manhattan <- function(data, p_col, max_points = PLOT_MAX_POINTS, keep_p = SIGNIF_P_SUGGESTIVE) {
  pvals <- data[[p_col]]
  keep_sig <- is.finite(pvals) & pvals <= keep_p
  sig_idx <- which(keep_sig)
  nonsig_idx <- which(!keep_sig)

  if (length(sig_idx) >= max_points) {
    return(data[sig_idx, , drop = FALSE])
  }

  remaining <- max_points - length(sig_idx)
  if (remaining <= 0 || length(nonsig_idx) == 0) {
    return(data[c(sig_idx, nonsig_idx), , drop = FALSE])
  }

  if (length(nonsig_idx) > remaining) {
    set.seed(42)
    nonsig_idx <- sample(nonsig_idx, remaining)
  }

  data[c(sig_idx, nonsig_idx), , drop = FALSE]
}

dir.create(dirname(output_merged), recursive = TRUE, showWarnings = FALSE)
write.table(df, file = output_merged, sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.na(output_png) && nzchar(output_png) && !is.na(output_pdf) && nzchar(output_pdf)) {
  df_plot <- downsample_for_manhattan(df, "p_gwas")
  p <- ggplot(df_plot, aes(x = BP_cum, y = -log10(p_gwas), color = CHROM)) +
    geom_point(size = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = -log10(SIGNIF_P_SUGGESTIVE), linetype = "dashed", color = "#ff7f0e") +
    geom_hline(yintercept = -log10(SIGNIF_P_MAIN), linetype = "dashed", color = "#d62728") +
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
}

if (!is.na(output_qq_png) && nzchar(output_qq_png) && !is.na(output_qq_pdf) && nzchar(output_qq_pdf)) {
  pvals <- df$p_gwas
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) > QQ_MAX_POINTS) {
    set.seed(42)
    pvals <- sample(pvals, QQ_MAX_POINTS)
  }

  if (length(pvals) > 0) {
    observed <- -log10(sort(pvals, decreasing = FALSE))
    expected <- -log10((seq_len(length(pvals)) - 0.5) / length(pvals))
    qq_df <- data.frame(expected = expected, observed = observed)

    qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
      geom_point(size = 0.8, alpha = 0.7, color = "#1f78b4") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#d62728") +
      labs(
        title = paste0("QQ Plot: Tractor Pairwise GWAS: ", pair_name),
        x = "Expected -log10(p)",
        y = "Observed -log10(p)"
      ) +
      theme_bw()
  } else {
    qq_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No valid p-values for QQ plot", size = 5) +
      labs(
        title = paste0("QQ Plot: Tractor Pairwise GWAS: ", pair_name),
        x = "Expected -log10(p)",
        y = "Observed -log10(p)"
      ) +
      theme_void()
  }

  dir.create(dirname(output_qq_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_qq_png, plot = qq_plot, width = 6, height = 6, dpi = 300)
  ggsave(filename = output_qq_pdf, plot = qq_plot, width = 6, height = 6)
}

build_manhattan <- function(data, p_col, title_text) {
  data_plot <- downsample_for_manhattan(data, p_col)
  ggplot(data_plot, aes_string(x = "BP_cum", y = paste0("-log10(", p_col, ")"), color = "CHROM")) +
    geom_point(size = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = -log10(SIGNIF_P_SUGGESTIVE), linetype = "dashed", color = "#ff7f0e") +
    geom_hline(yintercept = -log10(SIGNIF_P_MAIN), linetype = "dashed", color = "#d62728") +
    scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$BP_cum) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(
      title = title_text,
      x = "Chromosome",
      y = "-log10(p)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

build_qq <- function(data, p_col, title_text) {
  pvals <- data[[p_col]]
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) > QQ_MAX_POINTS) {
    set.seed(42)
    pvals <- sample(pvals, QQ_MAX_POINTS)
  }

  if (length(pvals) > 0) {
    observed <- -log10(sort(pvals, decreasing = FALSE))
    expected <- -log10((seq_len(length(pvals)) - 0.5) / length(pvals))
    qq_df <- data.frame(expected = expected, observed = observed)

    return(
      ggplot(qq_df, aes(x = expected, y = observed)) +
        geom_point(size = 0.8, alpha = 0.7, color = "#1f78b4") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#d62728") +
        labs(
          title = title_text,
          x = "Expected -log10(p)",
          y = "Observed -log10(p)"
        ) +
        theme_bw()
    )
  }

  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No valid p-values for QQ plot", size = 5) +
    labs(
      title = title_text,
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_void()
}

write_plot_set <- function(data, p_col, title_stub, manhattan_png, manhattan_pdf, qq_png, qq_pdf) {
  if (any(is.na(c(manhattan_png, manhattan_pdf, qq_png, qq_pdf))) ||
      !nzchar(manhattan_png) || !nzchar(manhattan_pdf) || !nzchar(qq_png) || !nzchar(qq_pdf)) {
    return(invisible(NULL))
  }

  manhattan_plot <- build_manhattan(data, p_col, paste0("Tractor Pairwise GWAS: ", pair_name, " (", title_stub, ")"))
  qq_plot <- build_qq(data, p_col, paste0("QQ Plot: Tractor Pairwise GWAS: ", pair_name, " (", title_stub, ")"))

  dir.create(dirname(manhattan_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = manhattan_png, plot = manhattan_plot, width = 12, height = 5, dpi = 300)
  ggsave(filename = manhattan_pdf, plot = manhattan_plot, width = 12, height = 5)
  ggsave(filename = qq_png, plot = qq_plot, width = 6, height = 6, dpi = 300)
  ggsave(filename = qq_pdf, plot = qq_plot, width = 6, height = 6)
}

write_plot_set(
  df,
  "p_dosage_i",
  "Dosage I",
  output_dosage_i_manhattan_png,
  output_dosage_i_manhattan_pdf,
  output_dosage_i_qq_png,
  output_dosage_i_qq_pdf
)

write_plot_set(
  df,
  "p_dosage_j",
  "Dosage J",
  output_dosage_j_manhattan_png,
  output_dosage_j_manhattan_pdf,
  output_dosage_j_qq_png,
  output_dosage_j_qq_pdf
)

write_plot_set(
  df,
  "p_local_ancestry",
  "Local Ancestry",
  output_local_manhattan_png,
  output_local_manhattan_pdf,
  output_local_qq_png,
  output_local_qq_pdf
)
