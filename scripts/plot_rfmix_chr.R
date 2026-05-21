#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: plot_rfmix_chr.R --msp <file> --chr <chr> --png <png> --pdf <pdf> --label0 <label> --label1 <label> --label2 <label> --color0 <color> --color1 <color> --color2 <color>")
}

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) != 1 || idx == length(args)) {
    stop(paste("Missing argument", flag))
  }
  args[idx + 1]
}

msp_file <- get_arg("--msp")
chr_id <- get_arg("--chr")
png_out <- get_arg("--png")
pdf_out <- get_arg("--pdf")
label0 <- get_arg("--label0")
label1 <- get_arg("--label1")
label2 <- get_arg("--label2")
color0 <- get_arg("--color0")
color1 <- get_arg("--color1")
color2 <- get_arg("--color2")

if (!file.exists(msp_file)) {
  stop(paste("MSP file not found:", msp_file))
}

dt <- fread(msp_file)
if (!all(c("#chm", "spos", "epos") %in% names(dt))) {
  stop("MSP file missing required columns: #chm, spos, epos")
}

sample_cols <- grep("\\.(0|1)$", names(dt), value = TRUE)
if (length(sample_cols) == 0) {
  stop("MSP file has no haplotype ancestry columns ending in .0/.1")
}

long_dt <- melt(
  dt,
  id.vars = c("#chm", "spos", "epos"),
  measure.vars = sample_cols,
  variable.name = "haplotype",
  value.name = "ancestry_code"
)

summary_dt <- long_dt[, .N, by = .(spos, epos, ancestry_code)]
summary_dt[, width := epos - spos]
summary_dt[, prop := N / sum(N), by = .(spos, epos)]
summary_dt <- summary_dt[order(spos, ancestry_code)]
summary_dt[, ymax := cumsum(prop), by = .(spos, epos)]
summary_dt[, ymin := shift(ymax, fill = 0), by = .(spos, epos)]
summary_dt[, ancestry_code := as.character(ancestry_code)]

label_map <- c("0" = label0, "1" = label1, "2" = label2)
color_map <- c("0" = color0, "1" = color1, "2" = color2)

summary_dt[, ancestry_label := ifelse(
  ancestry_code %in% names(label_map),
  label_map[ancestry_code],
  paste0("anc", ancestry_code)
)]

plot_df <- summary_dt[, .(
  x_start = spos,
  x_end = epos,
  y_min = ymin,
  y_max = ymax,
  ancestry_label = factor(ancestry_label, levels = unique(label_map))
)]

named_colors <- c(label0 = color0, label1 = color1, label2 = color2)
names(named_colors) <- c(label0, label1, label2)

p <- ggplot(plot_df) +
  geom_rect(aes(xmin = x_start, xmax = x_end, ymin = y_min, ymax = y_max, fill = ancestry_label), color = NA) +
  scale_fill_manual(values = named_colors, drop = FALSE) +
  labs(
    title = paste0("RFMix Local Ancestry Proportions - chr", chr_id),
    x = "Genomic position (bp)",
    y = "Ancestry proportion",
    fill = "Ancestry"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(filename = png_out, plot = p, width = 12, height = 4, dpi = 300)
ggsave(filename = pdf_out, plot = p, width = 12, height = 4)
