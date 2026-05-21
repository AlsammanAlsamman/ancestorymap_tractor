#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) != 1 || idx == length(args)) {
    stop(paste("Missing argument", flag))
  }
  args[idx + 1]
}

get_arg_vector <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) != 1) stop(paste("Missing argument", flag))
  start <- idx + 1
  if (start > length(args)) stop(paste("Missing values for", flag))
  end <- start
  while (end <= length(args) && !grepl("^--", args[end])) {
    end <- end + 1
  }
  args[start:(end - 1)]
}

input_tsv <- get_arg_vector("--input-tsv")
png_out <- get_arg("--png")
pdf_out <- get_arg("--pdf")
csv_out <- get_arg("--csv")

if (length(input_tsv) == 0) stop("No per-chromosome TSV files provided")

all_dt <- rbindlist(lapply(input_tsv, function(path) {
  if (!file.exists(path)) stop(paste("Missing input TSV:", path))
  fread(path)
}), fill = TRUE)

required_cols <- c("chr", "ancestry_code", "ancestry_label", "r2")
if (!all(required_cols %in% names(all_dt))) {
  stop("Input TSV files missing required columns: chr, ancestry_code, ancestry_label, r2")
}

all_dt[, chr_num := suppressWarnings(as.integer(as.character(chr)))]
if (all(is.na(all_dt$chr_num))) {
  all_dt[, chr_factor := factor(as.character(chr), levels = unique(as.character(chr)))]
} else {
  chr_levels <- all_dt[order(chr_num), unique(as.character(chr))]
  all_dt[, chr_factor := factor(as.character(chr), levels = chr_levels)]
}

all_dt[, ancestry_label := factor(as.character(ancestry_label), levels = unique(as.character(ancestry_label)))]

plot_dt <- all_dt[, .(chr = as.character(chr), chr_factor, ancestry_label, r2)]
plot_dt[, r2_label := ifelse(is.na(r2), "NA", sprintf("%.3f", r2))]

p <- ggplot(plot_dt, aes(x = ancestry_label, y = chr_factor, fill = r2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = r2_label), size = 3) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08306b",
    limits = c(0, 1),
    na.value = "grey90",
    name = expression(R^2)
  ) +
  labs(
    title = "RFMix ancestry vs PRS correlation matrix",
    x = "Ancestry",
    y = "Chromosome"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

wide_dt <- dcast(
  plot_dt,
  chr ~ ancestry_label,
  value.var = "r2"
)

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(csv_out), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = png_out, plot = p, width = 7, height = 5.5, dpi = 300)
ggsave(filename = pdf_out, plot = p, width = 7, height = 5.5)
fwrite(wide_dt, file = csv_out)
