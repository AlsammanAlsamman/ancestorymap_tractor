#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing required argument", flag), call. = FALSE)
    return(default)
  }
  if (idx[length(idx)] == length(args)) {
    stop(paste("Missing value for argument", flag), call. = FALSE)
  }
  args[idx[length(idx)] + 1]
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

sanitize_name <- function(x) {
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("(^_+|_+$)", "", out)
  ifelse(nzchar(out), out, "unknown_locus")
}

save_png <- function(plot_obj, filename, width = 2400, height = 900, res = 160) {
  if (capabilities("cairo")) {
    png(filename, width = width, height = height, res = res, type = "cairo")
  } else {
    png(filename, width = width, height = height, res = res)
  }
  print(plot_obj)
  dev.off()
}

detect_ancestries <- function(dt) {
  cols <- grep("^p_dosage_", names(dt), value = TRUE)
  sub("^p_dosage_", "", cols)
}

detect_local_pairs <- function(dt) {
  cols <- grep("^p_local_ancestry_", names(dt), value = TRUE)
  sub("^p_local_ancestry_", "", cols)
}

required_base_columns <- c("chr", "pos", "rsid", "locus_name", "gene")

input_tsv <- get_arg("--input", required = TRUE)
output_dir <- get_arg("--output-dir", required = TRUE)
local_logp_min <- safe_num(get_arg("--local-logp-min", default = "1"))
plot_width <- as.integer(safe_num(get_arg("--width", default = "2400")))
plot_height <- as.integer(safe_num(get_arg("--height", default = "900")))
top_label_nudge <- safe_num(get_arg("--top-label-nudge", default = "0.35"))

if (!file.exists(input_tsv)) stop(paste("Input TSV not found:", input_tsv), call. = FALSE)
if (!is.finite(local_logp_min)) stop("--local-logp-min must be numeric", call. = FALSE)
if (!is.finite(plot_width) || plot_width <= 0) stop("--width must be a positive integer", call. = FALSE)
if (!is.finite(plot_height) || plot_height <= 0) stop("--height must be a positive integer", call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_tsv, sep = "\t", data.table = TRUE, check.names = FALSE)

missing_base <- setdiff(required_base_columns, names(dt))
if (length(missing_base)) {
  stop(sprintf("Input TSV is missing required columns: %s", paste(missing_base, collapse = ", ")), call. = FALSE)
}

ancestries <- detect_ancestries(dt)
if (length(ancestries) != 3) {
  stop(sprintf("Expected exactly 3 ancestry dosage columns, found %d: %s", length(ancestries), paste(ancestries, collapse = ", ")), call. = FALSE)
}

local_pairs <- detect_local_pairs(dt)
if (!length(local_pairs)) {
  stop("No local ancestry p-value columns found in input TSV", call. = FALSE)
}

required_dyn <- c(
  unlist(lapply(ancestries, function(anc) c(paste0("p_dosage_", anc), paste0("OR_dosage_", anc)))),
  paste0("p_local_ancestry_", local_pairs)
)
missing_dyn <- setdiff(required_dyn, names(dt))
if (length(missing_dyn)) {
  stop(sprintf("Input TSV is missing required ancestry columns: %s", paste(missing_dyn, collapse = ", ")), call. = FALSE)
}

dt[, chr := as.character(chr)]
dt[, pos := safe_num(pos)]
for (nm in required_dyn) dt[, (nm) := safe_num(get(nm))]

local_cols <- paste0("p_local_ancestry_", local_pairs)
dt[, best_local_p := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = local_cols]
dt[!is.finite(best_local_p), best_local_p := NA_real_]
dt[, best_local_logp := -log10(pmax(best_local_p, 1e-300))]

filtered <- dt[is.finite(best_local_logp) & best_local_logp >= local_logp_min]
if (!nrow(filtered)) {
  stop(sprintf("No SNPs passed the local ancestry -log10(p) threshold of %s", local_logp_min), call. = FALSE)
}

threshold_dt <- data.table(
  threshold_label = c("5e-8", "5e-5"),
  threshold_y = c(-log10(5e-8), -log10(5e-5))
)

ancestry_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")
names(ancestry_colors) <- ancestries

summary_rows <- list()

for (locus in unique(filtered$locus_name)) {
  locus_dt <- filtered[locus_name == locus & is.finite(pos)]
  if (!nrow(locus_dt)) next

  long_dt <- rbindlist(lapply(ancestries, function(anc) {
    p_col <- paste0("p_dosage_", anc)
    or_col <- paste0("OR_dosage_", anc)
    data.table(
      chr = locus_dt$chr,
      pos = locus_dt$pos,
      rsid = locus_dt$rsid,
      locus_name = locus_dt$locus_name,
      gene = locus_dt$gene,
      ancestry = anc,
      p_value = locus_dt[[p_col]],
      odds_ratio = locus_dt[[or_col]]
    )
  }), use.names = TRUE)

  long_dt <- long_dt[is.finite(pos) & is.finite(p_value) & p_value > 0]
  if (!nrow(long_dt)) next
  long_dt[, logp := -log10(pmax(p_value, 1e-300))]
  long_dt[, ancestry := factor(ancestry, levels = ancestries)]

  top_idx <- which.min(long_dt$p_value)
  top_row <- long_dt[top_idx]
  locus_chr <- unique(as.character(locus_dt$chr))[1]
  locus_start <- min(locus_dt$pos, na.rm = TRUE)
  locus_end <- max(locus_dt$pos, na.rm = TRUE)
  top_gene <- unique(na.omit(trimws(as.character(top_row$gene))))
  top_gene <- if (length(top_gene)) top_gene[1] else "NA"

  label_dt <- copy(top_row)
  label_dt[, label := as.character(rsid)]

  title_text <- sprintf(
    "%s | chr%s:%s-%s | top SNP: %s | gene: %s",
    as.character(locus),
    locus_chr,
    format(as.integer(locus_start), big.mark = ",", scientific = FALSE),
    format(as.integer(locus_end), big.mark = ",", scientific = FALSE),
    as.character(top_row$rsid),
    top_gene
  )

  subtitle_text <- sprintf(
    "Filtered on max local ancestry -log10(p) >= %s across %s",
    format(local_logp_min, trim = TRUE),
    paste(local_pairs, collapse = ", ")
  )

  y_max <- max(c(long_dt$logp, threshold_dt$threshold_y), na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- 8
  y_max <- min(max(8, y_max + 0.8), 60)

  plot_obj <- ggplot(long_dt, aes(x = pos, y = logp, color = ancestry)) +
    geom_hline(data = threshold_dt, aes(yintercept = threshold_y, linetype = threshold_label), color = "#666666", linewidth = 0.45, inherit.aes = FALSE) +
    geom_point(size = 1.8, alpha = 0.9) +
    geom_text(
      data = label_dt,
      aes(label = label),
      inherit.aes = TRUE,
      nudge_y = top_label_nudge,
      size = 3.2,
      fontface = "bold",
      show.legend = FALSE
    ) +
    facet_wrap(~ ancestry, nrow = 1, scales = "fixed") +
    scale_color_manual(values = ancestry_colors, drop = FALSE) +
    scale_linetype_manual(values = c("5e-8" = "solid", "5e-5" = "dashed"), drop = FALSE) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Genomic position",
      y = expression(-log[10](p[dose]))
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  out_file <- file.path(output_dir, paste0(sanitize_name(locus), ".local_ancestry_filtered_dosage_panels.png"))
  save_png(plot_obj, out_file, width = plot_width, height = plot_height)

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    locus_name = as.character(locus),
    chr = locus_chr,
    start_pos = as.integer(locus_start),
    end_pos = as.integer(locus_end),
    top_snp = as.character(top_row$rsid),
    top_gene = top_gene,
    top_ancestry = as.character(top_row$ancestry),
    top_p_dosage = as.numeric(top_row$p_value),
    max_local_logp = max(locus_dt$best_local_logp, na.rm = TRUE),
    n_snps_plotted = nrow(unique(long_dt[, .(pos, rsid)])),
    plot_file = out_file
  )
}

if (!length(summary_rows)) {
  stop("No locus plots were produced after filtering and ancestry conversion", call. = FALSE)
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
summary_file <- file.path(output_dir, "locus_plot_manifest.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

cat(sprintf(
  "Wrote %d locus plots to %s and manifest to %s\n",
  nrow(summary_dt), output_dir, summary_file
))