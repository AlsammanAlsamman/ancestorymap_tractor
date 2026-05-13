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

required_flags <- c(
  "--msp", "--chr", "--phenotype", "--phenotype-format",
  "--fam-case-code", "--fam-control-code",
  "--table-sample-column", "--table-status-column", "--table-case-value", "--table-control-value",
  "--genes", "--bin-size-bp", "--max-genes-per-bin-label",
  "--png", "--pdf", "--tsv",
  "--label0", "--label1", "--label2", "--color0", "--color1", "--color2"
)

for (f in required_flags) {
  if (!any(args == f)) stop(paste("Missing required flag", f))
}

msp_file <- get_arg("--msp")
chr_id <- get_arg("--chr")
phenotype_file <- get_arg("--phenotype")
phenotype_format <- tolower(get_arg("--phenotype-format"))
fam_case_code <- get_arg("--fam-case-code")
fam_control_code <- get_arg("--fam-control-code")
table_sample_column <- get_arg("--table-sample-column")
table_status_column <- get_arg("--table-status-column")
table_case_value <- tolower(get_arg("--table-case-value"))
table_control_value <- tolower(get_arg("--table-control-value"))
gene_file <- get_arg("--genes")
bin_size_bp <- as.integer(get_arg("--bin-size-bp"))
max_genes_per_bin_label <- as.integer(get_arg("--max-genes-per-bin-label"))
png_out <- get_arg("--png")
pdf_out <- get_arg("--pdf")
tsv_out <- get_arg("--tsv")

label0 <- get_arg("--label0")
label1 <- get_arg("--label1")
label2 <- get_arg("--label2")
color0 <- get_arg("--color0")
color1 <- get_arg("--color1")
color2 <- get_arg("--color2")

if (!file.exists(msp_file)) stop(paste("MSP file not found:", msp_file))
if (!file.exists(phenotype_file)) stop(paste("Phenotype file not found:", phenotype_file))
if (!file.exists(gene_file)) stop(paste("Gene annotation file not found:", gene_file))
if (is.na(bin_size_bp) || bin_size_bp <= 0) stop("--bin-size-bp must be a positive integer")
if (is.na(max_genes_per_bin_label) || max_genes_per_bin_label <= 0) stop("--max-genes-per-bin-label must be a positive integer")

normalize_chr <- function(x) {
  sub("^chr", "", tolower(as.character(x)))
}

extract_sample_id <- function(hap_col) {
  x <- sub("\\.[01]$", "", hap_col)
  us <- gregexpr("_", x, fixed = TRUE)[[1]]
  if (length(us) > 0 && us[1] != -1) {
    for (p in us) {
      left <- substr(x, 1, p - 1)
      right <- substr(x, p + 1, nchar(x))
      if (left == right) {
        return(left)
      }
    }
  }
  x
}

load_phenotype <- function(path, format_mode) {
  if (format_mode == "fam") {
    ph <- fread(path, header = FALSE)
    if (ncol(ph) < 6) stop("FAM phenotype file must have at least 6 columns")
    out <- data.table(
      sample_id = as.character(ph[[2]]),
      status_raw = as.character(ph[[6]])
    )
    out[status_raw == fam_case_code, group := "case"]
    out[status_raw == fam_control_code, group := "control"]
    return(out[!is.na(group), .(sample_id, group)])
  }

  if (format_mode == "table") {
    ph <- fread(path)
    if (!(table_sample_column %in% names(ph))) stop(paste("Missing table sample column:", table_sample_column))
    if (!(table_status_column %in% names(ph))) stop(paste("Missing table status column:", table_status_column))

    out <- data.table(
      sample_id = as.character(ph[[table_sample_column]]),
      status_raw = tolower(as.character(ph[[table_status_column]]))
    )
    out[status_raw == table_case_value, group := "case"]
    out[status_raw == table_control_value, group := "control"]
    return(out[!is.na(group), .(sample_id, group)])
  }

  stop("--phenotype-format must be 'fam' or 'table'")
}

dt <- fread(msp_file, skip = "#chm")
required_cols <- c("#chm", "spos", "epos")
if (!all(required_cols %in% names(dt))) {
  stop("MSP file missing required columns: #chm, spos, epos")
}

dt[, chr_norm := normalize_chr(`#chm`)]
target_chr_norm <- normalize_chr(chr_id)
dt <- dt[which(dt[["chr_norm"]] == target_chr_norm)]
if (nrow(dt) == 0) {
  stop(paste("No rows in MSP for chromosome", chr_id))
}

meta_cols <- c("#chm", "spos", "epos", "sgpos", "egpos", "n snps", "chr_norm")
hap_cols <- setdiff(names(dt), meta_cols)
if (length(hap_cols) == 0) stop("MSP file has no haplotype ancestry columns")

phenotype_dt <- unique(load_phenotype(phenotype_file, phenotype_format))
if (nrow(phenotype_dt) == 0) stop("No usable case/control samples found in phenotype file")

hap_map <- data.table(hap_col = hap_cols)
hap_map[, sample_id := vapply(hap_col, extract_sample_id, FUN.VALUE = character(1))]
hap_map <- merge(hap_map, phenotype_dt, by = "sample_id", all.x = FALSE, all.y = FALSE)
hap_map <- unique(hap_map[, .(hap_col, sample_id, group)])
if (nrow(hap_map) == 0) {
  stop("No MSP haplotype columns matched phenotype samples")
}

long_dt <- melt(
  dt,
  id.vars = c("spos", "epos"),
  measure.vars = hap_map$hap_col,
  variable.name = "hap_col",
  value.name = "ancestry_code"
)
long_dt[, ancestry_code := as.character(ancestry_code)]

hap_group_map <- unique(hap_map[, .(hap_col, group)])
setkey(hap_group_map, hap_col)
setkey(long_dt, hap_col)
long_dt <- hap_group_map[long_dt]
long_dt <- long_dt[!is.na(group)]

if (nrow(long_dt) == 0) {
  stop("No long-format ancestry rows after phenotype matching")
}

chr_start <- min(dt$spos)
chr_end <- max(dt$epos)
bin_starts <- seq(chr_start, chr_end, by = bin_size_bp)
bin_dt <- data.table(
  bin_id = seq_along(bin_starts),
  bin_start = as.integer(bin_starts),
  bin_end = as.integer(pmin(bin_starts + bin_size_bp - 1, chr_end))
)

long_dt[, seg_start := as.integer(spos)]
long_dt[, seg_end := as.integer(epos)]
setkey(long_dt, seg_start, seg_end)
bin_overlap_dt <- copy(bin_dt)
setnames(bin_overlap_dt, c("bin_start", "bin_end"), c("ov_bin_start", "ov_bin_end"))
setkey(bin_overlap_dt, ov_bin_start, ov_bin_end)

ov <- foverlaps(
  long_dt[, .(seg_start, seg_end, group, ancestry_code)],
  bin_overlap_dt,
  by.x = c("seg_start", "seg_end"),
  by.y = c("ov_bin_start", "ov_bin_end"),
  nomatch = 0L
)

ov[, overlap_bp := pmax(0L, pmin(seg_end, ov_bin_end) - pmax(seg_start, ov_bin_start) + 1L)]
ov <- ov[overlap_bp > 0]

agg <- ov[, .(bp = sum(overlap_bp)), by = .(bin_id, group, ancestry_code)]

groups_present <- c("case", "control")
ancestries_present <- sort(unique(agg$ancestry_code))
complete_grid <- CJ(
  bin_id = unique(bin_dt$bin_id),
  group = groups_present,
  ancestry_code = ancestries_present,
  unique = TRUE
)

agg <- merge(
  complete_grid,
  agg,
  by = c("bin_id", "group", "ancestry_code"),
  all.x = TRUE
)

agg <- merge(agg, bin_dt, by = "bin_id", all.x = TRUE)
agg[is.na(bp), bp := 0]
agg[, total_bp := sum(bp), by = .(bin_id, group)]
agg[, prop := fifelse(total_bp > 0, bp / total_bp, 0)]

gene_dt <- fread(gene_file, header = FALSE)
if (ncol(gene_dt) < 6) stop("Gene annotation file must have at least 6 columns")

gene_dt <- gene_dt[, .(
  chr_norm = normalize_chr(V2),
  gene_start = as.integer(V3),
  gene_end = as.integer(V4),
  gene_name = as.character(V6)
)]
gene_dt <- gene_dt[which(gene_dt[["chr_norm"]] == target_chr_norm)]
gene_dt <- gene_dt[!is.na(gene_start) & !is.na(gene_end) & gene_start <= gene_end]

bin_gene_labels <- data.table(bin_id = bin_dt$bin_id, gene_label = "-")
if (nrow(gene_dt) > 0) {
  setkey(gene_dt, gene_start, gene_end)
  bin_gene_dt <- copy(bin_dt)
  setnames(bin_gene_dt, c("bin_start", "bin_end"), c("gene_bin_start", "gene_bin_end"))
  setkey(bin_gene_dt, gene_bin_start, gene_bin_end)
  gene_ov <- foverlaps(
    gene_dt[, .(gene_start, gene_end, gene_name)],
    bin_gene_dt,
    by.x = c("gene_start", "gene_end"),
    by.y = c("gene_bin_start", "gene_bin_end"),
    nomatch = 0L
  )

  if (nrow(gene_ov) > 0) {
    gene_labels <- gene_ov[
      , {
          g <- unique(gene_name)
          n_g <- length(g)
          if (n_g > max_genes_per_bin_label) {
            g <- c(g[1:max_genes_per_bin_label], "...")
          }
          .(gene_label = paste(g, collapse = ","))
        },
      by = .(bin_id)
    ]
    bin_gene_labels <- merge(bin_gene_labels, gene_labels, by = "bin_id", all.x = TRUE, suffixes = c("", ".new"))
    bin_gene_labels[!is.na(gene_label.new), gene_label := gene_label.new]
    bin_gene_labels[, gene_label.new := NULL]
  }
}

agg <- merge(agg, bin_gene_labels, by = "bin_id", all.x = TRUE)
agg[is.na(gene_label), gene_label := "-"]

label_map <- c("0" = label0, "1" = label1, "2" = label2)
agg[, ancestry_label := ifelse(
  ancestry_code %in% names(label_map),
  label_map[ancestry_code],
  paste0("anc", ancestry_code)
)]

group_levels <- c("case", "control")
agg[, group := factor(group, levels = group_levels)]

known_ancestries <- c(label0, label1, label2)
other_ancestries <- setdiff(unique(agg$ancestry_label), known_ancestries)
agg[, ancestry_label := factor(ancestry_label, levels = c(known_ancestries, sort(other_ancestries)))]

bin_label_dt <- unique(agg[, .(bin_id, gene_label)])
bin_label_dt[, bin_text := paste0("Bin", bin_id, "\n", gene_label)]
agg <- merge(agg, bin_label_dt[, .(bin_id, bin_text)], by = "bin_id", all.x = TRUE)
agg[, x_label := paste0(bin_text, "\n", as.character(group))]

ordered_x <- unique(agg[order(bin_id, group), x_label])
agg[, x_label := factor(x_label, levels = ordered_x)]

named_colors <- c(label0 = color0, label1 = color1, label2 = color2)
names(named_colors) <- c(label0, label1, label2)
missing_ancestry <- setdiff(levels(agg$ancestry_label), names(named_colors))
if (length(missing_ancestry) > 0) {
  missing_cols <- rep("#999999", length(missing_ancestry))
  names(missing_cols) <- missing_ancestry
  named_colors <- c(named_colors, missing_cols)
}

plot_title <- paste0(
  "RFMix binned ancestry by case/control - chr", chr_id,
  " (bin=", format(bin_size_bp, scientific = FALSE), " bp)"
)

p <- ggplot(agg, aes(x = x_label, y = prop, fill = ancestry_label)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = named_colors, drop = FALSE) +
  labs(
    title = plot_title,
    x = "Bin / Group (with intersecting genes)",
    y = "Ancestry proportion",
    fill = "Ancestry"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

summary_dt <- agg[, .(
  chr = chr_id,
  bin_id,
  bin_start,
  bin_end,
  gene_label,
  group = as.character(group),
  ancestry_code,
  ancestry_label = as.character(ancestry_label),
  bp,
  prop
)]

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(tsv_out), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = png_out, plot = p, width = 16, height = 7, dpi = 300)
ggsave(filename = pdf_out, plot = p, width = 16, height = 7)
fwrite(summary_dt, file = tsv_out, sep = "\t")
