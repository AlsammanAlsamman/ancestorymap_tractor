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
  "--q", "--chr", "--phenotype", "--phenotype-format", "--fam-case-code", "--fam-control-code", "--prs", "--prs-iid-column", "--prs-score-column",
  "--prs-in-regression-column", "--prs-in-regression-keep-value",
  "--ancestry-bins",
  "--png", "--pdf", "--tsv", "--cc-png", "--cc-pdf", "--cc-tsv",
  "--label0", "--label1", "--label2", "--color0", "--color1", "--color2"
)
for (f in required_flags) {
  if (!any(args == f)) stop(paste("Missing required flag", f))
}

q_file <- get_arg("--q")
chr_id <- get_arg("--chr")
phenotype_file <- get_arg("--phenotype")
phenotype_format <- tolower(get_arg("--phenotype-format"))
fam_case_code <- get_arg("--fam-case-code")
fam_control_code <- get_arg("--fam-control-code")
prs_file <- get_arg("--prs")
prs_iid_column <- get_arg("--prs-iid-column")
prs_score_column <- get_arg("--prs-score-column")
prs_in_reg_col <- get_arg("--prs-in-regression-column")
prs_in_reg_keep <- get_arg("--prs-in-regression-keep-value")
ancestry_bins <- as.integer(get_arg("--ancestry-bins"))
png_out <- get_arg("--png")
pdf_out <- get_arg("--pdf")
tsv_out <- get_arg("--tsv")
cc_png_out <- get_arg("--cc-png")
cc_pdf_out <- get_arg("--cc-pdf")
cc_tsv_out <- get_arg("--cc-tsv")

label0 <- get_arg("--label0")
label1 <- get_arg("--label1")
label2 <- get_arg("--label2")
color0 <- get_arg("--color0")
color1 <- get_arg("--color1")
color2 <- get_arg("--color2")

if (!file.exists(q_file)) stop(paste("Q file not found:", q_file))
if (!file.exists(phenotype_file)) stop(paste("Phenotype file not found:", phenotype_file))
if (!file.exists(prs_file)) stop(paste("PRS file not found:", prs_file))
if (is.na(ancestry_bins) || ancestry_bins < 2) stop("--ancestry-bins must be an integer >= 2")
if (phenotype_format != "fam") stop("--phenotype-format currently supports only 'fam'")

extract_sample_id <- function(q_sample) {
  x <- as.character(q_sample)
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

q_dt <- fread(q_file, skip = "#sample")
if (ncol(q_dt) < 4) {
  stop("Q file must have at least 4 columns: sample + three ancestry columns")
}

setnames(q_dt, old = names(q_dt)[1:4], new = c("sample_raw", "a0", "a1", "a2"))
q_dt <- q_dt[, .(sample_raw, a0 = as.numeric(a0), a1 = as.numeric(a1), a2 = as.numeric(a2))]
q_dt[, sample_id := vapply(sample_raw, extract_sample_id, FUN.VALUE = character(1))]

prs_dt <- fread(prs_file)
if (!(prs_iid_column %in% names(prs_dt))) stop(paste("PRS missing IID column:", prs_iid_column))
if (!(prs_score_column %in% names(prs_dt))) stop(paste("PRS missing score column:", prs_score_column))

if (prs_in_reg_col %in% names(prs_dt)) {
  prs_dt <- prs_dt[tolower(as.character(get(prs_in_reg_col))) == tolower(prs_in_reg_keep)]
}

prs_dt <- prs_dt[, .(
  sample_id = as.character(get(prs_iid_column)),
  prs = as.numeric(get(prs_score_column))
)]
prs_dt <- prs_dt[!is.na(prs)]

ph_dt <- fread(phenotype_file, header = FALSE)
if (ncol(ph_dt) < 6) stop("Phenotype FAM file must have at least 6 columns")
ph_dt <- ph_dt[, .(
  sample_id = as.character(V2),
  pheno = as.character(V6)
)]
ph_dt[pheno == fam_case_code, group := "case"]
ph_dt[pheno == fam_control_code, group := "control"]
ph_dt <- ph_dt[!is.na(group), .(sample_id, group)]
ph_dt <- unique(ph_dt)

merged <- merge(
  q_dt[, .(sample_id, a0, a1, a2)],
  prs_dt,
  by = "sample_id",
  all = FALSE
)
merged <- merge(merged, ph_dt, by = "sample_id", all = FALSE)

if (nrow(merged) < 3) {
  stop("Not enough matched samples between Q/PRS/phenotype files (need >=3)")
}

long_dt <- melt(
  merged,
  id.vars = c("sample_id", "prs", "group"),
  measure.vars = c("a0", "a1", "a2"),
  variable.name = "ancestry_code_col",
  value.name = "ancestry_prop"
)

code_map <- c("a0" = "0", "a1" = "1", "a2" = "2")
label_map <- c("0" = label0, "1" = label1, "2" = label2)
long_dt[, ancestry_code := code_map[ancestry_code_col]]
long_dt[, ancestry_label := label_map[ancestry_code]]
long_dt[, group := factor(group, levels = c("case", "control"))]

long_dt[
  , ancestry_bin := cut(
      ancestry_prop,
      breaks = seq(0, 1, length.out = ancestry_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
]

bin_dt <- long_dt[
  , .(
      ancestry_prop = mean(ancestry_prop, na.rm = TRUE),
      prs = mean(prs, na.rm = TRUE),
      n_samples = .N
    ),
  by = .(ancestry_code, ancestry_label, ancestry_bin)
]

corr_dt <- bin_dt[
  , {
      n <- .N
      if (n < 3 || sd(ancestry_prop, na.rm = TRUE) == 0 || sd(prs, na.rm = TRUE) == 0) {
        .(n_samples = n, r = NA_real_, r2 = NA_real_, p_value = NA_real_)
      } else {
        ct <- suppressWarnings(cor.test(ancestry_prop, prs, method = "pearson"))
        r_val <- unname(ct$estimate)
        .(n_samples = n, r = r_val, r2 = r_val^2, p_value = ct$p.value)
      }
    },
  by = .(ancestry_code, ancestry_label)
]

bin_cc_dt <- long_dt[
  , .(
      ancestry_prop = mean(ancestry_prop, na.rm = TRUE),
      prs = mean(prs, na.rm = TRUE),
      n_samples = .N
    ),
  by = .(group, ancestry_code, ancestry_label, ancestry_bin)
]

corr_cc_dt <- bin_cc_dt[
  , {
      n <- .N
      if (n < 3 || sd(ancestry_prop, na.rm = TRUE) == 0 || sd(prs, na.rm = TRUE) == 0) {
        .(n_samples = n, r = NA_real_, r2 = NA_real_, p_value = NA_real_)
      } else {
        ct <- suppressWarnings(cor.test(ancestry_prop, prs, method = "pearson"))
        r_val <- unname(ct$estimate)
        .(n_samples = n, r = r_val, r2 = r_val^2, p_value = ct$p.value)
      }
    },
  by = .(group, ancestry_code, ancestry_label)
]

cc_n_individuals_dt <- long_dt[
  , .(n_individuals = uniqueN(sample_id)),
  by = .(group, ancestry_code, ancestry_label)
]

corr_dt[, chr := as.character(chr_id)]
setcolorder(corr_dt, c("chr", "ancestry_code", "ancestry_label", "n_samples", "r", "r2", "p_value"))

corr_cc_dt[, chr := as.character(chr_id)]
setcolorder(corr_cc_dt, c("chr", "group", "ancestry_code", "ancestry_label", "n_samples", "r", "r2", "p_value"))

cc_annot_dt <- merge(
  corr_cc_dt,
  cc_n_individuals_dt,
  by = c("group", "ancestry_code", "ancestry_label"),
  all.x = TRUE
)
cc_annot_dt[, annot_label := paste0(
  "R2=", ifelse(is.na(r2), "NA", sprintf("%.3f", r2)),
  "\nN=", ifelse(is.na(n_individuals), "NA", as.character(n_individuals))
)]

named_colors <- c(label0 = color0, label1 = color1, label2 = color2)
names(named_colors) <- c(label0, label1, label2)

subtitle_text <- corr_dt[
  , paste0(ancestry_label, ": R2=", ifelse(is.na(r2), "NA", sprintf("%.3f", r2))),
  by = ancestry_label
][, paste(V1, collapse = " | ")]

long_dt[, ancestry_label := factor(ancestry_label, levels = c(label0, label1, label2))]

bin_dt[, ancestry_label := factor(ancestry_label, levels = c(label0, label1, label2))]
bin_cc_dt[, ancestry_label := factor(ancestry_label, levels = c(label0, label1, label2))]
bin_cc_dt[, group := factor(group, levels = c("case", "control"))]

p <- ggplot(bin_dt, aes(x = ancestry_prop, y = prs, color = ancestry_label)) +
  geom_point(aes(size = n_samples), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = named_colors, drop = FALSE) +
  scale_size_continuous(name = "Samples/bin") +
  facet_wrap(~ ancestry_label, nrow = 1, scales = "free_x") +
  labs(
    title = paste0("PRS vs binned RFMix ancestry contribution - chr", chr_id),
    subtitle = subtitle_text,
    x = paste0("Ancestry contribution (.Q), binned (", ancestry_bins, " bins)"),
    y = "PRS",
    color = "Ancestry"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

p_cc <- ggplot(bin_cc_dt, aes(x = ancestry_prop, y = prs, color = ancestry_label)) +
  geom_point(aes(size = n_samples), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  geom_text(
    data = cc_annot_dt,
    aes(x = -Inf, y = Inf, label = annot_label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.1,
    size = 3,
    color = "black"
  ) +
  scale_color_manual(values = named_colors, drop = FALSE) +
  scale_size_continuous(name = "Samples/bin") +
  facet_grid(group ~ ancestry_label, scales = "free_x") +
  labs(
    title = paste0("PRS vs binned RFMix ancestry contribution by case/control - chr", chr_id),
    x = paste0("Ancestry contribution (.Q), binned (", ancestry_bins, " bins)"),
    y = "PRS",
    color = "Ancestry"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(tsv_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cc_png_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cc_pdf_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cc_tsv_out), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = png_out, plot = p, width = 14, height = 4.2, dpi = 300)
ggsave(filename = pdf_out, plot = p, width = 14, height = 4.2)
fwrite(corr_dt, file = tsv_out, sep = "\t")

ggsave(filename = cc_png_out, plot = p_cc, width = 14, height = 7.5, dpi = 300)
ggsave(filename = cc_pdf_out, plot = p_cc, width = 14, height = 7.5)
fwrite(corr_cc_dt, file = cc_tsv_out, sep = "\t")
