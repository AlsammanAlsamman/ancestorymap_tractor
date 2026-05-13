#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
    return(default)
  }
  if (idx[length(idx)] == length(args)) {
    stop(sprintf("Missing value for argument: %s", flag), call. = FALSE)
  }
  args[idx[length(idx)] + 1]
}

parse_bool <- function(value) {
  tolower(trimws(as.character(value))) %in% c("1", "true", "yes", "y")
}

extract_base_iid <- function(sample_id) {
  x <- as.character(sample_id)
  us <- gregexpr("_", x, fixed = TRUE)[[1]]
  if (length(us) > 0 && us[1] != -1) {
    for (p in us) {
      left <- substr(x, 1, p - 1)
      right <- substr(x, p + 1, nchar(x))
      if (left == right) return(left)
    }
  }
  x
}

run_cmd <- function(cmd, cmd_args) {
  res <- system2(cmd, cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste0(
        "Command failed: ", cmd, " ", paste(cmd_args, collapse = " "), "\n",
        paste(res, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  invisible(res)
}

chr_sort_key <- function(ch) {
  s <- gsub("^chr", "", as.character(ch))
  if (grepl("^[0-9]+$", s)) return(c(0, as.integer(s)))
  c(1, match(s, unique(s)))
}

main <- function() {
  ancestry <- get_arg("--ancestry", required = TRUE)
  chromosomes <- strsplit(get_arg("--chromosomes", required = TRUE), ",", fixed = TRUE)[[1]]
  chromosomes <- chromosomes[nzchar(chromosomes)]
  input_dir <- get_arg("--input-dir", required = TRUE)
  output_dir <- get_arg("--output-dir", required = TRUE)
  merged_prefix <- get_arg("--merged-prefix", required = TRUE)
  gwas_prefix <- get_arg("--gwas-prefix", required = TRUE)
  gwas_tsv <- get_arg("--gwas-tsv", required = TRUE)
  manhattan_png <- get_arg("--manhattan-png", required = TRUE)
  phenotype <- get_arg("--phenotype", required = TRUE)
  phenotype_format <- get_arg("--phenotype-format", required = TRUE)
  table_iid_column <- get_arg("--table-iid-column", default = "IID")
  table_phenotype_column <- get_arg("--table-phenotype-column", default = "PHENO")
  covariates_file <- get_arg("--covariates-file", required = TRUE)
  covariates_format <- get_arg("--covariates-format", default = "plink_eigenvec")
  covariates_has_header <- parse_bool(get_arg("--covariates-has-header", default = "false"))
  covariates_fid_column <- get_arg("--covariates-fid-column", default = "FID")
  covariates_iid_column <- get_arg("--covariates-iid-column", default = "IID")
  n_pcs <- as.integer(get_arg("--n-pcs", default = "5"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pfile_prefixes <- character()
  for (chr_name in chromosomes) {
    prefix <- file.path(input_dir, sprintf("chr%s.dosage", chr_name))
    if (file.exists(paste0(prefix, ".pgen")) && file.exists(paste0(prefix, ".pvar")) && file.exists(paste0(prefix, ".psam"))) {
      pfile_prefixes <- c(pfile_prefixes, prefix)
    }
  }
  if (length(pfile_prefixes) == 0) {
    stop(sprintf("No chromosome PLINK dosage files found in %s", input_dir), call. = FALSE)
  }

  if (length(pfile_prefixes) == 1) {
    run_cmd("plink2", c("--pfile", pfile_prefixes[1], "--allow-extra-chr", "--make-pgen", "--out", merged_prefix))
  } else {
    merge_list <- file.path(output_dir, "merge_list.txt")
    merge_lines <- vapply(
      pfile_prefixes[-1],
      function(pref) paste0(pref, ".pgen ", pref, ".pvar ", pref, ".psam"),
      FUN.VALUE = character(1)
    )
    writeLines(merge_lines, merge_list)
    run_cmd(
      "plink2",
      c("--pfile", pfile_prefixes[1], "--pmerge-list", merge_list, "--allow-extra-chr", "--make-pgen", "--out", merged_prefix)
    )
  }

  merged_psam <- paste0(merged_prefix, ".psam")
  psam <- fread(merged_psam, data.table = TRUE)
  iid_col <- if ("IID" %in% names(psam)) "IID" else if ("#IID" %in% names(psam)) "#IID" else NA_character_
  if (is.na(iid_col)) stop(sprintf("Could not find IID column in %s", merged_psam), call. = FALSE)
  merged_iids <- as.character(psam[[iid_col]])

  sample_map <- setNames(merged_iids, vapply(merged_iids, extract_base_iid, FUN.VALUE = character(1)))

  if (tolower(phenotype_format) == "fam") {
    fam <- fread(phenotype, header = FALSE, data.table = TRUE)
    if (ncol(fam) < 6) stop(sprintf("FAM file must have at least 6 columns: %s", phenotype), call. = FALSE)
    pheno <- data.table(IID_BASE = as.character(fam[[2]]), PHENO = as.character(fam[[6]]))
  } else {
    ptab <- fread(phenotype, data.table = TRUE)
    if (!(table_iid_column %in% names(ptab)) || !(table_phenotype_column %in% names(ptab))) {
      stop("Phenotype table missing configured IID/phenotype columns", call. = FALSE)
    }
    pheno <- data.table(IID_BASE = as.character(ptab[[table_iid_column]]), PHENO = as.character(ptab[[table_phenotype_column]]))
  }
  pheno[, IID := sample_map[IID_BASE]]
  pheno <- pheno[!is.na(IID) & nzchar(PHENO)]
  pheno[, FID := IID]
  pheno <- unique(pheno[, .(FID, IID, PHENO)], by = "IID")
  if (nrow(pheno) == 0) stop("No phenotype rows matched merged PLINK sample IDs", call. = FALSE)
  pheno_file <- file.path(output_dir, "phenotype_for_plink.tsv")
  fwrite(pheno, pheno_file, sep = "\t")

  if (tolower(covariates_format) != "plink_eigenvec") {
    stop(sprintf("Unsupported covariates format: %s", covariates_format), call. = FALSE)
  }
  if (covariates_has_header) {
    cov <- fread(covariates_file, data.table = TRUE)
  } else {
    cov <- fread(covariates_file, header = FALSE, data.table = TRUE)
    if (ncol(cov) < 3) stop("Covariates file must have at least FID IID PC1", call. = FALSE)
    setnames(
      cov,
      old = names(cov),
      new = c(covariates_fid_column, covariates_iid_column, paste0("PC", seq_len(ncol(cov) - 2)))
    )
  }
  pc_cols <- names(cov)[grepl("^PC", names(cov), ignore.case = TRUE)]
  pc_cols <- head(pc_cols, n_pcs)
  if (length(pc_cols) < n_pcs) {
    stop(sprintf("Requested %d PCs; found %d", n_pcs, length(pc_cols)), call. = FALSE)
  }
  cov_out <- data.table(IID_BASE = as.character(cov[[covariates_iid_column]]))
  cov_out[, IID := sample_map[IID_BASE]]
  cov_out <- cov_out[!is.na(IID)]
  cov_out[, FID := IID]
  for (pc in pc_cols) { 
    cov_tmp <- cov[, c(covariates_iid_column, pc), with = FALSE] 
    names(cov_tmp) <- c("IID_BASE", pc) 
    cov_out <- merge(cov_out, cov_tmp, by = "IID_BASE", all.x = TRUE) 
  }
  cov_out <- unique(cov_out[, c("FID", "IID", pc_cols), with = FALSE], by = "IID")
  if (nrow(cov_out) == 0) stop("No covariate rows matched merged PLINK sample IDs", call. = FALSE)
  cov_file <- file.path(output_dir, "covariates_for_plink.tsv")
  fwrite(cov_out, cov_file, sep = "\t")

  run_cmd(
    "plink2",
    c(
      "--pfile", merged_prefix,
      "--allow-extra-chr",
      "--pheno", pheno_file,
      "--pheno-name", "PHENO",
      "--covar", cov_file,
      "--covar-name", paste(pc_cols, collapse = ","),
      "--glm", "hide-covar",
      "--out", gwas_prefix
    )
  )

  glm_candidates <- Sys.glob(paste0(gwas_prefix, ".PHENO.glm.*"))
  if (length(glm_candidates) == 0) {
    stop(sprintf("No PLINK GLM output found for prefix: %s", gwas_prefix), call. = FALSE)
  }
  logistic_hits <- glm_candidates[grepl("\\.logistic", glm_candidates)]
  glm_file <- if (length(logistic_hits) > 0) logistic_hits[1] else glm_candidates[1]

  g <- fread(glm_file, data.table = TRUE)
  if ("TEST" %in% names(g)) g <- g[TEST == "ADD"]
  chr_col <- if ("#CHROM" %in% names(g)) "#CHROM" else if ("CHROM" %in% names(g)) "CHROM" else NA_character_
  if (is.na(chr_col) || !("POS" %in% names(g)) || !("ID" %in% names(g)) || !("P" %in% names(g))) {
    stop(sprintf("Unexpected GLM columns in %s", glm_file), call. = FALSE)
  }

  g[, P := as.numeric(P)]
  g <- g[!is.na(P) & P > 0]
  g[, CHR := gsub("^chr", "", as.character(get(chr_col)))]
  g[, POS := as.numeric(POS)]
  g <- g[!is.na(POS)]

  keep <- c("CHR", "POS", "ID", "P")
  extras <- c("A1", "A1_FREQ", "BETA", "SE", "OR", "Z_STAT", "T_STAT", "ERRCODE")
  keep <- c(keep, extras[extras %in% names(g)])
  g_out <- g[, ..keep]
  fwrite(g_out, gwas_tsv, sep = "\t")

  if (nrow(g_out) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No GWAS points to plot", size = 5) +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      ggtitle(sprintf("%s ancestry dosage GWAS", ancestry))
    ggsave(manhattan_png, p, width = 14, height = 6, dpi = 180)
    return(invisible(NULL))
  }

  chr_levels <- unique(g_out$CHR)
  chr_levels <- chr_levels[order(vapply(chr_levels, function(ch) chr_sort_key(ch)[1], numeric(1)),
                                 vapply(chr_levels, function(ch) chr_sort_key(ch)[2], numeric(1)))]

  offsets <- data.table(CHR = chr_levels)
  offsets[, chr_max := g_out[CHR == .BY$CHR, max(POS, na.rm = TRUE)], by = CHR]
  offsets[, offset := shift(cumsum(chr_max + 1e6), fill = 0)]
  g_plot <- merge(g_out, offsets[, .(CHR, offset, chr_max)], by = "CHR", all.x = TRUE)
  g_plot[, x := POS + offset]
  g_plot[, y := -log10(P)]
  g_plot[, chr_idx := as.integer(factor(CHR, levels = chr_levels))]

  centers <- offsets[, .(CHR, center = offset + chr_max / 2)]

  p <- ggplot(g_plot, aes(x = x, y = y, color = factor(chr_idx %% 2))) +
    geom_point(size = 0.6, alpha = 0.8) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#cc2f2f", linewidth = 0.4) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e"), guide = "none") +
    scale_x_continuous(breaks = centers$center, labels = centers$CHR) +
    labs(
      title = sprintf("%s ancestry dosage GWAS", ancestry),
      x = "Chromosome",
      y = "-log10(P)"
    ) +
    theme_bw(base_size = 12)

  ggsave(manhattan_png, p, width = 14, height = 6, dpi = 180)
}

main()
