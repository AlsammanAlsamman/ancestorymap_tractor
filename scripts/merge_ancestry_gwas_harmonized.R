#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
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

parse_csv <- function(x) {
  parts <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  trimws(parts[nzchar(trimws(parts))])
}

parse_map <- function(x) {
  m <- list()
  for (p in parse_csv(x)) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    m[[trimws(kv[1])]] <- trimws(kv[2])
  }
  m
}

as_chr <- function(x) gsub("^chr", "", as.character(x))

collapse_unique <- function(x, sep = ";") {
  y <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (length(y) == 0) return(NA_character_)
  paste(y, collapse = sep)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

harmonize_one <- function(dt, anc, ref_keys) {
  keys <- c("CHR", "POS", "ID")
  setkeyv(dt, keys)
  work <- merge(ref_keys, dt, by = keys, all.x = TRUE)

  work[, P_raw := safe_numeric(P)]
  if ("BETA" %in% names(work)) {
    work[, BETA_raw := safe_numeric(BETA)]
  } else {
    work[, BETA_raw := NA_real_]
  }
  if ("OR" %in% names(work)) {
    work[, OR_raw := safe_numeric(OR)]
  } else {
    work[, OR_raw := NA_real_]
  }
  work[is.na(BETA_raw) & !is.na(OR_raw) & OR_raw > 0, BETA_raw := log(OR_raw)]

  work[, ALT_FREQ_CASE := safe_numeric(ALT_FREQ_CASE)]
  work[, REF_FREQ_CASE := safe_numeric(REF_FREQ_CASE)]
  work[, ALT_FREQ_CTRL := safe_numeric(ALT_FREQ_CTRL)]
  work[, REF_FREQ_CTRL := safe_numeric(REF_FREQ_CTRL)]

  work[, aligned := !is.na(A1) & !is.na(EA) & as.character(A1) == as.character(EA)]

  work[, OR_h := fifelse(
    is.na(OR_raw),
    NA_real_,
    fifelse(aligned, OR_raw, fifelse(OR_raw > 0, 1 / OR_raw, NA_real_))
  )]
  work[, BETA_h := fifelse(is.na(BETA_raw), NA_real_, fifelse(aligned, BETA_raw, -BETA_raw))]

  work[, EA_FREQ_CASE_h := fifelse(
    is.na(EA),
    NA_real_,
    fifelse(as.character(EA) == as.character(ALT), ALT_FREQ_CASE,
            fifelse(as.character(EA) == as.character(REF), REF_FREQ_CASE, NA_real_))
  )]
  work[, EA_FREQ_CTRL_h := fifelse(
    is.na(EA),
    NA_real_,
    fifelse(as.character(EA) == as.character(ALT), ALT_FREQ_CTRL,
            fifelse(as.character(EA) == as.character(REF), REF_FREQ_CTRL, NA_real_))
  )]

  out <- work[, .(
    CHR, POS, ID,
    p = P_raw,
    or = OR_h,
    beta = BETA_h,
    ea_freq_case = EA_FREQ_CASE_h,
    ea_freq_ctrl = EA_FREQ_CTRL_h,
    n_case = safe_numeric(OBS_CT_CASE),
    n_ctrl = safe_numeric(OBS_CT_CTRL)
  )]

  setnames(out,
           old = c("p", "or", "beta", "ea_freq_case", "ea_freq_ctrl", "n_case", "n_ctrl"),
           new = c(
             paste0("P_", anc),
             paste0("OR_", anc),
             paste0("BETA_", anc),
             paste0("EA_FREQ_CASE_", anc),
             paste0("EA_FREQ_CTRL_", anc),
             paste0("N_CASE_", anc),
             paste0("N_CTRL_", anc)
           ))
  out
}

main <- function() {
  ancestries <- parse_csv(get_arg("--ancestries", required = TRUE))
  if (length(ancestries) == 0) stop("No ancestries provided", call. = FALSE)

  gwas_map <- parse_map(get_arg("--gwas-with-freq-map", required = TRUE))
  pvar_map <- parse_map(get_arg("--pvar-map", required = TRUE))
  summary_tsv <- get_arg("--summary-tsv", required = TRUE)
  reference_ancestry <- get_arg("--reference-ancestry", default = ancestries[1])
  out_tsv <- get_arg("--output-tsv", required = TRUE)
  out_xlsx <- get_arg("--output-xlsx", required = TRUE)

  for (a in ancestries) {
    if (is.null(gwas_map[[a]]) || is.null(pvar_map[[a]])) {
      stop(sprintf("Missing file mapping for ancestry: %s", a), call. = FALSE)
    }
  }

  if (!(reference_ancestry %in% ancestries)) {
    warning(sprintf("Reference ancestry %s not in ancestry list; using %s", reference_ancestry, ancestries[1]))
    reference_ancestry <- ancestries[1]
  }

  dts <- list()
  ref_keys <- NULL

  for (a in ancestries) {
    g <- fread(gwas_map[[a]], data.table = TRUE)
    pv <- fread(pvar_map[[a]], data.table = TRUE)

    chr_col <- if ("#CHROM" %in% names(pv)) "#CHROM" else if ("CHROM" %in% names(pv)) "CHROM" else NA_character_
    if (is.na(chr_col)) stop(sprintf("Could not find chromosome column in pvar: %s", pvar_map[[a]]), call. = FALSE)
    if (!("POS" %in% names(pv)) || !("ID" %in% names(pv)) || !("REF" %in% names(pv)) || !("ALT" %in% names(pv))) {
      stop(sprintf("Unexpected pvar columns in %s", pvar_map[[a]]), call. = FALSE)
    }

    pv <- pv[, .(
      CHR = as_chr(get(chr_col)),
      POS = safe_numeric(POS),
      ID = as.character(ID),
      REF = as.character(REF),
      ALT = as.character(ALT)
    )]

    if (!("CHR" %in% names(g)) || !("POS" %in% names(g)) || !("ID" %in% names(g))) {
      stop(sprintf("Unexpected GWAS columns in %s", gwas_map[[a]]), call. = FALSE)
    }

    g[, CHR := as_chr(CHR)]
    g[, POS := safe_numeric(POS)]
    g[, ID := as.character(ID)]

    dt <- merge(g, pv, by = c("CHR", "POS", "ID"), all.x = TRUE)
    dts[[a]] <- dt

    if (a == reference_ancestry) {
      ref_keys <- dt[, .(CHR, POS, ID, EA = as.character(A1), NEA = fifelse(as.character(A1) == as.character(ALT), as.character(REF), as.character(ALT)))]
      ref_keys <- unique(ref_keys, by = c("CHR", "POS", "ID"))
      ref_keys[is.na(EA) | !nzchar(EA), EA := ALT]
      ref_keys[is.na(NEA) | !nzchar(NEA), NEA := REF]
    }
  }

  if (is.null(ref_keys) || nrow(ref_keys) == 0) {
    stop("Could not determine reference allele keys", call. = FALSE)
  }

  merged <- copy(ref_keys)
  for (a in ancestries) {
    part <- harmonize_one(dts[[a]], a, ref_keys)
    merged <- merge(merged, part, by = c("CHR", "POS", "ID"), all.x = TRUE)
  }

  summary <- fread(summary_tsv, data.table = TRUE)
  if (!("chr" %in% names(summary)) || !("pos" %in% names(summary))) {
    stop(sprintf("Summary table missing chr/pos columns: %s", summary_tsv), call. = FALSE)
  }
  summary[, CHR := as_chr(chr)]
  summary[, POS := safe_numeric(pos)]
  summary[, nearest_gene := if ("gene" %in% names(summary)) as.character(gene) else NA_character_]
  summary[, locus_id := if ("locus_name" %in% names(summary)) as.character(locus_name) else NA_character_]
  summary_ann <- unique(summary[, .(CHR, POS, nearest_gene, locus_id)], by = c("CHR", "POS"))

  # Build per-locus genomic intervals from summary to map unlabeled SNPs by proximity.
  valid_locus_rows <- summary[!is.na(locus_id) & nzchar(locus_id)]
  locus_bounds <- valid_locus_rows[, .(
    locus_start = min(POS, na.rm = TRUE),
    locus_end = max(POS, na.rm = TRUE)
  ), by = .(CHR, locus_id)]
  locus_gene_map <- valid_locus_rows[, .(
    locus_nearest_genes = collapse_unique(nearest_gene, sep = ",")
  ), by = .(CHR, locus_id)]

  merged <- merge(merged, summary_ann, by = c("CHR", "POS"), all.x = TRUE)

  # For SNPs without direct locus_name, assign nearest locus if within 250 Kbp to locus start/end.
  missing_idx <- which(is.na(merged$locus_id) | !nzchar(merged$locus_id))
  if (length(missing_idx) > 0 && nrow(locus_bounds) > 0) {
    for (i in missing_idx) {
      chr_i <- as.character(merged$CHR[i])
      pos_i <- safe_numeric(merged$POS[i])
      if (is.na(pos_i) || !nzchar(chr_i)) next

      cand <- locus_bounds[CHR == chr_i]
      if (nrow(cand) == 0) next

      cand[, dist_bp := fifelse(
        pos_i < locus_start,
        locus_start - pos_i,
        fifelse(pos_i > locus_end, pos_i - locus_end, 0)
      )]

      min_dist <- suppressWarnings(min(cand$dist_bp, na.rm = TRUE))
      if (is.finite(min_dist) && !is.na(min_dist) && min_dist <= 250000) {
        nearest_row <- cand[which.min(dist_bp)]
        merged$locus_id[i] <- nearest_row$locus_id[1]

        if (is.na(merged$nearest_gene[i]) || !nzchar(merged$nearest_gene[i])) {
          g <- locus_gene_map[CHR == chr_i & locus_id == nearest_row$locus_id[1], locus_nearest_genes]
          if (length(g) > 0 && !is.na(g[1]) && nzchar(g[1])) merged$nearest_gene[i] <- g[1]
        }
      }
    }
  }

  merged[is.na(locus_id) | !nzchar(locus_id), locus_id := paste0("chr", CHR, ":", POS)]

  p_cols <- paste0("P_", ancestries)
  present_p <- p_cols[p_cols %in% names(merged)]
  merged[, P_MIN := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = present_p]
  merged[is.infinite(P_MIN), P_MIN := NA_real_]

  setcolorder(merged, c(
    "locus_id", "nearest_gene", "CHR", "POS", "EA", "NEA", "ID",
    setdiff(names(merged), c("locus_id", "nearest_gene", "CHR", "POS", "EA", "NEA", "ID"))
  ))

  fwrite(merged, out_tsv, sep = "\t")

  has_sig <- function(dt, threshold) {
    if (length(present_p) == 0) return(rep(FALSE, nrow(dt)))
    mat <- as.matrix(dt[, ..present_p])
    apply(mat, 1, function(r) any(!is.na(r) & r < threshold))
  }

  any_5e8 <- merged[has_sig(merged, 5e-8)]
  any_5e5 <- merged[has_sig(merged, 5e-5)]
  any_5e3 <- merged[has_sig(merged, 5e-3)]

  # Locus-level ancestry significance lists at each threshold
  locus_sig <- merged[, {
    anc_5e8 <- ancestries[vapply(ancestries, function(a) {
      col <- paste0("P_", a)
      col %in% names(.SD) && any(!is.na(.SD[[col]]) & .SD[[col]] < 5e-8)
    }, logical(1))]
    anc_5e5 <- ancestries[vapply(ancestries, function(a) {
      col <- paste0("P_", a)
      col %in% names(.SD) && any(!is.na(.SD[[col]]) & .SD[[col]] < 5e-5)
    }, logical(1))]
    anc_5e3 <- ancestries[vapply(ancestries, function(a) {
      col <- paste0("P_", a)
      col %in% names(.SD) && any(!is.na(.SD[[col]]) & .SD[[col]] < 5e-3)
    }, logical(1))]

    # Gene at the top SNP across ancestries: SNP with minimum P_MIN in locus.
    min_idx <- if (all(is.na(P_MIN))) integer(0) else which(P_MIN == min(P_MIN, na.rm = TRUE))
    top_snp_gene <- if (length(min_idx) == 0) {
      NA_character_
    } else {
      collapse_unique(nearest_gene[min_idx], sep = ",")
    }
    all_genes <- collapse_unique(nearest_gene, sep = ",")

    list(
      significant_at_5e8 = collapse_unique(anc_5e8, sep = ","),
      significant_at_5e5 = collapse_unique(anc_5e5, sep = ","),
      significant_at_5e3 = collapse_unique(anc_5e3, sep = ","),
      all_genes_in_locus = all_genes,
      top_snp_gene_across_ancestries = top_snp_gene
    )
  }, by = .(locus_id)]

  # Loci significant at 5e-8 in at least one ancestry but not all ancestries
  locus_5e8_matrix <- merged[, {
    vals <- vapply(ancestries, function(a) {
      col <- paste0("P_", a)
      col %in% names(.SD) && any(!is.na(.SD[[col]]) & .SD[[col]] < 5e-8)
    }, logical(1))
    as.list(vals)
  }, by = .(locus_id)]
  setnames(locus_5e8_matrix, old = names(locus_5e8_matrix)[-1], new = paste0("sig5e8_", ancestries))

  sig_cols <- grep("^sig5e8_", names(locus_5e8_matrix), value = TRUE)
  sig_mat <- as.matrix(locus_5e8_matrix[, ..sig_cols])
  sig_count <- rowSums(sig_mat, na.rm = TRUE)
  loci_specific_5e8 <- copy(locus_5e8_matrix[sig_count >= 1 & sig_count < length(ancestries)])
  setcolorder(loci_specific_5e8, c("locus_id", setdiff(names(loci_specific_5e8), "locus_id")))
  loci_specific_5e8 <- merge(
    loci_specific_5e8,
    unique(locus_sig[, .(locus_id, all_genes_in_locus, top_snp_gene_across_ancestries)], by = "locus_id"),
    by = "locus_id",
    all.x = TRUE
  )

  wb <- createWorkbook()
  addWorksheet(wb, "snps_p_lt_5e8_any")
  if (!is.null(any_5e8) && nrow(any_5e8) > 0) writeData(wb, "snps_p_lt_5e8_any", any_5e8)

  addWorksheet(wb, "snps_p_lt_5e5_any")
  if (!is.null(any_5e5) && nrow(any_5e5) > 0) writeData(wb, "snps_p_lt_5e5_any", any_5e5)

  addWorksheet(wb, "snps_p_lt_5e3_any")
  if (!is.null(any_5e3) && nrow(any_5e3) > 0) writeData(wb, "snps_p_lt_5e3_any", any_5e3)

  addWorksheet(wb, "loci_5e8_not_all")
  if (!is.null(loci_specific_5e8) && nrow(loci_specific_5e8) > 0) writeData(wb, "loci_5e8_not_all", loci_specific_5e8)

  addWorksheet(wb, "locus_summary")
  if (!is.null(locus_sig) && nrow(locus_sig) > 0) writeData(wb, "locus_summary", locus_sig)

  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
}

main()
