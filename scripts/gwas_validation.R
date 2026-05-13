#!/usr/bin/env Rscript
# GWAS Ancestry Validation
# Classifies loci by ancestry significance, matches against external validation
# GWAS files, and produces a formatted Excel report.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# ─── Argument parsing ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
    return(default)
  }
  if (idx[length(idx)] == length(args))
    stop(sprintf("Missing value for argument: %s", flag), call. = FALSE)
  args[idx[length(idx)] + 1]
}

parse_csv <- function(x) {
  parts <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  trimws(parts[nzchar(trimws(parts))])
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

normalize_snp_id <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[is.na(y)] <- ""

  # If IDs come as chr:pos:rsid, extract the rs token.
  has_colon <- grepl(":", y, fixed = TRUE)
  if (any(has_colon)) {
    y[has_colon] <- vapply(strsplit(y[has_colon], ":", fixed = TRUE), function(parts) {
      rs_parts <- parts[grepl("^rs[0-9]+", parts)]
      if (length(rs_parts) > 0) rs_parts[length(rs_parts)] else parts[length(parts)]
    }, character(1))
  }

  # Trim common rsid decorations, e.g. rs12345.1 -> rs12345.
  has_rs <- grepl("rs[0-9]+", y)
  y[has_rs] <- sub("^.*?(rs[0-9]+).*$", "\\1", y[has_rs])
  y
}

first_present_col <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

harmonized_tsv <- get_arg("--harmonized-tsv",   required = TRUE)
out_xlsx       <- get_arg("--output-excel",      required = TRUE)
ancestries     <- parse_csv(get_arg("--ancestries",     required = TRUE))
val_ancestries <- parse_csv(get_arg("--val-ancestries", required = TRUE))
val_files      <- parse_csv(get_arg("--val-files",      required = TRUE))

if (length(val_ancestries) != length(val_files))
  stop("--val-ancestries and --val-files must have the same number of entries", call. = FALSE)

# ─── 1. Load harmonized data ──────────────────────────────────────────────────
message("Loading harmonized TSV: ", harmonized_tsv)
dt <- fread(harmonized_tsv, data.table = TRUE)

p_cols <- paste0("P_", ancestries)
p_cols <- p_cols[p_cols %in% names(dt)]
if (length(p_cols) == 0)
  stop("No P-value columns found for the specified ancestries in harmonized TSV", call. = FALSE)

dt[, (p_cols) := lapply(.SD, safe_num), .SDcols = p_cols]

# ─── 2. Classify loci ─────────────────────────────────────────────────────────
# Step A: for each locus, which ancestries have any SNP at P<5E-8?
# Step B: if none at 5E-8, which have P<5E-5?
# Step C: if none at 5E-5, which have P<5E-4?
# Label = sorted ancestries joined by "+"

classify_loci <- function(threshold) {
  dt[, lapply(
    setNames(p_cols, ancestries),
    function(pc) any(!is.na(get(pc)) & get(pc) < threshold)
  ), by = locus_id]
}

locus_5e8 <- classify_loci(5e-8)
locus_5e5 <- classify_loci(5e-5)
locus_5e4 <- classify_loci(5e-4)

locus_class <- data.table(locus_id = locus_5e8$locus_id)
keep_anc_cols <- intersect(ancestries, names(locus_5e8))  # ancestries present as cols

locus_class[, has_5e8 := apply(as.matrix(locus_5e8[, keep_anc_cols, with = FALSE]), 1, any)]
locus_class[, has_5e5 := apply(as.matrix(locus_5e5[, keep_anc_cols, with = FALSE]), 1, any)]
locus_class[, has_5e4 := apply(as.matrix(locus_5e4[, keep_anc_cols, with = FALSE]), 1, any)]

locus_class[, locus_ancestry_label := {
  mapply(function(lid, use5e8, use5e5, use5e4) {
    if (use5e8) {
      sigs <- keep_anc_cols[as.logical(unlist(locus_5e8[locus_id == lid, keep_anc_cols, with = FALSE]))]
    } else if (use5e5) {
      sigs <- keep_anc_cols[as.logical(unlist(locus_5e5[locus_id == lid, keep_anc_cols, with = FALSE]))]
    } else if (use5e4) {
      sigs <- keep_anc_cols[as.logical(unlist(locus_5e4[locus_id == lid, keep_anc_cols, with = FALSE]))]
    } else {
      sigs <- character(0)
    }
    if (length(sigs) == 0) "NONE" else paste(sort(sigs), collapse = "+")
  }, locus_id, has_5e8, has_5e5, has_5e4)
}]

locus_class[, classification_threshold := fifelse(
  has_5e8, "5E-8",
  fifelse(has_5e5, "5E-5", fifelse(has_5e4, "5E-4", "NONE"))
)]

sig_locus_class <- locus_class[locus_ancestry_label != "NONE"]
message(sprintf("Significant loci: %d (at 5E-8: %d, 5E-5 only: %d, 5E-4 only: %d)",
  nrow(sig_locus_class),
  sum(sig_locus_class$has_5e8),
  sum(!sig_locus_class$has_5e8 & sig_locus_class$has_5e5),
  sum(!sig_locus_class$has_5e8 & !sig_locus_class$has_5e5 & sig_locus_class$has_5e4)
))

# ─── 3. Load validation GWAS files ────────────────────────────────────────────
val_data <- list()
for (i in seq_along(val_ancestries)) {
  anc   <- val_ancestries[i]
  fpath <- val_files[i]
  if (!file.exists(fpath)) {
    warning(sprintf("Validation file not found for %s: %s — skipping", anc, fpath))
    next
  }
  message("Loading validation GWAS for ", anc, ": ", fpath)
  vd <- fread(fpath, data.table = TRUE)
  setnames(vd, tolower(trimws(names(vd))))   # normalise column names to lowercase

  rsid_col <- first_present_col(names(vd), c("rsid", "snpid", "snp", "markername", "rs_number", "rs"))
  pval_col <- first_present_col(names(vd), c("pvalue", "p", "pval", "p_value"))
  beta_col <- first_present_col(names(vd), c("beta", "b", "effect", "effectsize"))
  or_col   <- first_present_col(names(vd), c("or", "oddsratio", "odds_ratio"))

  if (is.na(rsid_col)) {
    stop(sprintf("Validation file for %s is missing SNP id column (checked: rsid/snpid/snp/markername)", anc), call. = FALSE)
  }
  if (is.na(pval_col)) {
    stop(sprintf("Validation file for %s is missing p-value column (checked: pvalue/p/pval/p_value)", anc), call. = FALSE)
  }

  vd[, rsid_norm := as.character(get(rsid_col))]
  vd[, id_match := normalize_snp_id(rsid_norm)]
  vd[, pvalue_norm := safe_num(get(pval_col))]
  if (!is.na(or_col)) {
    vd[, OR_val := safe_num(get(or_col))]
  } else if (!is.na(beta_col)) {
    vd[, OR_val := exp(safe_num(get(beta_col)))]
  } else {
    vd[, OR_val := NA_real_]
  }

  # Keep one row per normalized ID: the smallest p-value.
  vd <- vd[!is.na(id_match) & nzchar(id_match)]
  setorderv(vd, c("id_match", "pvalue_norm"), c(1, 1), na.last = TRUE)
  vd <- vd[!duplicated(id_match)]
  val_data[[anc]] <- vd[, .(rsid = rsid_norm, id_match, pvalue = pvalue_norm, OR_val)]
}

# ─── 4. Build main sheet: all SNPs from sig loci + validation columns ─────────
sig_loci_ids <- sig_locus_class$locus_id
main_dt <- dt[locus_id %in% sig_loci_ids]
main_dt[, id_match := normalize_snp_id(ID)]

# Attach locus classification columns
main_dt <- merge(
  main_dt,
  sig_locus_class[, .(locus_id, locus_ancestry_label, classification_threshold)],
  by = "locus_id", all.x = TRUE
)

# Attach validation P / OR / top-SNP columns
for (anc in val_ancestries) {
  p_v   <- paste0("P_",   anc, "_val")
  or_v  <- paste0("OR_",  anc, "_val")
  top_v <- paste0("top_", anc)

  main_dt[[p_v]]   <- NA_real_
  main_dt[[or_v]]  <- NA_real_
  main_dt[[top_v]] <- FALSE

  if (!(anc %in% names(val_data))) next
  vd <- val_data[[anc]]

  m_idx <- match(main_dt$id_match, vd$id_match)
  main_dt[[p_v]]  <- vd$pvalue[m_idx]
  main_dt[[or_v]] <- vd$OR_val[m_idx]

  # One top SNP per locus (lowest pvalue among matched SNPs in validation GWAS)
  for (lid in sig_loci_ids) {
    rows <- which(main_dt$locus_id == lid & !is.na(main_dt[[p_v]]))
    if (length(rows) == 0) next
    best <- rows[which.min(main_dt[[p_v]][rows])]
    main_dt[[top_v]][best] <- TRUE
  }
}

# ─── 5. Compute validation_pass per locus ─────────────────────────────────────
# Expected sig  = ancestries in the label that have a validation file → must have ≥1 SNP P<5E-8
# Unexpected    = validation ancestries NOT in the label              → must have 0 SNPs P<5E-8
# AFR (no val file): ignored, NA if label contains ONLY AFR

compute_val_pass <- function(lid, label) {
  label_ancs   <- strsplit(label, "+", fixed = TRUE)[[1]]
  expected     <- intersect(label_ancs,    val_ancestries)
  unexpected   <- setdiff(val_ancestries,  label_ancs)

  if (length(expected) == 0) return(NA)   # label only contains ancestries without val files

  locus_rows <- main_dt[locus_id == lid]

  for (anc in expected) {
    p_v <- paste0("P_", anc, "_val")
    if (!(p_v %in% names(locus_rows))) return(NA)
    if (!any(!is.na(locus_rows[[p_v]]) & locus_rows[[p_v]] < 5e-8)) return(FALSE)
  }
  for (anc in unexpected) {
    p_v <- paste0("P_", anc, "_val")
    if (!(p_v %in% names(locus_rows))) next
    if (any(!is.na(locus_rows[[p_v]]) & locus_rows[[p_v]] < 5e-8)) return(FALSE)
  }
  TRUE
}

val_pass_map <- sig_locus_class[, .(
  locus_id,
  validation_pass = mapply(compute_val_pass, locus_id, locus_ancestry_label)
)]
main_dt <- merge(main_dt, val_pass_map, by = "locus_id", all.x = TRUE)

message(sprintf("Loci passing validation: %d / %d",
  sum(val_pass_map$validation_pass %in% TRUE),
  nrow(val_pass_map)
))

# ─── 6. Build summary sheet ───────────────────────────────────────────────────
locus_genes <- dt[locus_id %in% sig_loci_ids,
  .(genes = paste(sort(unique(nearest_gene[!is.na(nearest_gene) & nzchar(nearest_gene)])),
                  collapse = ",")),
  by = locus_id]

summary_dt <- merge(sig_locus_class[, .(locus_id, locus_ancestry_label, classification_threshold)],
                    locus_genes, by = "locus_id", all.x = TRUE)
summary_dt <- merge(summary_dt, val_pass_map, by = "locus_id", all.x = TRUE)

for (anc in val_ancestries) {
  p_v    <- paste0("P_",   anc, "_val")
  col_nm <- paste0("sig_", anc, "_in_val_gwas")
  summary_dt[[col_nm]] <- vapply(summary_dt$locus_id, function(lid) {
    rows <- main_dt[locus_id == lid]
    if (!(p_v %in% names(rows))) return(NA)
    any(!is.na(rows[[p_v]]) & rows[[p_v]] < 5e-8)
  }, logical(1))
}

# ─── 7. Write Excel ───────────────────────────────────────────────────────────
dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

write_sheet <- function(ws_name, data_tbl) {
  addWorksheet(wb, ws_name)
  writeData(wb, ws_name, data_tbl)

  hdr_style <- createStyle(textRotation = 90, textDecoration = "bold",
                            halign = "center", valign = "center")
  addStyle(wb, ws_name, hdr_style,
           rows = 1, cols = seq_len(ncol(data_tbl)), gridExpand = TRUE)
  setRowHeights(wb, ws_name, rows = 1, heights = 120)

  if (nrow(data_tbl) == 0) return(invisible(NULL))
  data_rows <- 2:(nrow(data_tbl) + 1)

  # P-value columns → scientific notation + light red when < 5E-8
  all_p_nm <- grep("^P_", names(data_tbl), value = TRUE)
  p_idx    <- which(names(data_tbl) %in% all_p_nm)
  if (length(p_idx) > 0) {
    addStyle(wb, ws_name, createStyle(numFmt = "0.00E+00"),
             rows = data_rows, cols = p_idx, gridExpand = TRUE, stack = TRUE)
    red_sty <- createStyle(fgFill = "#FADBD8")
    for (ci in p_idx) {
      conditionalFormatting(wb, ws_name, cols = ci, rows = data_rows,
                            style = red_sty,
                            rule = sprintf("< %g", 5e-8), type = "expression")
    }
  }

  # OR columns → 2 decimal places
  or_idx <- grep("^OR_", names(data_tbl))
  if (length(or_idx) > 0) {
    addStyle(wb, ws_name, createStyle(numFmt = "0.00"),
             rows = data_rows, cols = or_idx, gridExpand = TRUE, stack = TRUE)
  }
}

write_sheet("gwas_validation", main_dt)
write_sheet("summary",         summary_dt)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx)
