#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
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

main <- function() {
  gwas_path <- get_arg("--gwas-tsv", required = TRUE)
  case_path <- get_arg("--case-afreq", required = TRUE)
  ctrl_path <- get_arg("--ctrl-afreq", required = TRUE)
  out_path <- get_arg("--output-tsv", required = TRUE)

  gwas <- fread(gwas_path, data.table = TRUE)
  if (!("ID" %in% names(gwas))) {
    stop(sprintf("GWAS table missing ID column: %s", gwas_path), call. = FALSE)
  }

  cases <- fread(case_path, data.table = TRUE)
  ctrls <- fread(ctrl_path, data.table = TRUE)

  required_cols <- c("ID", "ALT", "ALT_FREQS", "OBS_CT")
  for (nm in required_cols) {
    if (!(nm %in% names(cases))) stop(sprintf("Cases afreq missing %s: %s", nm, case_path), call. = FALSE)
    if (!(nm %in% names(ctrls))) stop(sprintf("Controls afreq missing %s: %s", nm, ctrl_path), call. = FALSE)
  }

  cases <- cases[, .(
    ID = as.character(ID),
    ALT_CASE = as.character(ALT),
    ALT_FREQ_CASE = as.numeric(ALT_FREQS),
    REF_FREQ_CASE = 1 - as.numeric(ALT_FREQS),
    OBS_CT_CASE = as.numeric(OBS_CT)
  )]

  ctrls <- ctrls[, .(
    ID = as.character(ID),
    ALT_CTRL = as.character(ALT),
    ALT_FREQ_CTRL = as.numeric(ALT_FREQS),
    REF_FREQ_CTRL = 1 - as.numeric(ALT_FREQS),
    OBS_CT_CTRL = as.numeric(OBS_CT)
  )]

  out <- merge(gwas, cases, by = "ID", all.x = TRUE)
  out <- merge(out, ctrls, by = "ID", all.x = TRUE)

  if ("A1" %in% names(out)) {
    out[, A1_MATCHES_ALT_CASE := fifelse(!is.na(ALT_CASE), as.character(A1) == ALT_CASE, NA)]
    out[, A1_MATCHES_ALT_CTRL := fifelse(!is.na(ALT_CTRL), as.character(A1) == ALT_CTRL, NA)]
  }

  setcolorder(out, c(
    "CHR", "POS", "ID",
    intersect(c("A1", "OR", "BETA", "SE", "P", "Z_STAT", "T_STAT", "ERRCODE"), names(out)),
    "ALT_CASE", "ALT_FREQ_CASE", "REF_FREQ_CASE", "OBS_CT_CASE",
    "ALT_CTRL", "ALT_FREQ_CTRL", "REF_FREQ_CTRL", "OBS_CT_CTRL",
    intersect(c("A1_MATCHES_ALT_CASE", "A1_MATCHES_ALT_CTRL"), names(out))
  ))

  fwrite(out, out_path, sep = "\t")
}

main()
