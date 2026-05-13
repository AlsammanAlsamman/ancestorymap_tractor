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

parse_ancestry_map <- function(value) {
  mapping <- list()
  parts <- strsplit(as.character(value), ",", fixed = TRUE)[[1]]
  for (part in parts) {
    if (!nzchar(part)) next
    kv <- strsplit(part, ":", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    mapping[[trimws(kv[1])]] <- trimws(kv[2])
  }
  mapping
}

add_sample_id_aliases_multi <- function(dt, iid_col, value_columns, fid_col = NULL) {
  value_columns <- as.character(value_columns)
  base <- copy(dt[, c(iid_col, value_columns), with = FALSE])
  setnames(base, iid_col, "sample_id")
  base[, sample_id := as.character(sample_id)]
  alias_list <- list(base)

  iid_iid <- copy(base)
  iid_iid[, sample_id := paste0(sample_id, "_", sample_id)]
  alias_list[[length(alias_list) + 1]] <- iid_iid

  if (!is.null(fid_col) && fid_col %in% names(dt)) {
    fid_iid <- copy(dt[, c(fid_col, iid_col, value_columns), with = FALSE])
    fid_iid[, sample_id := paste0(as.character(get(fid_col)), "_", as.character(get(iid_col)))]
    fid_iid[, c(fid_col, iid_col) := NULL]
    setcolorder(fid_iid, c("sample_id", value_columns))
    alias_list[[length(alias_list) + 1]] <- fid_iid
  }

  alias_dt <- unique(rbindlist(alias_list, fill = TRUE), by = "sample_id")
  alias_dt
}

load_phenotype <- function(args_list) {
  if (tolower(args_list$phenotype_format) == "fam") {
    fam <- fread(args_list$phenotype, header = FALSE, data.table = TRUE)
    if (ncol(fam) < 6) {
      stop(sprintf("FAM file must have at least 6 columns; found %d in %s", ncol(fam), args_list$phenotype), call. = FALSE)
    }

    phe <- data.table(
      FID = as.character(fam[[1]]),
      IID = as.character(fam[[2]]),
      SEX = as.character(fam[[5]]),
      PHENO = as.character(fam[[6]])
    )
    phe[, y := fifelse(PHENO == args_list$fam_case_code, 1.0,
                       fifelse(PHENO == args_list$fam_control_code, 0.0, NA_real_))]
    if (parse_bool(args_list$include_sex)) {
      phe[, sex_female := fifelse(SEX == "2", 1.0,
                                  fifelse(SEX == "1", 0.0, NA_real_))]
    }

    keep_cols <- c("y")
    if (parse_bool(args_list$include_sex) && "sex_female" %in% names(phe)) {
      keep_cols <- c(keep_cols, "sex_female")
    }

    alias <- add_sample_id_aliases_multi(phe, iid_col = "IID", fid_col = "FID", value_columns = keep_cols)
    alias[, IID := vapply(sample_id, extract_base_iid, FUN.VALUE = character(1))]
    return(alias)
  }

  phe <- fread(args_list$phenotype, data.table = TRUE)
  if (!(args_list$table_iid_column %in% names(phe))) {
    stop(sprintf("IID column '%s' not found in %s", args_list$table_iid_column, args_list$phenotype), call. = FALSE)
  }
  if (!(args_list$table_phenotype_column %in% names(phe))) {
    stop(sprintf("Phenotype column '%s' not found in %s", args_list$table_phenotype_column, args_list$phenotype), call. = FALSE)
  }

  phe[, y := as.numeric(get(args_list$table_phenotype_column))]
  keep_cols <- c("y")
  if (parse_bool(args_list$include_sex) && "SEX" %in% names(phe)) {
    phe[, sex_female := fifelse(as.character(SEX) == "2", 1.0,
                                fifelse(as.character(SEX) == "1", 0.0, NA_real_))]
    keep_cols <- c(keep_cols, "sex_female")
  }
  alias <- add_sample_id_aliases_multi(
    phe,
    iid_col = args_list$table_iid_column,
    fid_col = if ("FID" %in% names(phe)) "FID" else NULL,
    value_columns = keep_cols
  )
  alias[, IID := vapply(sample_id, extract_base_iid, FUN.VALUE = character(1))]
  alias
}

load_covariates <- function(args_list) {
  if (!nzchar(args_list$covariates_file)) {
    return(list(data = NULL, pc_columns = character()))
  }

  if (tolower(args_list$covariates_format) != "plink_eigenvec") {
    stop(sprintf("Unsupported covariates format: %s", args_list$covariates_format), call. = FALSE)
  }

  has_header <- parse_bool(args_list$covariates_has_header)
  cov_dt <- fread(args_list$covariates_file, header = has_header, data.table = TRUE)
  if (!has_header) {
    if (ncol(cov_dt) < 3) {
      stop(sprintf(
        "Covariates file must contain at least FID, IID, and one PC column; found %d in %s",
        ncol(cov_dt), args_list$covariates_file
      ), call. = FALSE)
    }
    setnames(
      cov_dt,
      old = names(cov_dt),
      new = c(args_list$covariates_fid_column, args_list$covariates_iid_column,
              paste0("PC", seq_len(ncol(cov_dt) - 2)))
    )
  }

  pc_columns <- names(cov_dt)[grepl("^PC", names(cov_dt), ignore.case = TRUE)]
  pc_columns <- pc_columns[seq_len(min(length(pc_columns), args_list$n_pcs))]
  if (length(pc_columns) < args_list$n_pcs) {
    stop(sprintf(
      "Requested %d PCs, but only found %d in %s",
      args_list$n_pcs, length(pc_columns), args_list$covariates_file
    ), call. = FALSE)
  }

  for (col in pc_columns) cov_dt[, (col) := as.numeric(get(col))]

  alias <- add_sample_id_aliases_multi(
    cov_dt,
    iid_col = args_list$covariates_iid_column,
    fid_col = if (args_list$covariates_fid_column %in% names(cov_dt)) args_list$covariates_fid_column else NULL,
    value_columns = pc_columns
  )
  alias <- unique(alias, by = "sample_id")
  list(data = alias, pc_columns = pc_columns)
}

load_age <- function(args_list) {
  if (!parse_bool(args_list$include_age) || !nzchar(args_list$age_file)) {
    return(NULL)
  }
  if (!file.exists(args_list$age_file)) {
    message(sprintf("Age file was requested but not found: %s; age will be skipped.", args_list$age_file))
    return(NULL)
  }

  age_dt <- fread(args_list$age_file, data.table = TRUE)
  if (!(args_list$age_iid_column %in% names(age_dt)) || !(args_list$age_column %in% names(age_dt))) {
    stop(sprintf(
      "Age file must contain columns '%s' and '%s'.",
      args_list$age_iid_column, args_list$age_column
    ), call. = FALSE)
  }
  age_dt[, age := as.numeric(get(args_list$age_column))]
  alias <- add_sample_id_aliases_multi(
    age_dt,
    iid_col = args_list$age_iid_column,
    fid_col = if ("FID" %in% names(age_dt)) "FID" else NULL,
    value_columns = "age"
  )
  unique(alias, by = "sample_id")
}

load_global_ancestry <- function(args_list, ancestry_map, chromosomes) {
  frames <- list()
  for (chr_name in chromosomes) {
    path <- gsub("{chr}", chr_name, args_list$rfmix_q_pattern, fixed = TRUE)
    header_lines <- readLines(path, n = 2, warn = FALSE)
    if (length(header_lines) < 2) {
      stop(sprintf("RFMix Q file %s does not contain the expected two-line header", path), call. = FALSE)
    }

    header_fields <- strsplit(sub("^#", "", trimws(header_lines[2])), "[[:space:]]+")[[1]]
    q_dt <- fread(path, skip = 2, header = FALSE, data.table = TRUE)
    if (ncol(q_dt) != length(header_fields)) {
      stop(sprintf(
        "RFMix Q file %s header/data column mismatch: header has %d fields but data has %d columns",
        path, length(header_fields), ncol(q_dt)
      ), call. = FALSE)
    }

    setnames(q_dt, header_fields)
    sample_col <- if ("sample" %in% names(q_dt)) "sample" else names(q_dt)[1]
    setnames(q_dt, old = sample_col, new = "sample_id")

    keep_cols <- c("sample_id")
    rename_old <- character()
    rename_new <- character()
    for (idx in names(ancestry_map)) {
      if (idx %in% names(q_dt)) {
        keep_cols <- c(keep_cols, idx)
        rename_old <- c(rename_old, idx)
        rename_new <- c(rename_new, paste0("global_", ancestry_map[[idx]], "_prop"))
      }
    }

    if (length(rename_old) == 0) {
      stop(sprintf("Failed to parse any ancestry-proportion columns from %s", path), call. = FALSE)
    }

    q_dt <- q_dt[, ..keep_cols]
    setnames(q_dt, rename_old, rename_new)
    for (col in setdiff(names(q_dt), "sample_id")) q_dt[, (col) := as.numeric(get(col))]
    frames[[length(frames) + 1]] <- q_dt
  }

  global_dt <- rbindlist(frames, fill = TRUE)
  value_cols <- setdiff(names(global_dt), "sample_id")
  if (length(value_cols) == 0) {
    stop("No global ancestry columns were retained from the RFMix .Q files.", call. = FALSE)
  }
  global_dt <- global_dt[, lapply(.SD, mean, na.rm = TRUE), by = sample_id, .SDcols = value_cols]
  global_dt
}

select_snps <- function(summary_tsv, ancestry_labels, p_threshold, max_snps_per_ancestry, use_top_snp_per_locus) {
  summary_dt <- fread(summary_tsv, data.table = TRUE)
  summary_dt[, chr := as.character(chr)]
  summary_dt[, pos := as.integer(pos)]

  if (!("locus_name" %in% names(summary_dt))) summary_dt[, locus_name := NA_character_]
  if (!("gene" %in% names(summary_dt))) summary_dt[, gene := NA_character_]
  summary_dt[, locus_name := as.character(locus_name)]
  summary_dt[, gene := as.character(gene)]
  summary_dt[is.na(locus_name) | !nzchar(trimws(locus_name)), locus_name := gene]
  summary_dt[is.na(locus_name) | !nzchar(trimws(locus_name)), locus_name := paste(chr, pos, rsid, sep = ":")]

  out <- list()
  for (label in ancestry_labels) {
    p_col <- paste0("p_dosage_", label)
    or_col <- paste0("OR_dosage_", label)
    if (!(p_col %in% names(summary_dt)) || !(or_col %in% names(summary_dt))) {
      stop(sprintf("Required columns missing from summary table: %s and/or %s", p_col, or_col), call. = FALSE)
    }

    sub_dt <- summary_dt[, .(
      chr,
      pos,
      rsid,
      locus_name = if ("locus_name" %in% names(summary_dt)) as.character(locus_name) else NA_character_,
      gene = if ("gene" %in% names(summary_dt)) as.character(gene) else NA_character_,
      p_value = as.numeric(get(p_col)),
      odds_ratio = as.numeric(get(or_col))
    )]
    sub_dt <- sub_dt[!is.na(pos) & !is.na(p_value) & !is.na(odds_ratio) & is.finite(odds_ratio) & odds_ratio > 0]
    setorder(sub_dt, p_value, pos)

    filtered <- sub_dt[p_value <= p_threshold]
    selection_mode <- if (nrow(filtered) > 0) "threshold" else "fallback_top_hits"
    if (nrow(filtered) == 0) filtered <- copy(sub_dt)

    if (use_top_snp_per_locus && "locus_name" %in% names(filtered)) {
      filtered <- filtered[!duplicated(locus_name)]
    } else {
      filtered <- filtered[!duplicated(paste(chr, pos, rsid, sep = ":"))]
    }

    filtered <- head(filtered, max_snps_per_ancestry)
    filtered[, ancestry := label]
    filtered[, selection_mode := selection_mode]
    filtered[, log_or_weight := log(odds_ratio)]
    filtered[, variant_key := paste(chr, pos, rsid, sep = ":")]
    filtered[, found_in_dosage := FALSE]
    out[[length(out) + 1]] <- filtered
  }

  rbindlist(out, fill = TRUE)
}

load_sample_ids_from_first_dosage <- function(args_list, ancestry_map, chromosomes) {
  first_anc <- names(ancestry_map)[order(as.integer(names(ancestry_map)))][1]
  first_chr <- chromosomes[1]
  path <- gsub("{chr}", first_chr, args_list$tractor_dosage_pattern, fixed = TRUE)
  path <- gsub("{anc}", first_anc, path, fixed = TRUE)
  header <- fread(path, nrows = 0, data.table = TRUE, check.names = FALSE)
  names(header)[-(1:5)]
}

compute_dosage_scores <- function(args_list, selected_dt, ancestry_map, chromosomes) {
  sample_ids <- load_sample_ids_from_first_dosage(args_list, ancestry_map, chromosomes)
  result_dt <- data.table(sample_id = sample_ids)
  found_keys <- character()
  ancestry_by_label <- setNames(names(ancestry_map), unlist(ancestry_map))

  for (label in names(ancestry_by_label)) {
    feature_name <- paste0("hap_dosage_", label)
    accum <- rep(0, length(sample_ids))
    label_selected <- selected_dt[ancestry == label]
    if (nrow(label_selected) == 0) {
      result_dt[, (feature_name) := accum]
      next
    }

    anc_idx <- ancestry_by_label[[label]]
    for (chr_name in chromosomes) {
      chr_selected <- label_selected[chr == as.character(chr_name)]
      if (nrow(chr_selected) == 0) next

      dosage_path <- gsub("{chr}", chr_name, args_list$tractor_dosage_pattern, fixed = TRUE)
      dosage_path <- gsub("{anc}", anc_idx, dosage_path, fixed = TRUE)

      ids_file <- tempfile(fileext = ".txt")
      writeLines(unique(as.character(chr_selected$rsid)), ids_file)
      cmd <- sprintf("awk 'NR==FNR {ids[$1]=1; next} FNR==1 || ($3 in ids)' %s %s",
                     shQuote(ids_file), shQuote(dosage_path))
      dosage_dt <- tryCatch(
        fread(cmd = cmd, data.table = TRUE, check.names = FALSE, showProgress = FALSE),
        error = function(e) NULL
      )
      unlink(ids_file)

      if (is.null(dosage_dt) || nrow(dosage_dt) == 0) next
      if (ncol(dosage_dt) <= 5) next

      sample_cols <- names(dosage_dt)[-(1:5)]
      if (!identical(sample_cols, sample_ids)) {
        stop(sprintf("Sample order mismatch detected in dosage file: %s", dosage_path), call. = FALSE)
      }

      weight_dt <- unique(chr_selected[, .(ID = as.character(rsid), log_or_weight = as.numeric(log_or_weight))], by = "ID")
      dosage_dt[, ID := as.character(ID)]
      dosage_dt <- merge(dosage_dt, weight_dt, by = "ID", all = FALSE)
      if (nrow(dosage_dt) == 0) next

      numeric_mat <- as.matrix(dosage_dt[, ..sample_cols])
      storage.mode(numeric_mat) <- "numeric"
      weights <- dosage_dt$log_or_weight
      accum <- accum + colSums(sweep(numeric_mat, 1, weights, `*`), na.rm = TRUE)
      found_keys <- c(found_keys, paste(dosage_dt$CHROM, dosage_dt$POS, dosage_dt$ID, label, sep = ":"))
    }

    result_dt[, (feature_name) := accum]
  }

  selected_dt[, found_in_dosage := paste(chr, pos, rsid, ancestry, sep = ":") %in% found_keys]
  list(dosage = result_dt, selected = selected_dt)
}

build_feature_table <- function(phenotype_dt, global_dt, covariates_dt, age_dt, dosage_dt) {
  feature_dt <- merge(phenotype_dt, global_dt, by = "sample_id", all = FALSE)
  feature_dt <- merge(feature_dt, dosage_dt, by = "sample_id", all.x = TRUE)
  if (!is.null(covariates_dt)) feature_dt <- merge(feature_dt, covariates_dt, by = "sample_id", all.x = TRUE)
  if (!is.null(age_dt)) feature_dt <- merge(feature_dt, age_dt, by = "sample_id", all.x = TRUE)
  feature_dt <- unique(feature_dt, by = "sample_id")
  feature_dt
}

compute_auc <- function(y_true, y_prob) {
  y_true <- as.integer(y_true)
  y_prob <- as.numeric(y_prob)
  keep <- is.finite(y_true) & is.finite(y_prob)
  y_true <- y_true[keep]
  y_prob <- y_prob[keep]
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(y_prob, ties.method = "average")
  (sum(ranks[y_true == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

compute_pr_auc <- function(y_true, y_prob) {
  y_true <- as.integer(y_true)
  y_prob <- as.numeric(y_prob)
  ord <- order(-y_prob)
  y_sorted <- y_true[ord]
  if (sum(y_sorted == 1) == 0) return(NA_real_)
  tp <- cumsum(y_sorted == 1)
  fp <- cumsum(y_sorted == 0)
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / sum(y_sorted == 1)
  recall_prev <- c(0, recall[-length(recall)])
  sum((recall - recall_prev) * precision, na.rm = TRUE)
}

make_roc_df <- function(predictions_dt, metrics_dt) {
  out <- list()
  for (model_name in unique(predictions_dt$model)) {
    sub <- predictions_dt[model == model_name]
    y_true <- as.integer(sub$y)
    y_prob <- as.numeric(sub$predicted_probability)
    thresholds <- sort(unique(y_prob), decreasing = TRUE)
    thresholds <- c(Inf, thresholds, -Inf)
    roc_dt <- rbindlist(lapply(thresholds, function(thr) {
      pred <- ifelse(y_prob >= thr, 1L, 0L)
      tp <- sum(pred == 1 & y_true == 1)
      fp <- sum(pred == 1 & y_true == 0)
      fn <- sum(pred == 0 & y_true == 1)
      tn <- sum(pred == 0 & y_true == 0)
      data.table(
        model = model_name,
        threshold = thr,
        tpr = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
        fpr = ifelse(fp + tn > 0, fp / (fp + tn), NA_real_)
      )
    }), fill = TRUE)
    auc_val <- metrics_dt[model == model_name, auc][1]
    roc_dt[, legend_label := sprintf("%s (AUC=%.3f)", model_name, auc_val)]
    out[[length(out) + 1]] <- roc_dt
  }
  rbindlist(out, fill = TRUE)
}

make_calibration_df <- function(predictions_dt) {
  out <- list()
  for (model_name in unique(predictions_dt$model)) {
    sub <- predictions_dt[model == model_name & is.finite(predicted_probability)]
    if (nrow(sub) == 0) next
    probs <- sub$predicted_probability
    breaks <- unique(quantile(probs, probs = seq(0, 1, length.out = 11), na.rm = TRUE, type = 8))
    if (length(breaks) < 3) {
      breaks <- seq(0, 1, length.out = 11)
    }
    sub[, bin := cut(predicted_probability, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
    calib <- sub[, .(
      mean_predicted = mean(predicted_probability, na.rm = TRUE),
      observed_fraction = mean(y, na.rm = TRUE),
      n = .N
    ), by = .(model, bin)]
    out[[length(out) + 1]] <- calib
  }
  rbindlist(out, fill = TRUE)
}

standardize_train_test <- function(train_dt, test_dt, feature_cols) {
  for (col in feature_cols) {
    med <- suppressWarnings(median(train_dt[[col]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    train_dt[is.na(get(col)), (col) := med]
    test_dt[is.na(get(col)), (col) := med]

    mu <- suppressWarnings(mean(train_dt[[col]], na.rm = TRUE))
    if (!is.finite(mu)) mu <- 0
    sigma <- suppressWarnings(sd(train_dt[[col]], na.rm = TRUE))
    if (!is.finite(sigma) || sigma == 0) sigma <- 1

    train_dt[, (col) := (get(col) - mu) / sigma]
    test_dt[, (col) := (get(col) - mu) / sigma]
  }
  list(train = train_dt, test = test_dt)
}

build_formula <- function(feature_cols) {
  if (length(feature_cols) == 0) return(y ~ 1)
  as.formula(paste("y ~", paste(feature_cols, collapse = " + ")))
}

evaluate_model <- function(feature_dt, model_name, feature_cols, cv_folds, random_state) {
  model_cols <- unique(c("sample_id", "IID", "y", feature_cols))
  model_dt <- copy(feature_dt[, ..model_cols])
  for (col in feature_cols) model_dt[, (col) := as.numeric(get(col))]
  model_dt <- model_dt[!is.na(y)]

  y <- as.integer(model_dt$y)
  n_cases <- sum(y == 1L)
  n_controls <- sum(y == 0L)
  formula_obj <- build_formula(feature_cols)

  fallback_predictions <- function(mode, note) {
    prevalence <- if (length(y) > 0) mean(y == 1L, na.rm = TRUE) else 0.5
    if (!is.finite(prevalence)) prevalence <- 0.5
    probs <- rep(prevalence, nrow(model_dt))
    preds <- ifelse(probs >= 0.5, 1L, 0L)

    tp <- sum(preds == 1L & y == 1L)
    tn <- sum(preds == 0L & y == 0L)
    fp <- sum(preds == 1L & y == 0L)
    fn <- sum(preds == 0L & y == 1L)
    sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
    specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)

    metrics_dt <- data.table(
      model = model_name,
      evaluation_mode = mode,
      note = note,
      n_samples = nrow(model_dt),
      n_cases = n_cases,
      n_controls = n_controls,
      n_features = length(feature_cols),
      cv_folds = ifelse(mode == "cross_validation", cv_folds, 1L),
      auc = compute_auc(y, probs),
      pr_auc = compute_pr_auc(y, probs),
      brier_score = mean((probs - y)^2, na.rm = TRUE),
      accuracy = mean(preds == y, na.rm = TRUE),
      balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
      sensitivity = sensitivity,
      specificity = specificity
    )

    predictions_dt <- model_dt[, .(sample_id, IID, y)]
    predictions_dt[, model := model_name]
    predictions_dt[, fold := 1L]
    predictions_dt[, predicted_probability := probs]
    predictions_dt[, predicted_label := preds]

    coef_names <- c("intercept", feature_cols)
    coefficients_dt <- data.table(
      model = model_name,
      feature = coef_names,
      coefficient = NA_real_
    )

    list(predictions = predictions_dt, metrics = metrics_dt, coefficients = coefficients_dt)
  }

  if (nrow(model_dt) == 0) {
    return(fallback_predictions("unavailable", "No non-missing phenotype rows after merging features."))
  }

  if (n_cases == 0L || n_controls == 0L) {
    note <- sprintf("Only one phenotype class remained after merging features (cases=%d, controls=%d).", n_cases, n_controls)
    return(fallback_predictions("single_class", note))
  }

  n_splits <- min(cv_folds, n_cases, n_controls)
  if (n_splits < 2) {
    note <- sprintf("Too few samples per class for cross-validation (cases=%d, controls=%d); using fallback predictions.", n_cases, n_controls)
    return(fallback_predictions("insufficient_cv", note))
  }

  set.seed(random_state)
  fold <- integer(nrow(model_dt))
  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  fold[case_idx] <- sample(rep(seq_len(n_splits), length.out = length(case_idx)))
  fold[control_idx] <- sample(rep(seq_len(n_splits), length.out = length(control_idx)))

  probs <- rep(NA_real_, nrow(model_dt))

  for (fold_idx in seq_len(n_splits)) {
    train_idx <- which(fold != fold_idx)
    test_idx <- which(fold == fold_idx)

    train_dt <- copy(model_dt[train_idx, c("y", feature_cols), with = FALSE])
    test_dt <- copy(model_dt[test_idx, c("y", feature_cols), with = FALSE])
    scaled <- standardize_train_test(train_dt, test_dt, feature_cols)
    train_dt <- scaled$train
    test_dt <- scaled$test

    fit <- suppressWarnings(glm(formula_obj, data = as.data.frame(train_dt), family = binomial()))
    probs[test_idx] <- suppressWarnings(as.numeric(predict(fit, newdata = as.data.frame(test_dt), type = "response")))
  }

  if (all(!is.finite(probs))) probs[] <- mean(y == 1L, na.rm = TRUE)
  if (any(!is.finite(probs))) probs[!is.finite(probs)] <- median(probs[is.finite(probs)], na.rm = TRUE)
  preds <- ifelse(probs >= 0.5, 1L, 0L)

  tp <- sum(preds == 1L & y == 1L)
  tn <- sum(preds == 0L & y == 0L)
  fp <- sum(preds == 1L & y == 0L)
  fn <- sum(preds == 0L & y == 1L)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)

  metrics_dt <- data.table(
    model = model_name,
    evaluation_mode = "cross_validation",
    note = "ok",
    n_samples = nrow(model_dt),
    n_cases = n_cases,
    n_controls = n_controls,
    n_features = length(feature_cols),
    cv_folds = n_splits,
    auc = compute_auc(y, probs),
    pr_auc = compute_pr_auc(y, probs),
    brier_score = mean((probs - y)^2, na.rm = TRUE),
    accuracy = mean(preds == y, na.rm = TRUE),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity
  )

  predictions_dt <- model_dt[, .(sample_id, IID, y)]
  predictions_dt[, model := model_name]
  predictions_dt[, fold := fold]
  predictions_dt[, predicted_probability := probs]
  predictions_dt[, predicted_label := preds]

  full_dt <- copy(model_dt[, c("y", feature_cols), with = FALSE])
  scaled_full <- standardize_train_test(copy(full_dt), copy(full_dt), feature_cols)
  full_fit <- suppressWarnings(glm(formula_obj, data = as.data.frame(scaled_full$train), family = binomial()))
  coef_vec <- coef(full_fit)
  coefficients_dt <- data.table(
    model = model_name,
    feature = names(coef_vec),
    coefficient = as.numeric(coef_vec)
  )
  coefficients_dt[feature == "(Intercept)", feature := "intercept"]

  list(predictions = predictions_dt, metrics = metrics_dt, coefficients = coefficients_dt)
}

evaluate_model_repeated_holdout <- function(feature_dt, model_name, feature_cols, n_repeats, test_size, random_state) {
  model_cols <- unique(c("sample_id", "IID", "y", feature_cols))
  model_dt <- copy(feature_dt[, ..model_cols])
  for (col in feature_cols) model_dt[, (col) := as.numeric(get(col))]
  model_dt <- model_dt[!is.na(y)]

  y <- as.integer(model_dt$y)
  n_cases <- sum(y == 1L)
  n_controls <- sum(y == 0L)
  formula_obj <- build_formula(feature_cols)

  fallback_predictions <- function(mode, note) {
    prevalence <- if (length(y) > 0) mean(y == 1L, na.rm = TRUE) else 0.5
    if (!is.finite(prevalence)) prevalence <- 0.5
    probs <- rep(prevalence, nrow(model_dt))
    preds <- ifelse(probs >= 0.5, 1L, 0L)

    tp <- sum(preds == 1L & y == 1L)
    tn <- sum(preds == 0L & y == 0L)
    fp <- sum(preds == 1L & y == 0L)
    fn <- sum(preds == 0L & y == 1L)
    sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
    specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)

    metrics_dt <- data.table(
      model = model_name,
      evaluation_mode = mode,
      note = note,
      n_samples = nrow(model_dt),
      n_cases = n_cases,
      n_controls = n_controls,
      n_features = length(feature_cols),
      cv_folds = 1L,
      auc = compute_auc(y, probs),
      pr_auc = compute_pr_auc(y, probs),
      brier_score = mean((probs - y)^2, na.rm = TRUE),
      accuracy = mean(preds == y, na.rm = TRUE),
      balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
      sensitivity = sensitivity,
      specificity = specificity,
      n_repeats = n_repeats,
      holdout_test_size = test_size,
      n_test_predictions = nrow(model_dt),
      n_unique_test_samples = uniqueN(model_dt$sample_id)
    )

    predictions_dt <- model_dt[, .(sample_id, IID, y)]
    predictions_dt[, model := model_name]
    predictions_dt[, fold := 1L]
    predictions_dt[, repeat_id := 1L]
    predictions_dt[, predicted_probability := probs]
    predictions_dt[, predicted_label := preds]

    coef_names <- c("intercept", feature_cols)
    coefficients_dt <- data.table(
      model = model_name,
      feature = coef_names,
      coefficient = NA_real_
    )

    list(predictions = predictions_dt, metrics = metrics_dt, coefficients = coefficients_dt)
  }

  if (nrow(model_dt) == 0) {
    return(fallback_predictions("unavailable", "No non-missing phenotype rows after merging features."))
  }

  if (n_cases == 0L || n_controls == 0L) {
    note <- sprintf("Only one phenotype class remained after merging features (cases=%d, controls=%d).", n_cases, n_controls)
    return(fallback_predictions("single_class", note))
  }

  if (!is.finite(test_size) || test_size <= 0 || test_size >= 1) {
    stop(sprintf("test_size must be in (0,1), got %s", as.character(test_size)), call. = FALSE)
  }

  test_cases <- max(1L, floor(n_cases * test_size))
  test_controls <- max(1L, floor(n_controls * test_size))
  if (test_cases >= n_cases || test_controls >= n_controls) {
    note <- sprintf("Holdout split too aggressive for class counts (cases=%d, controls=%d).", n_cases, n_controls)
    return(fallback_predictions("insufficient_holdout", note))
  }

  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  prediction_frames <- list()

  for (repeat_idx in seq_len(n_repeats)) {
    set.seed(random_state + repeat_idx - 1L)

    test_case_idx <- sample(case_idx, size = test_cases, replace = FALSE)
    test_control_idx <- sample(control_idx, size = test_controls, replace = FALSE)
    test_idx <- sort(c(test_case_idx, test_control_idx))
    train_idx <- setdiff(seq_len(nrow(model_dt)), test_idx)

    train_dt <- copy(model_dt[train_idx, c("y", feature_cols), with = FALSE])
    test_dt <- copy(model_dt[test_idx, c("y", feature_cols), with = FALSE])
    scaled <- standardize_train_test(train_dt, test_dt, feature_cols)
    train_dt <- scaled$train
    test_dt <- scaled$test

    fit <- suppressWarnings(glm(formula_obj, data = as.data.frame(train_dt), family = binomial()))
    probs <- suppressWarnings(as.numeric(predict(fit, newdata = as.data.frame(test_dt), type = "response")))
    if (all(!is.finite(probs))) {
      probs <- rep(mean(train_dt$y == 1L, na.rm = TRUE), length(test_idx))
    }
    if (any(!is.finite(probs))) {
      probs[!is.finite(probs)] <- median(probs[is.finite(probs)], na.rm = TRUE)
    }

    preds <- ifelse(probs >= 0.5, 1L, 0L)
    rep_dt <- model_dt[test_idx, .(sample_id, IID, y)]
    rep_dt[, model := model_name]
    rep_dt[, fold := repeat_idx]
    rep_dt[, repeat_id := repeat_idx]
    rep_dt[, predicted_probability := probs]
    rep_dt[, predicted_label := preds]
    prediction_frames[[length(prediction_frames) + 1]] <- rep_dt
  }

  predictions_dt <- rbindlist(prediction_frames, fill = TRUE)

  y_eval <- as.integer(predictions_dt$y)
  probs_eval <- as.numeric(predictions_dt$predicted_probability)
  preds_eval <- as.integer(predictions_dt$predicted_label)

  tp <- sum(preds_eval == 1L & y_eval == 1L)
  tn <- sum(preds_eval == 0L & y_eval == 0L)
  fp <- sum(preds_eval == 1L & y_eval == 0L)
  fn <- sum(preds_eval == 0L & y_eval == 1L)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)

  metrics_dt <- data.table(
    model = model_name,
    evaluation_mode = "repeated_stratified_holdout_80_20",
    note = sprintf("%d stratified random repeats with %.0f%% train / %.0f%% test", n_repeats, (1 - test_size) * 100, test_size * 100),
    n_samples = nrow(model_dt),
    n_cases = n_cases,
    n_controls = n_controls,
    n_features = length(feature_cols),
    cv_folds = n_repeats,
    auc = compute_auc(y_eval, probs_eval),
    pr_auc = compute_pr_auc(y_eval, probs_eval),
    brier_score = mean((probs_eval - y_eval)^2, na.rm = TRUE),
    accuracy = mean(preds_eval == y_eval, na.rm = TRUE),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    n_repeats = n_repeats,
    holdout_test_size = test_size,
    n_test_predictions = nrow(predictions_dt),
    n_unique_test_samples = uniqueN(predictions_dt$sample_id)
  )

  full_dt <- copy(model_dt[, c("y", feature_cols), with = FALSE])
  scaled_full <- standardize_train_test(copy(full_dt), copy(full_dt), feature_cols)
  full_fit <- suppressWarnings(glm(formula_obj, data = as.data.frame(scaled_full$train), family = binomial()))
  coef_vec <- coef(full_fit)
  coefficients_dt <- data.table(
    model = model_name,
    feature = names(coef_vec),
    coefficient = as.numeric(coef_vec)
  )
  coefficients_dt[feature == "(Intercept)", feature := "intercept"]

  list(predictions = predictions_dt, metrics = metrics_dt, coefficients = coefficients_dt)
}

make_plots <- function(predictions_dt, metrics_dt, roc_png, calibration_png) {
  roc_dt <- make_roc_df(predictions_dt, metrics_dt)
  if (nrow(roc_dt) == 0 || !all(c("fpr", "tpr") %in% names(roc_dt))) {
    p_roc <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "ROC unavailable\n(no evaluable case/control split)", size = 5) +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = "Admixture-validation ROC curves", x = "False positive rate", y = "True positive rate") +
      theme_bw()
  } else {
    p_roc <- ggplot(roc_dt, aes(x = fpr, y = tpr, color = legend_label)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      labs(
        title = "Admixture-validation ROC curves",
        x = "False positive rate",
        y = "True positive rate",
        color = NULL
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  ggsave(roc_png, p_roc, width = 7, height = 6, dpi = 200)

  calib_dt <- make_calibration_df(predictions_dt)
  if (nrow(calib_dt) == 0 || !all(c("mean_predicted", "observed_fraction") %in% names(calib_dt))) {
    p_cal <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Calibration unavailable\n(no evaluable case/control split)", size = 5) +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = "Admixture-validation calibration", x = "Mean predicted probability", y = "Observed case fraction") +
      theme_bw()
  } else {
    p_cal <- ggplot(calib_dt, aes(x = mean_predicted, y = observed_fraction, color = model)) +
      geom_point(size = 2, na.rm = TRUE) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      labs(
        title = "Admixture-validation calibration",
        x = "Mean predicted probability",
        y = "Observed case fraction",
        color = NULL
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  ggsave(calibration_png, p_cal, width = 7, height = 6, dpi = 200)
}

main <- function() {
  args_list <- list(
    summary_tsv = get_arg("--summary-tsv", required = TRUE),
    rfmix_q_pattern = get_arg("--rfmix-q-pattern", required = TRUE),
    tractor_dosage_pattern = get_arg("--tractor-dosage-pattern", required = TRUE),
    chromosomes = get_arg("--chromosomes", required = TRUE),
    ancestry_map = get_arg("--ancestry-map", required = TRUE),
    phenotype = get_arg("--phenotype", required = TRUE),
    phenotype_format = get_arg("--phenotype-format", required = TRUE),
    fam_case_code = get_arg("--fam-case-code", default = "2"),
    fam_control_code = get_arg("--fam-control-code", default = "1"),
    table_iid_column = get_arg("--table-iid-column", default = "IID"),
    table_phenotype_column = get_arg("--table-phenotype-column", default = "PHENO"),
    covariates_file = get_arg("--covariates-file", default = ""),
    covariates_format = get_arg("--covariates-format", default = "plink_eigenvec"),
    covariates_has_header = get_arg("--covariates-has-header", default = "false"),
    covariates_fid_column = get_arg("--covariates-fid-column", default = "FID"),
    covariates_iid_column = get_arg("--covariates-iid-column", default = "IID"),
    n_pcs = as.integer(get_arg("--n-pcs", default = "5")),
    include_sex = get_arg("--include-sex", default = "true"),
    include_age = get_arg("--include-age", default = "false"),
    age_file = get_arg("--age-file", default = ""),
    age_iid_column = get_arg("--age-iid-column", default = "IID"),
    age_column = get_arg("--age-column", default = "AGE"),
    global_reference_ancestry = get_arg("--global-reference-ancestry", default = "AMR"),
    snp_pvalue_threshold = as.numeric(get_arg("--snp-pvalue-threshold", default = "0.01")),
    max_snps_per_ancestry = as.integer(get_arg("--max-snps-per-ancestry", default = "50")),
    use_top_snp_per_locus = get_arg("--use-top-snp-per-locus", default = "true"),
    cv_folds = as.integer(get_arg("--cv-folds", default = "5")),
    validation_design = get_arg("--validation-design", default = "cross_validation"),
    n_repeats = as.integer(get_arg("--n-repeats", default = "10")),
    holdout_test_size = as.numeric(get_arg("--holdout-test-size", default = "0.2")),
    random_state = as.integer(get_arg("--random-state", default = "42")),
    sample_features_out = get_arg("--sample-features-out", required = TRUE),
    selected_snps_out = get_arg("--selected-snps-out", required = TRUE),
    cv_predictions_out = get_arg("--cv-predictions-out", required = TRUE),
    metrics_out = get_arg("--metrics-out", required = TRUE),
    coefficients_out = get_arg("--coefficients-out", required = TRUE),
    roc_png = get_arg("--roc-png", required = TRUE),
    calibration_png = get_arg("--calibration-png", required = TRUE)
  )

  chromosomes <- Filter(nzchar, strsplit(as.character(args_list$chromosomes), ",", fixed = TRUE)[[1]])
  ancestry_map <- parse_ancestry_map(args_list$ancestry_map)
  ancestry_labels <- unname(unlist(ancestry_map[order(as.integer(names(ancestry_map)))]))

  phenotype_dt <- load_phenotype(args_list)
  covariate_info <- load_covariates(args_list)
  covariates_dt <- covariate_info$data
  pc_columns <- covariate_info$pc_columns
  age_dt <- load_age(args_list)
  global_dt <- load_global_ancestry(args_list, ancestry_map, chromosomes)

  selected_dt <- select_snps(
    args_list$summary_tsv,
    ancestry_labels = ancestry_labels,
    p_threshold = args_list$snp_pvalue_threshold,
    max_snps_per_ancestry = args_list$max_snps_per_ancestry,
    use_top_snp_per_locus = parse_bool(args_list$use_top_snp_per_locus)
  )

  dosage_info <- compute_dosage_scores(args_list, selected_dt, ancestry_map, chromosomes)
  dosage_dt <- dosage_info$dosage
  selected_dt <- dosage_info$selected

  feature_dt <- build_feature_table(phenotype_dt, global_dt, covariates_dt, age_dt, dosage_dt)
  message(sprintf(
    "Merged validation feature table: n=%d, cases=%d, controls=%d",
    nrow(feature_dt),
    sum(feature_dt$y == 1, na.rm = TRUE),
    sum(feature_dt$y == 0, na.rm = TRUE)
  ))

  global_cols <- paste0("global_", ancestry_labels[ancestry_labels != args_list$global_reference_ancestry], "_prop")
  dosage_cols <- paste0("hap_dosage_", ancestry_labels)
  optional_cols <- character()
  if (parse_bool(args_list$include_sex) && "sex_female" %in% names(feature_dt)) optional_cols <- c(optional_cols, "sex_female")
  if (parse_bool(args_list$include_age) && "age" %in% names(feature_dt) && any(!is.na(feature_dt$age))) optional_cols <- c(optional_cols, "age")

  model_definitions <- list(
    pcs_only = unique(c(pc_columns, optional_cols)),
    pcs_plus_global = unique(c(pc_columns, optional_cols, global_cols)),
    full_admixture = unique(c(pc_columns, optional_cols, global_cols, dosage_cols))
  )

  prediction_frames <- list()
  metrics_frames <- list()
  coefficient_frames <- list()

  for (model_name in names(model_definitions)) {
    feature_cols <- model_definitions[[model_name]]
    feature_cols <- feature_cols[feature_cols %in% names(feature_dt)]
    if (tolower(args_list$validation_design) == "repeated_stratified_holdout_80_20") {
      res <- evaluate_model_repeated_holdout(
        feature_dt = feature_dt,
        model_name = model_name,
        feature_cols = feature_cols,
        n_repeats = args_list$n_repeats,
        test_size = args_list$holdout_test_size,
        random_state = args_list$random_state
      )
    } else {
      res <- evaluate_model(
        feature_dt = feature_dt,
        model_name = model_name,
        feature_cols = feature_cols,
        cv_folds = args_list$cv_folds,
        random_state = args_list$random_state
      )
    }
    prediction_frames[[length(prediction_frames) + 1]] <- res$predictions
    metrics_frames[[length(metrics_frames) + 1]] <- res$metrics
    coefficient_frames[[length(coefficient_frames) + 1]] <- res$coefficients
    message(sprintf("Evaluated %s with %d features: AUC=%.4f", model_name, length(feature_cols), res$metrics$auc[1]))
  }

  predictions_dt <- rbindlist(prediction_frames, fill = TRUE)
  metrics_dt <- rbindlist(metrics_frames, fill = TRUE)
  coefficients_dt <- rbindlist(coefficient_frames, fill = TRUE)

  dir.create(dirname(args_list$sample_features_out), recursive = TRUE, showWarnings = FALSE)
  fwrite(feature_dt, args_list$sample_features_out, sep = "\t")
  fwrite(selected_dt, args_list$selected_snps_out, sep = "\t")
  fwrite(predictions_dt, args_list$cv_predictions_out, sep = "\t")
  fwrite(metrics_dt, args_list$metrics_out, sep = "\t")
  fwrite(coefficients_dt, args_list$coefficients_out, sep = "\t")
  make_plots(predictions_dt, metrics_dt, args_list$roc_png, args_list$calibration_png)

  message(sprintf(
    "Wrote sample features to %s, metrics to %s, and plots to %s / %s",
    args_list$sample_features_out,
    args_list$metrics_out,
    args_list$roc_png,
    args_list$calibration_png
  ))
}

main()
