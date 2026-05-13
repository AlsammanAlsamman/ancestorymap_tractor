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
  if (idx[length(idx)] == length(args)) stop(paste("Missing value for argument", flag), call. = FALSE)
  args[idx[length(idx)] + 1]
}

parse_ancestry_map <- function(value) {
  mapping <- list()
  parts <- strsplit(as.character(value), ",", fixed = TRUE)[[1]]
  for (part in parts) {
    if (!nzchar(part)) next
    kv <- strsplit(part, ":", fixed = TRUE)[[1]]
    if (length(kv) == 2) mapping[[trimws(kv[1])]] <- trimws(kv[2])
  }
  mapping
}

sanitize_name <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("(^_+|_+$)", "", out)
  ifelse(nzchar(out), out, "unknown_locus")
}

make_feature_name <- function(locus_name, ancestry_label) {
  paste0("dosage__", ancestry_label, "__", sanitize_name(locus_name))
}

compute_auc <- function(y_true, scores) {
  y_true <- as.integer(y_true)
  keep <- is.finite(scores) & !is.na(y_true)
  y_true <- y_true[keep]
  scores <- scores[keep]
  n1 <- sum(y_true == 1L)
  n0 <- sum(y_true == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(scores, ties.method = "average")
  (sum(r[y_true == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

save_png <- function(plot_obj, filename, width = 1400, height = 1000, res = 150) {
  if (capabilities("cairo")) {
    png(filename, width = width, height = height, res = res, type = "cairo")
  } else {
    png(filename, width = width, height = height, res = res)
  }
  print(plot_obj)
  dev.off()
}

load_sample_ids_from_first_dosage <- function(pattern, ancestry_map, chromosomes) {
  first_idx <- names(ancestry_map)[order(as.integer(names(ancestry_map)))][1]
  first_chr <- chromosomes[1]
  path <- gsub("{chr}", first_chr, pattern, fixed = TRUE)
  path <- gsub("{anc}", first_idx, path, fixed = TRUE)
  header <- fread(path, nrows = 0, data.table = TRUE, check.names = FALSE)
  names(header)[-(1:5)]
}

compute_locus_scores <- function(selected_dt, dosage_pattern, ancestry_map, chromosomes) {
  sample_ids <- load_sample_ids_from_first_dosage(dosage_pattern, ancestry_map, chromosomes)
  scores_dt <- data.table(sample_id = sample_ids)
  meta_list <- list()
  ancestry_by_label <- setNames(names(ancestry_map), unlist(ancestry_map))

  selected_dt[, chr := as.character(chr)]
  selected_dt[, locus_name := fifelse(is.na(locus_name) | !nzchar(locus_name), as.character(gene), as.character(locus_name))]

  for (label in names(ancestry_by_label)) {
    anc_idx <- ancestry_by_label[[label]]
    label_selected <- selected_dt[ancestry == label]
    if (nrow(label_selected) == 0) next

    for (chr_name in chromosomes) {
      chr_selected <- label_selected[chr == as.character(chr_name)]
      if (nrow(chr_selected) == 0) next

      dosage_path <- gsub("{chr}", chr_name, dosage_pattern, fixed = TRUE)
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
      if (is.null(dosage_dt) || nrow(dosage_dt) == 0 || ncol(dosage_dt) <= 5) next

      sample_cols <- names(dosage_dt)[-(1:5)]
      if (!identical(sample_cols, sample_ids)) {
        stop(sprintf("Sample order mismatch detected in dosage file: %s", dosage_path), call. = FALSE)
      }

      dosage_dt[, ID := as.character(ID)]
      weight_dt <- unique(chr_selected[, .(
        ID = as.character(rsid),
        locus_name = as.character(locus_name),
        log_or_weight = as.numeric(log_or_weight)
      )], by = c("ID", "locus_name"))
      dosage_dt <- merge(dosage_dt, weight_dt, by = "ID", all = FALSE)
      if (nrow(dosage_dt) == 0) next

      for (loc in unique(dosage_dt$locus_name)) {
        loc_dt <- dosage_dt[locus_name == loc]
        if (nrow(loc_dt) == 0) next
        mat <- as.matrix(loc_dt[, ..sample_cols])
        storage.mode(mat) <- "numeric"
        weights <- as.numeric(loc_dt$log_or_weight)
        weights[!is.finite(weights)] <- 0
        score <- colSums(sweep(mat, 1, weights, `*`), na.rm = TRUE)
        feature_name <- make_feature_name(loc, label)

        if (feature_name %in% names(scores_dt)) {
          scores_dt[, (feature_name) := get(feature_name) + score]
        } else {
          scores_dt[, (feature_name) := score]
        }

        meta_list[[length(meta_list) + 1]] <- data.table(
          feature = feature_name,
          locus_name = as.character(loc),
          ancestry = as.character(label),
          chr = as.character(chr_name),
          n_snps = uniqueN(loc_dt$ID)
        )
      }
    }
  }

  meta_dt <- unique(rbindlist(meta_list, fill = TRUE), by = c("feature", "locus_name", "ancestry"))
  list(scores = scores_dt, meta = meta_dt)
}

prepare_model_data <- function(sample_features_dt, locus_scores_dt, global_ref_ancestry) {
  model_dt <- merge(sample_features_dt, locus_scores_dt, by = "sample_id", all.x = TRUE)
  locus_feature_cols <- grep("^dosage__", names(model_dt), value = TRUE)
  for (col in locus_feature_cols) set(model_dt, which(is.na(model_dt[[col]])), col, 0)

  base_cols <- c(grep("^PC", names(model_dt), value = TRUE), "sex_female")
  base_cols <- intersect(base_cols, names(model_dt))
  global_cols <- grep("^global_.*_prop$", names(model_dt), value = TRUE)
  if (nzchar(global_ref_ancestry)) {
    drop_col <- paste0("global_", global_ref_ancestry, "_prop")
    global_cols <- setdiff(global_cols, drop_col)
  }
  feature_cols <- unique(c(base_cols, global_cols, locus_feature_cols))

  model_dt <- model_dt[!is.na(y)]
  for (col in feature_cols) {
    vals <- as.numeric(model_dt[[col]])
    med <- suppressWarnings(median(vals[is.finite(vals)], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    vals[!is.finite(vals)] <- med
    sdv <- sd(vals)
    if (is.na(sdv) || sdv == 0) {
      model_dt[, (col) := 0]
    } else {
      model_dt[, (col) := as.numeric(scale(vals))]
    }
  }

  non_zero_cols <- feature_cols[vapply(feature_cols, function(col) sd(model_dt[[col]]) > 0, logical(1))]
  list(data = model_dt, feature_cols = non_zero_cols, locus_feature_cols = intersect(locus_feature_cols, non_zero_cols))
}

fit_importance_model <- function(model_dt, feature_cols) {
  if (length(feature_cols) == 0) stop("No non-zero features available for the importance model.", call. = FALSE)
  formula_txt <- paste("y ~", paste(feature_cols, collapse = " + "))
  fit <- suppressWarnings(glm(as.formula(formula_txt), data = model_dt, family = binomial(), control = list(maxit = 100)))
  probs <- suppressWarnings(predict(fit, newdata = model_dt, type = "response"))
  preds <- as.integer(probs >= 0.5)
  y <- as.integer(model_dt$y)

  tp <- sum(preds == 1 & y == 1)
  tn <- sum(preds == 0 & y == 0)
  fp <- sum(preds == 1 & y == 0)
  fn <- sum(preds == 0 & y == 1)

  metrics_dt <- data.table(
    model = "locus_dosage_importance",
    n_samples = nrow(model_dt),
    n_cases = sum(y == 1),
    n_controls = sum(y == 0),
    n_features = length(feature_cols),
    auc = compute_auc(y, probs),
    accuracy = mean(preds == y),
    sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  )

  coef_mat <- summary(fit)$coefficients
  coef_dt <- data.table(
    feature = rownames(coef_mat),
    coefficient = coef_mat[, 1],
    std_error = coef_mat[, 2],
    z_value = coef_mat[, 3],
    p_value = coef_mat[, 4]
  )
  rownames(coef_mat) <- NULL
  list(metrics = metrics_dt, coefficients = coef_dt)
}

args_list <- list(
  sample_features_tsv = get_arg("--sample-features-tsv", required = TRUE),
  selected_snps_tsv = get_arg("--selected-snps-tsv", required = TRUE),
  tractor_dosage_pattern = get_arg("--tractor-dosage-pattern", required = TRUE),
  chromosomes = strsplit(get_arg("--chromosomes", required = TRUE), ",", fixed = TRUE)[[1]],
  ancestry_map = parse_ancestry_map(get_arg("--ancestry-map", required = TRUE)),
  global_reference_ancestry = get_arg("--global-reference-ancestry", default = "AMR"),
  top_n_loci = as.integer(get_arg("--top-n-loci", default = "20")),
  locus_scores_out = get_arg("--locus-scores-out", required = TRUE),
  locus_importance_out = get_arg("--locus-importance-out", required = TRUE),
  locus_overall_out = get_arg("--locus-overall-out", required = TRUE),
  ancestry_importance_out = get_arg("--ancestry-importance-out", required = TRUE),
  metrics_out = get_arg("--metrics-out", required = TRUE),
  top_loci_png = get_arg("--top-loci-png", required = TRUE),
  ancestry_png = get_arg("--ancestry-png", required = TRUE),
  heatmap_png = get_arg("--heatmap-png", required = TRUE)
)

sample_features_dt <- fread(args_list$sample_features_tsv, data.table = TRUE)
selected_dt <- fread(args_list$selected_snps_tsv, data.table = TRUE)
selected_dt[, log_or_weight := as.numeric(log_or_weight)]
if (!("locus_name" %in% names(selected_dt))) selected_dt[, locus_name := NA_character_]
if (!("gene" %in% names(selected_dt))) selected_dt[, gene := NA_character_]
selected_dt[, locus_name := as.character(locus_name)]
selected_dt[, gene := as.character(gene)]
selected_dt[is.na(locus_name) | !nzchar(trimws(locus_name)), locus_name := gene]
selected_dt[is.na(locus_name) | !nzchar(trimws(locus_name)), locus_name := paste(chr, pos, rsid, sep = ":")]

locus_res <- compute_locus_scores(
  selected_dt = selected_dt,
  dosage_pattern = args_list$tractor_dosage_pattern,
  ancestry_map = args_list$ancestry_map,
  chromosomes = args_list$chromosomes
)

prepared <- prepare_model_data(
  sample_features_dt = sample_features_dt,
  locus_scores_dt = locus_res$scores,
  global_ref_ancestry = args_list$global_reference_ancestry
)
fit_res <- fit_importance_model(prepared$data, prepared$feature_cols)

importance_dt <- merge(locus_res$meta, fit_res$coefficients, by = "feature", all.x = TRUE)
importance_dt <- importance_dt[!is.na(coefficient)]
importance_dt[, abs_coefficient := abs(coefficient)]
importance_dt[, direction := fifelse(coefficient >= 0, "risk_increasing", "protective")]
setorder(importance_dt, -abs_coefficient)

locus_overall_dt <- importance_dt[, .(
  total_importance = sum(abs_coefficient, na.rm = TRUE),
  signed_effect_sum = sum(coefficient, na.rm = TRUE),
  dominant_ancestry = ancestry[which.max(abs_coefficient)],
  n_ancestry_features = .N
), by = locus_name][order(-total_importance)]

ancestry_importance_dt <- importance_dt[, .(
  total_importance = sum(abs_coefficient, na.rm = TRUE),
  mean_importance = mean(abs_coefficient, na.rm = TRUE),
  n_locus_features = .N
), by = ancestry][order(-total_importance)]

fwrite(locus_res$scores, args_list$locus_scores_out, sep = "\t")
fwrite(importance_dt, args_list$locus_importance_out, sep = "\t")
fwrite(locus_overall_dt, args_list$locus_overall_out, sep = "\t")
fwrite(ancestry_importance_dt, args_list$ancestry_importance_out, sep = "\t")
fwrite(fit_res$metrics, args_list$metrics_out, sep = "\t")

top_n <- min(args_list$top_n_loci, nrow(locus_overall_dt))
top_loci <- head(locus_overall_dt$locus_name, top_n)
plot_dt <- importance_dt[locus_name %in% top_loci]
plot_dt[, locus_name := factor(locus_name, levels = rev(locus_overall_dt[locus_name %in% top_loci]$locus_name))]

p_top <- ggplot(plot_dt, aes(x = locus_name, y = abs_coefficient, fill = ancestry)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Top locus dosage importance in the prediction model",
    x = "Locus",
    y = "Absolute standardized coefficient"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")
save_png(p_top, args_list$top_loci_png)

p_anc <- ggplot(ancestry_importance_dt, aes(x = reorder(ancestry, total_importance), y = total_importance, fill = ancestry)) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    title = "Overall ancestry dosage importance",
    x = "Ancestry",
    y = "Total absolute dosage importance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
save_png(p_anc, args_list$ancestry_png, width = 1100, height = 800)

heatmap_dt <- importance_dt[locus_name %in% top_loci]
heatmap_dt[, locus_name := factor(locus_name, levels = rev(locus_overall_dt[locus_name %in% top_loci]$locus_name))]
p_heat <- ggplot(heatmap_dt, aes(x = ancestry, y = locus_name, fill = abs_coefficient)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  labs(
    title = "Locus-by-ancestry dosage importance heatmap",
    x = "Ancestry",
    y = "Locus",
    fill = "|Coefficient|"
  ) +
  theme_bw(base_size = 12)
save_png(p_heat, args_list$heatmap_png, width = 1200, height = 1000)

message(sprintf(
  "Wrote locus importance to %s, ancestry importance to %s, and plots to %s / %s / %s",
  args_list$locus_importance_out,
  args_list$ancestry_importance_out,
  args_list$top_loci_png,
  args_list$ancestry_png,
  args_list$heatmap_png
))
