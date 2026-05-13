#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
})

options(bitmapType = "cairo")

ANCESTRIES <- character()
CLASSIFICATION_LABELS <- c()
CLASSIFICATION_NOTES <- c()
CLASSIFICATION_COLORS <- c()

build_classification_metadata <- function(ancestries) {
  labels <- c()
  notes <- c()
  colors <- c()
  base_palette <- c("#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F")
  pair_palette <- c("#3CAEA3", "#A1B83A", "#B27652", "#7E9C6F", "#9D7A68")

  for (i in seq_along(ancestries)) {
    anc <- ancestries[i]
    key <- paste0(anc, "_specific")
    labels[key] <- paste0(anc, "-specific")
    notes[key] <- paste0("Signal passes the QC rules in ", anc, " only.")
    colors[key] <- base_palette[((i - 1) %% length(base_palette)) + 1]
  }

  pair_idx <- 1
  if (length(ancestries) >= 3) {
    for (k in seq.int(2, length(ancestries) - 1)) {
      cmb <- utils::combn(ancestries, k, simplify = FALSE)
      for (pair in cmb) {
        key <- paste0(paste(pair, collapse = "_"), "_shared")
        labels[key] <- paste0(paste(pair, collapse = " + "), " shared")
        notes[key] <- paste0("Signal is shared between ", paste(pair, collapse = ", "), ".")
        colors[key] <- pair_palette[((pair_idx - 1) %% length(pair_palette)) + 1]
        pair_idx <- pair_idx + 1
      }
    }
  }

  labels["shared_all"] <- paste0("Shared across all ", length(ancestries), " ancestries")
  notes["shared_all"] <- "Signal is present across all configured ancestries."
  colors["shared_all"] <- "#C76BA2"

  labels["no_clear_signal"] <- "No clear ancestry-specific pattern"
  notes["no_clear_signal"] <- "No ancestry reached the support threshold under the current QC rules."
  colors["no_clear_signal"] <- "#B8BEC8"

  list(labels = labels, notes = notes, colors = colors)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --name value pairs.")
  }
  out <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    key <- gsub("-", "_", key)
    out[[key]] <- args[[i + 1]]
  }
  out
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

not_blank <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

unique_join <- function(x, sep = ";") {
  vals <- trimws(as.character(x))
  vals <- vals[not_blank(vals)]
  if (!length(vals)) return("")
  paste(unique(vals), collapse = sep)
}

format_region <- function(chr, start, end) {
  if (is.na(start) || is.na(end)) return(paste0("chr", chr))
  paste0("chr", chr, ":", format(start, big.mark = ",", scientific = FALSE), "-", format(end, big.mark = ",", scientific = FALSE))
}

build_parameters_table <- function(args) {
  data.frame(
    parameter = c(
      "cohort_maf_min",
      "ancestry_maf_min",
      "local_ancestry_min_logp",
      "dosage_gwas_max_p",
      "min_supporting_snps",
      "clustering_method",
      "clustering_metric"
    ),
    value = c(
      args$cohort_maf_min,
      args$ancestry_maf_min,
      args$local_ancestry_min_logp,
      args$dosage_gwas_max_p,
      args$min_supporting_snps,
      args$clustering_method,
      args$clustering_metric
    ),
    description = c(
      "Keep SNPs with cohort MAF at or above this threshold.",
      "Require at least one ancestry-specific MAF to meet this minimum.",
      "Filter SNPs first by requiring -log10(local ancestry p-value) to be above this threshold.",
      "Among SNPs that pass the local ancestry filter, require ancestry-specific GWAS dosage p-value at or below this threshold.",
      "Minimum number of supporting SNPs to call an ancestry active for a locus.",
      "Hierarchical clustering linkage method used for the dendrogram.",
      "Distance metric used to compare loci in the dendrogram."
    ),
    stringsAsFactors = FALSE
  )
}

build_conclusions_table <- function(report) {
  if (!nrow(report)) {
    return(data.frame(
      section = c("summary", "summary"),
      item = c("total_loci", "message"),
      value = c("0", "No loci were available for clustering."),
      stringsAsFactors = FALSE
    ))
  }

  class_counts <- sort(table(report$classification), decreasing = TRUE)
  strongest <- report[order(safe_num(report$lead_best_p_dosage), report$chr, safe_num(report$start_pos)), , drop = FALSE]
  strongest <- head(strongest, 5)
  rows <- data.frame(
    section = c(
      "summary", "summary", "summary", "summary", "summary", "summary"
    ),
    item = c(
      "total_loci",
      "one_ancestry_loci",
      "shared_subset_loci",
      "shared_all_ancestries_loci",
      "no_clear_signal_loci",
      "most_common_pattern"
    ),
    value = c(
      nrow(report),
      sum(safe_num(report$n_active_ancestries) == 1, na.rm = TRUE),
      sum(safe_num(report$n_active_ancestries) >= 2 & safe_num(report$n_active_ancestries) < length(ANCESTRIES), na.rm = TRUE),
      sum(safe_num(report$n_active_ancestries) == length(ANCESTRIES), na.rm = TRUE),
      sum(safe_num(report$n_active_ancestries) == 0, na.rm = TRUE),
      unname(CLASSIFICATION_LABELS[names(class_counts)[1]])
    ),
    stringsAsFactors = FALSE
  )

  if (length(class_counts)) {
    rows <- rbind(
      rows,
      data.frame(
        section = "classification_counts",
        item = unname(ifelse(names(class_counts) %in% names(CLASSIFICATION_LABELS), CLASSIFICATION_LABELS[names(class_counts)], names(class_counts))),
        value = as.integer(class_counts),
        stringsAsFactors = FALSE
      )
    )
  }

  if (nrow(strongest)) {
    rows <- rbind(
      rows,
      data.frame(
        section = "top_loci",
        item = paste0(strongest$locus_name, " | ", strongest$lead_gene, " | ", strongest$region),
        value = unname(ifelse(strongest$classification %in% names(CLASSIFICATION_LABELS), CLASSIFICATION_LABELS[strongest$classification], strongest$classification)),
        stringsAsFactors = FALSE
      )
    )
  }

  rows
}

build_conclusion_lines <- function(report) {
  if (!nrow(report)) {
    return(c("No loci were available for clustering."))
  }
  class_counts <- sort(table(report$classification), decreasing = TRUE)
  top_loci <- report[order(safe_num(report$lead_best_p_dosage), report$chr, safe_num(report$start_pos)), , drop = FALSE]
  top_loci <- head(top_loci, 3)
  lines <- c(
    paste("Total loci clustered:", nrow(report)),
    paste("Single-ancestry loci:", sum(safe_num(report$n_active_ancestries) == 1, na.rm = TRUE)),
    paste("Shared across ancestry subsets:", sum(safe_num(report$n_active_ancestries) >= 2 & safe_num(report$n_active_ancestries) < length(ANCESTRIES), na.rm = TRUE)),
    paste("Shared across all ancestries:", sum(safe_num(report$n_active_ancestries) == length(ANCESTRIES), na.rm = TRUE)),
    paste("No-clear-signal loci:", sum(safe_num(report$n_active_ancestries) == 0, na.rm = TRUE)),
    paste("Most common pattern:", unname(CLASSIFICATION_LABELS[names(class_counts)[1]]))
  )
  if (nrow(top_loci)) {
    lines <- c(lines, paste("Top lead loci:", paste(paste0(top_loci$locus_name, " (", top_loci$lead_gene, ")"), collapse = ", ")))
  }
  lines
}

build_cluster_feature_matrix <- function(report) {
  mat <- data.frame(
    n_active_ancestries = safe_num(report$n_active_ancestries),
    n_snps_passing_filters = safe_num(report$n_snps_passing_filters),
    n_local_support_snps = safe_num(report$n_local_support_snps),
    stringsAsFactors = FALSE
  )

  total <- safe_num(report$n_snps_total)
  total[is.na(total) | total == 0] <- NA_real_
  mat$fraction_snps_passing <- safe_num(report$n_snps_passing_filters) / total

  for (anc in ANCESTRIES) {
    mat[[paste0("n_signal_", anc)]] <- safe_num(report[[paste0("n_signal_snps_", anc)]])
    pvals <- safe_num(report[[paste0("best_p_dosage_", anc)]])
    pvals[is.na(pvals) | pvals <= 0] <- 1
    ors <- safe_num(report[[paste0("best_OR_dosage_", anc)]])
    ors[is.na(ors) | ors <= 0] <- 1
    mat[[paste0("logp_", anc)]] <- -log10(pmax(pvals, 1e-300))
    mat[[paste0("abslogor_", anc)]] <- abs(log(pmax(ors, 1e-12)))
  }

  localp <- safe_num(report$best_p_local_ancestry)
  localp[is.na(localp) | localp <= 0] <- 1
  mat$logp_local <- -log10(pmax(localp, 1e-300))

  for (j in seq_len(ncol(mat))) {
    vals <- mat[[j]]
    vals[!is.finite(vals)] <- NA_real_
    med <- suppressWarnings(stats::median(vals, na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    vals[is.na(vals)] <- med
    mat[[j]] <- vals
  }

  scaled <- scale(mat)
  scaled[is.na(scaled)] <- 0
  scaled
}

apply_common_sheet_format <- function(wb, sheet, df) {
  header_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#1F4E78", halign = "center", textDecoration = "bold", border = "Bottom")
  body_style <- createStyle(valign = "top", wrapText = TRUE)

  addStyle(wb, sheet, header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
  if (nrow(df) > 0) {
    addStyle(wb, sheet, body_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
  }
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:ncol(df), widths = "auto")
}

apply_number_styles <- function(wb, sheet, df) {
  if (!nrow(df)) return(invisible(NULL))
  sci_style <- createStyle(numFmt = "0.00E+00")
  or_style <- createStyle(numFmt = "0.00")
  maf_style <- createStyle(numFmt = "0.000")
  int_style <- createStyle(numFmt = "#,##0")

  for (j in seq_along(names(df))) {
    nm <- names(df)[j]
    rows <- 2:(nrow(df) + 1)
    if (grepl("^(best_p|p_)", nm)) {
      addStyle(wb, sheet, sci_style, rows = rows, cols = j, gridExpand = TRUE, stack = TRUE)
    } else if (grepl("^(best_OR|OR_)", nm)) {
      addStyle(wb, sheet, or_style, rows = rows, cols = j, gridExpand = TRUE, stack = TRUE)
    } else if (nm == "cohort_maf" || grepl("_maf$", nm)) {
      addStyle(wb, sheet, maf_style, rows = rows, cols = j, gridExpand = TRUE, stack = TRUE)
    } else if (nm %in% c("start_pos", "end_pos", "boundary_start", "boundary_end", "lead_pos", "n_snps_total", "n_snps_passing_filters", "n_local_support_snps", "n_active_ancestries") || grepl("^n_signal_snps_", nm)) {
      addStyle(wb, sheet, int_style, rows = rows, cols = j, gridExpand = TRUE, stack = TRUE)
    }
  }
}

write_excel_report <- function(locus_report, filtered_snps, parameters, output_path) {
  wb <- createWorkbook(creator = "GitHub Copilot")
  conclusions <- build_conclusions_table(locus_report)

  addWorksheet(wb, "Classification")
  writeDataTable(wb, "Classification", locus_report, withFilter = TRUE, tableStyle = "TableStyleMedium2")
  apply_common_sheet_format(wb, "Classification", locus_report)
  apply_number_styles(wb, "Classification", locus_report)

  if (nrow(locus_report) && "classification" %in% names(locus_report)) {
    for (i in seq_len(nrow(locus_report))) {
      cls <- as.character(locus_report$classification[i])
      fill <- unname(CLASSIFICATION_COLORS[ifelse(cls %in% names(CLASSIFICATION_COLORS), cls, "no_clear_signal")])
      row_style <- createStyle(fgFill = fill)
      addStyle(wb, "Classification", row_style, rows = i + 1, cols = 1:ncol(locus_report), gridExpand = TRUE, stack = TRUE)
    }
  }

  addWorksheet(wb, "Parameters")
  writeDataTable(wb, "Parameters", parameters, withFilter = TRUE, tableStyle = "TableStyleMedium9")
  apply_common_sheet_format(wb, "Parameters", parameters)

  addWorksheet(wb, "Conclusions")
  writeDataTable(wb, "Conclusions", conclusions, withFilter = TRUE, tableStyle = "TableStyleMedium4")
  apply_common_sheet_format(wb, "Conclusions", conclusions)

  addWorksheet(wb, "Filtered_SNPs")
  writeDataTable(wb, "Filtered_SNPs", filtered_snps, withFilter = TRUE, tableStyle = "TableStyleLight9")
  apply_common_sheet_format(wb, "Filtered_SNPs", filtered_snps)
  apply_number_styles(wb, "Filtered_SNPs", filtered_snps)

  saveWorkbook(wb, output_path, overwrite = TRUE)
}

render_dendrogram_svg <- function(report, output_path, method = "ward", metric = "euclidean") {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  if (!nrow(report)) {
    grDevices::png(output_path, width = 1600, height = 520, res = 150, type = "cairo")
    par(mar = c(1, 1, 3, 1))
    plot.new()
    text(0.05, 0.7, "No loci were available for clustering.", adj = c(0, 0.5), cex = 1.3)
    dev.off()
    return(invisible(NULL))
  }

  report$lead_gene[!not_blank(report$lead_gene)] <- "NA"
  computed_region <- mapply(format_region, report$chr, safe_num(report$boundary_start), safe_num(report$boundary_end))
  report$region <- ifelse(not_blank(report$region), report$region, computed_region)
  labels <- paste0(report$locus_name, " | ", report$lead_gene)
  colors <- unname(CLASSIFICATION_COLORS[ifelse(report$classification %in% names(CLASSIFICATION_COLORS), report$classification, "no_clear_signal")])

  if (nrow(report) == 1) {
    grDevices::png(output_path, width = 1900, height = 650, res = 150, type = "cairo")
    par(mar = c(1, 1, 4, 1))
    plot.new()
    title(main = "Locus ancestry dendrogram", sub = "Only one locus passed into the clustering stage.", line = 1)
    points(0.08, 0.60, pch = 22, cex = 2.2, bg = colors[1], col = "#333333")
    text(0.13, 0.60, labels[1], adj = c(0, 0.5), cex = 1)
    text(0.08, 0.33, paste0("Conclusion: ", CLASSIFICATION_NOTES[[report$classification[1]]]), adj = c(0, 0.5), cex = 0.95)
    dev.off()
    return(invisible(NULL))
  }

  features <- build_cluster_feature_matrix(report)
  metric <- if (metric %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) metric else "euclidean"
  method_map <- c(ward = "ward.D2", complete = "complete", average = "average", single = "single", mcquitty = "mcquitty", median = "median", centroid = "centroid")
  hclust_method <- if (!is.null(method_map[[method]])) method_map[[method]] else method

  hc <- hclust(dist(features, method = metric), method = hclust_method)
  ord <- hc$order
  summary_lines <- build_conclusion_lines(report)

  fig_height <- max(8, 0.42 * nrow(report) + 4)
  grDevices::png(output_path, width = 2700, height = max(1200, as.integer(fig_height * 150)), res = 150, type = "cairo")
  layout(matrix(c(1, 2), nrow = 2), heights = c(4.6, 1.4))

  par(mar = c(12, 4.5, 7.2, 6.5), xpd = NA)
  plot(
    hc,
    labels = FALSE,
    hang = -1,
    main = "",
    sub = "",
    xlab = "",
    ylab = "Cluster height"
  )
  title(main = "Locus ancestry dendrogram", line = 4.8)
  mtext(
    paste("Colored labels indicate ancestry classification. Linkage:", hclust_method, "| distance:", metric),
    side = 3,
    line = 3.2,
    cex = 0.92
  )
  usr <- par("usr")
  label_y <- usr[3] - 0.05 * diff(usr[3:4])
  x_pos <- seq_along(ord)
  points(x_pos, rep(label_y, length(x_pos)), pch = 21, bg = colors[ord], col = "#333333", cex = 1.25)
  text(x_pos, rep(label_y, length(x_pos)), labels = labels[ord], srt = 90, adj = 1, cex = 0.82, font = 2, col = colors[ord])
  legend("topright", inset = c(0.015, 0.02), legend = unname(CLASSIFICATION_LABELS), fill = unname(CLASSIFICATION_COLORS), border = "#555555", bg = "white", cex = 0.78, title = "Ancestry class")

  par(mar = c(1.5, 1.5, 0.8, 1.5))
  plot.new()
  rect(0.01, 0.06, 0.99, 0.96, col = "#F8FAFD", border = "#CBD5E1", lwd = 1.2)
  text(0.03, 0.90, "Conclusions", adj = c(0, 0.5), font = 2, cex = 1.1)
  y_vals <- seq(0.78, 0.16, length.out = length(summary_lines))
  for (i in seq_along(summary_lines)) {
    text(0.035, y_vals[i], paste0("• ", summary_lines[i]), adj = c(0, 0.5), cex = 0.95)
  }

  dev.off()
}

render_locus_panel_tree <- function(report, filtered_snps, output_path, method = "ward", metric = "euclidean") {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  if (!nrow(report)) {
    grDevices::png(output_path, width = 1800, height = 900, res = 150, type = "cairo")
    par(mar = c(1, 1, 3, 1))
    plot.new()
    text(0.05, 0.7, "No loci were available for the combined Manhattan/tree plot.", adj = c(0, 0.5), cex = 1.3)
    dev.off()
    return(invisible(NULL))
  }

  report$lead_gene[!not_blank(report$lead_gene)] <- "NA"
  plot_ancestries <- ANCESTRIES
  ancestry_colors <- setNames(rep("#4C78A8", length(plot_ancestries)), plot_ancestries)
  for (anc in plot_ancestries) {
    key <- paste0(anc, "_specific")
    if (key %in% names(CLASSIFICATION_COLORS)) ancestry_colors[[anc]] <- CLASSIFICATION_COLORS[[key]]
  }
  threshold_levels <- c(`5e-8` = -log10(5e-8), `5e-5` = -log10(5e-5))
  threshold_colors <- c(`5e-8` = "#D62728", `5e-5` = "#F1B61B")
  method_map <- c(ward = "ward.D2", complete = "complete", average = "average", single = "single", mcquitty = "mcquitty", median = "median", centroid = "centroid")
  hclust_method <- if (!is.null(method_map[[method]])) method_map[[method]] else method
  distance_metric <- if (metric %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) metric else "euclidean"

  if (nrow(report) > 1) {
    features <- build_cluster_feature_matrix(report)
    hc <- hclust(dist(features, method = distance_metric), method = hclust_method)
    dend <- as.dendrogram(hc)
    ord <- rev(order.dendrogram(dend))
  } else {
    hc <- NULL
    dend <- NULL
    ord <- 1L
  }

  ordered_report <- report[ord, , drop = FALSE]
  row_labels <- paste0(ordered_report$locus_name, "\n", ordered_report$lead_gene)
  row_colors <- unname(CLASSIFICATION_COLORS[ifelse(ordered_report$classification %in% names(CLASSIFICATION_COLORS), ordered_report$classification, "no_clear_signal")])

  all_logp <- unname(threshold_levels)
  for (anc in plot_ancestries) {
    vals <- safe_num(filtered_snps[[paste0("p_dosage_", anc)]])
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals)) all_logp <- c(all_logp, -log10(pmax(vals, 1e-300)))
  }
  ymax <- max(all_logp, na.rm = TRUE)
  if (!is.finite(ymax)) ymax <- 8
  ymax <- max(8, min(ymax + 0.4, 35))

  n <- nrow(ordered_report)
  n_panels <- length(plot_ancestries)
  n_cols <- n_panels + 1
  layout_mat <- matrix(0, nrow = n + 2, ncol = n_cols)
  layout_mat[1, ] <- 1
  layout_mat[2, ] <- seq.int(2, n_cols + 1)
  next_id <- n_cols + 2
  for (i in seq_len(n)) {
    layout_mat[i + 2, seq_len(n_panels)] <- seq.int(next_id, next_id + n_panels - 1)
    next_id <- next_id + n_panels
  }
  tree_panel_id <- next_id
  layout_mat[3:(n + 2), n_cols] <- tree_panel_id

  img_height <- max(2000, 185 * (n + 3))
  grDevices::png(output_path, width = 3600, height = img_height, res = 170, type = "cairo")
  panel_widths <- c(rep(3.0, n_panels), 4.1)
  layout(layout_mat, widths = panel_widths, heights = c(0.95, 0.82, rep(1.0, n)))
  par(oma = c(1.2, 0.5, 1.2, 0.5))

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  text(0.5, 0.78, "Locus clustering with ancestry-specific Manhattan panels", cex = 1.30, font = 2)
  legend(
    x = 0.5,
    y = 0.28,
    legend = unname(CLASSIFICATION_LABELS),
    fill = unname(CLASSIFICATION_COLORS),
    horiz = TRUE,
    bty = "n",
    xjust = 0.5,
    yjust = 0.5,
    cex = 0.72,
    inset = 0
  )

  for (anc in plot_ancestries) {
    par(mar = if (anc == plot_ancestries[1]) c(0.2, 10.0, 0.6, 0.45) else c(0.2, 0.65, 0.6, 0.45))
    plot.new()
    rect(0.06, 0.16, 0.94, 0.84, col = "white", border = ancestry_colors[[anc]], lwd = 1.7)
    text(0.5, 0.50, anc, cex = 1.25, font = 2, col = ancestry_colors[[anc]])
  }
  par(mar = c(0.2, 0.8, 0.6, 0.8))
  plot.new()
  rect(0.06, 0.16, 0.94, 0.84, col = "white", border = "#666666", lwd = 1.5)
  text(0.5, 0.50, "Clustering tree", cex = 1.1, font = 2)

  for (i in seq_len(n)) {
    locus <- ordered_report$locus_name[i]
    locus_snps <- filtered_snps[filtered_snps$locus_name == locus, , drop = FALSE]
    xvals <- safe_num(locus_snps$pos)
    xvals <- xvals[is.finite(xvals)]
    if (!length(xvals)) xvals <- c(0, 1)
    xlim <- range(xvals, na.rm = TRUE)
    if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)

    for (anc in plot_ancestries) {
      par(mar = if (anc == plot_ancestries[1]) c(0.35, 10.0, 0.15, 0.45) else c(0.35, 0.65, 0.15, 0.45), xaxs = "i", yaxs = "i", xpd = FALSE)
      plot(xlim, c(0, ymax), type = "n", axes = FALSE, xlab = "", ylab = "")
      line_x1 <- xlim[1] + 0.01 * diff(xlim)
      line_x2 <- xlim[2] - 0.01 * diff(xlim)
      segments(line_x1, threshold_levels["5e-5"], line_x2, threshold_levels["5e-5"], col = threshold_colors["5e-5"], lty = 2, lwd = 1)
      segments(line_x1, threshold_levels["5e-8"], line_x2, threshold_levels["5e-8"], col = threshold_colors["5e-8"], lty = 2, lwd = 1)

      pvals <- safe_num(locus_snps[[paste0("p_dosage_", anc)]])
      xpos <- safe_num(locus_snps$pos)
      ypos <- -log10(pmax(pvals, 1e-300))
      keep <- is.finite(xpos) & is.finite(ypos)
      if (sum(keep)) {
        points(xpos[keep], pmin(ypos[keep], ymax), pch = 16, cex = 0.54, col = ancestry_colors[[anc]])
      }
      box(col = "#C9CED6", lwd = 0.9)

      usr <- par("usr")
      rect(usr[2] - 0.018 * diff(usr[1:2]), usr[3], usr[2], usr[4], col = row_colors[i], border = NA)

      if (anc == plot_ancestries[1]) {
        text(xlim[1] - 0.060 * diff(xlim), ymax / 2, row_labels[i], xpd = NA, srt = 90, adj = 0.5, cex = 1.06, font = 2, col = "black")
      }
    }
  }

  par(mar = c(2.5, 0.8, 0.15, 0.8), xpd = NA)
  if (!is.null(dend) && n > 1) {
    tree_limit <- max(hc$height, na.rm = TRUE) * 1.03
    plot(dend, horiz = TRUE, leaflab = "none", yaxt = "n", xaxt = "n", axes = FALSE, main = "", xlab = "", xlim = c(tree_limit, 0))
    axis(1, at = pretty(c(0, tree_limit)), labels = rev(pretty(c(0, tree_limit))), cex.axis = 0.72)
    box(col = "#777777")
  } else {
    plot.new()
    plot.window(xlim = c(1, 0), ylim = c(0.5, 1.5), xaxs = "i", yaxs = "i")
    segments(1, 1, 0.35, 1, lwd = 1.4, col = "#444444")
    points(1, 1, pch = 15, cex = 1.0, col = row_colors[1])
    axis(1, at = c(1, 0.5, 0), labels = c("0", "0.5", "1"), cex.axis = 0.72)
    box(col = "#777777")
  }
  title(xlab = "Height", line = 1.0)
}

main <- function() {
  args <- parse_args()
  required <- c("input_locus", "input_snp", "output_excel", "output_dendrogram", "output_panel_tree")
  missing <- required[!required %in% names(args)]
  if (length(missing)) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  locus_report <- read.delim(args$input_locus, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  filtered_snps <- read.delim(args$input_snp, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  detected_ancestries <- sub("^n_signal_snps_", "", grep("^n_signal_snps_", names(locus_report), value = TRUE))
  if (!length(detected_ancestries)) {
    detected_ancestries <- sub("^best_p_dosage_", "", grep("^best_p_dosage_", names(locus_report), value = TRUE))
  }
  if (!length(detected_ancestries)) {
    detected_ancestries <- sub("_maf$", "", grep("_maf$", names(filtered_snps), value = TRUE))
    detected_ancestries <- setdiff(detected_ancestries, "cohort")
  }
  if (length(detected_ancestries)) {
    ANCESTRIES <<- unique(detected_ancestries)
  }
  meta <- build_classification_metadata(ANCESTRIES)
  CLASSIFICATION_LABELS <<- meta$labels
  CLASSIFICATION_NOTES <<- meta$notes
  CLASSIFICATION_COLORS <<- meta$colors

  if (!"classification_label" %in% names(locus_report) && "classification" %in% names(locus_report)) {
    locus_report$classification_label <- unname(ifelse(locus_report$classification %in% names(CLASSIFICATION_LABELS), CLASSIFICATION_LABELS[locus_report$classification], locus_report$classification))
  }
  if (!"classification_note" %in% names(locus_report) && "classification" %in% names(locus_report)) {
    locus_report$classification_note <- unname(ifelse(locus_report$classification %in% names(CLASSIFICATION_NOTES), CLASSIFICATION_NOTES[locus_report$classification], locus_report$classification))
  }

  if (!"genes" %in% names(locus_report) && "lead_gene" %in% names(locus_report)) {
    locus_report$genes <- locus_report$lead_gene
  }

  if (!"boundary_start" %in% names(locus_report)) {
    locus_report$boundary_start <- if ("start_pos" %in% names(locus_report)) safe_num(locus_report$start_pos) else NA_real_
  }
  if (!"boundary_end" %in% names(locus_report)) {
    locus_report$boundary_end <- if ("end_pos" %in% names(locus_report)) safe_num(locus_report$end_pos) else NA_real_
  }

  if (!is.null(args$loci_file) && file.exists(args$loci_file)) {
    loci_meta <- read.delim(args$loci_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    names(loci_meta)[names(loci_meta) == "locus"] <- "locus_name"
    if (all(c("locus_name", "chr", "start", "end") %in% names(loci_meta)) && all(c("locus_name", "chr") %in% names(locus_report))) {
      loci_meta$chr <- as.character(loci_meta$chr)
      locus_report$chr <- as.character(locus_report$chr)
      merged <- merge(locus_report, loci_meta[, c("locus_name", "chr", "start", "end")], by = c("locus_name", "chr"), all.x = TRUE, sort = FALSE, suffixes = c("", "_loci"))
      merged$boundary_start <- ifelse(!is.na(safe_num(merged$start)), safe_num(merged$start), safe_num(merged$boundary_start))
      merged$boundary_end <- ifelse(!is.na(safe_num(merged$end)), safe_num(merged$end), safe_num(merged$boundary_end))
      merged$start <- NULL
      merged$end <- NULL
      locus_report <- merged
    }
  }

  if (!"region" %in% names(locus_report)) {
    locus_report$region <- mapply(format_region, locus_report$chr, safe_num(locus_report$boundary_start), safe_num(locus_report$boundary_end))
  }

  dynamic_anc_cols <- unlist(lapply(ANCESTRIES, function(anc) c(
    paste0("n_signal_snps_", anc),
    paste0("lead_rsid_", anc),
    paste0("best_p_dosage_", anc),
    paste0("best_OR_dosage_", anc)
  )))

  preferred_order <- c(
    "chr", "locus_name", "region", "boundary_start", "boundary_end", "genes",
    "n_snps_total", "n_snps_passing_filters", "n_local_support_snps",
    "active_ancestries", "n_active_ancestries", "classification", "classification_label", "classification_note",
    "lead_rsid", "lead_gene", "lead_pos", "lead_best_ancestry", "lead_best_p_dosage",
    "best_p_local_ancestry"
  )
  preferred_order <- c(preferred_order[preferred_order != "best_p_local_ancestry"], dynamic_anc_cols, "best_p_local_ancestry")
  keep_cols <- preferred_order[preferred_order %in% names(locus_report)]
  other_cols <- setdiff(names(locus_report), keep_cols)
  locus_report <- locus_report[, c(keep_cols, other_cols), drop = FALSE]

  parameters <- build_parameters_table(args)

  dir.create(dirname(args$output_excel), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$output_dendrogram), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$output_panel_tree), recursive = TRUE, showWarnings = FALSE)
  write_excel_report(locus_report, filtered_snps, parameters, args$output_excel)
  render_dendrogram_svg(locus_report, args$output_dendrogram, args$clustering_method %||% "ward", args$clustering_metric %||% "euclidean")
  render_locus_panel_tree(locus_report, filtered_snps, args$output_panel_tree, args$clustering_method %||% "ward", args$clustering_metric %||% "euclidean")

  cat(sprintf("Wrote R-based Excel report to %s, dendrogram to %s, and combined Manhattan/tree plot to %s\n", args$output_excel, args$output_dendrogram, args$output_panel_tree))
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) && nzchar(x)) x else y

main()
