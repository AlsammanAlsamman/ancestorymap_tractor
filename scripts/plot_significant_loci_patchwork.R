#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(locuszoomr)
  library(AnnotationFilter)
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

safe_num <- function(x) suppressWarnings(as.numeric(x))

normalize_build <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", as.character(x)))
}

ensdb_pkg_for_build <- function(genome_build) {
  b <- normalize_build(genome_build)
  if (b %in% c("grch37", "hg19", "b37")) return("EnsDb.Hsapiens.v75")
  if (b %in% c("grch38", "hg38", "b38")) return("EnsDb.Hsapiens.v86")
  stop(sprintf("Unsupported genome build '%s'. Use one of: GRCh37/hg19 or GRCh38/hg38.", genome_build), call. = FALSE)
}

make_locuszoom_gg <- function(dt, p_col, title, ensdb_pkg, ylim_vals) {
  # Clean column set for locuszoomr
  # Requires: chrom, pos, p, rsid (optional)
  plot_dt <- dt[, .(chrom = CHR, pos = POS, p = get(p_col), rsid = ID)]
  
  # Initialize locus object
  loc <- tryCatch({
    locus(
      data = plot_dt,
      ens_db = ensdb_pkg,
      seqname = unique(dt$CHR),
      xrange = c(min(dt$POS), max(dt$POS))
    )
  }, error = function(e) {
    message("locus() failed for ", title, ": ", e$message)
    return(NULL)
  })

  if (is.null(loc)) return(ggplot() + labs(title = paste("Error:", title)) + theme_void())

  # Plot with locuszoomr, then convert to ggplot
  # Note: gg_locus() does the GWAS panel
  p <- gg_locus(loc) + 
       labs(title = title) +
       ylim(ylim_vals) +
       theme_bw() +
       theme(panel.grid.minor = element_blank())
  
  return(p)
}

main <- function() {
  harmonized_tsv <- get_arg("--harmonized-tsv", required = TRUE)
  ancestries     <- parse_csv(get_arg("--ancestries", "AFR,AMR,EAS,EUR"))
  genome_build   <- get_arg("--genome-build", "GRCh37")
  out_plot_dir   <- get_arg("--output-plot-dir", "results/loci_plots")
  out_excel      <- get_arg("--output-excel", "results/loci_report.xlsx")

  dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)
  ensdb_pkg <- ensdb_pkg_for_build(genome_build)
  
  message("Loading data: ", harmonized_tsv)
  dt <- fread(harmonized_tsv)
  
  # Identify significant loci (P < 5e-8 in any ancestry)
  p_cols <- paste0("P_", ancestries)
  dt[, min_p := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = p_cols]
  
  sig_snps <- dt[min_p < 5e-8]
  if (nrow(sig_snps) == 0) {
    message("No significant associations found.")
    file.create(out_excel)
    return()
  }

  # Cluster SNPs into loci (e.g. within 500kb)
  sig_snps <- sig_snps[order(CHR, POS)]
  sig_snps[, locus_id := cumsum(c(1, (diff(POS) > 500000) | (diff(CHR) != 0)))]
  
  # Prepare report
  report_list <- list()
  
  # Loop each locus
  loci_ids <- unique(sig_snps$locus_id)
  for (lid in loci_ids) {
    locus_snps <- sig_snps[locus_id == lid]
    lead_snp   <- locus_snps[min_p == min(min_p)][1]
    
    chrom     <- lead_snp$CHR
    start_pos <- max(0, lead_snp$POS - 500000)
    end_pos   <- lead_snp$POS + 500000
    loc_id    <- sprintf("chr%s_%s", chrom, lead_snp$POS)
    
    # Get all SNPs in region from main DT
    loc_dt <- dt[CHR == chrom & POS >= start_pos & POS <= end_pos]
    
    # Add to report
    report_list[[as.character(lid)]] <- locus_snps
    
    # Plotting
    message("Plotting locus: ", loc_id)
    pos_label <- sprintf("%s:%s-%s", chrom, start_pos, end_pos)
    
    # Global Y-axis limit for GWAS panels
    p_vals <- unlist(loc_dt[, ..p_cols])
    max_logp <- max(-log10(p_vals[p_vals > 0]), na.rm = TRUE)
    ylim_vals <- c(0, ceiling(max_logp + 1))

    # Single locuszoomr object for gene track
    loc_for_genes <- tryCatch(
      locus(
        data = loc_dt[, .(chrom=CHR, pos=POS, p=get(p_cols[1]), rsid=ID)],
        seqname = chrom,
        xrange = c(start_pos, end_pos),
        ens_db = ensdb_pkg
      ),
      error = function(e) NULL
    )

    # Panels
    panels <- lapply(ancestries, function(anc) {
      panel_title <- sprintf("Ancestry: %s  |  Lead: %s  |  %s", anc, lead_snp$ID, pos_label)
      make_locuszoom_gg(loc_dt, paste0("P_", anc), panel_title, ensdb_pkg, ylim_vals)
    })

    # Combined layout
    final_plot <- wrap_plots(panels, ncol = 1)
    
    if (!is.null(loc_for_genes)) {
      gene_panel <- tryCatch(gg_addgenes(loc_for_genes), error = function(e) NULL)
      if (!is.null(gene_panel)) {
        # Using a very safe height allocation
        final_plot <- final_plot / gene_panel + plot_layout(heights = c(rep(1, length(panels)), 0.5))
      }
    }

    out_file <- file.path(out_plot_dir, paste0(gsub(":", "_", loc_id), ".pdf"))
    
    # Try multiple dimensions if it fails
    tryCatch({
      ggsave(out_file, plot = final_plot, width = 12, height = 4 * length(panels) + 3, device = "pdf")
    }, error = function(e) {
      message("First ggsave failed, trying larger dimensions...")
      ggsave(out_file, plot = final_plot, width = 18, height = 6 * length(panels) + 10, device = "pdf")
    })
  }
  
  # Save Excel
  write.xlsx(report_list, file = out_excel)
  message("Done.")
}

main()
