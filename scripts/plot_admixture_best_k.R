#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

option_list <- list(
  make_option(c("--best-k-file"), type = "character", help = "Path to best_k.txt"),
  make_option(c("--fam-file"), type = "character", help = "Merged .fam file"),
  make_option(c("--sample-groups"), type = "character", help = "Sample group TSV"),
  make_option(c("--q-dir"), type = "character", help = "Directory with admixture K Q files"),
  make_option(c("--cohort-label"), type = "character", default = "Cohort", help = "Cohort label"),
  make_option(c("--group-colors"), type = "character", default = "", help = "Comma-separated label=#hex"),
  make_option(c("--cohort-color"), type = "character", default = "#7f7f7f", help = "Cohort color"),
  make_option(c("--output-png"), type = "character", help = "Output PNG path")
)

opt <- parse_args(OptionParser(option_list = option_list))

best_k <- as.integer(readLines(opt$`best-k-file`, n = 1, warn = FALSE))
if (is.na(best_k)) {
  stop("Could not parse best K")
}

fam <- read.table(opt$`fam-file`, stringsAsFactors = FALSE)
if (ncol(fam) < 2) {
  stop("FAM file must contain at least two columns")
}
iids <- fam[[2]]

q_file <- file.path(opt$`q-dir`, paste0("k", best_k, ".Q"))
q_mat <- as.matrix(read.table(q_file, stringsAsFactors = FALSE))
if (nrow(q_mat) != length(iids)) {
  stop("Q rows do not match number of samples in FAM")
}

q_df <- as.data.frame(q_mat)
colnames(q_df) <- paste0("Ancestry_", seq_len(ncol(q_df)))
q_df$IID <- iids
q_df$sample_index <- seq_len(nrow(q_df))

groups <- read_tsv(opt$`sample-groups`, show_col_types = FALSE)
plot_df <- q_df %>% left_join(groups, by = "IID")
plot_df$GROUP[is.na(plot_df$GROUP)] <- opt$`cohort-label`

plot_df <- plot_df %>% arrange(GROUP, sample_index)
plot_df$plot_index <- seq_len(nrow(plot_df))

long_df <- plot_df %>%
  select(IID, GROUP, plot_index, starts_with("Ancestry_")) %>%
  pivot_longer(cols = starts_with("Ancestry_"), names_to = "Component", values_to = "Proportion")

group_color_map <- c()
if (nchar(opt$`group-colors`) > 0) {
  tokens <- strsplit(opt$`group-colors`, ",")[[1]]
  for (tok in tokens) {
    pair <- strsplit(tok, "=", fixed = TRUE)[[1]]
    if (length(pair) == 2) {
      group_color_map[pair[1]] <- pair[2]
    }
  }
}
if (!(opt$`cohort-label` %in% names(group_color_map))) {
  group_color_map[opt$`cohort-label`] <- opt$`cohort-color`
}

component_colors <- scales::hue_pal()(length(unique(long_df$Component)))
names(component_colors) <- sort(unique(long_df$Component))

p <- ggplot(long_df, aes(x = plot_index, y = Proportion, fill = Component)) +
  geom_col(width = 1.0) +
  geom_point(
    data = plot_df,
    aes(x = plot_index, y = -0.03, color = GROUP),
    inherit.aes = FALSE,
    shape = 15,
    size = 0.7
  ) +
  scale_fill_manual(values = component_colors) +
  scale_color_manual(values = group_color_map, drop = FALSE) +
  scale_y_continuous(limits = c(-0.06, 1.0), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(
    title = paste0("Unsupervised ADMIXTURE (Best K = ", best_k, ")"),
    x = "Samples (ordered by group)",
    y = "Ancestry proportion",
    fill = "Components",
    color = "Sample group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(opt$`output-png`, plot = p, width = 16, height = 7, dpi = 300)
