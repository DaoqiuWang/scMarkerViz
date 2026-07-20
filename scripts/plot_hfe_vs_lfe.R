#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else "D:/Users/Downloads/scRNAtoolVis-master/HFE_vs_LFE_overall_all_DEGs.csv"
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(getwd(), "output")
package_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "dplyr", "ggrepel", "rlang")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

pkgload_available <- requireNamespace("pkgload", quietly = TRUE)
if (pkgload_available) {
  pkgload::load_all(package_dir, quiet = TRUE)
} else {
  source(file.path(package_dir, "R", "data-core.R"))
  source(file.path(package_dir, "R", "layout-labels.R"))
  source(file.path(package_dir, "R", "marker-effect-plot.R"))
}

markers <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("gene", "celltype", "avg_log2FC", "p_val_adj")
missing_columns <- setdiff(required_columns, names(markers))
if (length(missing_columns)) {
  stop("Input is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

message("Loaded ", nrow(markers), " rows across ", dplyr::n_distinct(markers$celltype), " cell types.")
message("Adjusted P-value missing rows: ", sum(is.na(markers$p_val_adj)))

common_options <- list(
  data = markers,
  gene = "gene",
  group = "celltype",
  effect = "avg_log2FC",
  p_value = "p_val",
  p_adjust = "p_val_adj",
  effect_cutoff = 0.25,
  significance_cutoff = 0.05,
  significance_by = "p_adjust",
  show = "all",
  color_by = "status",
  label = "top",
  label_n = 5,
  label_by = "effect",
  label_from = "significant",
  point_layout = "spread",
  spread_by = "all",
  spread_width = 0.8,
  point_size = 0.8,
  group_label_angle = 70,
  base_size = 12
)

vertical <- do.call(marker_effect_plot, c(common_options, list(layout = "vertical"))) +
  ggplot2::labs(
    title = "Differential expression across cell types",
    subtitle = "HFE vs LFE — top 5 significant genes per direction",
    y = "Average log2 fold change"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0.5)
  )

horizontal <- do.call(marker_effect_plot, c(common_options, list(layout = "horizontal"))) +
  ggplot2::labs(
    title = "Differential expression across cell types",
    subtitle = "HFE vs LFE — horizontal layout",
    y = "Average log2 fold change"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0.5)
  )

n_groups <- dplyr::n_distinct(markers$celltype)
vertical_width <- max(12, n_groups * 0.85)
horizontal_height <- max(10, n_groups * 0.65)

vertical_pdf <- file.path(output_dir, "HFE_vs_LFE_marker_effect_vertical.pdf")
vertical_png <- file.path(output_dir, "HFE_vs_LFE_marker_effect_vertical.png")
horizontal_pdf <- file.path(output_dir, "HFE_vs_LFE_marker_effect_horizontal.pdf")

ggplot2::ggsave(vertical_pdf, vertical, width = vertical_width, height = 10, units = "in")
ggplot2::ggsave(vertical_png, vertical, width = vertical_width, height = 10, units = "in", dpi = 300)
ggplot2::ggsave(horizontal_pdf, horizontal, width = 12, height = horizontal_height, units = "in", limitsize = FALSE)

message("Saved: ", normalizePath(vertical_pdf, winslash = "/", mustWork = FALSE))
message("Saved: ", normalizePath(vertical_png, winslash = "/", mustWork = FALSE))
message("Saved: ", normalizePath(horizontal_pdf, winslash = "/", mustWork = FALSE))
