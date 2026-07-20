.smv_abort <- function(message, class = NULL) {
  rlang::abort(message, class = c(class, "scMarkerViz_error"))
}

.smv_match_column <- function(data, supplied, candidates, argument, required = TRUE) {
  if (!is.null(supplied)) {
    if (!rlang::is_string(supplied) || !supplied %in% names(data)) {
      .smv_abort(sprintf("Column `%s` selected by `%s` was not found.", supplied, argument), "scMarkerViz_missing_column")
    }
    return(supplied)
  }

  found <- intersect(candidates, names(data))
  if (length(found) == 1L) {
    return(found)
  }
  if (length(found) > 1L) {
    .smv_abort(sprintf(
      "Multiple candidate columns were found for `%s`: %s. Select one explicitly.",
      argument, paste(sprintf("`%s`", found), collapse = ", ")
    ), "scMarkerViz_ambiguous_column")
  }
  if (required) {
    .smv_abort(sprintf(
      "Could not identify `%s`. Supply its column name explicitly. Available columns: %s.",
      argument, paste(sprintf("`%s`", names(data)), collapse = ", ")
    ), "scMarkerViz_missing_column")
  }
  NULL
}

normalize_marker_data <- function(data, gene, group, effect, p_value, p_adjust,
                                  significance_by) {
  if (!is.data.frame(data)) {
    .smv_abort("`data` must be a data frame.", "scMarkerViz_invalid_data")
  }
  if (nrow(data) == 0L) {
    .smv_abort("`data` must contain at least one row.", "scMarkerViz_empty_data")
  }

  gene <- .smv_match_column(data, gene, c("gene", "feature", "symbol", "names"), "gene")
  group <- .smv_match_column(data, group, c("cluster", "celltype", "cell_type", "group"), "group")
  effect <- .smv_match_column(data, effect, c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "logfoldchanges"), "effect")
  p_value <- .smv_match_column(data, p_value, c("p_val", "pvalue", "p_value", "pvals"), "p_value", required = significance_by == "p_value")
  p_adjust <- .smv_match_column(data, p_adjust, c("p_val_adj", "padj", "p_adjust", "FDR", "pvals_adj"), "p_adjust", required = significance_by == "p_adjust")

  result <- data
  result$.smv_gene <- as.character(data[[gene]])
  result$.smv_group <- as.character(data[[group]])
  result$.smv_effect <- data[[effect]]
  result$.smv_p_value <- if (is.null(p_value)) NA_real_ else data[[p_value]]
  result$.smv_p_adjust <- if (is.null(p_adjust)) NA_real_ else data[[p_adjust]]
  result$.smv_row_id <- seq_len(nrow(data))

  if (!is.numeric(result$.smv_effect)) {
    .smv_abort(sprintf("Effect column `%s` must be numeric.", effect), "scMarkerViz_invalid_column")
  }
  for (column in c(".smv_p_value", ".smv_p_adjust")) {
    values <- result[[column]]
    if (!is.numeric(values)) {
      .smv_abort(sprintf("Column represented by `%s` must be numeric.", column), "scMarkerViz_invalid_column")
    }
    invalid <- is.finite(values) & (values < 0 | values > 1)
    if (any(invalid, na.rm = TRUE)) {
      .smv_abort("P values must lie between 0 and 1.", "scMarkerViz_invalid_p_value")
    }
  }

  invalid_key <- is.na(result$.smv_gene) | !nzchar(result$.smv_gene) |
    is.na(result$.smv_group) | !nzchar(result$.smv_group)
  invalid_effect <- !is.finite(result$.smv_effect)
  dropped <- sum(invalid_key | invalid_effect)
  if (dropped > 0L) {
    rlang::warn(sprintf("Removed %d row(s) with missing gene/group values or non-finite effects.", dropped))
    result <- result[!(invalid_key | invalid_effect), , drop = FALSE]
  }
  if (nrow(result) == 0L) {
    .smv_abort("No valid rows remain after data validation.", "scMarkerViz_empty_data")
  }
  result
}

classify_markers <- function(data, effect_cutoff, significance_cutoff, significance_by) {
  selected_p <- switch(significance_by,
    p_adjust = data$.smv_p_adjust,
    p_value = data$.smv_p_value,
    none = rep(0, nrow(data))
  )
  passes_p <- if (significance_by == "none") rep(TRUE, nrow(data)) else !is.na(selected_p) & selected_p <= significance_cutoff
  significant <- abs(data$.smv_effect) >= effect_cutoff & passes_p
  direction <- ifelse(data$.smv_effect > 0, "up", ifelse(data$.smv_effect < 0, "down", "zero"))
  status <- ifelse(significant & direction == "up", "up", ifelse(significant & direction == "down", "down", "not_significant"))

  data$.smv_direction <- factor(direction, levels = c("down", "zero", "up"))
  data$.smv_significant <- significant
  data$.smv_status <- factor(status, levels = c("down", "not_significant", "up"))
  data
}

validate_group_order <- function(data, group_order) {
  groups <- unique(data$.smv_group)
  if (is.null(group_order)) return(groups)
  group_order <- as.character(group_order)
  if (anyDuplicated(group_order)) .smv_abort("`group_order` must not contain duplicates.", "scMarkerViz_invalid_group_order")
  missing <- setdiff(groups, group_order)
  unknown <- setdiff(group_order, groups)
  if (length(missing) || length(unknown)) {
    .smv_abort(sprintf(
      "`group_order` must contain every observed group exactly once. Missing: %s. Unknown: %s.",
      if (length(missing)) paste(missing, collapse = ", ") else "none",
      if (length(unknown)) paste(unknown, collapse = ", ") else "none"
    ), "scMarkerViz_invalid_group_order")
  }
  group_order
}
