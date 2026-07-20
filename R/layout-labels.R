compute_point_positions <- function(data, group_order, point_layout, spread_by,
                                    spread_width, jitter_width, seed) {
  data$.smv_group <- factor(data$.smv_group, levels = group_order)
  data$.smv_group_position <- match(as.character(data$.smv_group), group_order)

  if (point_layout == "jitter") {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    data$.smv_x <- data$.smv_group_position + stats::runif(nrow(data), -jitter_width, jitter_width)
    return(data)
  }

  spread_variables <- if (spread_by == "direction") c(".smv_group", ".smv_direction") else ".smv_group"
  data <- dplyr::group_by(data, dplyr::across(dplyr::all_of(spread_variables)))
  data <- dplyr::arrange(data, .data$.smv_gene, .data$.smv_row_id, .by_group = TRUE)
  data <- dplyr::mutate(
    data,
    .smv_spread_rank = dplyr::row_number(),
    .smv_spread_n = dplyr::n(),
    .smv_offset = ifelse(
      .data$.smv_spread_n == 1L,
      0,
      ((.data$.smv_spread_rank - 1) / (.data$.smv_spread_n - 1) - 0.5) * spread_width
    ),
    .smv_x = .data$.smv_group_position + .data$.smv_offset
  )
  dplyr::ungroup(data)
}

select_marker_labels <- function(data, shown_data, label, label_n, label_by,
                                 label_from, label_direction, label_genes,
                                 missing_labels, significance_by) {
  if (label == "none" || label_n == 0L && is.null(label_genes)) return(data[0, , drop = FALSE])

  if (label == "custom") {
    if (is.null(label_genes)) .smv_abort("`label_genes` is required when `label = \"custom\"`.", "scMarkerViz_missing_labels")
    requested <- normalize_requested_labels(label_genes)
    if (is.null(requested$group)) {
      selected <- data[data$.smv_gene %in% requested$gene, , drop = FALSE]
      found <- unique(selected$.smv_gene)
      missing <- setdiff(unique(requested$gene), found)
    } else {
      keys <- paste(data$.smv_group, data$.smv_gene, sep = "\r")
      requested_keys <- paste(requested$group, requested$gene, sep = "\r")
      selected <- data[keys %in% requested_keys, , drop = FALSE]
      missing_keys <- setdiff(unique(requested_keys), unique(keys))
      missing <- gsub("\r", "/", missing_keys, fixed = TRUE)
    }
    handle_missing_labels(missing, missing_labels)
    return(selected[!duplicated(selected[c(".smv_group", ".smv_gene")]), , drop = FALSE])
  }

  candidates <- switch(label_from,
    significant = data[data$.smv_significant, , drop = FALSE],
    shown = shown_data,
    all = data
  )
  if (label_direction != "both") {
    candidates <- candidates[as.character(candidates$.smv_direction) == label_direction, , drop = FALSE]
  } else {
    candidates <- candidates[as.character(candidates$.smv_direction) %in% c("up", "down"), , drop = FALSE]
  }
  if (nrow(candidates) == 0L) return(candidates)

  candidates$.smv_label_direction <- as.character(candidates$.smv_direction)
  if (label_by == "significance" && significance_by == "none") {
    .smv_abort("`label_by = \"significance\"` requires a P-value significance mode.", "scMarkerViz_invalid_option")
  }
  selected_p <- if (significance_by == "p_adjust") candidates$.smv_p_adjust else candidates$.smv_p_value
  candidates$.smv_label_score <- if (label_by == "effect") abs(candidates$.smv_effect) else -selected_p
  candidates <- dplyr::group_by(candidates, .data$.smv_group, .data$.smv_label_direction)
  candidates <- dplyr::arrange(candidates, dplyr::desc(.data$.smv_label_score), .data$.smv_gene, .data$.smv_row_id, .by_group = TRUE)
  candidates <- dplyr::slice_head(candidates, n = label_n)
  candidates <- dplyr::ungroup(candidates)
  candidates[!duplicated(candidates[c(".smv_group", ".smv_gene")]), , drop = FALSE]
}

normalize_requested_labels <- function(label_genes) {
  if (is.character(label_genes)) return(list(group = NULL, gene = label_genes))
  if (is.list(label_genes) && !is.data.frame(label_genes)) {
    if (is.null(names(label_genes)) || any(!nzchar(names(label_genes)))) {
      .smv_abort("A list supplied to `label_genes` must be named by group.", "scMarkerViz_invalid_labels")
    }
    return(list(
      group = rep(names(label_genes), lengths(label_genes)),
      gene = unlist(label_genes, use.names = FALSE)
    ))
  }
  if (is.data.frame(label_genes) && all(c("group", "gene") %in% names(label_genes))) {
    return(list(group = as.character(label_genes$group), gene = as.character(label_genes$gene)))
  }
  .smv_abort("`label_genes` must be a character vector, named list, or data frame with `group` and `gene` columns.", "scMarkerViz_invalid_labels")
}

handle_missing_labels <- function(missing, action) {
  if (!length(missing) || action == "ignore") return(invisible(NULL))
  message <- sprintf("Requested marker(s) were not found: %s.", paste(missing, collapse = ", "))
  if (action == "error") .smv_abort(message, "scMarkerViz_missing_labels")
  rlang::warn(message)
}

build_group_layout <- function(data, group_order, group_colors, effect_cutoff) {
  groups <- data.frame(
    .smv_group = factor(group_order, levels = group_order),
    .smv_group_position = seq_along(group_order),
    stringsAsFactors = FALSE
  )
  ranges <- dplyr::group_by(data, .data$.smv_group)
  ranges <- dplyr::summarise(
    ranges,
    .smv_ymin = min(.data$.smv_effect),
    .smv_ymax = max(.data$.smv_effect),
    .groups = "drop"
  )
  groups <- dplyr::left_join(groups, ranges, by = ".smv_group")
  groups$.smv_xmin <- groups$.smv_group_position - 0.5
  groups$.smv_xmax <- groups$.smv_group_position + 0.5
  groups$.smv_color <- unname(group_colors[as.character(groups$.smv_group)])
  groups$.smv_band_min <- -effect_cutoff
  groups$.smv_band_max <- effect_cutoff
  groups
}
