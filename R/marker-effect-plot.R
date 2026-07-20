#' Plot marker effects across cell groups
#'
#' Creates a cluster-wise marker effect plot from a differential expression
#' table. Unlike a classical volcano plot, the x axis represents groups and
#' the y axis represents an effect such as average log2 fold-change.
#'
#' @param data A data frame containing marker results.
#' @param gene,group,effect,p_value,p_adjust Column names. Common Seurat column
#'   names are detected when these are `NULL`.
#' @param effect_cutoff Non-negative absolute effect threshold.
#' @param significance_cutoff P-value threshold in `[0, 1]`.
#' @param significance_by Use adjusted P values, raw P values, or no P-value
#'   criterion.
#' @param show Show all valid markers or significant markers only.
#' @param color_by Color points by status, direction, significance, group, or
#'   continuous effect.
#' @param label Label top markers, custom markers, or no markers.
#' @param label_n Number of labels selected per group and direction.
#' @param label_by Rank automatic labels by absolute effect or significance.
#' @param label_from Candidate population for automatic labels.
#' @param label_direction Directions eligible for automatic labels.
#' @param label_genes Character vector, named list, or data frame with `group`
#'   and `gene` columns.
#' @param missing_labels How missing custom markers are handled.
#' @param group_order Optional complete ordering of observed groups.
#' @param layout Vertical or horizontal layout.
#' @param point_layout Deterministic spread or seeded jitter.
#' @param spread_by Spread all points together or separately by direction.
#' @param spread_width,jitter_width Width of point displacement within a group.
#' @param seed Seed used only for jitter placement.
#' @param point_size Point size.
#' @param point_alpha Point alpha, or `NULL` for an automatic value.
#' @param palette Named point palette. For continuous effect coloring, a vector
#'   of two colors may be supplied.
#' @param group_palette Named group colors, an unnamed color vector, or `NULL`.
#' @param group_background,group_band Draw group regions and the zero-centered
#'   group band.
#' @param effect_limits Optional length-two visible effect range. Applied with
#'   Cartesian zoom so rows are not removed.
#' @param group_label_angle Rotation angle in degrees for cell-group labels.
#' @param group_label_size Text size for cell-group labels.
#' @param base_size Base font size.
#' @param ... Additional arguments passed to [ggrepel::geom_text_repel()].
#'
#' @return A `ggplot` object.
#' @export
marker_effect_plot <- function(
    data,
    gene = NULL,
    group = NULL,
    effect = NULL,
    p_value = NULL,
    p_adjust = NULL,
    effect_cutoff = 0.25,
    significance_cutoff = 0.05,
    significance_by = c("p_adjust", "p_value", "none"),
    show = c("all", "significant"),
    color_by = c("status", "direction", "significance", "group", "effect"),
    label = c("top", "custom", "none"),
    label_n = 5,
    label_by = c("effect", "significance"),
    label_from = c("significant", "shown", "all"),
    label_direction = c("both", "up", "down"),
    label_genes = NULL,
    missing_labels = c("warn", "ignore", "error"),
    group_order = NULL,
    layout = c("vertical", "horizontal"),
    point_layout = c("spread", "jitter"),
    spread_by = c("all", "direction"),
    spread_width = 0.8,
    jitter_width = 0.18,
    seed = 123,
    point_size = 1,
    point_alpha = NULL,
    palette = NULL,
    group_palette = NULL,
    group_background = TRUE,
    group_band = TRUE,
    effect_limits = NULL,
    group_label_angle = 0,
    group_label_size = 3,
    base_size = 12,
    ...) {
  significance_by <- rlang::arg_match(significance_by)
  show <- rlang::arg_match(show)
  color_by <- rlang::arg_match(color_by)
  label <- rlang::arg_match(label)
  label_by <- rlang::arg_match(label_by)
  label_from <- rlang::arg_match(label_from)
  label_direction <- rlang::arg_match(label_direction)
  missing_labels <- rlang::arg_match(missing_labels)
  layout <- rlang::arg_match(layout)
  point_layout <- rlang::arg_match(point_layout)
  spread_by <- rlang::arg_match(spread_by)
  validate_plot_options(effect_cutoff, significance_cutoff, label_n, spread_width,
    jitter_width, seed, point_size, point_alpha, effect_limits, group_label_angle,
    group_label_size)

  marker_data <- normalize_marker_data(data, gene, group, effect, p_value, p_adjust, significance_by)
  marker_data <- classify_markers(marker_data, effect_cutoff, significance_cutoff, significance_by)
  group_order <- validate_group_order(marker_data, group_order)
  marker_data <- compute_point_positions(marker_data, group_order, point_layout, spread_by, spread_width, jitter_width, seed)
  shown_data <- if (show == "significant") marker_data[marker_data$.smv_significant, , drop = FALSE] else marker_data
  if (nrow(shown_data) == 0L) .smv_abort("No markers satisfy the requested display criteria.", "scMarkerViz_empty_plot")

  label_data <- select_marker_labels(marker_data, shown_data, label, label_n, label_by,
    label_from, label_direction, label_genes, missing_labels, significance_by)
  group_colors <- resolve_group_palette(group_order, group_palette)
  group_data <- build_group_layout(shown_data, group_order, group_colors, effect_cutoff)
  group_effect_min <- stats::aggregate(
    shown_data$.smv_effect,
    list(.smv_group = shown_data$.smv_group),
    min
  )
  group_effect_max <- stats::aggregate(
    shown_data$.smv_effect,
    list(.smv_group = shown_data$.smv_group),
    max
  )
  min_index <- match(group_data$.smv_group, group_effect_min$.smv_group)
  max_index <- match(group_data$.smv_group, group_effect_max$.smv_group)
  group_data$.smv_xmin_background <- group_data$.smv_group_position - 0.425
  group_data$.smv_xmax_background <- group_data$.smv_group_position + 0.425
  group_data$.smv_ymin_background <- group_effect_min$x[min_index] - 0.3
  group_data$.smv_ymax_background <- group_effect_max$x[max_index] + 0.3
  point_alpha <- point_alpha %||% automatic_point_alpha(shown_data)

  plot <- ggplot2::ggplot(shown_data, ggplot2::aes(x = .data$.smv_x, y = .data$.smv_effect))
  if (group_background) {
    plot <- plot + ggplot2::geom_rect(
      data = group_data,
      ggplot2::aes(
        xmin = .data$.smv_xmin_background,
        xmax = .data$.smv_xmax_background,
        ymin = .data$.smv_ymin_background,
        ymax = .data$.smv_ymax_background
      ),
      inherit.aes = FALSE, fill = "#F2F2F2", color = NA, alpha = 1
    )
  }
  plot <- add_point_layer(plot, color_by, palette, group_colors, point_size, point_alpha)
  if (group_band) {
    plot <- plot +
      ggplot2::geom_rect(
        data = group_data,
        ggplot2::aes(xmin = .data$.smv_xmin, xmax = .data$.smv_xmax, ymin = .data$.smv_band_min, ymax = .data$.smv_band_max, fill = .data$.smv_group),
        inherit.aes = FALSE, color = "grey20", linewidth = 0.3, alpha = 1
      ) +
      ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
      ggplot2::geom_text(
        data = group_data,
        ggplot2::aes(x = .data$.smv_group_position, y = 0, label = .data$.smv_group),
        inherit.aes = FALSE, fontface = "bold", size = group_label_size,
        angle = group_label_angle
      )
  }
  if (nrow(label_data)) {
    plot <- plot + ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(x = .data$.smv_x, y = .data$.smv_effect, label = .data$.smv_gene),
      inherit.aes = FALSE, size = 3, max.overlaps = Inf, seed = seed, ...
    )
  }

  plot <- plot +
    ggplot2::scale_x_continuous(
      breaks = seq_along(group_order), labels = if (group_band) NULL else group_order,
      limits = c(0.5, length(group_order) + 0.5), expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Cell group", y = "Effect", color = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = if (group_band) ggplot2::element_blank() else ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  if (layout == "horizontal") {
    plot <- plot +
      ggplot2::coord_flip(ylim = effect_limits, clip = "off") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(),
        axis.ticks.x = ggplot2::element_line(color = "grey20")
      )
  } else if (!is.null(effect_limits)) {
    plot <- plot + ggplot2::coord_cartesian(ylim = effect_limits, clip = "off")
  }
  plot
}

validate_plot_options <- function(effect_cutoff, significance_cutoff, label_n,
                                  spread_width, jitter_width, seed, point_size,
                                  point_alpha, effect_limits, group_label_angle,
                                  group_label_size) {
  scalar_number <- function(x, nonnegative = TRUE) is.numeric(x) && length(x) == 1L && is.finite(x) && (!nonnegative || x >= 0)
  if (!scalar_number(effect_cutoff)) .smv_abort("`effect_cutoff` must be one finite non-negative number.", "scMarkerViz_invalid_option")
  if (!scalar_number(significance_cutoff) || significance_cutoff > 1) .smv_abort("`significance_cutoff` must be between 0 and 1.", "scMarkerViz_invalid_option")
  if (!scalar_number(label_n) || label_n != as.integer(label_n)) .smv_abort("`label_n` must be a non-negative integer.", "scMarkerViz_invalid_option")
  if (!scalar_number(spread_width) || spread_width > 1) .smv_abort("`spread_width` must be between 0 and 1.", "scMarkerViz_invalid_option")
  if (!scalar_number(jitter_width) || jitter_width > 0.5) .smv_abort("`jitter_width` must be between 0 and 0.5.", "scMarkerViz_invalid_option")
  if (!scalar_number(seed, FALSE)) .smv_abort("`seed` must be one finite number.", "scMarkerViz_invalid_option")
  if (!scalar_number(point_size)) .smv_abort("`point_size` must be one finite non-negative number.", "scMarkerViz_invalid_option")
  if (!is.null(point_alpha) && (!scalar_number(point_alpha) || point_alpha > 1)) .smv_abort("`point_alpha` must be `NULL` or between 0 and 1.", "scMarkerViz_invalid_option")
  if (!is.null(effect_limits) && (!is.numeric(effect_limits) || length(effect_limits) != 2L || any(!is.finite(effect_limits)) || effect_limits[1] >= effect_limits[2])) {
    .smv_abort("`effect_limits` must be two increasing finite numbers.", "scMarkerViz_invalid_option")
  }
  if (!scalar_number(group_label_angle, FALSE)) {
    .smv_abort("`group_label_angle` must be one finite number.", "scMarkerViz_invalid_option")
  }
  if (!scalar_number(group_label_size)) {
    .smv_abort("`group_label_size` must be one finite non-negative number.", "scMarkerViz_invalid_option")
  }
}

resolve_group_palette <- function(groups, palette) {
  n <- length(groups)
  if (is.null(palette)) palette <- grDevices::hcl.colors(n, palette = "Dynamic")
  if (!is.character(palette) || length(palette) < n) .smv_abort(sprintf("`group_palette` must provide at least %d colors.", n), "scMarkerViz_invalid_palette")
  if (!is.null(names(palette))) {
    missing <- setdiff(groups, names(palette))
    if (length(missing)) .smv_abort(sprintf("`group_palette` is missing colors for: %s.", paste(missing, collapse = ", ")), "scMarkerViz_invalid_palette")
    return(palette[groups])
  }
  stats::setNames(palette[seq_len(n)], groups)
}

automatic_point_alpha <- function(data) {
  max_group <- max(table(data$.smv_group))
  if (max_group <= 100) 0.8 else if (max_group <= 500) 0.5 else 0.3
}

add_point_layer <- function(plot, color_by, palette, group_colors, point_size, point_alpha) {
  default_status <- c(down = "#47C68D", not_significant = "#C9CDD3", up = "#7F3BD8")
  if (color_by == "effect") {
    colors <- palette %||% c("#47C68D", "#7F3BD8")
    if (length(colors) < 2L) .smv_abort("Continuous `palette` must contain at least two colors.", "scMarkerViz_invalid_palette")
    return(plot +
      ggplot2::geom_point(ggplot2::aes(color = .data$.smv_effect), size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_gradient2(low = colors[1], mid = "grey85", high = colors[length(colors)], midpoint = 0))
  }
  if (color_by == "group") {
    return(plot + ggplot2::geom_point(ggplot2::aes(color = .data$.smv_group), size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_manual(values = group_colors))
  }
  if (color_by == "direction") {
    colors <- palette %||% c(down = "#47C68D", zero = "#C9CDD3", up = "#7F3BD8")
    mapping <- ggplot2::aes(color = .data$.smv_direction)
  } else if (color_by == "significance") {
    colors <- palette %||% c(`FALSE` = "#C9CDD3", `TRUE` = "#7F3BD8")
    mapping <- ggplot2::aes(color = .data$.smv_significant)
  } else {
    colors <- palette %||% default_status
    mapping <- ggplot2::aes(color = .data$.smv_status)
  }
  plot + ggplot2::geom_point(mapping, size = point_size, alpha = point_alpha) + ggplot2::scale_color_manual(values = colors, drop = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
