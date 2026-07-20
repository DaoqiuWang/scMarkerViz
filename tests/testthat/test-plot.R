test_that("vertical and horizontal plots return ggplot objects", {
  vertical <- marker_effect_plot(example_markers, label_n = 1)
  horizontal <- marker_effect_plot(example_markers, label_n = 1, layout = "horizontal")
  expect_s3_class(vertical, "ggplot")
  expect_s3_class(horizontal, "ggplot")
  expect_silent(ggplot2::ggplot_build(vertical))
  expect_silent(ggplot2::ggplot_build(horizontal))
})

test_that("cell-group labels support rotation", {
  plot <- marker_effect_plot(example_markers, label_n = 1, group_label_angle = 30)
  built <- ggplot2::ggplot_build(plot)
  expect_true(any(vapply(built$data, function(layer) {
    "angle" %in% names(layer) && any(layer$angle == 30)
  }, logical(1))))
  expect_error(
    marker_effect_plot(example_markers, group_label_angle = NA_real_),
    class = "scMarkerViz_invalid_option"
  )
})

test_that("custom labels support vectors and named lists", {
  vector_plot <- marker_effect_plot(
    example_markers, label = "custom", label_genes = c("gene1", "gene8")
  )
  list_plot <- marker_effect_plot(
    example_markers, label = "custom", label_genes = list(T = "gene1", B = "gene8")
  )
  expect_s3_class(vector_plot, "ggplot")
  expect_s3_class(list_plot, "ggplot")
})

test_that("significant-only display errors clearly when empty", {
  expect_error(
    marker_effect_plot(example_markers, effect_cutoff = 100, show = "significant"),
    class = "scMarkerViz_empty_plot"
  )
})

test_that("more than nine groups receive sufficient colors", {
  many <- example_markers[rep(seq_len(nrow(example_markers)), 4), ]
  many$cluster <- rep(paste0("group", seq_len(12)), each = 6)
  plot <- marker_effect_plot(many, label = "none")
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})
