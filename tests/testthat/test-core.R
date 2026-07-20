test_that("common Seurat columns are normalized without modifying input", {
  original <- example_markers
  normalized <- scMarkerViz:::normalize_marker_data(
    example_markers, NULL, NULL, NULL, NULL, NULL, "p_adjust"
  )
  expect_equal(example_markers, original)
  expect_true(all(c(".smv_gene", ".smv_group", ".smv_effect") %in% names(normalized)))
})

test_that("classification respects effect and adjusted P-value boundaries", {
  normalized <- scMarkerViz:::normalize_marker_data(
    example_markers, NULL, NULL, NULL, NULL, NULL, "p_adjust"
  )
  classified <- scMarkerViz:::classify_markers(normalized, 0.25, 0.05, "p_adjust")
  expect_equal(sum(classified$.smv_significant), 12)
  expect_true(all(as.character(classified$.smv_status[classified$.smv_effect > 0.25 & classified$.smv_p_adjust <= 0.05]) == "up"))
})

test_that("group order must be complete and unique", {
  normalized <- scMarkerViz:::normalize_marker_data(
    example_markers, NULL, NULL, NULL, NULL, NULL, "p_adjust"
  )
  expect_equal(scMarkerViz:::validate_group_order(normalized, c("Mono", "T", "B")), c("Mono", "T", "B"))
  expect_error(scMarkerViz:::validate_group_order(normalized, c("T", "B")), class = "scMarkerViz_invalid_group_order")
  expect_error(scMarkerViz:::validate_group_order(normalized, c("T", "B", "B", "Mono")), class = "scMarkerViz_invalid_group_order")
})

test_that("spread stays inside each group and is deterministic", {
  normalized <- scMarkerViz:::normalize_marker_data(
    example_markers, NULL, NULL, NULL, NULL, NULL, "p_adjust"
  )
  classified <- scMarkerViz:::classify_markers(normalized, 0.25, 0.05, "p_adjust")
  positioned <- scMarkerViz:::compute_point_positions(classified, c("T", "B", "Mono"), "spread", "all", 0.8, 0.18, 123)
  expect_true(all(abs(positioned$.smv_x - positioned$.smv_group_position) <= 0.4 + 1e-12))
  positioned2 <- scMarkerViz:::compute_point_positions(classified, c("T", "B", "Mono"), "spread", "all", 0.8, 0.18, 999)
  expect_equal(positioned$.smv_x, positioned2$.smv_x)
})

test_that("automatic labels select each group and direction", {
  normalized <- scMarkerViz:::normalize_marker_data(
    example_markers, NULL, NULL, NULL, NULL, NULL, "p_adjust"
  )
  classified <- scMarkerViz:::classify_markers(normalized, 0.25, 0.05, "p_adjust")
  positioned <- scMarkerViz:::compute_point_positions(classified, c("T", "B", "Mono"), "spread", "all", 0.8, 0.18, 123)
  labels <- scMarkerViz:::select_marker_labels(positioned, positioned, "top", 1, "effect", "significant", "both", NULL, "warn", "p_adjust")
  expect_equal(nrow(labels), 6)
  expect_equal(nrow(unique(labels[c(".smv_group", ".smv_gene")])), 6)
})
