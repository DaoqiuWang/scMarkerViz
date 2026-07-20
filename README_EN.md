# scMarkerViz

[中文说明](README.md)

`scMarkerViz` creates reproducible cluster-wise marker effect plots from single-cell differential expression results. The group or cell type is displayed along one axis, while an effect such as average log2 fold change is displayed along the other. Upregulated and downregulated genes appear on opposite sides of zero, and representative genes can be labeled automatically.

The package accepts a regular data frame, does not require a Seurat object, and returns a standard `ggplot` object that can be customized with additional `ggplot2` layers and themes.

## Installation

### Install from a local source directory

```r
install.packages("remotes")
remotes::install_local("scMarkerViz")
```

### Install from a source archive

```r
install.packages(
  "scMarkerViz_0.1.1.tar.gz",
  repos = NULL,
  type = "source"
)
```

### Install from GitHub

Replace `YOUR_USERNAME` with the GitHub account that hosts the repository:

```r
install.packages("remotes")
remotes::install_github("YOUR_USERNAME/scMarkerViz")
```

Load the package and verify its version:

```r
library(scMarkerViz)
packageVersion("scMarkerViz")
```

## Requirements

- R >= 4.1.0
- `dplyr`
- `ggplot2`
- `ggrepel`
- `rlang`

The dependencies are normally installed automatically. To install them manually:

```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "rlang"))
```

Optional development dependencies:

- `testthat >= 3.0.0` for tests
- `pkgload` for loading the package from source
- `remotes` for local or GitHub installation

## Input Data

The input must be a non-empty `data.frame`. Each row normally represents one gene-level differential expression result for one cell type or cluster.

### Required information

| Information | Requirement | Automatically detected column names |
|---|---|---|
| Gene | Non-empty gene name | `gene`, `feature`, `symbol`, `names` |
| Group | Non-empty cell type, cluster, or group name | `cluster`, `celltype`, `cell_type`, `group` |
| Effect | Finite numeric value, typically log2 fold change | `avg_log2FC`, `avg_logFC`, `log2FoldChange`, `logFC`, `logfoldchanges` |
| Raw P value | Numeric values between 0 and 1 | `p_val`, `pvalue`, `p_value`, `pvals` |
| Adjusted P value | Numeric values between 0 and 1 | `p_val_adj`, `padj`, `p_adjust`, `FDR`, `pvals_adj` |

An adjusted P-value column is required when `significance_by = "p_adjust"`. A raw P-value column is required when `significance_by = "p_value"`. P-value columns are optional when `significance_by = "none"`.

If multiple candidate columns are present for the same role, specify the desired column explicitly.

### Example data

```r
markers <- data.frame(
  gene = c("CD3D", "IL7R", "MS4A1", "CD79A"),
  celltype = c("T cell", "T cell", "B cell", "B cell"),
  avg_log2FC = c(1.2, 0.8, 1.5, 1.1),
  p_val = c(1e-10, 1e-6, 1e-12, 1e-8),
  p_val_adj = c(1e-8, 1e-4, 1e-10, 1e-6)
)
```

### Explicit column mapping

```r
p <- marker_effect_plot(
  data = results,
  gene = "symbol",
  group = "cell_type",
  effect = "logFC",
  p_value = "PValue",
  p_adjust = "FDR"
)
```

### Data validation

- Rows with missing or empty gene/group names are removed with a warning.
- Rows with `NA`, `Inf`, or `-Inf` effect values are removed with a warning.
- P-value columns must be numeric, and finite values must be between 0 and 1.
- A missing P value does not pass the significance criterion, but the row can still be shown with `show = "all"`.
- The function stops if no valid rows remain.

## Quick Start

```r
library(scMarkerViz)

p <- marker_effect_plot(
  data = markers,
  effect_cutoff = 0.25,
  significance_cutoff = 0.05,
  label_n = 5,
  group_label_angle = 30
)

p
```

Save the result:

```r
ggplot2::ggsave(
  "marker_effect.png",
  p,
  width = 12,
  height = 10,
  dpi = 300
)

ggplot2::ggsave(
  "marker_effect.pdf",
  p,
  width = 12,
  height = 10
)
```

## Function Reference

```r
marker_effect_plot(
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
  base_size = 12,
  ...
)
```

## Parameters

### Data and column mapping

- `data`: Differential expression results as a data frame.
- `gene`: Gene column name. If `NULL`, a common column name is detected.
- `group`: Cell-type or cluster column name. If `NULL`, a common column name is detected.
- `effect`: Numeric effect column, such as `avg_log2FC`.
- `p_value`: Raw P-value column.
- `p_adjust`: Adjusted P-value column.

### Significance and display

- `effect_cutoff`: Non-negative absolute effect threshold. Default: `0.25`.
- `significance_cutoff`: P-value threshold between 0 and 1. Default: `0.05`.
- `significance_by`: Use adjusted P values (`"p_adjust"`), raw P values (`"p_value"`), or no P-value criterion (`"none"`).
- `show`: Show all valid markers (`"all"`) or significant markers only (`"significant"`).

A marker is significant when its absolute effect reaches `effect_cutoff` and its selected P value is no greater than `significance_cutoff`. With `significance_by = "none"`, only the effect threshold is used.

### Point colors

- `color_by = "status"`: Colors significant upregulated, significant downregulated, and non-significant markers separately.
- `color_by = "direction"`: Colors markers by positive, negative, or zero effect.
- `color_by = "significance"`: Colors markers by significance status.
- `color_by = "group"`: Colors markers by cell type or cluster.
- `color_by = "effect"`: Uses a continuous effect gradient.

Custom status colors:

```r
palette <- c(
  down = "#47C68D",
  not_significant = "#C9CDD3",
  up = "#7F3BD8"
)

p <- marker_effect_plot(
  markers,
  color_by = "status",
  palette = palette
)
```

Continuous effect colors:

```r
p <- marker_effect_plot(
  markers,
  color_by = "effect",
  palette = c("#47C68D", "#7F3BD8")
)
```

Group colors may be supplied as an unnamed vector or a vector named by group:

```r
group_colors <- c(
  "T cell" = "#E64B35",
  "B cell" = "#4DBBD5",
  "Myeloid" = "#00A087"
)

p <- marker_effect_plot(
  markers,
  group_palette = group_colors,
  group_order = names(group_colors)
)
```

A named `group_palette` must contain a color for every observed group.

### Gene labels

- `label = "top"`: Automatically label top markers.
- `label = "custom"`: Label genes supplied through `label_genes`.
- `label = "none"`: Disable gene labels.
- `label_n`: Number of automatic labels per group and direction. Default: `5`.
- `label_by = "effect"`: Rank candidates by absolute effect.
- `label_by = "significance"`: Rank candidates by P value. This requires a P-value significance mode.
- `label_from`: Select candidates from significant, shown, or all markers.
- `label_direction`: Label both, upregulated, or downregulated markers.
- `missing_labels`: Warn, ignore, or stop when requested custom genes are absent.

Custom labels can be supplied as a character vector:

```r
p <- marker_effect_plot(
  markers,
  label = "custom",
  label_genes = c("CD3D", "MS4A1")
)
```

As a named list:

```r
p <- marker_effect_plot(
  markers,
  label = "custom",
  label_genes = list(
    "T cell" = c("CD3D", "IL7R"),
    "B cell" = c("MS4A1", "CD79A")
  )
)
```

Or as a data frame with `group` and `gene` columns:

```r
requested_labels <- data.frame(
  group = c("T cell", "B cell"),
  gene = c("CD3D", "MS4A1")
)

p <- marker_effect_plot(
  markers,
  label = "custom",
  label_genes = requested_labels
)
```

Additional arguments in `...` are passed to `ggrepel::geom_text_repel()`:

```r
p <- marker_effect_plot(
  markers,
  box.padding = 0.6,
  point.padding = 0.4,
  segment.color = "grey50",
  force = 2,
  min.segment.length = 0
)
```

### Group order and orientation

Set an explicit complete group order:

```r
p <- marker_effect_plot(
  markers,
  group_order = c("Astrocytes", "Microglia", "Oligodendrocytes")
)
```

`group_order` must contain every observed group exactly once.

Use a horizontal layout:

```r
p <- marker_effect_plot(
  markers,
  layout = "horizontal"
)
```

Rotate the text inside the group label bands:

```r
p <- marker_effect_plot(
  markers,
  group_label_angle = 30
)
```

`group_label_angle` is measured in degrees and accepts positive or negative finite values. It is supported in both vertical and horizontal layouts.

### Point placement

- `point_layout = "spread"`: Deterministic placement sorted by gene and source row. Default.
- `point_layout = "jitter"`: Seeded random displacement.
- `spread_by = "all"`: Spread all markers in a group together.
- `spread_by = "direction"`: Spread upregulated, downregulated, and zero-effect markers separately.
- `spread_width`: Total deterministic spread width from 0 to 1. Default: `0.8`.
- `jitter_width`: One-sided random displacement from 0 to 0.5. Default: `0.18`.
- `seed`: Seed used for jitter and label placement. Default: `123`.

Example:

```r
p <- marker_effect_plot(
  markers,
  point_layout = "jitter",
  jitter_width = 0.2,
  seed = 2026
)
```

### Point appearance

- `point_size`: Non-negative point size. Default: `1`.
- `point_alpha`: Point opacity from 0 to 1.

When `point_alpha = NULL`, opacity is selected automatically from the largest group size:

- Up to 100 markers: `0.8`
- 101–500 markers: `0.5`
- More than 500 markers: `0.3`

### Group backgrounds and label bands

`group_background = TRUE` draws an independent light-gray background rectangle for each group. Each rectangle:

- has a fixed width of `0.85`;
- is centered on the group position;
- begins at the group's minimum shown effect minus `0.3`;
- ends at the group's maximum shown effect plus `0.3`;
- is rendered below points, group bands, and text.

Disable these rectangles with:

```r
group_background = FALSE
```

`group_band = TRUE` draws an opaque colored band around zero and prints the group name inside it. The band spans from `-effect_cutoff` to `effect_cutoff`. When disabled, group names are shown as regular axis labels.

```r
p <- marker_effect_plot(
  markers,
  group_band = FALSE
)
```

### Effect limits and text size

Set the visible effect range without dropping rows:

```r
p <- marker_effect_plot(
  markers,
  effect_limits = c(-4, 4)
)
```

`effect_limits` must contain two increasing finite numbers. In the horizontal layout it still controls the effect-value range.

`base_size` controls the base font size used by `theme_minimal()`. Default: `12`.

## Seurat Example

`marker_effect_plot()` does not require Seurat, but it can use the output of `FindAllMarkers()`:

```r
markers <- Seurat::FindAllMarkers(seurat_object)

# Add a gene column if genes are stored in row names.
if (!"gene" %in% names(markers)) {
  markers$gene <- rownames(markers)
}

p <- marker_effect_plot(
  data = markers,
  gene = "gene",
  group = "cluster",
  effect = "avg_log2FC",
  p_value = "p_val",
  p_adjust = "p_val_adj",
  group_label_angle = 30
)
```

## Complete Example

```r
library(scMarkerViz)

markers <- read.csv(
  "HFE_vs_LFE_overall_all_DEGs.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

p <- marker_effect_plot(
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
  label_direction = "both",
  point_layout = "spread",
  spread_by = "all",
  spread_width = 0.8,
  point_size = 0.8,
  group_background = TRUE,
  group_band = TRUE,
  group_label_angle = 30,
  base_size = 12
) +
  ggplot2::labs(
    title = "Differential expression across cell types",
    subtitle = "HFE vs LFE — top 5 significant genes per direction",
    y = "Average log2 fold change"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0.5)
  )

ggplot2::ggsave(
  "HFE_vs_LFE_marker_effect_vertical.png",
  p,
  width = 12,
  height = 10,
  dpi = 300
)
```

## Customizing the Returned Plot

The result is a standard `ggplot` object:

```r
p <- marker_effect_plot(markers, group_label_angle = 30)

p +
  ggplot2::labs(
    title = "Differential expression across cell types",
    subtitle = "Top markers per direction",
    y = "Average log2 fold change"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
```

## Troubleshooting

### A required column cannot be detected

Specify all column mappings explicitly:

```r
marker_effect_plot(
  results,
  gene = "symbol",
  group = "annotation",
  effect = "logFC",
  p_adjust = "FDR"
)
```

### Only raw P values are available

```r
marker_effect_plot(
  results,
  significance_by = "p_value",
  p_value = "pvalue"
)
```

### No P values are available

```r
marker_effect_plot(
  results,
  significance_by = "none",
  label_by = "effect"
)
```

### No markers remain with `show = "significant"`

Reduce `effect_cutoff`, increase `significance_cutoff`, or use `show = "all"` to inspect the data first.

### Gene labels overlap

Reduce `label_n` or adjust `ggrepel` options:

```r
marker_effect_plot(
  markers,
  label_n = 3,
  box.padding = 0.8,
  force = 2
)
```

### Group names are long

Rotate group labels and increase the output dimensions:

```r
marker_effect_plot(
  markers,
  group_label_angle = 30
)
```

## Help

After installation, open the R help page with:

```r
?marker_effect_plot
help("marker_effect_plot", package = "scMarkerViz")
```
