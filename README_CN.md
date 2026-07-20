# scMarkerViz 使用说明

[English documentation](README.md)

`scMarkerViz` 用于将单细胞差异表达结果绘制为按细胞群/cluster 排列的 marker effect 图。横轴表示细胞类型或 cluster，纵轴表示 `avg_log2FC` 等效应值；上调和下调基因分布在零点两侧，并可自动标注每组的代表性基因。

该包直接接收差异表达结果数据框，不依赖 Seurat 对象，返回标准 `ggplot` 对象，可以继续叠加 `ggplot2` 图层或主题。

## 1. 安装

### 1.1 从源码目录安装

```r
install.packages("remotes")
remotes::install_local("scMarkerViz")
```

也可以使用 R 自带的安装命令：

```r
install.packages("scMarkerViz_0.1.1.tar.gz", repos = NULL, type = "source")
```

安装后加载并检查版本：

```r
library(scMarkerViz)
packageVersion("scMarkerViz")
```

### 1.2 开发环境中加载

在包源码目录中开发时，可以使用：

```r
install.packages("pkgload")
pkgload::load_all("scMarkerViz")
```

## 2. 依赖

### 2.1 R 版本

- R >= 4.1.0

### 2.2 必需依赖

安装本包时通常会自动安装：

- `ggplot2`：绘图
- `dplyr`：数据整理
- `ggrepel`：基因标签避让
- `rlang`：参数匹配、错误和警告处理

手动安装方式：

```r
install.packages(c("ggplot2", "dplyr", "ggrepel", "rlang"))
```

### 2.3 可选依赖

- `testthat >= 3.0.0`：运行包测试时需要
- `remotes`：从本地源码目录安装时可用
- `pkgload`：开发阶段直接加载源码时可用

## 3. 输入数据要求

输入必须是至少包含一行的 `data.frame`。每一行通常代表某个基因在某个细胞类型或 cluster 中的一条差异表达结果。

### 3.1 必需信息

| 信息 | 内容要求 | 可自动识别的列名 |
|---|---|---|
| 基因 | 非空基因名称，可转换为字符型 | `gene`, `feature`, `symbol`, `names` |
| 分组 | 非空的细胞类型、cluster 或其他分组名称 | `cluster`, `celltype`, `cell_type`, `group` |
| 效应值 | 有限数值，通常为 log2 fold change | `avg_log2FC`, `avg_logFC`, `log2FoldChange`, `logFC`, `logfoldchanges` |

当 `significance_by = "p_adjust"` 时，还必须提供调整后 P 值；当 `significance_by = "p_value"` 时，必须提供原始 P 值。

| 显著性信息 | 内容要求 | 可自动识别的列名 |
|---|---|---|
| 原始 P 值 | 数值型，范围 0–1 | `p_val`, `pvalue`, `p_value`, `pvals` |
| 调整后 P 值 | 数值型，范围 0–1 | `p_val_adj`, `padj`, `p_adjust`, `FDR`, `pvals_adj` |

如果同一种信息存在多个候选列，包不会自行猜测，需要通过 `gene`、`group`、`effect`、`p_value` 或 `p_adjust` 明确指定。

### 3.2 推荐数据格式

```r
markers <- data.frame(
  gene = c("CD3D", "IL7R", "MS4A1", "CD79A"),
  celltype = c("T", "T", "B", "B"),
  avg_log2FC = c(1.2, 0.8, 1.5, 1.1),
  p_val = c(1e-10, 1e-6, 1e-12, 1e-8),
  p_val_adj = c(1e-8, 1e-4, 1e-10, 1e-6)
)
```

### 3.3 自定义列名

```r
results <- data.frame(
  symbol = c("CD3D", "MS4A1"),
  cell_type = c("T cell", "B cell"),
  logFC = c(1.2, -1.1),
  PValue = c(1e-8, 1e-7),
  FDR = c(1e-6, 1e-5)
)

p <- marker_effect_plot(
  data = results,
  gene = "symbol",
  group = "cell_type",
  effect = "logFC",
  p_value = "PValue",
  p_adjust = "FDR"
)
```

### 3.4 数据清理规则

- 基因名或分组名缺失、为空的行会被移除并给出警告。
- 效应值为 `NA`、`Inf` 或 `-Inf` 的行会被移除并给出警告。
- P 值必须为数值型；有限值必须位于 0–1。
- 当前显著性模式所需的 P 值为 `NA` 时，该行不会被判定为显著，但在 `show = "all"` 时仍可显示。
- 清理后没有有效数据时会停止运行。

## 4. 基本用法

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

保存图片：

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

## 5. 完整函数形式

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

对于具有多个可选值的参数，只需传入其中一个字符串，例如 `layout = "horizontal"`。

## 6. 所有参数说明

### 6.1 数据与列映射

#### `data`

必需。差异表达结果 `data.frame`。

#### `gene`

基因列名，单个字符串。默认 `NULL` 时自动识别。常用设置：

```r
gene = "gene"
```

#### `group`

细胞类型或 cluster 列名，单个字符串。默认 `NULL` 时自动识别。

```r
group = "celltype"
```

#### `effect`

效应值列名，单个字符串，例如 `avg_log2FC`。

```r
effect = "avg_log2FC"
```

#### `p_value`

原始 P 值列名。使用 `significance_by = "p_value"` 时必需。

#### `p_adjust`

调整后 P 值列名。使用默认的 `significance_by = "p_adjust"` 时必需。

### 6.2 显著性判定与数据显示

#### `effect_cutoff`

非负数，默认 `0.25`。基因需满足：

```r
abs(effect) >= effect_cutoff
```

才可能被判定为显著。该值同时决定细胞类型标签色带在效应轴上的范围，即 `-effect_cutoff` 到 `effect_cutoff`。

#### `significance_cutoff`

0–1 之间的数，默认 `0.05`。与 `significance_by` 指定的 P 值比较。

#### `significance_by`

显著性所依据的统计量：

- `"p_adjust"`：使用调整后 P 值，默认且通常推荐。
- `"p_value"`：使用原始 P 值。
- `"none"`：忽略 P 值，只根据绝对效应阈值判断。

最终显著条件为绝对效应达到阈值，并且所选 P 值不大于显著性阈值。

#### `show`

控制绘制哪些点：

- `"all"`：显示所有有效数据，默认。
- `"significant"`：仅显示被判定为显著的数据。

如果选择 `"significant"` 后没有任何数据，函数会报错。

### 6.3 点的颜色

#### `color_by`

控制点的着色方式：

- `"status"`：默认。显著上调、显著下调和不显著分别着色。
- `"direction"`：按效应方向 `up`、`down`、`zero` 着色，不考虑显著性。
- `"significance"`：按是否显著着色。
- `"group"`：按细胞类型或 cluster 着色。
- `"effect"`：按连续效应值使用渐变色。

#### `palette`

点颜色配置，格式取决于 `color_by`。

`color_by = "status"` 时建议使用具名颜色：

```r
palette = c(
  down = "#47C68D",
  not_significant = "#C9CDD3",
  up = "#7F3BD8"
)
```

`color_by = "direction"`：

```r
palette = c(
  down = "#47C68D",
  zero = "#C9CDD3",
  up = "#7F3BD8"
)
```

`color_by = "significance"`：

```r
palette = c(
  `FALSE` = "#C9CDD3",
  `TRUE` = "#7F3BD8"
)
```

`color_by = "effect"` 时至少提供两个颜色，第一个用于负效应，最后一个用于正效应，中点为灰色：

```r
palette = c("#47C68D", "#7F3BD8")
```

`color_by = "group"` 时点颜色来自 `group_palette`，而不是 `palette`。

#### `group_palette`

细胞类型标签框及 `color_by = "group"` 时的颜色。默认使用 `hcl.colors(..., palette = "Dynamic")`。

可以提供未命名颜色向量，其长度至少等于分组数：

```r
group_palette = c("#E64B35", "#4DBBD5", "#00A087")
```

也可以提供按分组命名的颜色；此时必须覆盖所有分组：

```r
group_palette = c(
  "T cell" = "#E64B35",
  "B cell" = "#4DBBD5",
  "Myeloid" = "#00A087"
)
```

### 6.4 基因标签

#### `label`

标签选择方式：

- `"top"`：自动选择每个组的代表基因，默认。
- `"custom"`：使用 `label_genes` 指定基因。
- `"none"`：不显示基因标签。

#### `label_n`

每个分组、每个方向选择的标签数，默认 `5`。例如默认 `label_direction = "both"` 时，每组最多标注 5 个上调和 5 个下调基因。

必须是非负整数；设为 `0` 可关闭自动标签。

#### `label_by`

自动标签排序依据：

- `"effect"`：按绝对效应值排序，默认。
- `"significance"`：优先选择 P 值更小的基因。不能与 `significance_by = "none"` 同时使用。

#### `label_from`

自动标签的候选数据范围：

- `"significant"`：只从显著基因中选择，默认。
- `"shown"`：从当前实际显示的数据中选择。
- `"all"`：从全部有效数据中选择，即使某些数据未显示。

#### `label_direction`

自动标签允许的方向：

- `"both"`：上调和下调都标注，默认。
- `"up"`：只标注正效应基因。
- `"down"`：只标注负效应基因。

零效应基因不会进入自动标签。

#### `label_genes`

在 `label = "custom"` 时指定要标注的基因，支持三种格式。

字符向量：在所有组中查找这些基因：

```r
label_genes = c("CD3D", "MS4A1")
```

按组命名的列表：

```r
label_genes = list(
  "T cell" = c("CD3D", "IL7R"),
  "B cell" = c("MS4A1", "CD79A")
)
```

包含 `group` 和 `gene` 两列的数据框：

```r
label_genes = data.frame(
  group = c("T cell", "B cell"),
  gene = c("CD3D", "MS4A1")
)
```

#### `missing_labels`

自定义标签未在数据中找到时的处理方式：

- `"warn"`：给出警告并继续，默认。
- `"ignore"`：忽略缺失标签。
- `"error"`：停止并报错。

#### `...`

传递给 `ggrepel::geom_text_repel()` 的其他参数，例如：

```r
marker_effect_plot(
  markers,
  box.padding = 0.6,
  point.padding = 0.4,
  segment.color = "grey50",
  segment.linewidth = 0.4,
  force = 2,
  min.segment.length = 0
)
```

函数内部已经设置 `size = 3`、`max.overlaps = Inf` 和 `seed = seed`，不要在 `...` 中重复传入这些参数。

### 6.5 分组顺序与布局

#### `group_order`

完整的分组显示顺序。必须包含数据中所有分组，不能重复，也不能包含未知分组。

```r
group_order = c("Astrocytes", "Microglia", "Oligodendrocytes")
```

默认 `NULL`，使用分组在数据中首次出现的顺序。

#### `layout`

图形方向：

- `"vertical"`：默认。细胞类型沿水平轴排列，效应值沿垂直轴排列。
- `"horizontal"`：翻转坐标，细胞类型沿垂直方向排列，效应值沿水平方向排列。

```r
marker_effect_plot(markers, layout = "horizontal")
```

#### `group_label_angle`

细胞类型标签的旋转角度，单位为度，默认 `0`。支持正数和负数。

```r
group_label_angle = 30
group_label_angle = 45
group_label_angle = 90
group_label_angle = -30
```

此参数旋转彩色标签框内部的分组文字；vertical 和 horizontal 布局均可使用。

### 6.6 点的位置

#### `point_layout`

组内点的水平排列方式：

- `"spread"`：默认。按基因名称和原始行号进行确定性排序，再均匀展开；相同数据每次得到相同位置。
- `"jitter"`：使用带种子的随机抖动。

#### `spread_by`

仅在 `point_layout = "spread"` 时使用：

- `"all"`：每个组的所有点一起展开，默认。
- `"direction"`：每个组内按上调、下调和零效应分别展开。

#### `spread_width`

确定性展开的总宽度，范围 0–1，默认 `0.8`。值越大，点在每个分组内分布得越宽。

```r
spread_width = 0.6
```

#### `jitter_width`

随机抖动的单侧宽度，范围 0–0.5，默认 `0.18`。点的位置在分组中心加减该值之间。

#### `seed`

`point_layout = "jitter"` 以及 `ggrepel` 标签布局所用随机种子，默认 `123`。函数会恢复调用前的全局随机数状态。

### 6.7 点的外观

#### `point_size`

散点大小，默认 `1`，必须为非负有限数。

#### `point_alpha`

散点透明度，范围 0–1。默认 `NULL` 时根据单个分组中的最大点数自动选择：

- 不超过 100 个点：`0.8`
- 101–500 个点：`0.5`
- 超过 500 个点：`0.3`

也可明确设置：

```r
point_alpha = 0.7
```

### 6.8 分组背景和标签框

#### `group_background`

逻辑值，默认 `TRUE`。控制是否绘制每个 cluster 独立的浅灰背景矩形。

当前背景矩形规则为：

- 宽度固定为 `0.85`；
- 水平中心为对应分组中心；
- 下边界为该分组显示数据的最小效应值减 `0.3`；
- 上边界为该分组显示数据的最大效应值加 `0.3`；
- 背景位于散点、标签框和文字下方。

关闭方式：

```r
group_background = FALSE
```

#### `group_band`

逻辑值，默认 `TRUE`。控制是否在零点附近绘制彩色细胞类型标签框，并在框内显示分组名称。

标签框垂直范围是：

```r
-effect_cutoff 到 effect_cutoff
```

标签框使用不透明的分组颜色。设置为 `FALSE` 时，不绘制彩色标签框，分组名会改为普通坐标轴标签。

### 6.9 坐标范围与字号

#### `effect_limits`

长度为 2 的递增有限数值向量，用于设置可见效应范围。例如：

```r
effect_limits = c(-4, 4)
```

使用 Cartesian zoom，不会因缩放而先删除数据行。horizontal 布局中同样表示效应值范围。

#### `base_size`

`theme_minimal()` 的基础字号，默认 `12`。

## 7. 常用示例

### 7.1 Seurat 差异表达结果

如果基因名保存在行名中，应先转换为普通列：

```r
markers <- Seurat::FindAllMarkers(seurat_object)
markers$gene <- rownames(markers)

p <- marker_effect_plot(
  markers,
  gene = "gene",
  group = "cluster",
  effect = "avg_log2FC",
  p_value = "p_val",
  p_adjust = "p_val_adj"
)
```

如果 `FindAllMarkers()` 已经返回 `gene` 列，则不需要从行名复制。

### 7.2 仅显示显著基因

```r
p <- marker_effect_plot(
  markers,
  effect_cutoff = 0.5,
  significance_cutoff = 0.01,
  significance_by = "p_adjust",
  show = "significant"
)
```

### 7.3 忽略 P 值，只按 fold change

```r
p <- marker_effect_plot(
  markers,
  significance_by = "none",
  effect_cutoff = 0.5,
  label_by = "effect"
)
```

### 7.4 自定义颜色

```r
p <- marker_effect_plot(
  markers,
  color_by = "status",
  palette = c(
    down = "#2CA25F",
    not_significant = "#D9D9D9",
    up = "#756BB1"
  )
)
```

### 7.5 按细胞类型着色

```r
cell_colors <- c(
  "Astrocytes" = "#E64B35",
  "Microglia" = "#4DBBD5",
  "Oligodendrocytes" = "#00A087"
)

p <- marker_effect_plot(
  markers,
  color_by = "group",
  group_palette = cell_colors,
  group_order = names(cell_colors)
)
```

### 7.6 连续效应渐变

```r
p <- marker_effect_plot(
  markers,
  color_by = "effect",
  palette = c("#47C68D", "#7F3BD8")
)
```

### 7.7 自定义基因标签

```r
p <- marker_effect_plot(
  markers,
  label = "custom",
  label_genes = list(
    "T cells" = c("CD3D", "IL7R"),
    "B cells" = c("MS4A1", "CD79A")
  ),
  missing_labels = "warn",
  group_label_angle = 30
)
```

### 7.8 横向布局

```r
p <- marker_effect_plot(
  markers,
  layout = "horizontal",
  group_label_angle = 30,
  effect_limits = c(-4, 4)
)
```

### 7.9 随机抖动布局

```r
p <- marker_effect_plot(
  markers,
  point_layout = "jitter",
  jitter_width = 0.2,
  seed = 2026
)
```

### 7.10 继续使用 ggplot2 定制

函数返回标准 `ggplot` 对象：

```r
p <- marker_effect_plot(markers, group_label_angle = 30)

p <- p +
  ggplot2::labs(
    title = "Differential expression across cell types",
    subtitle = "Top five significant genes per direction",
    y = "Average log2 fold change"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0.5),
    legend.position = "right"
  )
```

## 8. 完整 HFE 与 LFE 示例

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

## 9. 返回值

`marker_effect_plot()` 返回一个 `ggplot` 对象。可以直接打印、保存，或继续添加 `ggplot2::labs()`、`ggplot2::theme()` 等图层。

```r
p <- marker_effect_plot(markers)
class(p)
```

## 10. 常见问题

### 找不到基因、分组或效应列

显式指定列名：

```r
marker_effect_plot(
  results,
  gene = "symbol",
  group = "annotation",
  effect = "logFC",
  p_adjust = "FDR"
)
```

### 数据只有原始 P 值

```r
marker_effect_plot(
  results,
  significance_by = "p_value",
  p_value = "pvalue"
)
```

### 数据没有任何 P 值

```r
marker_effect_plot(
  results,
  significance_by = "none",
  label_by = "effect"
)
```

### 显著基因模式没有点

降低 `effect_cutoff` 或提高 `significance_cutoff`，也可以先使用 `show = "all"` 检查数据。

### 标签过多或重叠

减少 `label_n`，或通过 `...` 调整 `ggrepel`：

```r
marker_effect_plot(
  markers,
  label_n = 3,
  box.padding = 0.8,
  force = 2
)
```

### 分组名称较长

使用标签角度参数：

```r
marker_effect_plot(
  markers,
  group_label_angle = 30
)
```

也可以增加输出图的宽度或高度。

### 修改图例位置

```r
marker_effect_plot(markers) +
  ggplot2::theme(legend.position = "right")
```

## 11. 查看 R 帮助

安装并加载包后运行：

```r
?marker_effect_plot
help("marker_effect_plot", package = "scMarkerViz")
```
