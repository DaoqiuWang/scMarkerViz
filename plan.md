# scMarkerViz 开发计划

## 1. 项目目标

`scMarkerViz` 是一个面向单细胞差异表达结果的轻量 R 绘图包。首个版本聚焦 cluster-wise marker effect plot：按细胞群展示基因效应量分布，并提供清晰、可复现的显著性编码和 marker 标注。

本包不把这种图误称为经典火山图。经典火山图应以 log fold-change 为横轴、`-log10(P)` 为纵轴；本包首版实现的是 marker effect plot。

## 2. 首版范围

### 包含

- 普通 `data.frame` 输入，不绑定 Seurat 对象。
- 可配置 gene、group、effect、原始 P 值及校正 P 值列。
- 自动识别常见 Seurat 列名。
- 基于 effect cutoff 与显著性 cutoff 的稳定分类。
- 显示全部基因或只显示显著基因。
- 按综合状态、方向、显著性、group 或连续 effect 着色。
- 每组、每个方向自动选择 top markers。
- 支持全局 gene 向量、按 group 命名的 gene list，以及 group/gene 两列表格。
- 确定性 spread 和可复现 jitter 两种点布局。
- 纵向和横向布局。
- 每组只绘制一次的背景、零点色带和 group 标签。
- 返回标准 `ggplot` 对象。
- 单元测试与基础包检查。

### 暂不包含

- 极坐标布局。
- 在包内运行差异表达分析。
- 强制依赖 Seurat。
- 自动保存文件。
- Shiny 界面。
- 经典 volcano plot（计划在后续版本作为独立函数实现）。

## 3. API 设计

主函数：

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
  base_size = 12,
  ...
)
```

## 4. 内部数据协议

输入标准化后使用以下内部列：

- `.smv_gene`
- `.smv_group`
- `.smv_effect`
- `.smv_p_value`
- `.smv_p_adjust`
- `.smv_row_id`
- `.smv_direction`
- `.smv_significant`
- `.smv_status`
- `.smv_x`

内部列使用包前缀，避免覆盖用户原始列。函数不修改用户输入对象。

## 5. 统计规则

- `direction` 只描述 effect 的符号：`up`、`down`、`zero`。
- `significant` 同时要求 `abs(effect) >= effect_cutoff`，并在启用 P 值判定时要求所选 P 值 `<= significance_cutoff`。
- `status` 为 `up`、`down` 或 `not_significant`。
- 缺失的显著性值默认视为不显著，不静默删除整行。
- 自动标签默认只从显著基因中选择。
- top 标签按 group 和方向分别选择，明确不保留 ties，并使用 gene 与原始行号保证稳定排序。

## 6. 架构

1. `normalize_marker_data()`：解析列映射并建立内部标准表。
2. `validate_marker_options()`：验证 cutoff、枚举、顺序、颜色等参数。
3. `classify_markers()`：生成方向、显著性及综合状态。
4. `compute_point_positions()`：生成确定性 spread 或带 seed 的 jitter 位置。
5. `select_marker_labels()`：选择自动或自定义标签。
6. `build_group_layout()`：生成每组一行的背景、色带及标签数据。
7. `marker_effect_plot()`：组合图层并应用纵向或横向布局。

## 7. 依赖策略

核心 Imports：

- `ggplot2`
- `dplyr`
- `rlang`
- `ggrepel`
- `scales`

不使用 `Depends: tidyverse`，不将 Seurat 放入 Imports。

测试 Suggests：

- `testthat`

## 8. 错误处理

- 缺少列时指出参数名、目标列和可用列。
- 非数值 effect/P 值立即报错。
- P 值超出 `[0, 1]` 报错。
- 空数据和筛选后无可显示数据给出明确错误。
- `group_order` 重复、遗漏或包含未知组时给出明确错误。
- 自定义 marker 缺失时按 `missing_labels` 处理。
- 非法枚举通过 `rlang::arg_match()` 在入口处理。

## 9. 测试计划

- 列自动识别与显式映射。
- 输入不被修改。
- 显著性边界与 NA 行为。
- 每组每方向 top N、去重和稳定排序。
- 三种自定义 marker 输入。
- group 顺序验证。
- spread 位置始终位于本组区域。
- jitter 固定 seed 后可复现。
- group layout 每组恰好一行。
- 纵向和横向均返回 `ggplot`。
- 超过九个 group 时颜色充足。
- 非法参数与空结果错误信息。

## 10. 里程碑

### M1：包骨架与数据核心

完成元数据、标准化、验证、分类和标签选择。

### M2：首个绘图版本

完成纵向、横向、spread、jitter、背景区域、group band、图例和主题。

### M3：质量保障

补充示例、文档、单元测试，运行 `testthat` 与 `R CMD check`。

### M4：后续版本

在首版 API 稳定后，再独立设计极坐标布局和经典 volcano plot。

## 11. 首版验收标准

- 常见 Seurat `FindAllMarkers()` 表格可一行调用。
- 任意合法列名可显式映射。
- 默认展示所有点，并明确区分不显著、上调和下调。
- group 背景、色带和名称不会按基因重复绘制。
- 默认 spread 完全确定且可复现。
- 自动标签数量和排序可预测。
- 横向布局不依赖与纵向冲突的布尔参数。
- 返回值可继续叠加任意 ggplot2 图层。
- 测试通过，包检查不存在 error 或 warning。
