example_markers <- data.frame(
  gene = paste0("gene", seq_len(18)),
  cluster = rep(c("T", "B", "Mono"), each = 6),
  avg_log2FC = rep(c(-1.2, -0.7, -0.1, 0.2, 0.8, 1.4), 3),
  p_val = rep(c(0.001, 0.01, 0.4, 0.2, 0.02, 0.0001), 3),
  p_val_adj = rep(c(0.005, 0.04, 0.8, 0.5, 0.03, 0.001), 3),
  stringsAsFactors = FALSE
)
