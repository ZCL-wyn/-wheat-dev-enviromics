library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(readr)
library(grid)

setwd("/mnt/7t_storage/zhangcl/TKW/TKW_GS_Complete_NoPCA")

geno5_data  <- read_csv("Geno5_fold/基因型5折交叉验证_最终结果.csv")
envloo_data <- read_csv("EnvLOO/环境留一交叉验证_最终结果.csv")
combo_data  <- read_csv("Geno5_EnvLOO/环境留一基因型5折组合交叉验证_最终结果.csv")

preprocess_data <- function(data, scenario_name) {
  data %>%
    filter(框架 == "BGLR", 预测类型 == "截距+斜率") %>%
    filter(方法类型 %in% c("MAS", "GBLUP", "MAS+GBLUP")) %>%
    mutate(
      场景   = scenario_name,
      实际值 = as.numeric(实际值),
      预测值 = as.numeric(预测值),
      环境ID = as.factor(环境ID)
    )
}

geno5_processed  <- preprocess_data(geno5_data,  "基因型5折交叉验证")
envloo_processed <- preprocess_data(envloo_data, "环境留一交叉验证")
combo_processed  <- preprocess_data(combo_data,  "组合交叉验证")

all_data <- bind_rows(geno5_processed, envloo_processed, combo_processed)

all_data$方法类型 <- factor(
  all_data$方法类型,
  levels = c("MAS", "GBLUP", "MAS+GBLUP"),
  labels = c("Marker-Assisted\nSelection (MAS)",
             "Genomic Best Linear\nUnbiased Prediction (GBLUP)",
             "MAS+GBLUP\nCombined")
)

all_data$场景 <- factor(
  all_data$场景,
  levels = c("环境留一交叉验证", "基因型5折交叉验证", "组合交叉验证"),
  labels = c("环境留一验证\n(已知基因型 × 新环境)",
             "基因型5折验证\n(新基因型 × 已知环境)",
             "组合交叉验证\n(新基因型 × 新环境)")
)

get_color_palette <- function(n) {
  base_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
  )
  if (n <= length(base_colors)) base_colors[1:n]
  else hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
}

env_list_all <- sort(unique(all_data$环境ID))
n_env_all <- length(env_list_all)
global_color_palette <- get_color_palette(n_env_all)
names(global_color_palette) <- as.character(env_list_all)

top8_env <- all_data %>%
  count(环境ID, sort = TRUE, name = "n") %>%
  slice_head(n = 8) %>%
  pull(环境ID) %>%
  as.character()

top8_color_mapping <- global_color_palette[names(global_color_palette) %in% top8_env]
top8_color_mapping <- top8_color_mapping[top8_env]

env_r_data <- all_data %>%
  filter(as.character(环境ID) %in% top8_env) %>%
  group_by(方法类型, 场景, 环境ID) %>%
  summarise(
    r_env = cor(实际值, 预测值, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(env_label = sprintf("%s (r = %.2f)", as.character(环境ID), r_env))

all_data <- all_data %>%
  left_join(env_r_data, by = c("方法类型", "场景", "环境ID"))

label_map <- all_data %>%
  filter(!is.na(env_label)) %>%
  distinct(env_label, 环境ID) %>%
  mutate(环境ID = as.character(环境ID))

label_color_mapping <- top8_color_mapping[as.character(label_map$环境ID)]
names(label_color_mapping) <- label_map$env_label

correlation_data <- all_data %>%
  group_by(方法类型, 场景) %>%
  summarise(
    r    = round(cor(实际值, 预测值, use = "complete.obs"), 3),
    rmse = round(sqrt(mean((实际值 - 预测值)^2, na.rm = TRUE)), 2),
    .groups = "drop"
  )

expand_range <- function(rng, frac = 0.30) {
  span <- diff(rng)
  c(rng[1] - frac * span, rng[2] + frac * span)
}

x_range0 <- range(all_data$实际值, na.rm = TRUE)
y_range0 <- range(all_data$预测值, na.rm = TRUE)

x_range <- expand_range(x_range0, 0.30)
y_range <- expand_range(y_range0, 0.30)

top_center_x <- x_range[1] + 0.50 * diff(x_range)
top_y_1 <- y_range[2] - 0.03 * diff(y_range)
top_y_2 <- y_range[2] - 0.11 * diff(y_range)

make_inset_grob <- function(scenario_label) {
  arrow_target <- "2"
  if (grepl("环境留一验证", scenario_label)) arrow_target <- "3"
  if (grepl("组合交叉验证", scenario_label)) arrow_target <- "4"

  pos <- list(
    "1" = c(0.25, 0.75),
    "2" = c(0.75, 0.75),
    "3" = c(0.25, 0.25),
    "4" = c(0.75, 0.25)
  )
  p1 <- pos[["1"]]
  pT <- pos[[arrow_target]]

  # 箭头更长一点：frac_keep = 0.80（更接近目标格，但仍避免压住数字）
  shrink_arrow <- function(a, b, frac_keep = 0.80) {
    c(a[1] + frac_keep * (b[1] - a[1]), a[2] + frac_keep * (b[2] - a[2]))
  }
  pT2 <- shrink_arrow(p1, pT, frac_keep = 0.80)

  grobTree(
    rectGrob(x = 0.5, y = 0.5, width = 1, height = 1,
             gp = gpar(col = "black", fill = NA, lwd = 0.8)),
    segmentsGrob(x0 = 0.5, x1 = 0.5, y0 = 0, y1 = 1, gp = gpar(lwd = 0.8)),
    segmentsGrob(x0 = 0, x1 = 1, y0 = 0.5, y1 = 0.5, gp = gpar(lwd = 0.8)),

    textGrob("1", x = 0.25, y = 0.75, gp = gpar(fontsize = 9, fontface = "bold", col = "black")),
    textGrob("2", x = 0.75, y = 0.75, gp = gpar(fontsize = 9, col = ifelse(arrow_target=="2","black","grey70"))),
    textGrob("3", x = 0.25, y = 0.25, gp = gpar(fontsize = 9, col = ifelse(arrow_target=="3","black","grey70"))),
    textGrob("4", x = 0.75, y = 0.25, gp = gpar(fontsize = 9, col = ifelse(arrow_target=="4","black","grey70"))),

    segmentsGrob(
      x0 = p1[1], y0 = p1[2],
      x1 = pT2[1], y1 = pT2[2],
      arrow = arrow(type = "closed", length = unit(0.08, "inches")),
      gp = gpar(lwd = 1.0, col = "black")
    ),

    textGrob("Environment", x = 0.5, y = 1.15, gp = gpar(fontsize = 8)),
    textGrob("Tested",   x = 0.25, y = 1.05, gp = gpar(fontsize = 7)),
    textGrob("Untested", x = 0.75, y = 1.05, gp = gpar(fontsize = 7)),

    textGrob("Genotype", x = -0.15, y = 0.5, rot = 90, gp = gpar(fontsize = 8)),
    textGrob("Tested",   x = -0.05, y = 0.75, rot = 90, gp = gpar(fontsize = 7)),
    textGrob("Untested", x = -0.05, y = 0.25, rot = 90, gp = gpar(fontsize = 7))
  )
}

create_scatter_plot <- function(data, method_type, scenario_name) {
  plot_data <- data %>% filter(方法类型 == method_type, 场景 == scenario_name)
  plot_data_top8  <- plot_data %>% filter(!is.na(env_label))
  plot_data_other <- plot_data %>% filter(is.na(env_label))

  stats_info <- correlation_data %>%
    filter(方法类型 == method_type, 场景 == scenario_name)

  inset_grob <- make_inset_grob(as.character(scenario_name))

  ggplot() +
    geom_point(
      data = plot_data_other,
      aes(x = 实际值, y = 预测值),
      alpha = 0.35, size = 1.3, color = "grey70"
    ) +
    geom_point(
      data = plot_data_top8,
      aes(x = 实际值, y = 预测值, color = env_label),
      alpha = 0.80, size = 1.5
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    xlim(x_range[1], x_range[2]) +
    ylim(y_range[1], y_range[2]) +
    labs(x = "Observed TKW (g)", y = "Predicted TKW (g)", color = "Environment") +

    annotate("text", x = top_center_x, y = top_y_1,
             label = sprintf("r = %.3f", stats_info$r),
             hjust = 0.5, vjust = 1, size = 3.6) +
    annotate("text", x = top_center_x, y = top_y_2,
             label = sprintf("RMSE = %.2f", stats_info$rmse),
             hjust = 0.5, vjust = 1, size = 3.6) +

    annotation_custom(
      grob = inset_grob,
      xmin = x_range[1] + 0.02 * diff(x_range),
      xmax = x_range[1] + 0.30 * diff(x_range),
      ymin = y_range[2] - 0.30 * diff(y_range),
      ymax = y_range[2] - 0.02 * diff(y_range)
    ) +

    scale_color_manual(values = label_color_mapping) +

    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = unit(c(0.25, 0.2, 0.2, 0.2), "cm"),

      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text  = element_text(size = 7),
      legend.key.size = unit(0.28, "cm")
    ) +
    guides(color = guide_legend(
      title.position = "top",
      title.hjust = 0,
      override.aes = list(alpha = 1, size = 2)
    ))
}

method_types <- levels(all_data$方法类型)
scenarios    <- levels(all_data$场景)

plots <- list()
for (i in seq_along(method_types)) {
  for (j in seq_along(scenarios)) {
    plots[[(i - 1) * 3 + j]] <- create_scatter_plot(all_data, method_types[i], scenarios[j])
  }
}

# -------------------------
# ✅ 修改后的标签添加函数：使用npc坐标系统，固定在绘图区域左上角
# -------------------------
add_outer_panel_label <- function(plot, label) {
  g <- ggplotGrob(plot)

  # 使用npc坐标系统，精确控制在绘图区域左上角
  label_grob <- textGrob(
    label,
    x = unit(0.02, "npc"),   # 距离左边缘2%
    y = unit(0.98, "npc"),   # 距离上边缘2%
    just = c("left", "top"),  # 左对齐，顶部对齐
    gp = gpar(fontsize = 16, fontface = "bold", col = "black")
  )

  grobTree(g, label_grob)
}

panel_labels <- c("a","b","c","d","e","f","g","h","k")
plots_labeled <- list()
for (k in seq_along(plots)) {
  plots_labeled[[k]] <- add_outer_panel_label(plots[[k]], panel_labels[[k]])
}

scenario_labels <- lapply(scenarios, function(x) {
  textGrob(x, gp = gpar(fontsize = 11, fontface = "bold"))
})

method_labels <- lapply(method_types, function(x) {
  textGrob(x, rot = 90, gp = gpar(fontsize = 11, fontface = "bold"))
})

layout_matrix <- rbind(
  c(NA, 10, 11, 12),
  c(13, 1,  2,  3 ),
  c(14, 4,  5,  6 ),
  c(15, 7,  8,  9 )
)

all_grobs <- c(plots_labeled, scenario_labels, method_labels)

final_plot <- grid.arrange(
  grobs = all_grobs,
  layout_matrix = layout_matrix,
  heights = c(0.5, 3, 3, 3),
  widths  = c(0.5, 3, 3, 3),
  top = NULL,
  padding = unit(0.5, "cm")
)

ggsave("TKW_GS_Performance_Top8Env_withInset.png", final_plot, width = 16, height = 14, dpi = 300)
ggsave("TKW_GS_Performance_Top8Env_withInset.pdf", final_plot, width = 16, height = 14)

cat("\n已保存：TKW_GS_Performance_Top8Env_withInset.png / .pdf\n")
cat("\nTop 8 环境ID（按样本量排序）：\n")
print(top8_env)

cat("\n全局环境颜色映射（用于整篇文章一致性）：\n")
print(global_color_palette)

cat("\nTop8 环境颜色映射（用于本图）：\n")
print(top8_color_mapping)