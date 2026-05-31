# ============================================================
# 风格：散点(每环境均值点) + 线性回归线 + 线尾QTL标签(含PVE)
# X轴：PAR_TEMP&GS67真实数值（非等距），环境名称标注在X轴上方
# Y轴：QTL Effect (g)
#
# 修复1：回归线基于PAR_TEMP&GS67真实值建模，保证线性
# 修复2：X轴向右延申 + 关闭裁剪 + 增大右边距，避免QTL标签显示不全
# 修复3：环境名旋转45度，并用ggrepel在X轴上方避免相互挤压
# 修复4：X轴刻度保留三位有效数字（signif=3）
# 修复5：QTL标签改为：QTLname (Phenotype_Variance_Explained(%))
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggrepel)
  library(scales)
})

# ==================== 1. 路径设置 ====================
base_dir   <- "/mnt/7t_storage/zhangcl/TKW/"
output_dir <- file.path(base_dir, "QTL_plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

slope_qtl_path <- file.path(base_dir, "QTL4.csv")
base_env_path  <- file.path(base_dir, "GAPIT_BLINK_Results/Single_Environment")
ecs_path       <- file.path(base_dir, "ECs_results.csv")

env_dirs <- c("2024NY", "2024PY", "2024YL", "2024ZMD",
              "2025YLW", "2025YLZ", "2025ZMDW", "2025ZMDZ")

# ==================== 2. 读取数据 ====================
cat("读取QTL数据...\n")
slope_qtl <- read_csv(slope_qtl_path, show_col_types = FALSE)

cat("读取环境协变量数据...\n")
ecs_data <- read_csv(ecs_path, show_col_types = FALSE)
if ("env" %in% colnames(ecs_data)) ecs_data <- ecs_data %>% rename(env_code = env)

# ==================== 2.1 提取每个QTL的PVE，并做成标签 ====================
# 兼容列名：Phenotype_Variance_Explained(%) 可能有大小写/空格差异
pve_col <- "Phenotype_Variance_Explained(%)"
if (!(pve_col %in% colnames(slope_qtl))) {
  hit <- grep("Phenotype_Variance_Explained", colnames(slope_qtl), ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) {
    message(sprintf("未找到列 '%s'，使用 '%s' 代替", pve_col, hit[1]))
    pve_col <- hit[1]
  } else {
    stop("QTL.csv 中未找到 Phenotype_Variance_Explained(%) 列")
  }
}

# 每个QTLname只取一个PVE（如果重复，取第一个非NA；你也可以改成mean/median）
qtl_pve <- slope_qtl %>%
  select(QTLname, PVE = all_of(pve_col)) %>%
  filter(!is.na(QTLname)) %>%
  group_by(QTLname) %>%
  summarise(
    PVE = dplyr::first(PVE[!is.na(PVE)]),
    .groups = "drop"
  ) %>%
  mutate(
    # 格式化：保留2位小数（可改）
    PVE_fmt = ifelse(is.na(PVE), "NA", sprintf("%.2f%%", as.numeric(PVE))),
    QTL_label = paste0(QTLname, " (", PVE_fmt, ")")
  )

# ==================== 3. 提取每个环境的效应值 ====================
all_effects <- slope_qtl %>%
  select(QTLname, SNP, Chromosome, Position) %>%
  distinct()

cat("\n正在从各个环境文件中提取效应值...\n")
for (env in env_dirs) {
  file_path <- file.path(base_env_path, env,
                         "GAPIT.Association.GWAS_Results.TKW.BLINK.Trait(NYC).csv")
  if (file.exists(file_path)) {
    env_data <- read_csv(file_path, show_col_types = FALSE)
    env_effects <- env_data %>% select(SNP, Effect) %>% rename(!!env := Effect)
    all_effects <- all_effects %>% left_join(env_effects, by = "SNP")

    matched <- sum(!is.na(all_effects[[env]]))
    total   <- nrow(all_effects)
    cat(sprintf("环境 %s: 匹配 %d/%d (%.1f%%)\n", env, matched, total, 100*matched/total))
  } else {
    warning(paste("文件不存在:", file_path))
    all_effects[[env]] <- NA
  }
}

# ==================== 4. 长格式 + 合并 PAR_TEMP&GS67 ====================
long_format <- all_effects %>%
  pivot_longer(cols = all_of(env_dirs),
               names_to = "env_code",
               values_to = "Effect") %>%
  filter(!is.na(Effect)) %>%
  select(QTLname, SNP, Chromosome, Position, env_code, Effect)

x_var_name <- "PAR_TEMP&GS67"
if (!(x_var_name %in% colnames(ecs_data))) {
  hit <- grep("PAR_TEMP&GS67", colnames(ecs_data), ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) {
    message(sprintf("未找到列 '%s'，使用 '%s' 代替", x_var_name, hit[1]))
    x_var_name <- hit[1]
  } else {
    stop("ECs_results.csv 中未找到 PAR_TEMP&GS67 列")
  }
}

long_format <- long_format %>%
  left_join(ecs_data %>% select(env_code, Xvar = all_of(x_var_name)), by = "env_code")

data_clean <- long_format %>%
  filter(!is.na(Effect), !is.na(Xvar))

# ==================== 5. 每QTL×每环境：均值点 ====================
qtl_summary <- data_clean %>%
  group_by(QTLname, env_code, Xvar) %>%
  summarise(
    mean_effect = mean(Effect, na.rm = TRUE),
    n_snps = n(),
    .groups = "drop"
  )

# ==================== 6. 选择要绘制的QTL ====================
qtl_list_from_image <- c("QTKW-4A-2", "TaGW2L-7D", "QTKW-1A",
                         "QTKW-2D", "QTKW-7D-1", "QTKW-4A-1")

available_qtl <- unique(qtl_summary$QTLname)
selected_qtl <- intersect(qtl_list_from_image, available_qtl)
if (length(selected_qtl) == 0) selected_qtl <- available_qtl

plot_df <- qtl_summary %>%
  filter(QTLname %in% selected_qtl) %>%
  mutate(QTLname = factor(QTLname, levels = selected_qtl))

# 合并PVE标签到数据中
plot_df <- plot_df %>%
  left_join(qtl_pve %>% select(QTLname, QTL_label), by = "QTLname") %>%
  mutate(QTL_label = ifelse(is.na(QTL_label), as.character(QTLname), QTL_label))

# ==================== 7. 环境顺序（按真实Xvar） ====================
env_order <- plot_df %>%
  select(env_code, Xvar) %>%
  distinct() %>%
  arrange(Xvar)

cat("\n环境按PAR_TEMP&GS67排序 (真实值):\n")
print(env_order)

# ==================== 8. 回归线数据（基于真实Xvar） ====================
lm_data <- data.frame()
for (qtl in selected_qtl) {
  qtl_data <- plot_df %>% filter(QTLname == qtl)
  if (nrow(qtl_data) >= 2) {
    lm_model <- lm(mean_effect ~ Xvar, data = qtl_data)
    x_range <- range(env_order$Xvar, na.rm = TRUE)
    x_seq <- seq(x_range[1], x_range[2], length.out = 250)
    y_pred <- predict(lm_model, newdata = data.frame(Xvar = x_seq))
    lm_data <- rbind(lm_data, data.frame(
      QTLname = qtl,
      Xvar = x_seq,
      mean_effect = y_pred
    ))
  }
}

# ==================== 9. 颜色 ====================
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
            "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
if (length(selected_qtl) > length(colors)) {
  if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
  colors <- viridis::viridis(length(selected_qtl))
} else {
  colors <- colors[1:length(selected_qtl)]
}
names(colors) <- selected_qtl

# ==================== 10. QTL标签位置（回归线末端） ====================
label_positions <- lm_data %>%
  group_by(QTLname) %>%
  filter(Xvar == max(Xvar)) %>%
  ungroup()

# 把PVE标签带上
label_positions <- label_positions %>%
  left_join(plot_df %>% distinct(QTLname, QTL_label), by = "QTLname")

# 标签X偏移，让文字离开末端点，避免压线
label_positions <- label_positions %>%
  mutate(label_x = Xvar + (max(env_order$Xvar) - min(env_order$Xvar)) * 0.02)

# ==================== 11. X轴刻度：三位有效数字（并避免挤压） ====================
x_min <- min(env_order$Xvar, na.rm = TRUE)
x_max <- max(env_order$Xvar, na.rm = TRUE)
x_pad <- (x_max - x_min) * 0.20

# 默认：每个环境一个刻度（若挤压，你可改成最多显示6个刻度）
x_breaks <- env_order$Xvar
x_labels <- format(signif(x_breaks, 3), trim = TRUE, scientific = FALSE)

# ==================== 12. 环境名放到上方（ggrepel） ====================
y_range <- range(plot_df$mean_effect, na.rm = TRUE)
y_padding <- (y_range[2] - y_range[1]) * 0.15
y_top <- y_range[2] + y_padding + (y_range[2] - y_range[1]) * 0.05

env_label_df <- env_order %>%
  mutate(y_top = y_top)

# ==================== 13. 主题 ====================
custom_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 25)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 12)),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 5)),
    axis.text.y  = element_text(colour = "black", size = 10),
    legend.position = "none",
    plot.margin = margin(90, 140, 25, 25), # 上给环境名，右给QTL+PVE标签
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# ==================== 14. 作图 ====================
p <- ggplot() +
  geom_point(
    data = plot_df,
    aes(x = Xvar, y = mean_effect, color = QTLname),
    size = 2.5, alpha = 0.7
  ) +
  geom_line(
    data = lm_data,
    aes(x = Xvar, y = mean_effect, color = QTLname, group = QTLname),
    linewidth = 1.0
  ) +
  # QTL线尾标签：QTLname (PVE%)
  geom_text_repel(
    data = label_positions,
    aes(x = label_x, y = mean_effect, label = QTL_label, color = QTLname),
    hjust = 0,
    vjust = 0.5,
    size = 3.8,
    fontface = "italic",
    direction = "y",
    nudge_x = 0.5,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.size = 0.3,
    min.segment.length = 0,
    seed = 123
  ) +
  scale_x_continuous(
    name = "PAR_TEMP&GS67",
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  # 环境名：45度，上方，防挤压
  geom_text_repel(
    data = env_label_df,
    aes(x = Xvar, y = y_top, label = env_code),
    inherit.aes = FALSE,
    angle = 45,
    direction = "x",
    nudge_y = 0.5,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2,
    min.segment.length = 0,
    seed = 456,
    size = 3.2,
    color = "black",
    hjust = 0.5,
    vjust = 0
  ) +
  labs(y = "QTL Effect (g)") +
  scale_color_manual(values = colors) +
  custom_theme +
  annotate(
    "text",
    x = x_min - (x_max - x_min) * 0.02,
    y = y_range[2] + y_padding * 0.8,
    label = "(b)",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4.5,
    color = "black"
  ) +
  coord_cartesian(
    xlim = c(x_min, x_max + x_pad),
    ylim = c(y_range[1] - y_padding, y_top + (y_top - y_range[1]) * 0.15),
    clip = "off"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

# ==================== 15. 保存 ====================
png_file <- file.path(output_dir, "QTL_by_PAR_TEMP_actual_values_withPVE.png")
pdf_file <- file.path(output_dir, "QTL_by_PAR_TEMP_actual_values_withPVE.pdf")

ggsave(png_file, p, width = 13, height = 7.2, dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 13, height = 7.2, bg = "white")

cat("\n图形已保存:\n")
cat("  - ", png_file, "\n", sep = "")
cat("  - ", pdf_file, "\n", sep = "")

# ==================== 16. 回归分析（基于真实Xvar） ====================
lm_results <- data.frame()
for (qtl in selected_qtl) {
  qtl_data <- plot_df %>% filter(QTLname == qtl)
  if (nrow(qtl_data) >= 2) {
    m <- lm(mean_effect ~ Xvar, data = qtl_data)
    s <- summary(m)
    lm_results <- rbind(lm_results, data.frame(
      QTLname = qtl,
      斜率 = round(coef(m)[2], 4),
      截距 = round(coef(m)[1], 4),
      R平方 = round(s$r.squared, 4),
      P值 = ifelse(s$coefficients[2, 4] < 0.001, "< 0.001",
                   round(s$coefficients[2, 4], 6)),
      显著性 = ifelse(s$coefficients[2, 4] < 0.001, "***",
                      ifelse(s$coefficients[2, 4] < 0.01, "**",
                             ifelse(s$coefficients[2, 4] < 0.05, "*", "ns"))),
      环境数 = nrow(qtl_data)
    ))
  }
}

if (nrow(lm_results) > 0) {
  regression_file <- file.path(output_dir, "QTL_regression_by_PAR_TEMP.csv")
  write_csv(lm_results, regression_file)
  cat("\n回归结果已保存到: ", regression_file, "\n", sep = "")
  cat("\n回归分析结果 (基于PAR_TEMP&GS67真实值):\n")
  print(lm_results)
}

# ==================== 17. 数据摘要 ====================
summary_stats <- plot_df %>%
  group_by(QTLname) %>%
  summarise(
    环境数 = n_distinct(env_code),
    环境列表 = paste(sort(unique(env_code)), collapse = ", "),
    PAR_TEMP_GS67范围 = sprintf("%.3f - %.3f", min(Xvar), max(Xvar)),
    效应值范围 = sprintf("%.3f - %.3f", min(mean_effect), max(mean_effect)),
    平均效应值 = sprintf("%.3f ± %.3f", mean(mean_effect), sd(mean_effect)),
    .groups = "drop"
  ) %>%
  left_join(qtl_pve %>% select(QTLname, PVE, PVE_fmt), by = "QTLname")

summary_file <- file.path(output_dir, "QTL_summary_by_PAR_TEMP_withPVE.csv")
write_csv(summary_stats, summary_file)
cat("\n数据摘要已保存到: ", summary_file, "\n", sep = "")
cat("\n数据摘要:\n")
print(summary_stats)

cat("\n完成。\n")
