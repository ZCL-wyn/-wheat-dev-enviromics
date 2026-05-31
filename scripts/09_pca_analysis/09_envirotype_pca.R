# 清理环境
rm(list = ls())

# 设置工作目录
setwd("C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\01缺失的遗传力\\")

# 加载所需的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(gridExtra)
library(viridis)

# 检查并安装缺失的包
required_packages <- c("ggplot2", "dplyr", "tidyr", "cowplot", "viridis", "grid", "gridExtra")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 读取数据
planting <- read.csv("planting.csv")
ecs_results <- read.csv("ECs_results.csv")

# 查看数据结构
cat("=== 数据文件信息 ===\n")
cat("planting.csv 维度:", dim(planting), "\n")
cat("planting.csv 列名:", paste(colnames(planting), collapse = ", "), "\n")
cat("\n")
cat("ECs_results.csv 维度:", dim(ecs_results), "\n")
cat("ECs_results.csv 列名:", paste(colnames(ecs_results), collapse = ", "), "\n")
cat("\n")

# 显示前几行数据
cat("planting.csv 前3行:\n")
print(head(planting, 3))
cat("\n")
cat("ECs_results.csv 前3行:\n")
print(head(ecs_results, 3))
cat("\n")

# 检查是否包含必要的列
required_planting_cols <- c("env", "latitude", "longitude", "state")
missing_planting_cols <- setdiff(required_planting_cols, colnames(planting))
if (length(missing_planting_cols) > 0) {
  stop(paste("planting.csv 缺少以下必要的列:", paste(missing_planting_cols, collapse = ", ")))
}

if (!"env" %in% colnames(ecs_results)) {
  stop("ECs_results.csv 必须包含 'env' 列")
}

# 合并数据
cat("正在合并数据...\n")
merged_data <- merge(planting, ecs_results, by = "env", all.x = TRUE)
cat("合并后数据维度:", dim(merged_data), "\n")
cat("合并后列名:", paste(colnames(merged_data), collapse = ", "), "\n")
cat("\n")

# 检查是否有缺失值
na_count <- sum(is.na(merged_data))
cat("合并数据中的缺失值总数:", na_count, "\n")
if (na_count > 0) {
  na_by_col <- colSums(is.na(merged_data))
  na_cols <- names(na_by_col[na_by_col > 0])
  cat("有缺失值的列数量:", length(na_cols), "\n")
  cat("缺失值最多的前10列:\n")
  na_by_col_sorted <- sort(na_by_col[na_by_col > 0], decreasing = TRUE)
  print(head(na_by_col_sorted, 10))
}
cat("\n")

# 准备环境因子数据用于 PCA
# 获取环境因子列（排除 env 列）
env_factor_cols <- setdiff(colnames(ecs_results), "env")
cat("环境因子总数:", length(env_factor_cols), "\n")
cat("前5个环境因子:", paste(head(env_factor_cols, 5), collapse = ", "), "\n")
cat("\n")

# 创建环境因子矩阵
env_matrix <- as.matrix(ecs_results[, env_factor_cols])
rownames(env_matrix) <- ecs_results$env

# 分析环境因子矩阵
cat("=== 环境因子矩阵分析 ===\n")
cat("原始环境因子矩阵维度:", dim(env_matrix), "\n")

# 1. 检查每列的缺失值
missing_counts <- colSums(is.na(env_matrix))
cat("各环境因子缺失值数量统计:\n")
print(summary(missing_counts))

# 找出有缺失值的列
cols_with_na <- names(missing_counts[missing_counts > 0])
cat("有缺失值的环境因子数量:", length(cols_with_na), "\n")
if (length(cols_with_na) > 0) {
  cat("前10个有缺失值的环境因子:\n")
  print(head(cols_with_na, 10))
}

# 找出缺失值比例高的列（比如超过50%缺失）
missing_prop <- missing_counts / nrow(env_matrix)
high_missing_cols <- names(missing_prop[missing_prop > 0.5])
if (length(high_missing_cols) > 0) {
  cat("缺失值比例超过50%的环境因子数量:", length(high_missing_cols), "\n")
} else {
  cat("没有环境因子缺失值比例超过50%\n")
}
cat("\n")

# 2. 检查每列的标准差（忽略缺失值）
col_sd <- apply(env_matrix, 2, function(x) {
  x_non_missing <- x[!is.na(x)]
  if (length(x_non_missing) < 2) return(0)  # 少于2个非缺失值，标准差为0
  sd(x_non_missing, na.rm = TRUE)
})

cat("各环境因子标准差统计:\n")
print(summary(col_sd))

# 找出标准差为0或接近0的列
zero_sd_cols <- names(col_sd[col_sd < 1e-10])  # 使用小阈值而不是精确为0
if (length(zero_sd_cols) > 0) {
  cat("标准差为零或接近零的环境因子数量:", length(zero_sd_cols), "\n")
} else {
  cat("没有标准差为零的环境因子\n")
}
cat("\n")

# 3. 确定要移除的环境因子
cols_to_remove <- unique(c(cols_with_na, high_missing_cols, zero_sd_cols))
cat("总共要移除的环境因子数量:", length(cols_to_remove), "\n")
cat("移除原因:\n")
cat("  - 有缺失值:", length(cols_with_na), "\n")
cat("  - 缺失值比例>50%:", length(high_missing_cols), "\n")
cat("  - 标准差为零:", length(zero_sd_cols), "\n")

# 4. 创建处理后的环境因子矩阵（移除有问题的列）
cols_to_keep <- setdiff(env_factor_cols, cols_to_remove)
cat("要保留的环境因子数量:", length(cols_to_keep), "\n")

if (length(cols_to_keep) == 0) {
  stop("错误：没有环境因子可以用于PCA分析！请检查数据。")
}

env_matrix_processed <- env_matrix[, cols_to_keep, drop = FALSE]
cat("处理后的环境因子矩阵维度:", dim(env_matrix_processed), "\n")

# 5. 再次检查处理后的矩阵是否有缺失值
cat("处理后矩阵缺失值数量:", sum(is.na(env_matrix_processed)), "\n")
if (sum(is.na(env_matrix_processed)) > 0) {
  cat("错误：处理后矩阵仍有缺失值！\n")
  na_by_col_processed <- colSums(is.na(env_matrix_processed))
  na_cols_processed <- names(na_by_col_processed[na_by_col_processed > 0])
  cat("仍有缺失值的列:", paste(na_cols_processed, collapse = ", "), "\n")
  
  # 移除这些列
  cat("移除这些列...\n")
  cols_to_keep <- setdiff(cols_to_keep, na_cols_processed)
  cols_to_remove <- c(cols_to_remove, na_cols_processed)
  env_matrix_processed <- env_matrix_processed[, cols_to_keep, drop = FALSE]
  cat("再次处理后矩阵维度:", dim(env_matrix_processed), "\n")
  cat("再次处理后缺失值数量:", sum(is.na(env_matrix_processed)), "\n")
}

# 6. 检查并处理无穷值
cat("处理后矩阵无穷值数量:", sum(is.infinite(env_matrix_processed)), "\n")
if (sum(is.infinite(env_matrix_processed)) > 0) {
  cat("警告：处理后的矩阵包含无穷值！\n")
  inf_cols <- apply(env_matrix_processed, 2, function(x) any(is.infinite(x)))
  if (any(inf_cols)) {
    cat("包含无穷值的列数量:", sum(inf_cols), "\n")
    # 移除包含无穷值的列
    env_matrix_processed <- env_matrix_processed[, !inf_cols, drop = FALSE]
    cols_to_keep <- cols_to_keep[!inf_cols]
    cat("移除无穷值列后矩阵维度:", dim(env_matrix_processed), "\n")
  }
}
cat("\n")

# 7. 标准化环境因子数据（PCA前通常需要标准化）
cat("正在标准化环境因子数据...\n")
env_matrix_scaled <- scale(env_matrix_processed, center = TRUE, scale = TRUE)

# 检查标准化后的矩阵
cat("标准化后矩阵维度:", dim(env_matrix_scaled), "\n")

# 8. 检查标准化后的矩阵是否有缺失值或无穷值
cat("标准化后矩阵缺失值数量:", sum(is.na(env_matrix_scaled)), "\n")
if (sum(is.na(env_matrix_scaled)) > 0) {
  cat("警告：标准化后矩阵有缺失值！\n")
  
  # 找出包含缺失值的列
  na_cols_scaled <- apply(env_matrix_scaled, 2, function(x) any(is.na(x)))
  if (any(na_cols_scaled)) {
    cat("包含缺失值的列数量:", sum(na_cols_scaled), "\n")
    cat("这些列将被移除:\n")
    na_col_names <- colnames(env_matrix_scaled)[na_cols_scaled]
    print(na_col_names)
    
    # 移除这些列
    env_matrix_scaled <- env_matrix_scaled[, !na_cols_scaled, drop = FALSE]
    cols_to_keep <- cols_to_keep[!na_cols_scaled]
    cols_to_remove <- c(cols_to_remove, na_col_names)
    cat("移除缺失值列后矩阵维度:", dim(env_matrix_scaled), "\n")
  }
}

cat("标准化后矩阵无穷值数量:", sum(is.infinite(env_matrix_scaled)), "\n")
if (sum(is.infinite(env_matrix_scaled)) > 0) {
  cat("警告：标准化后矩阵有无穷值！\n")
  
  # 找出包含无穷值的列
  inf_cols_scaled <- apply(env_matrix_scaled, 2, function(x) any(is.infinite(x)))
  if (any(inf_cols_scaled)) {
    cat("包含无穷值的列数量:", sum(inf_cols_scaled), "\n")
    
    # 移除这些列
    env_matrix_scaled <- env_matrix_scaled[, !inf_cols_scaled, drop = FALSE]
    cols_to_keep <- cols_to_keep[!inf_cols_scaled]
    cat("移除无穷值列后矩阵维度:", dim(env_matrix_scaled), "\n")
  }
}

# 9. 最终检查矩阵完整性
cat("\n最终检查矩阵完整性...\n")
if (any(is.na(env_matrix_scaled))) {
  cat("错误：矩阵中仍有缺失值！详细统计:\n")
  na_count_by_col <- colSums(is.na(env_matrix_scaled))
  print(na_count_by_col[na_count_by_col > 0])
  stop("无法进行PCA分析，请检查数据。")
}

if (any(is.infinite(env_matrix_scaled))) {
  cat("错误：矩阵中仍有无穷值！详细统计:\n")
  inf_count_by_col <- colSums(is.infinite(env_matrix_scaled))
  print(inf_count_by_col[inf_count_by_col > 0])
  stop("无法进行PCA分析，请检查数据。")
}

cat("矩阵完整性检查通过！\n")

# 检查标准化后的统计
cat("标准化后列均值（应该接近0）:\n")
col_means_scaled <- colMeans(env_matrix_scaled)
print(summary(col_means_scaled))
cat("标准化后列标准差（应该为1）:\n")
col_sd_scaled <- apply(env_matrix_scaled, 2, sd)
print(summary(col_sd_scaled))
cat("\n")

# 进行 PCA 分析
cat("正在进行 PCA 分析...\n")
pca_result <- prcomp(env_matrix_scaled, center = FALSE, scale. = FALSE)  # 已经标准化过了

# 提取 PCA 结果
eigenvalues <- pca_result$sdev^2
total_variance <- sum(eigenvalues)
variance_proportion <- eigenvalues / total_variance
cumulative_variance <- cumsum(variance_proportion)

# 获取主成分得分
pc_scores <- pca_result$x

# 创建输出目录
outfolder <- "2_PC_ecov"
if (!dir.exists(outfolder)) {
  dir.create(outfolder, recursive = TRUE)
}

# 保存 PCA 结果
save(pca_result, eigenvalues, variance_proportion, cumulative_variance, 
     file = paste0(outfolder, "/pca_results.RData"))
cat("PCA 结果已保存到:", paste0(outfolder, "/pca_results.RData"), "\n")

# 打印 PCA 摘要
cat("\n=== PCA 分析结果摘要 ===\n")
cat("总方差:", round(total_variance, 2), "\n")
cat("主成分总数:", length(eigenvalues), "\n")
cat("前10个主成分解释的方差比例:\n")
for (i in 1:min(10, length(variance_proportion))) {
  cat(sprintf("  PC%d: %.3f (%.2f%%, 累计: %.2f%%)\n", 
              i, eigenvalues[i], variance_proportion[i] * 100, cumulative_variance[i] * 100))
}

# 创建完整的 rotate_PC 函数
cat("\n正在创建坐标轴旋转函数...\n")
rotate_PC <- function(PC, latitude, longitude, signPC = c(1, 1)) {
  # 确保输入正确
  if (ncol(PC) < 2) {
    stop("PC 矩阵必须至少有两列")
  }
  if (length(latitude) != nrow(PC) || length(longitude) != nrow(PC)) {
    stop("latitude 和 longitude 的长度必须与 PC 的行数相同")
  }
  
  # 应用符号调整
  PC_adjusted <- sweep(PC[, 1:2], 2, signPC, FUN = "*")
  
  # 计算最优旋转角度（使 PC1 与经度、PC2 与纬度最大相关）
  angles <- seq(0, 360, by = 1)  # 1度间隔搜索
  best_correlation <- -Inf
  best_angle <- 0
  best_PC_rotated <- NULL
  
  for (angle in angles) {
    # 将角度转换为弧度
    theta <- angle * pi / 180
    
    # 创建旋转矩阵
    rotation_matrix <- matrix(c(cos(theta), -sin(theta), 
                                sin(theta), cos(theta)), 
                              nrow = 2, byrow = TRUE)
    
    # 旋转 PC
    PC_rotated <- PC_adjusted %*% rotation_matrix
    
    # 计算与经纬度的相关性
    cor_longitude <- cor(PC_rotated[, 1], longitude, use = "complete.obs")
    cor_latitude <- cor(PC_rotated[, 2], latitude, use = "complete.obs")
    total_correlation <- abs(cor_longitude) + abs(cor_latitude)
    
    if (total_correlation > best_correlation) {
      best_correlation = total_correlation
      best_angle = angle
      best_PC_rotated = PC_rotated
    }
  }
  
  cat(sprintf("最优旋转角度: %.1f 度\n", best_angle))
  cat(sprintf("旋转后 PC1 与经度的相关性: %.3f\n", cor(best_PC_rotated[, 1], longitude, use = "complete.obs")))
  cat(sprintf("旋转后 PC2 与纬度的相关性: %.3f\n", cor(best_PC_rotated[, 2], latitude, use = "complete.obs")))
  
  # 应用最优旋转
  theta <- best_angle * pi / 180
  rotation_matrix <- matrix(c(cos(theta), -sin(theta), 
                              sin(theta), cos(theta)), 
                            nrow = 2, byrow = TRUE)
  PC_rotated <- PC_adjusted %*% rotation_matrix
  
  # 计算旋转后的轴方向
  x_axis <- c(1, 0) %*% rotation_matrix
  y_axis <- c(0, 1) %*% rotation_matrix
  
  # 归一化到单位长度
  x_axis <- x_axis / sqrt(sum(x_axis^2))
  y_axis <- y_axis / sqrt(sum(y_axis^2))
  
  # 创建结果列表
  result <- list(
    PC_rotated = PC_rotated,
    rotation_angle = best_angle,
    x_axis = x_axis,
    y_axis = y_axis,
    PC_original = PC_adjusted,
    correlations = list(
      longitude_cor = cor(PC_rotated[, 1], longitude, use = "complete.obs"),
      latitude_cor = cor(PC_rotated[, 2], latitude, use = "complete.obs")
    )
  )
  
  return(result)
}

# 准备用于旋转的数据
cat("\n准备用于坐标轴旋转的数据...\n")
# 获取有环境因子数据的 env
env_with_data <- rownames(env_matrix_scaled)

# 获取这些 env 的位置信息
env_locations <- merged_data %>% 
  filter(env %in% env_with_data) %>%
  select(env, latitude, longitude, state)

# 检查是否有重复的 env
duplicate_env <- env_locations$env[duplicated(env_locations$env)]
if (length(duplicate_env) > 0) {
  cat("警告: 发现重复的 env, 使用第一个记录:\n")
  print(duplicate_env)
  env_locations <- env_locations[!duplicated(env_locations$env), ]
}

# 确保顺序匹配
rownames(pc_scores) <- ecs_results$env
pc_scores_filtered <- pc_scores[env_locations$env, ]

# 检查数据一致性
if (!all(env_locations$env == rownames(pc_scores_filtered))) {
  stop("环境因子数据和位置数据的顺序不匹配")
}

# 执行坐标轴旋转
cat("正在执行坐标轴旋转...\n")
rotation_result <- rotate_PC(
  PC = pc_scores_filtered,
  latitude = env_locations$latitude,
  longitude = env_locations$longitude,
  signPC = c(1, 1)
)

# 创建用于绘图的数据框
plot_data <- data.frame(
  env = env_locations$env,
  PC1_original = pc_scores_filtered[, 1],
  PC2_original = pc_scores_filtered[, 2],
  PC1_rotated = rotation_result$PC_rotated[, 1],
  PC2_rotated = rotation_result$PC_rotated[, 2],
  latitude = env_locations$latitude,
  longitude = env_locations$longitude,
  state = env_locations$state,
  stringsAsFactors = FALSE
)

# 计算每个 env 的显示标签（如果 env 太长，可以截取）
plot_data$env_display <- plot_data$env

# 计算方差解释比例
var_pc1 <- variance_proportion[1] * 100
var_pc2 <- variance_proportion[2] * 100

# 创建轴标签（包含方差解释比例）
x_label <- paste0("PC1 (", round(var_pc1, 1), "%)")
y_label <- paste0("PC2 (", round(var_pc2, 1), "%)")

# 创建旋转轴箭头数据（包含反向延长线）
# 确定坐标轴的范围
pc1_range <- range(plot_data$PC1_rotated, na.rm = TRUE)
pc2_range <- range(plot_data$PC2_rotated, na.rm = TRUE)

# 计算轴的长度（使用数据范围的70%）
axis_length_pc1 <- diff(pc1_range) * 0.7 / 2
axis_length_pc2 <- diff(pc2_range) * 0.7 / 2
axis_length <- max(axis_length_pc1, axis_length_pc2)

# PC1轴：正方向和反方向
x_axis_arrow_forward <- data.frame(
  x = 0, y = 0,
  xend = rotation_result$x_axis[1] * axis_length,
  yend = rotation_result$x_axis[2] * axis_length
)

x_axis_arrow_backward <- data.frame(
  x = 0, y = 0,
  xend = -rotation_result$x_axis[1] * axis_length,
  yend = -rotation_result$x_axis[2] * axis_length
)

# PC2轴：正方向和反方向
y_axis_arrow_forward <- data.frame(
  x = 0, y = 0,
  xend = rotation_result$y_axis[1] * axis_length,
  yend = rotation_result$y_axis[2] * axis_length
)

y_axis_arrow_backward <- data.frame(
  x = 0, y = 0,
  xend = -rotation_result$y_axis[1] * axis_length,
  yend = -rotation_result$y_axis[2] * axis_length
)

# 创建PC1 vs PC2 散点图，使用渐变色表示纬度，env标签在图中
cat("\n正在创建 Panel A: 旋转后的 PC1 vs PC2 图...\n")

# 计算纬度的范围
lat_range <- range(plot_data$latitude, na.rm = TRUE)
cat("纬度范围: ", lat_range[1], " 到 ", lat_range[2], "\n")

# 计算每个env点的中点，用于标签放置
# 这里我们简化，直接用每个点的坐标作为标签位置

# 创建PC1 vs PC2散点图
pca_plot <- ggplot(plot_data, aes(x = PC1_rotated, y = PC2_rotated)) +
  # 添加PC1轴的反向延长线（无箭头）
  geom_segment(data = x_axis_arrow_backward,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "red", size = 0.8, alpha = 0.7, linetype = "solid") +
  # 添加PC1轴的正向箭头
  geom_segment(data = x_axis_arrow_forward, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red", size = 0.8, alpha = 0.7) +
  # 添加PC2轴的反向延长线（无箭头）
  geom_segment(data = y_axis_arrow_backward,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "blue", size = 0.8, alpha = 0.7, linetype = "solid") +
  # 添加PC2轴的正向箭头
  geom_segment(data = y_axis_arrow_forward,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "blue", size = 0.8, alpha = 0.7) +
  # 添加数据点，使用渐变色表示纬度，所有点都涂上颜色
  geom_point(aes(fill = latitude), 
             size = 3, shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  # 添加env标签在点上
  geom_text(aes(label = env), size = 2.5, color = "black", vjust = 1.5) +
  # 添加坐标轴标签
  annotate("text", 
           x = x_axis_arrow_forward$xend * 1.1, 
           y = x_axis_arrow_forward$yend * 1.1,
           label = "PC1", color = "red", size = 4.5, fontface = "bold") +
  annotate("text",
           x = y_axis_arrow_forward$xend * 1.1,
           y = y_axis_arrow_forward$yend * 1.1,
           label = "PC2", color = "blue", size = 4.5, fontface = "bold") +
  # 设置坐标轴
  labs(
    x = x_label,
    y = y_label,
    fill = "Latitude"
  ) +
  # 应用渐变色标尺 - 使用viridis的mako选项，避免黄色不明显的问题
  scale_fill_viridis_c(
    option = "mako",  # mako色系没有黄色，更适合表示纬度
    direction = -1,   # 反向颜色，使高纬度显示为亮色
    breaks = seq(floor(lat_range[1]), ceiling(lat_range[2]), by = 2),
    guide = guide_colorbar(
      title = "Latitude",
      title.position = "top",
      barwidth = unit(0.5, "cm"),
      barheight = unit(2, "cm"),
      ticks = FALSE
    )
  ) +
  # 设置主题 - 参考您提供的代码风格
  theme_minimal() +
  theme(
    legend.position = c(0.92, 0.15),  # 将颜色图例放在图形内部右上角
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    plot.title = element_blank(),  # 移除标题
    panel.grid = element_blank(),  # 移除所有网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 添加坐标轴线
  ) +
  # 设置坐标轴范围
  scale_x_continuous(
    limits = c(
      min(plot_data$PC1_rotated, -axis_length * 1.2, na.rm = TRUE),
      max(plot_data$PC1_rotated, axis_length * 1.2, na.rm = TRUE)
    )
  ) +
  scale_y_continuous(
    limits = c(
      min(plot_data$PC2_rotated, -axis_length * 1.2, na.rm = TRUE),
      max(plot_data$PC2_rotated, axis_length * 1.2, na.rm = TRUE)
    )
  )

#==================================================
# Panel B: 累计贡献率图 - 改进版本（不要柱状图）
#==================================================

cat("\n正在创建 Panel B: 累计贡献率图...\n")

# 创建 scree plot 数据
scree_data <- data.frame(
  PC = 1:length(variance_proportion),
  Eigenvalue = eigenvalues,
  Variance = variance_proportion * 100,
  Cumulative = cumulative_variance * 100
)

# 计算95%和50%阈值的位置
threshold_95 <- which(cumulative_variance * 100 >= 95)[1]
threshold_50 <- which(cumulative_variance * 100 >= 50)[1]

# 如果没有达到阈值，设为最后一个成分
if (is.na(threshold_95)) threshold_95 <- length(cumulative_variance)
if (is.na(threshold_50)) threshold_50 <- length(cumulative_variance)

# 确定要绘制的主成分数量：取前200个或所有主成分（如果少于200）
n_components_to_plot <- min(200, length(cumulative_variance))

cat("将绘制前", n_components_to_plot, "个主成分的累计方差图\n")
cat("累计方差达到95%需要", threshold_95, "个主成分\n")
cat("累计方差达到50%需要", threshold_50, "个主成分\n")

# 将累计方差转换为比例（0-1）
variance_df_subset <- scree_data[1:n_components_to_plot, ]
variance_df_subset$Cumulative_Proportion <- variance_df_subset$Cumulative / 100

# 获取50%和95%阈值线的y值
y_50 <- cumulative_variance[threshold_50] / 100
y_95 <- cumulative_variance[threshold_95] / 100

# 创建累计贡献率图
cumulative_plot <- ggplot(variance_df_subset, 
                          aes(x = PC, y = Cumulative_Proportion)) +
  # 绘制折线
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange", size = 1) +
  
  # 添加95%阈值线和标注
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = threshold_95, linetype = "dashed", color = "blue", linewidth = 0.8) +
  
  # 添加50%阈值线和标注
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = threshold_50, linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 在交点处添加标签
  annotate("point", 
           x = threshold_95, 
           y = y_95, 
           color = "blue", 
           size = 3) +
  annotate("text", 
           x = threshold_95, 
           y = y_95 + 0.03, 
           label = paste0(threshold_95, " (95%)"), 
           color = "blue", 
           size = 4,
           hjust = 0.5,
           vjust = 0) +
  
  annotate("point", 
           x = threshold_50, 
           y = y_50, 
           color = "red", 
           size = 3) +
  annotate("text", 
           x = threshold_50, 
           y = y_50 + 0.03, 
           label = paste0(threshold_50, " (50%)"), 
           color = "red", 
           size = 4,
           hjust = 0.5,
           vjust = 0) +
  
  # 设置坐标轴和标题
  labs(
    x = "Number of eigenvectors",
    y = "Proportion of variance"
  ) +
  scale_x_continuous(
    breaks = seq(0, n_components_to_plot, by = 25),
    limits = c(0, n_components_to_plot),
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) {
      # 每25个显示一个标签，其他为空
      ifelse(x %% 25 == 0, as.character(x), "")
    }
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  # 设置主题
  theme_minimal() +
  theme(
    plot.title = element_blank(),  # 移除标题
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # 旋转x轴标签
    panel.grid = element_blank(),  # 移除所有网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 添加坐标轴线
  )

# 如果特征向量太多，调整x轴标签密度
if (n_components_to_plot > 100) {
  # 对于大量特征向量，每50个显示一个标签
  cumulative_plot <- cumulative_plot +
    scale_x_continuous(
      breaks = seq(0, n_components_to_plot, by = 50),
      limits = c(0, n_components_to_plot),
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) {
        ifelse(x %% 50 == 0, as.character(x), "")
      }
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

#==================================================
# 组合图形
#==================================================

cat("\n正在组合图形...\n")

# 使用cowplot包组合图形
combined_plot <- plot_grid(
  pca_plot + theme(plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm")),  # 右边增加边距给图例
  cumulative_plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
  ncol = 2,
  labels = c("A", "B"),  # 使用大写字母标签
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1.2, 1)  # 调整两图宽度比例，给左边图更多空间
)

# 保存图形
output_file <- paste0(outfolder, "/Figure_PCA_analysis_final.png")
cat("正在保存图形到:", output_file, "\n")

png(
  output_file,
  res = 300,
  units = "in",
  width = 18,  # 增加宽度以适应图例和标签
  height = 6,
  bg = "white"
)
print(combined_plot)
dev.off()

cat("图形已保存到:", output_file, "\n")

# 保存为PDF
ggsave(
  filename = paste0(outfolder, "/Figure_PCA_analysis_final.pdf"),
  plot = combined_plot,
  width = 18,
  height = 6
)

# 显示图形
cat("\n正在显示图形...\n")
print(combined_plot)

# 保存详细结果到文本文件
cat("\n正在保存详细结果到文本文件...\n")
sink(paste0(outfolder, "/pca_analysis_report.txt"))
cat("=== PCA 分析报告 ===\n")
cat("分析日期:", as.character(Sys.Date()), "\n")
cat("分析时间:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("1. 数据概览\n")
cat("planting.csv 行数:", nrow(planting), " 列数:", ncol(planting), "\n")
cat("ECs_results.csv 行数:", nrow(ecs_results), " 列数:", ncol(ecs_results), "\n")
cat("环境因子总数:", length(env_factor_cols), "\n")
cat("处理后保留的环境因子数量:", length(cols_to_keep), "\n")
cat("移除的环境因子数量:", length(cols_to_remove), "\n")
cat("用于 PCA 的环境数量:", nrow(env_matrix_scaled), "\n")
cat("不同环境的数量:", length(unique(plot_data$env)), "\n\n")

cat("2. 地理分布\n")
cat("纬度范围:", min(plot_data$latitude, na.rm = TRUE), "到", max(plot_data$latitude, na.rm = TRUE), "\n")
cat("经度范围:", min(plot_data$longitude, na.rm = TRUE), "到", max(plot_data$longitude, na.rm = TRUE), "\n")
cat("州/地区数量:", length(unique(plot_data$state)), "\n\n")

cat("3. 数据预处理\n")
cat("移除的环境因子原因:\n")
cat("  有缺失值:", length(cols_with_na), "个\n")
cat("  缺失值比例超过50%:", length(high_missing_cols), "个\n")
cat("  标准差为零或接近零:", length(zero_sd_cols), "个\n")
cat("处理后的环境因子矩阵维度:", dim(env_matrix_processed), "\n\n")

cat("4. PCA 结果摘要\n")
cat("总方差:", round(total_variance, 4), "\n")
cat("主成分总数:", length(eigenvalues), "\n\n")

cat("前10个主成分的特征值和方差解释比例:\n")
cat("PC\t特征值\t方差比例(%)\t累计方差(%)\n")
for (i in 1:min(10, length(eigenvalues))) {
  cat(sprintf("%d\t%.4f\t%.2f\t\t%.2f\n", 
              i, eigenvalues[i], variance_proportion[i]*100, cumulative_variance[i]*100))
}
cat("\n")

cat("5. 坐标轴旋转结果\n")
cat("最优旋转角度:", rotation_result$rotation_angle, "度\n")
cat("旋转后 PC1 与经度的相关性:", round(rotation_result$correlations$longitude_cor, 4), "\n")
cat("旋转后 PC2 与纬度的相关性:", round(rotation_result$correlations$latitude_cor, 4), "\n\n")

cat("6. 方差阈值\n")
cat(sprintf("达到 %.0f%% 累计方差所需的主成分数量: %d\n", 95, threshold_95))
cat(sprintf("达到 %.0f%% 累计方差所需的主成分数量: %d\n", 50, threshold_50))
cat(sprintf("前 %d 个主成分解释的累计方差: %.2f%%\n", 
            threshold_95, cumulative_variance[threshold_95]*100))
cat(sprintf("前 %d 个主成分解释的累计方差: %.2f%%\n", 
            threshold_50, cumulative_variance[threshold_50]*100))
cat("\n")

cat("7. 图形输出\n")
cat("图形文件:", output_file, "\n")
cat("图形尺寸: 18 x 6 英寸\n")
cat("图形分辨率: 300 DPI\n")
cat("颜色方案: viridis::mako (避免黄色不明显)\n")
cat("标签显示: 每个环境点都显示env标签\n")
cat("\n")

cat("8. 处理的环境因子详情\n")
cat("保留的环境因子数量:", length(cols_to_keep), "\n")
cat("移除的环境因子数量:", length(cols_to_remove), "\n")
if (length(cols_to_remove) > 0) {
  cat("移除的环境因子列表（前20个）:\n")
  for (i in 1:min(20, length(cols_to_remove))) {
    cat(sprintf("  %d. %s\n", i, cols_to_remove[i]))
  }
  if (length(cols_to_remove) > 20) {
    cat(sprintf("  ... 还有 %d 个环境因子\n", length(cols_to_remove) - 20))
  }
}
sink()

cat("详细报告已保存到:", paste0(outfolder, "/pca_analysis_report.txt"), "\n")

# 保存用于绘图的中间数据
save(plot_data, rotation_result, scree_data, 
     env_matrix_processed, env_matrix_scaled,
     file = paste0(outfolder, "/plot_data.RData"))
cat("绘图数据已保存到:", paste0(outfolder, "/plot_data.RData"), "\n")

# 保存移除的环境因子列表
if (length(cols_to_remove) > 0) {
  removed_factors <- data.frame(
    Factor = cols_to_remove,
    Reason = ifelse(cols_to_remove %in% cols_with_na, "Has missing values",
             ifelse(cols_to_remove %in% high_missing_cols, 
                    "High missing values (>50%)",
                    "Zero or near-zero standard deviation"))
  )
  write.csv(removed_factors, 
            file = paste0(outfolder, "/removed_environmental_factors.csv"),
            row.names = FALSE)
  cat("移除的环境因子列表已保存到:", 
      paste0(outfolder, "/removed_environmental_factors.csv"), "\n")
}

# 保存保留的环境因子列表
if (length(cols_to_keep) > 0) {
  retained_factors <- data.frame(
    Factor = cols_to_keep,
    Standard_Deviation_Before_Scaling = col_sd[cols_to_keep]
  )
  write.csv(retained_factors, 
            file = paste0(outfolder, "/retained_environmental_factors.csv"),
            row.names = FALSE)
  cat("保留的环境因子列表已保存到:", 
      paste0(outfolder, "/retained_environmental_factors.csv"), "\n")
}

# 保存环境信息
env_info <- plot_data %>%
  select(env, latitude, longitude, state, PC1_rotated, PC2_rotated) %>%
  arrange(state, env)
write.csv(env_info, 
          file = paste0(outfolder, "/environment_info.csv"),
          row.names = FALSE)
cat("环境信息已保存到:", paste0(outfolder, "/environment_info.csv"), "\n")

cat("\n=== 分析完成 ===\n")
cat("所有输出文件已保存到 '", outfolder, "' 目录中\n", sep = "")
cat("保留的环境因子:", length(cols_to_keep), "/", length(env_factor_cols), "\n")
cat("分析的环境数量:", length(unique(plot_data$env)), "\n")
cat("图形已生成:", output_file, "\n")