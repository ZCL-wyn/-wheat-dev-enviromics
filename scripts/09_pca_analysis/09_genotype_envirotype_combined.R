# 清理环境
rm(list = ls())

# 设置工作目录
setwd("C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\01缺失的遗传力\\")

# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 加载必要的库
required_packages <- c("ggplot2", "dplyr", "tidyr", "factoextra", 
                      "cowplot", "viridis", "RColorBrewer", "data.table",
                      "grid", "gridExtra", "matrixStats")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 创建输出目录
output_dir <- "Combined_PCA_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##############################################################################
# 第一部分：基因型PCA分析
##############################################################################

cat("=== 开始基因型PCA分析 ===\n")

# 读取基因型数据
geno_data <- read.csv("myGD2.csv", header = TRUE, row.names = 1)

# 查看数据结构和维度
cat("基因型数据维度:", dim(geno_data), "\n")
cat("样本数:", nrow(geno_data), "\n")
cat("SNP数:", ncol(geno_data), "\n")

# 转换为数值矩阵
geno_matrix <- as.matrix(geno_data)

# 处理缺失值 - 用列均值填补
for (j in 1:ncol(geno_matrix)) {
  if (any(is.na(geno_matrix[, j]))) {
    col_mean <- mean(geno_matrix[, j], na.rm = TRUE)
    if (!is.na(col_mean)) {
      geno_matrix[is.na(geno_matrix[, j]), j] <- col_mean
    } else {
      geno_matrix[is.na(geno_matrix[, j]), j] <- 0
    }
  }
}

# 检查是否有无限值并处理
if (any(is.infinite(geno_matrix))) {
  geno_matrix[is.infinite(geno_matrix)] <- NA
  for (j in 1:ncol(geno_matrix)) {
    if (any(is.na(geno_matrix[, j]))) {
      col_mean <- mean(geno_matrix[, j], na.rm = TRUE)
      if (!is.na(col_mean)) {
        geno_matrix[is.na(geno_matrix[, j]), j] <- col_mean
      } else {
        geno_matrix[is.na(geno_matrix[, j]), j] <- 0
      }
    }
  }
}

# 1. 计算G矩阵（基因组关系矩阵）
cat("\n=== 计算G矩阵 ===\n")

# 中心化（不减均值）
X <- scale(geno_matrix, center = TRUE, scale = FALSE)

# 计算G矩阵：G = XX^T
G <- tcrossprod(X)

# 标准化G矩阵：除以对角线元素的均值
diag_mean <- mean(diag(G))
G <- G / diag_mean

cat("G矩阵维度:", dim(G), "\n")

# 2. 对G矩阵进行特征值分解
cat("\n=== 对G矩阵进行特征值分解 ===\n")
G_EVD <- eigen(G)

# 查看特征值
cat("特征值前10个:", round(G_EVD$values[1:10], 4), "\n")

# 设置行名
rownames(G_EVD$vectors) <- rownames(G)

# 选择显著的特征值（大于1E-8）
index <- which(G_EVD$values > 1e-8)
cat("显著特征值数量:", length(index), "\n")

# 计算主成分：特征向量乘以特征值的平方根
G_PC <- sweep(G_EVD$vectors[, index], 2, sqrt(G_EVD$values[index]), FUN = '*')
colnames(G_PC) <- paste0("PC", 1:ncol(G_PC))

# 3. 计算解释方差
cat("\n=== 计算解释方差 ===\n")
explained_variance_geno <- (G_EVD$values[index] / sum(G_EVD$values[index])) * 100
cumulative_variance_geno <- cumsum(explained_variance_geno)

# 创建解释方差数据框
variance_df_geno <- data.frame(
  PC = 1:length(explained_variance_geno),
  Variance = explained_variance_geno,
  Cumulative = cumulative_variance_geno
)

# 4. 创建主成分数据框用于绘图
pca_df_geno <- as.data.frame(G_PC[, 1:3])  # 取前3个主成分
pca_df_geno$Sample <- rownames(pca_df_geno)

# 5. 计算PC1值的绝对距离（从中心0的距离）
pc1_values_geno <- G_PC[, 1]
pc1_abs_distance_geno <- abs(pc1_values_geno)

# 6. 绘制基因型PCA的A图: PC1 vs PC2 散点图
cat("\n绘制基因型PCA的A图: PC1 vs PC2...\n")

geno_pca_plot <- ggplot(pca_df_geno, aes(x = PC1, y = PC2, color = pc1_abs_distance_geno)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("blue", "green", "yellow", "red"),
    name = "",
    guide = guide_colorbar(
      title = NULL,
      barwidth = unit(0.5, "cm"),
      barheight = unit(2, "cm"),
      ticks = FALSE
    )
  ) +
  labs(
    x = paste0("PC1 (", round(explained_variance_geno[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance_geno[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.92, 0.15),  # 将图例放在图形内部右上角
    legend.title = element_blank(),
    legend.text = element_blank(),
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
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25))

# 7. 绘制基因型PCA的B图: 累计贡献率图
cat("\n绘制基因型PCA的B图: 累计贡献率图...\n")

# 计算95%和50%阈值的位置
threshold_95_geno <- which(cumulative_variance_geno >= 95)[1]
threshold_50_geno <- which(cumulative_variance_geno >= 50)[1]

# 如果没有达到阈值，设为最后一个成分
if (is.na(threshold_95_geno)) threshold_95_geno <- length(cumulative_variance_geno)
if (is.na(threshold_50_geno)) threshold_50_geno <- length(cumulative_variance_geno)

# 确定要绘制的主成分数量：取前4000个或所有主成分（如果少于4000）
n_components_to_plot_geno <- min(4000, length(cumulative_variance_geno))

cat("将绘制前", n_components_to_plot_geno, "个主成分的累计方差图\n")
cat("累计方差达到95%需要", threshold_95_geno, "个主成分\n")
cat("累计方差达到50%需要", threshold_50_geno, "个主成分\n")

# 将累计方差转换为比例（0-1）
variance_df_subset_geno <- variance_df_geno[1:n_components_to_plot_geno, ]
variance_df_subset_geno$Cumulative_Proportion <- variance_df_subset_geno$Cumulative / 100

# 获取50%和95%阈值线的y值
y_50_geno <- cumulative_variance_geno[threshold_50_geno] / 100
y_95_geno <- cumulative_variance_geno[threshold_95_geno] / 100

# 创建累计贡献率图
geno_cumulative_plot <- ggplot(variance_df_subset_geno, 
                               aes(x = PC, y = Cumulative_Proportion)) +
  # 绘制折线
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange", size = 1) +
  
  # 添加95%阈值线和标注
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = threshold_95_geno, linetype = "dashed", color = "blue", linewidth = 0.8) +
  
  # 添加50%阈值线和标注
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = threshold_50_geno, linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 在交点处添加标签
  annotate("point", 
           x = threshold_95_geno, 
           y = y_95_geno, 
           color = "blue", 
           size = 3) +
  annotate("text", 
           x = threshold_95_geno, 
           y = y_95_geno + 0.03, 
           label = paste0(threshold_95_geno, " (95%)"), 
           color = "blue", 
           size = 4,
           hjust = 0.5,
           vjust = 0) +
  
  annotate("point", 
           x = threshold_50_geno, 
           y = y_50_geno, 
           color = "red", 
           size = 3) +
  annotate("text", 
           x = threshold_50_geno, 
           y = y_50_geno + 0.03, 
           label = paste0(threshold_50_geno, " (50%)"), 
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
    breaks = seq(0, n_components_to_plot_geno, by = 50),
    limits = c(0, n_components_to_plot_geno),
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) {
      # 每50个显示一个标签，其他为空
      ifelse(x %% 50 == 0, as.character(x), "")
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
if (n_components_to_plot_geno > 1000) {
  # 对于大量特征向量，每100个显示一个标签
  geno_cumulative_plot <- geno_cumulative_plot +
    scale_x_continuous(
      breaks = seq(0, n_components_to_plot_geno, by = 100),
      limits = c(0, n_components_to_plot_geno),
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) {
        ifelse(x %% 100 == 0, as.character(x), "")
      }
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

cat("基因型PCA分析完成！\n\n")

##############################################################################
# 第二部分：环境PCA分析
##############################################################################

cat("=== 开始环境PCA分析 ===\n")

# 读取数据
planting <- read.csv("planting.csv")
ecs_results <- read.csv("ECs_results.csv")

# 查看数据结构
cat("planting.csv 维度:", dim(planting), "\n")
cat("planting.csv 列名:", paste(colnames(planting), collapse = ", "), "\n")
cat("\n")
cat("ECs_results.csv 维度:", dim(ecs_results), "\n")
cat("ECs_results.csv 列名:", paste(colnames(ecs_results), collapse = ", "), "\n")
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

# 准备环境因子数据用于 PCA
# 获取环境因子列（排除 env 列）
env_factor_cols <- setdiff(colnames(ecs_results), "env")
cat("环境因子总数:", length(env_factor_cols), "\n")

# 创建环境因子矩阵
env_matrix <- as.matrix(ecs_results[, env_factor_cols])
rownames(env_matrix) <- ecs_results$env

# 分析环境因子矩阵
cat("=== 环境因子矩阵分析 ===\n")
cat("原始环境因子矩阵维度:", dim(env_matrix), "\n")

# 1. 检查每列的缺失值
missing_counts <- colSums(is.na(env_matrix))
cols_with_na <- names(missing_counts[missing_counts > 0])
cat("有缺失值的环境因子数量:", length(cols_with_na), "\n")

# 2. 检查每列的标准差（忽略缺失值）
col_sd <- apply(env_matrix, 2, function(x) {
  x_non_missing <- x[!is.na(x)]
  if (length(x_non_missing) < 2) return(0)
  sd(x_non_missing, na.rm = TRUE)
})

# 找出标准差为0或接近0的列
zero_sd_cols <- names(col_sd[col_sd < 1e-10])
cat("标准差为零或接近零的环境因子数量:", length(zero_sd_cols), "\n")

# 3. 确定要移除的环境因子
cols_to_remove <- unique(c(cols_with_na, zero_sd_cols))
cat("总共要移除的环境因子数量:", length(cols_to_remove), "\n")

# 4. 创建处理后的环境因子矩阵（移除有问题的列）
cols_to_keep <- setdiff(env_factor_cols, cols_to_remove)
cat("要保留的环境因子数量:", length(cols_to_keep), "\n")

if (length(cols_to_keep) == 0) {
  stop("错误：没有环境因子可以用于PCA分析！请检查数据。")
}

env_matrix_processed <- env_matrix[, cols_to_keep, drop = FALSE]
cat("处理后的环境因子矩阵维度:", dim(env_matrix_processed), "\n")

# 5. 再次检查处理后的矩阵是否有缺失值
if (sum(is.na(env_matrix_processed)) > 0) {
  # 移除这些列
  na_by_col_processed <- colSums(is.na(env_matrix_processed))
  na_cols_processed <- names(na_by_col_processed[na_by_col_processed > 0])
  cols_to_keep <- setdiff(cols_to_keep, na_cols_processed)
  cols_to_remove <- c(cols_to_remove, na_cols_processed)
  env_matrix_processed <- env_matrix_processed[, cols_to_keep, drop = FALSE]
  cat("再次处理后矩阵维度:", dim(env_matrix_processed), "\n")
}

# 6. 检查并处理无穷值
if (sum(is.infinite(env_matrix_processed)) > 0) {
  inf_cols <- apply(env_matrix_processed, 2, function(x) any(is.infinite(x)))
  if (any(inf_cols)) {
    env_matrix_processed <- env_matrix_processed[, !inf_cols, drop = FALSE]
    cols_to_keep <- cols_to_keep[!inf_cols]
    cat("移除无穷值列后矩阵维度:", dim(env_matrix_processed), "\n")
  }
}

# 7. 标准化环境因子数据（PCA前通常需要标准化）
cat("正在标准化环境因子数据...\n")
env_matrix_scaled <- scale(env_matrix_processed, center = TRUE, scale = TRUE)

# 8. 检查标准化后的矩阵是否有缺失值或无穷值
if (sum(is.na(env_matrix_scaled)) > 0) {
  # 找出包含缺失值的列
  na_cols_scaled <- apply(env_matrix_scaled, 2, function(x) any(is.na(x)))
  if (any(na_cols_scaled)) {
    env_matrix_scaled <- env_matrix_scaled[, !na_cols_scaled, drop = FALSE]
    cols_to_keep <- cols_to_keep[!na_cols_scaled]
    cols_to_remove <- c(cols_to_remove, colnames(env_matrix_scaled)[na_cols_scaled])
  }
}

if (sum(is.infinite(env_matrix_scaled)) > 0) {
  # 找出包含无穷值的列
  inf_cols_scaled <- apply(env_matrix_scaled, 2, function(x) any(is.infinite(x)))
  if (any(inf_cols_scaled)) {
    env_matrix_scaled <- env_matrix_scaled[, !inf_cols_scaled, drop = FALSE]
    cols_to_keep <- cols_to_keep[!inf_cols_scaled]
  }
}

# 9. 最终检查矩阵完整性
if (any(is.na(env_matrix_scaled)) || any(is.infinite(env_matrix_scaled))) {
  stop("矩阵完整性检查未通过，无法进行PCA分析，请检查数据。")
}

cat("矩阵完整性检查通过！\n")

# 进行 PCA 分析
cat("正在进行 PCA 分析...\n")
pca_result <- prcomp(env_matrix_scaled, center = FALSE, scale. = FALSE)

# 提取 PCA 结果
eigenvalues_env <- pca_result$sdev^2
total_variance_env <- sum(eigenvalues_env)
variance_proportion_env <- eigenvalues_env / total_variance_env
cumulative_variance_env <- cumsum(variance_proportion_env)

# 获取主成分得分
pc_scores_env <- pca_result$x

# 打印 PCA 摘要
cat("\n=== PCA 分析结果摘要 ===\n")
cat("总方差:", round(total_variance_env, 2), "\n")
cat("主成分总数:", length(eigenvalues_env), "\n")
cat("前10个主成分解释的方差比例:\n")
for (i in 1:min(10, length(variance_proportion_env))) {
  cat(sprintf("  PC%d: %.3f (%.2f%%, 累计: %.2f%%)\n", 
              i, eigenvalues_env[i], variance_proportion_env[i] * 100, cumulative_variance_env[i] * 100))
}

# 创建 rotate_PC 函数
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
  angles <- seq(0, 360, by = 1)
  best_correlation <- -Inf
  best_angle <- 0
  best_PC_rotated <- NULL
  
  for (angle in angles) {
    theta <- angle * pi / 180
    rotation_matrix <- matrix(c(cos(theta), -sin(theta), 
                                sin(theta), cos(theta)), 
                              nrow = 2, byrow = TRUE)
    
    PC_rotated <- PC_adjusted %*% rotation_matrix
    
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
rownames(pc_scores_env) <- ecs_results$env
pc_scores_filtered <- pc_scores_env[env_locations$env, ]

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
plot_data_env <- data.frame(
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

# 计算每个 env 的显示标签
plot_data_env$env_display <- plot_data_env$env

# 计算方差解释比例
var_pc1_env <- variance_proportion_env[1] * 100
var_pc2_env <- variance_proportion_env[2] * 100

# 创建轴标签（包含方差解释比例）
x_label_env <- paste0("PC1 (", round(var_pc1_env, 1), "%)")
y_label_env <- paste0("PC2 (", round(var_pc2_env, 1), "%)")

# 创建旋转轴箭头数据
pc1_range <- range(plot_data_env$PC1_rotated, na.rm = TRUE)
pc2_range <- range(plot_data_env$PC2_rotated, na.rm = TRUE)

axis_length_pc1 <- diff(pc1_range) * 0.7 / 2
axis_length_pc2 <- diff(pc2_range) * 0.7 / 2
axis_length_env <- max(axis_length_pc1, axis_length_pc2)

# PC1轴：正方向和反方向
x_axis_arrow_forward_env <- data.frame(
  x = 0, y = 0,
  xend = rotation_result$x_axis[1] * axis_length_env,
  yend = rotation_result$x_axis[2] * axis_length_env
)

x_axis_arrow_backward_env <- data.frame(
  x = 0, y = 0,
  xend = -rotation_result$x_axis[1] * axis_length_env,
  yend = -rotation_result$x_axis[2] * axis_length_env
)

# PC2轴：正方向和反方向
y_axis_arrow_forward_env <- data.frame(
  x = 0, y = 0,
  xend = rotation_result$y_axis[1] * axis_length_env,
  yend = rotation_result$y_axis[2] * axis_length_env
)

y_axis_arrow_backward_env <- data.frame(
  x = 0, y = 0,
  xend = -rotation_result$y_axis[1] * axis_length_env,
  yend = -rotation_result$y_axis[2] * axis_length_env
)

# 计算纬度的范围
lat_range_env <- range(plot_data_env$latitude, na.rm = TRUE)
cat("纬度范围: ", lat_range_env[1], " 到 ", lat_range_env[2], "\n")

# 创建环境PCA的A图: 旋转后的 PC1 vs PC2 图
cat("\n正在创建环境PCA的A图: 旋转后的 PC1 vs PC2 图...\n")

env_pca_plot <- ggplot(plot_data_env, aes(x = PC1_rotated, y = PC2_rotated)) +
  # 添加PC1轴的反向延长线（无箭头）
  geom_segment(data = x_axis_arrow_backward_env,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "red", size = 0.8, alpha = 0.7, linetype = "solid") +
  # 添加PC1轴的正向箭头
  geom_segment(data = x_axis_arrow_forward_env, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red", size = 0.8, alpha = 0.7) +
  # 添加PC2轴的反向延长线（无箭头）
  geom_segment(data = y_axis_arrow_backward_env,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "blue", size = 0.8, alpha = 0.7, linetype = "solid") +
  # 添加PC2轴的正向箭头
  geom_segment(data = y_axis_arrow_forward_env,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "blue", size = 0.8, alpha = 0.7) +
  # 添加数据点，使用渐变色表示纬度
  geom_point(aes(fill = latitude), 
             size = 3, shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  # 添加env标签在点上
  geom_text(aes(label = env), size = 2.5, color = "black", vjust = 1.5) +
  # 添加坐标轴标签
  annotate("text", 
           x = x_axis_arrow_forward_env$xend * 1.1, 
           y = x_axis_arrow_forward_env$yend * 1.1,
           label = "PC1", color = "red", size = 4.5, fontface = "bold") +
  annotate("text",
           x = y_axis_arrow_forward_env$xend * 1.1,
           y = y_axis_arrow_forward_env$yend * 1.1,
           label = "PC2", color = "blue", size = 4.5, fontface = "bold") +
  # 设置坐标轴
  labs(
    x = x_label_env,
    y = y_label_env,
    fill = "Latitude"
  ) +
  # 应用渐变色标尺 - 使用viridis的mako选项
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    breaks = seq(floor(lat_range_env[1]), ceiling(lat_range_env[2]), by = 2),
    guide = guide_colorbar(
      title = "Latitude",
      title.position = "top",
      barwidth = unit(0.5, "cm"),
      barheight = unit(2, "cm"),
      ticks = FALSE
    )
  ) +
  # 设置主题 - 与基因型PCA一致
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
      min(plot_data_env$PC1_rotated, -axis_length_env * 1.2, na.rm = TRUE),
      max(plot_data_env$PC1_rotated, axis_length_env * 1.2, na.rm = TRUE)
    )
  ) +
  scale_y_continuous(
    limits = c(
      min(plot_data_env$PC2_rotated, -axis_length_env * 1.2, na.rm = TRUE),
      max(plot_data_env$PC2_rotated, axis_length_env * 1.2, na.rm = TRUE)
    )
  )

# 创建环境PCA的B图: 累计贡献率图（风格与基因型PCA的B图完全一致）
cat("\n正在创建环境PCA的B图: 累计贡献率图...\n")

# 创建 scree plot 数据
scree_data_env <- data.frame(
  PC = 1:length(variance_proportion_env),
  Eigenvalue = eigenvalues_env,
  Variance = variance_proportion_env * 100,
  Cumulative = cumulative_variance_env * 100
)

# 计算95%和50%阈值的位置
threshold_95_env <- which(cumulative_variance_env * 100 >= 95)[1]
threshold_50_env <- which(cumulative_variance_env * 100 >= 50)[1]

# 如果没有达到阈值，设为最后一个成分
if (is.na(threshold_95_env)) threshold_95_env <- length(cumulative_variance_env)
if (is.na(threshold_50_env)) threshold_50_env <- length(cumulative_variance_env)

# 确定要绘制的主成分数量：取前200个或所有主成分（如果少于200）
n_components_to_plot_env <- min(200, length(cumulative_variance_env))

cat("将绘制前", n_components_to_plot_env, "个主成分的累计方差图\n")
cat("累计方差达到95%需要", threshold_95_env, "个主成分\n")
cat("累计方差达到50%需要", threshold_50_env, "个主成分\n")

# 将累计方差转换为比例（0-1）
variance_df_subset_env <- scree_data_env[1:n_components_to_plot_env, ]
variance_df_subset_env$Cumulative_Proportion <- variance_df_subset_env$Cumulative / 100

# 获取50%和95%阈值线的y值
y_50_env <- cumulative_variance_env[threshold_50_env] / 100
y_95_env <- cumulative_variance_env[threshold_95_env] / 100

# 创建累计贡献率图 - 风格与基因型PCA完全一致
env_cumulative_plot <- ggplot(variance_df_subset_env, 
                              aes(x = PC, y = Cumulative_Proportion)) +
  # 绘制折线 - 与基因型PCA相同的颜色
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange", size = 1) +
  
  # 添加95%阈值线和标注 - 与基因型PCA相同
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = threshold_95_env, linetype = "dashed", color = "blue", linewidth = 0.8) +
  
  # 添加50%阈值线和标注 - 与基因型PCA相同
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = threshold_50_env, linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 在交点处添加标签 - 与基因型PCA相同
  annotate("point", 
           x = threshold_95_env, 
           y = y_95_env, 
           color = "blue", 
           size = 3) +
  annotate("text", 
           x = threshold_95_env, 
           y = y_95_env + 0.03, 
           label = paste0(threshold_95_env, " (95%)"), 
           color = "blue", 
           size = 4,
           hjust = 0.5,
           vjust = 0) +
  
  annotate("point", 
           x = threshold_50_env, 
           y = y_50_env, 
           color = "red", 
           size = 3) +
  annotate("text", 
           x = threshold_50_env, 
           y = y_50_env + 0.03, 
           label = paste0(threshold_50_env, " (50%)"), 
           color = "red", 
           size = 4,
           hjust = 0.5,
           vjust = 0) +
  
  # 设置坐标轴和标题 - 与基因型PCA相同
  labs(
    x = "Number of eigenvectors",
    y = "Proportion of variance"
  ) +
  scale_x_continuous(
    breaks = seq(0, n_components_to_plot_env, by = 25),
    limits = c(0, n_components_to_plot_env),
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
  
  # 设置主题 - 与基因型PCA完全相同
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
if (n_components_to_plot_env > 100) {
  # 对于大量特征向量，每50个显示一个标签
  env_cumulative_plot <- env_cumulative_plot +
    scale_x_continuous(
      breaks = seq(0, n_components_to_plot_env, by = 50),
      limits = c(0, n_components_to_plot_env),
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) {
        ifelse(x %% 50 == 0, as.character(x), "")
      }
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

cat("环境PCA分析完成！\n\n")

##############################################################################
# 第三部分：组合四个图形
##############################################################################

cat("=== 开始组合四个图形 ===\n")

# 调整图形边距
geno_pca_plot_adj <- geno_pca_plot + 
  theme(plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"))

geno_cumulative_plot_adj <- geno_cumulative_plot + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

env_pca_plot_adj <- env_pca_plot + 
  theme(plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"))

env_cumulative_plot_adj <- env_cumulative_plot + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# 创建第一行（基因型PCA）
row1 <- plot_grid(
  geno_pca_plot_adj,
  geno_cumulative_plot_adj,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  rel_widths = c(0.9, 1.1)
)

# 创建第二行（环境PCA）
row2 <- plot_grid(
  env_pca_plot_adj,
  env_cumulative_plot_adj,
  ncol = 2,
  labels = c("C", "D"),
  label_size = 18,
  label_fontface = "bold",
  rel_widths = c(0.9, 1.1)
)

# 组合两行
combined_four_panel <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(1, 1)
)

# 保存高分辨率图形
cat("\n=== 保存高分辨率图形 ===\n")

# 1. 保存为600ppi的JPG
jpeg_filename <- file.path(output_dir, "Four_Panel_PCA_Analysis_600dpi.jpg")
cat("保存为600ppi JPG:", jpeg_filename, "\n")

jpeg(
  filename = jpeg_filename,
  width = 16,
  height = 12,
  units = "in",
  res = 600,  # 600 DPI
  quality = 100  # 最高质量
)
print(combined_four_panel)
dev.off()

# 2. 保存为PDF
pdf_filename <- file.path(output_dir, "Four_Panel_PCA_Analysis.pdf")
cat("保存为PDF:", pdf_filename, "\n")

ggsave(
  filename = pdf_filename,
  plot = combined_four_panel,
  width = 16,
  height = 12,
  device = "pdf",
  dpi = 600
)

# 3. 保存为高分辨率PNG
png_filename <- file.path(output_dir, "Four_Panel_PCA_Analysis_600dpi.png")
cat("保存为600ppi PNG:", png_filename, "\n")

png(
  filename = png_filename,
  width = 16,
  height = 12,
  units = "in",
  res = 600,
  bg = "white"
)
print(combined_four_panel)
dev.off()

# 4. 保存为TIFF格式（适合出版）
tiff_filename <- file.path(output_dir, "Four_Panel_PCA_Analysis_600dpi.tiff")
cat("保存为600ppi TIFF:", tiff_filename, "\n")

tiff(
  filename = tiff_filename,
  width = 16,
  height = 12,
  units = "in",
  res = 600,
  compression = "lzw",  # 无损压缩
  bg = "white"
)
print(combined_four_panel)
dev.off()

# 5. 保存为EPS格式（矢量图）
eps_filename <- file.path(output_dir, "Four_Panel_PCA_Analysis.eps")
cat("保存为EPS:", eps_filename, "\n")

ggsave(
  filename = eps_filename,
  plot = combined_four_panel,
  width = 16,
  height = 12,
  device = "eps",
  dpi = 600
)

# 保存单个图形
cat("\n保存单个图形...\n")

# 基因型PCA单个图形
ggsave(file.path(output_dir, "Genotype_PCA_A.jpg"), 
       geno_pca_plot, width = 8, height = 6, dpi = 600)
ggsave(file.path(output_dir, "Genotype_PCA_B.jpg"), 
       geno_cumulative_plot, width = 8, height = 6, dpi = 600)

# 环境PCA单个图形
ggsave(file.path(output_dir, "Environment_PCA_C.jpg"), 
       env_pca_plot, width = 8, height = 6, dpi = 600)
ggsave(file.path(output_dir, "Environment_PCA_D.jpg"), 
       env_cumulative_plot, width = 8, height = 6, dpi = 600)

##############################################################################
# 第四部分：保存分析结果
##############################################################################

cat("\n=== 保存分析结果 ===\n")

# 1. 保存基因型PCA结果
write.csv(pca_df_geno, file.path(output_dir, "Genotype_PCA_Scores.csv"), row.names = FALSE)
write.csv(variance_df_geno, file.path(output_dir, "Genotype_PCA_Variance_Explained.csv"), row.names = FALSE)

# 保存特征值和特征向量
eigen_results_geno <- data.frame(
  Eigenvalue = G_EVD$values[index],
  Explained_Variance = explained_variance_geno,
  Cumulative_Variance = cumulative_variance_geno
)
write.csv(eigen_results_geno, file.path(output_dir, "Genotype_Eigen_Results.csv"), row.names = FALSE)

# 2. 保存环境PCA结果
write.csv(plot_data_env, file.path(output_dir, "Environment_PCA_Scores.csv"), row.names = FALSE)
write.csv(scree_data_env, file.path(output_dir, "Environment_PCA_Variance_Explained.csv"), row.names = FALSE)

# 保存环境信息
env_info <- plot_data_env %>%
  select(env, latitude, longitude, state, PC1_rotated, PC2_rotated) %>%
  arrange(state, env)
write.csv(env_info, file.path(output_dir, "Environment_Info.csv"), row.names = FALSE)

# 保存移除的环境因子列表
if (length(cols_to_remove) > 0) {
  removed_factors <- data.frame(
    Factor = cols_to_remove,
    Reason = ifelse(cols_to_remove %in% cols_with_na, "Has missing values",
             ifelse(cols_to_remove %in% zero_sd_cols, 
                    "Zero or near-zero standard deviation",
                    "Other reasons"))
  )
  write.csv(removed_factors, 
            file = file.path(output_dir, "Removed_Environmental_Factors.csv"),
            row.names = FALSE)
}

# 保存保留的环境因子列表
if (length(cols_to_keep) > 0) {
  retained_factors <- data.frame(
    Factor = cols_to_keep,
    Standard_Deviation_Before_Scaling = col_sd[cols_to_keep]
  )
  write.csv(retained_factors, 
            file = file.path(output_dir, "Retained_Environmental_Factors.csv"),
            row.names = FALSE)
}

##############################################################################
# 第五部分：输出摘要统计
##############################################################################

cat("\n=== 分析摘要 ===\n")

# 创建详细的报告文件
sink(file.path(output_dir, "PCA_Analysis_Report.txt"))

cat("========================================\n")
cat("           PCA 分析报告\n")
cat("========================================\n")
cat("分析日期:", as.character(Sys.Date()), "\n")
cat("分析时间:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("第一部分：基因型PCA分析\n")
cat("------------------------\n")
cat("样本数量:", nrow(geno_matrix), "\n")
cat("SNP数量:", ncol(geno_matrix), "\n")
cat("G矩阵维度:", dim(G), "\n")
cat("显著特征值数量:", length(index), "\n\n")

cat("基于G矩阵的PCA解释方差（前10个主成分）:\n")
for (i in 1:min(10, length(explained_variance_geno))) {
  cat(sprintf("PC%d: %.2f%% (累计: %.2f%%)\n", 
              i, explained_variance_geno[i], cumulative_variance_geno[i]))
}

cat(sprintf("\n解释95%%方差所需主成分数: %d\n", threshold_95_geno))
cat(sprintf("解释50%%方差所需主成分数: %d\n", threshold_50_geno))

# 如果累计方差最终没有达到95%，给出提示
if (max(cumulative_variance_geno) < 95) {
  cat(sprintf("\n注意：所有主成分的累计方差仅为%.2f%%，未达到95%%\n", max(cumulative_variance_geno)))
  cat("可能原因：数据维度高，特征值分布较均匀\n")
}

cat("\n\n第二部分：环境PCA分析\n")
cat("------------------------\n")
cat("环境数量:", nrow(env_matrix_scaled), "\n")
cat("保留的环境因子数量:", length(cols_to_keep), "\n")
cat("移除的环境因子数量:", length(cols_to_remove), "\n")
cat("总方差:", round(total_variance_env, 2), "\n")
cat("主成分总数:", length(eigenvalues_env), "\n\n")

cat("前10个主成分解释的方差比例:\n")
for (i in 1:min(10, length(variance_proportion_env))) {
  cat(sprintf("PC%d: %.3f (%.2f%%, 累计: %.2f%%)\n", 
              i, eigenvalues_env[i], variance_proportion_env[i] * 100, cumulative_variance_env[i] * 100))
}

cat("\n坐标轴旋转结果:\n")
cat("最优旋转角度:", rotation_result$rotation_angle, "度\n")
cat("旋转后 PC1 与经度的相关性:", round(rotation_result$correlations$longitude_cor, 4), "\n")
cat("旋转后 PC2 与纬度的相关性:", round(rotation_result$correlations$latitude_cor, 4), "\n")

cat(sprintf("\n解释95%%方差所需主成分数: %d\n", threshold_95_env))
cat(sprintf("解释50%%方差所需主成分数: %d\n", threshold_50_env))

cat("\n\n第三部分：图形输出\n")
cat("------------------------\n")
cat("四面板组合图形已保存为以下格式:\n")
cat("1. 600ppi JPG:", jpeg_filename, "\n")
cat("2. PDF:", pdf_filename, "\n")
cat("3. 600ppi PNG:", png_filename, "\n")
cat("4. 600ppi TIFF:", tiff_filename, "\n")
cat("5. EPS矢量图:", eps_filename, "\n")

cat("\n图形尺寸: 16 x 12 英寸\n")
cat("分辨率: 600 DPI\n")
cat("颜色方案:\n")
cat("  - 基因型PCA: 蓝-绿-黄-红渐变色\n")
cat("  - 环境PCA: viridis::mako色系\n")
cat("  - 累计方差图: 橙色折线，蓝红阈值线\n")

cat("\n\n第四部分：文件保存\n")
cat("------------------------\n")
cat("所有输出文件已保存到目录:", output_dir, "\n")
cat("包括:\n")
cat("  1. 图形文件（多种格式）\n")
cat("  2. PCA得分和方差解释文件\n")
cat("  3. 特征值和特征向量结果\n")
cat("  4. 环境因子信息\n")
cat("  5. 分析报告\n")

cat("\n========================================\n")
cat("          分析完成\n")
cat("========================================\n")

sink()

# 在控制台显示摘要
cat("\n=== 最终摘要 ===\n")
cat("分析完成！所有结果已保存到目录:", output_dir, "\n")
cat("生成的图形:\n")
cat("1. 四面板组合图 (A, B, C, D)\n")
cat("2. 单个图形 (A, B, C, D)\n")
cat("\n图形规格:\n")
cat("- 尺寸: 16 x 12 英寸\n")
cat("- 分辨率: 600 DPI\n")
cat("- 格式: JPG, PDF, PNG, TIFF, EPS\n")
cat("\n数据规格:\n")
cat("- 基因型样本:", nrow(geno_matrix), "\n")
cat("- 环境样本:", nrow(env_matrix_scaled), "\n")
cat("- 基因型SNP数:", ncol(geno_matrix), "\n")
cat("- 环境因子数:", length(cols_to_keep), "（保留）\n")
cat("\n主成分分析结果:\n")
cat("- 基因型PCA: PC1解释", round(explained_variance_geno[1], 1), "% 方差\n")
cat("- 环境PCA: PC1解释", round(var_pc1_env, 1), "% 方差\n")

# 显示图形
cat("\n正在显示组合图形...\n")
print(combined_four_panel)

cat("\n=== 全部分析完成！ ===\n")