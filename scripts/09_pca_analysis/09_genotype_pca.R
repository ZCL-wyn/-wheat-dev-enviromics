setwd("C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\01缺失的遗传力\\")

# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 加载必要的库
required_packages <- c("ggplot2", "dplyr", "tidyr", "factoextra", 
                      "cowplot", "viridis", "RColorBrewer", "data.table")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 读取基因型数据
geno_data <- read.csv("myGD2.csv", header = TRUE, row.names = 1)

# 查看数据结构和维度
cat("数据维度:", dim(geno_data), "\n")
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
explained_variance <- (G_EVD$values[index] / sum(G_EVD$values[index])) * 100
cumulative_variance <- cumsum(explained_variance)

# 创建解释方差数据框
variance_df <- data.frame(
  PC = 1:length(explained_variance),
  Variance = explained_variance,
  Cumulative = cumulative_variance
)

# 4. 创建主成分数据框用于绘图
pca_df <- as.data.frame(G_PC[, 1:3])  # 取前3个主成分
pca_df$Sample <- rownames(pca_df)

# 5. 绘制图形
cat("\n=== 开始绘制图形 ===\n")

# A图: PC1 vs PC2 散点图 - 根据PC1距离设置颜色
cat("绘制A图: PC1 vs PC2...\n")

# 计算PC1值的绝对距离（从中心0的距离）
pc1_values <- G_PC[, 1]
pc1_abs_distance <- abs(pc1_values)

# 创建PC1 vs PC2散点图 - 根据PC1的绝对距离设置颜色
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = pc1_abs_distance)) +
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
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
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

# B图: 累计贡献率图 - 修改部分
cat("绘制B图: 累计贡献率图...\n")

# 计算95%和50%阈值的位置
threshold_95 <- which(cumulative_variance >= 95)[1]
threshold_50 <- which(cumulative_variance >= 50)[1]

# 如果没有达到阈值，设为最后一个成分
if (is.na(threshold_95)) threshold_95 <- length(cumulative_variance)
if (is.na(threshold_50)) threshold_50 <- length(cumulative_variance)

# 确定要绘制的主成分数量：取前4000个或所有主成分（如果少于4000）
n_components_to_plot <- min(4000, length(cumulative_variance))

cat("将绘制前", n_components_to_plot, "个主成分的累计方差图\n")
cat("累计方差达到95%需要", threshold_95, "个主成分\n")
cat("累计方差达到50%需要", threshold_50, "个主成分\n")

# 将累计方差转换为比例（0-1）
variance_df_subset <- variance_df[1:n_components_to_plot, ]
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
    breaks = seq(0, n_components_to_plot, by = 50),
    limits = c(0, n_components_to_plot),
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
if (n_components_to_plot > 1000) {
  # 对于大量特征向量，每100个显示一个标签
  cumulative_plot <- cumulative_plot +
    scale_x_continuous(
      breaks = seq(0, n_components_to_plot, by = 100),
      limits = c(0, n_components_to_plot),
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) {
        ifelse(x %% 100 == 0, as.character(x), "")
      }
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

# 6. 组合图形
cat("\n=== 组合图形 ===\n")

# 使用cowplot包组合图形
combined_plot <- plot_grid(
  pca_plot + theme(plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm")),  # 右边增加边距给图例
  cumulative_plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
  ncol = 2,
  labels = c("a", "b"),  # 使用小写字母标签
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(0.9, 1.1)  # 调整两图宽度比例
)

# 7. 保存图形
cat("\n=== 保存图形 ===\n")

# 保存组合图形为PDF
ggsave(
  filename = "G_PCA_Analysis_Two_Panels_no_grid.pdf",
  plot = combined_plot,
  width = 16,
  height = 6,
  dpi = 300
)

# 保存组合图形为PNG
ggsave(
  filename = "G_PCA_Analysis_Two_Panels_no_grid.png",
  plot = combined_plot,
  width = 16,
  height = 6,
  dpi = 300,
  bg = "white"
)

# 保存单个图形
ggsave("G_PCA_PC1_vs_PC2_no_grid.pdf", pca_plot, width = 6, height = 6, dpi = 300)
ggsave("G_PCA_Cumulative_Variance_no_grid.pdf", cumulative_plot, width = 8, height = 6, dpi = 300)

# 8. 保存分析结果
cat("\n=== 保存分析结果 ===\n")

# 保存PCA结果
write.csv(pca_df, "G_PCA_Scores.csv", row.names = FALSE)
write.csv(variance_df, "G_PCA_Variance_Explained.csv", row.names = FALSE)

# 保存特征值和特征向量
eigen_results <- data.frame(
  Eigenvalue = G_EVD$values[index],
  Explained_Variance = explained_variance,
  Cumulative_Variance = cumulative_variance
)
write.csv(eigen_results, "G_Eigen_Results.csv", row.names = FALSE)

# 9. 输出摘要统计
cat("\n=== 分析摘要 ===\n")
cat("样本数量:", nrow(geno_matrix), "\n")
cat("SNP数量:", ncol(geno_matrix), "\n")
cat("G矩阵维度:", dim(G), "\n")
cat("\n基于G矩阵的PCA解释方差（前10个主成分）:\n")
for (i in 1:min(10, length(explained_variance))) {
  cat(sprintf("PC%d: %.2f%% (累计: %.2f%%)\n", 
              i, explained_variance[i], cumulative_variance[i]))
}

# 找出解释95%和50%方差所需的主成分数
cat(sprintf("\n解释95%%方差所需主成分数: %d\n", threshold_95))
cat(sprintf("解释50%%方差所需主成分数: %d\n", threshold_50))

# 如果累计方差最终没有达到95%，给出提示
if (max(cumulative_variance) < 95) {
  cat(sprintf("\n注意：所有主成分的累计方差仅为%.2f%%，未达到95%%\n", max(cumulative_variance)))
  cat("可能原因：数据维度高，特征值分布较均匀\n")
}

cat("\n分析完成！所有结果已保存到当前目录。\n")

# 10. 显示图形
print(combined_plot)