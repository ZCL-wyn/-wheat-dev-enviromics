# 1. 加载必要的包
library(dplyr)
library(tidyr)
library(moments)

# 2. 读取数据
data <- read.table("TKW.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 查看数据结构
cat("数据结构:\n")
str(data)
cat("\n前6行数据:\n")
head(data)

# 3. 计算各品种各环境的初始平均值
initial_means <- data %>%
  group_by(ID, ENV) %>%
  summarise(
    TKW_Mean = mean(TKW, na.rm = TRUE),
    TKW_SD = sd(TKW, na.rm = TRUE),
    REP_Count = n(),
    .groups = 'drop'
  )

cat("\n=== 初始平均值统计 ===\n")
print(head(initial_means))

# 4. 识别和删除离群值（使用IQR方法，更稳健）
# 计算每个环境的四分位距
env_iqr_stats <- data %>%
  group_by(ENV) %>%
  summarise(
    Q1 = quantile(TKW, 0.25, na.rm = TRUE),
    Q3 = quantile(TKW, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower_Bound = Q1 - 1.5 * IQR,
    Upper_Bound = Q3 + 1.5 * IQR,
    .groups = 'drop'
  )

# 合并环境统计信息到原始数据
data_with_bounds <- data %>%
  left_join(env_iqr_stats, by = "ENV")

# 标记离群值（超过1.5倍IQR范围）
data_clean <- data_with_bounds %>%
  mutate(
    is_outlier = TKW < Lower_Bound | TKW > Upper_Bound
  ) %>%
  filter(!is_outlier) %>%
  select(ID, ENV, REP, TKW)

cat("\n=== 离群值处理结果 ===\n")
cat("原始数据行数:", nrow(data), "\n")
cat("清洗后数据行数:", nrow(data_clean), "\n")
cat("删除的离群值数量:", nrow(data) - nrow(data_clean), "\n")
cat("离群值比例:", round((nrow(data) - nrow(data_clean))/nrow(data) * 100, 2), "%\n")

# 5. 使用清洗后的数据重新计算各品种各环境的平均值
clean_means <- data_clean %>%
  group_by(ID, ENV) %>%
  summarise(
    TKW_Mean = mean(TKW, na.rm = TRUE),
    TKW_SD = sd(TKW, na.rm = TRUE),
    TKW_CV = sd(TKW, na.rm = TRUE) / mean(TKW, na.rm = TRUE) * 100,
    REP_Count = n(),
    .groups = 'drop'
  )

cat("\n=== 清洗后各品种各环境平均值（完整表） ===\n")
print(clean_means)

# 6. 输出平均值到文件（ID, ENV, TKW格式）
# 创建平均值表
mean_table <- clean_means %>%
  select(ID, ENV, TKW_Mean) %>%
  rename(TKW = TKW_Mean)  # 重命名为TKW以符合要求

# 保存平均值到CSV文件
write.csv(mean_table, "TKW_mean_values.csv", row.names = FALSE)

cat("\n已保存平均值到文件: TKW_mean_values.csv\n")
cat("文件包含列: ID, ENV, TKW\n")

# 7. 将数据转换为宽格式（用于计算环境间相关系数）
wide_data <- clean_means %>%
  select(ID, ENV, TKW_Mean) %>%
  pivot_wider(
    names_from = ENV,
    values_from = TKW_Mean,
    values_fill = NA
  )

cat("\n=== 宽格式数据（前6行） ===\n")
print(head(wide_data))

# 8. 计算环境间的相关系数
correlation_matrix <- cor(wide_data[, -1], use = "pairwise.complete.obs")

cat("\n=== 环境间相关系数矩阵 ===\n")
print(correlation_matrix)

# 9. 计算各环境的综合统计量
env_summary <- wide_data[, -1] %>%  # 去掉ID列
  summarise_all(list(
    Mean = ~mean(., na.rm = TRUE),
    SD = ~sd(., na.rm = TRUE),
    CV = ~sd(., na.rm = TRUE)/mean(., na.rm = TRUE)*100,
    Min = ~min(., na.rm = TRUE),
    Max = ~max(., na.rm = TRUE),
    Skewness = ~skewness(., na.rm = TRUE),
    Kurtosis = ~kurtosis(., na.rm = TRUE),
    N = ~sum(!is.na(.))
  ))

# 重塑环境统计结果使其更易读
env_summary_long <- env_summary %>%
  pivot_longer(
    cols = everything(),
    names_to = "Stat_Env",
    values_to = "Value"
  ) %>%
  separate(Stat_Env, into = c("Env", "Statistic"), sep = "_") %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  )

cat("\n=== 各环境统计汇总 ===\n")
print(env_summary_long)

# 10. 计算各品种的综合统计量
genotype_summary <- wide_data %>%
  mutate(
    Mean = apply(wide_data[, -1], 1, mean, na.rm = TRUE),
    SD = apply(wide_data[, -1], 1, sd, na.rm = TRUE),
    CV = SD / Mean * 100,
    Min = apply(wide_data[, -1], 1, min, na.rm = TRUE),
    Max = apply(wide_data[, -1], 1, max, na.rm = TRUE),
    Range = Max - Min,
    Stability_Index = SD / Mean
  ) %>%
  select(ID, Mean, SD, CV, Min, Max, Range, Stability_Index)

cat("\n=== 各品种统计汇总（前10行） ===\n")
print(head(genotype_summary, 10))

# 11. 保存所有结果到文件
# 保存清洗后的原始数据
write.csv(data_clean, "TKW_cleaned_data.csv", row.names = FALSE)

# 保存各品种各环境平均值（带统计信息）
write.csv(clean_means, "TKW_means_with_stats.csv", row.names = FALSE)

# 保存平均值表（仅ID, ENV, TKW）
write.table(mean_table, "TKW_mean_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 保存宽格式数据
write.csv(wide_data, "TKW_wide_format.csv", row.names = FALSE)

# 保存相关系数矩阵
write.csv(correlation_matrix, "TKW_correlation_matrix.csv")

# 保存环境统计汇总
write.csv(env_summary_long, "TKW_environment_summary.csv", row.names = FALSE)

# 保存品种统计汇总
write.csv(genotype_summary, "TKW_genotype_summary.csv", row.names = FALSE)

# 12. 生成分析报告
sink("TKW_analysis_summary.txt")
cat("TKW数据分析报告\n")
cat("================\n\n")
cat("分析日期:", date(), "\n\n")

cat("1. 数据概览\n")
cat("   原始数据行数:", nrow(data), "\n")
cat("   清洗后数据行数:", nrow(data_clean), "\n")
cat("   删除的离群值数量:", nrow(data) - nrow(data_clean), "\n")
cat("   离群值比例:", round((nrow(data) - nrow(data_clean))/nrow(data) * 100, 2), "%\n")
cat("   品种数量:", n_distinct(data_clean$ID), "\n")
cat("   环境数量:", n_distinct(data_clean$ENV), "\n")
cat("   环境列表:", paste(unique(data_clean$ENV), collapse = ", "), "\n\n")

cat("2. 平均值表结构\n")
cat("   文件: TKW_mean_table.txt\n")
cat("   列名: ID, ENV, TKW\n")
cat("   行数:", nrow(mean_table), "\n")
cat("   平均值范围:", round(min(mean_table$TKW), 2), "-", round(max(mean_table$TKW), 2), "\n\n")

cat("3. 环境间相关性\n")
cat("   相关系数矩阵已保存到: TKW_correlation_matrix.csv\n")
cat("   平均相关系数:", round(mean(correlation_matrix[upper.tri(correlation_matrix)], na.rm = TRUE), 3), "\n")
cat("   最小相关系数:", round(min(correlation_matrix[upper.tri(correlation_matrix)], na.rm = TRUE), 3), "\n")
cat("   最大相关系数:", round(max(correlation_matrix[upper.tri(correlation_matrix)], na.rm = TRUE), 3), "\n\n")

cat("4. 环境统计量摘要\n")
for(i in 1:nrow(env_summary_long)) {
  cat("   环境", env_summary_long$Env[i], ":\n")
  cat("     平均值:", round(env_summary_long$Mean[i], 2), "\n")
  cat("     标准差:", round(env_summary_long$SD[i], 2), "\n")
  cat("     变异系数:", round(env_summary_long$CV[i], 2), "%\n")
  cat("     偏度:", round(env_summary_long$Skewness[i], 3), "\n")
  cat("     峰度:", round(env_summary_long$Kurtosis[i], 3), "\n")
  cat("     品种数:", env_summary_long$N[i], "\n")
  cat("\n")
}

cat("5. 生成的文件列表\n")
cat("   - TKW_mean_table.txt: 各品种各环境平均值表\n")
cat("   - TKW_means_with_stats.csv: 带统计信息的平均值表\n")
cat("   - TKW_cleaned_data.csv: 清洗后的原始数据\n")
cat("   - TKW_wide_format.csv: 宽格式数据\n")
cat("   - TKW_correlation_matrix.csv: 相关系数矩阵\n")
cat("   - TKW_environment_summary.csv: 环境统计汇总\n")
cat("   - TKW_genotype_summary.csv: 品种统计汇总\n")
cat("   - TKW_analysis_summary.txt: 本分析报告\n")

sink()

cat("\n=== 分析完成 ===\n")
cat("主要输出文件:\n")
cat("1. TKW_mean_table.txt - 包含ID, ENV, TKW三列的平均值表\n")
cat("2. TKW_means_with_stats.csv - 包含更多统计信息的平均值表\n")
cat("3. TKW_correlation_matrix.csv - 环境间相关系数矩阵\n")
cat("4. TKW_environment_summary.csv - 各环境统计量汇总\n")
cat("5. TKW_analysis_summary.txt - 详细分析报告\n")

# 13. 可选：简单可视化
# 安装ggplot2包（如果未安装）
if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

# 创建平均值分布图
ggplot(mean_table, aes(x = ENV, y = TKW)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "各环境TKW平均值分布", 
       x = "环境", 
       y = "TKW平均值") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("TKW_mean_distribution.png", width = 8, height = 6, dpi = 300)

# 创建相关系数热图
if(!require(reshape2)) {
  install.packages("reshape2")
  library(reshape2)
}

cor_melted <- melt(correlation_matrix)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "相关系数") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  labs(title = "环境间TKW平均值相关系数热图", 
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("TKW_correlation_heatmap.png", width = 8, height = 6, dpi = 300)

cat("\n已生成可视化图表:\n")
cat("- TKW_mean_distribution.png: 平均值分布箱线图\n")
cat("- TKW_correlation_heatmap.png: 相关系数热图\n")