# ==============================================================================
# 完整分析流程：从数据读取到贝叶斯方差组分分析
# ==============================================================================

# 第一部分：安装和加载必要的包
# ==============================================================================
cat("=== 第一部分：加载必要的包 ===\n")

# 检查并安装必要的包
required_packages <- c("BGLR", "data.table", "dplyr", "tidyr", "moments", "ggplot2", "reshape2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0) {
  cat("安装必要的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)
}

# 加载包
library(BGLR)
library(data.table)
library(dplyr)
library(tidyr)
library(moments)
library(ggplot2)
library(reshape2)

cat("所有必要的包已加载完成!\n\n")

# 第二部分：读取数据并计算平均值
# ==============================================================================
cat("=== 第二部分：读取数据并计算平均值 ===\n")

# 读取原始数据
cat("读取原始数据TKW.txt...\n")
data <- fread("TKW.txt", header = TRUE, sep = "\t", na.strings = "NA")

cat("原始数据结构:\n")
str(data)
cat("\n前6行数据:\n")
print(head(data))

# 计算各品种各环境的平均值
cat("\n计算各品种各环境平均值...\n")
mean_table <- data %>%
  group_by(ID, ENV) %>%
  summarise(
    TKW = mean(TKW, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(genotype = ID, env_code = ENV)

cat("平均值表结构:\n")
str(mean_table)
cat("\n平均值表前6行:\n")
print(head(mean_table))

# 保存平均值表
write.table(mean_table, "TKW_mean_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n已保存平均值表到: TKW_mean_table.txt\n")

# 第三部分：数据预处理
# ==============================================================================
cat("\n=== 第三部分：数据预处理 ===\n")

# 指定要删除的环境
env_to_remove <- c("2019TS", "2020XC")  # 可以根据需要修改
cat("指定删除的环境:", paste(env_to_remove, collapse = ", "), "\n")

# 删除指定环境
pheno_original <- nrow(mean_table)

# 使用dplyr的filter函数来删除指定环境
pheno <- mean_table %>% filter(!env_code %in% env_to_remove)

pheno_after <- nrow(pheno)
cat("删除环境后，观测值数量从", pheno_original, "减少到", pheno_after, "\n")

# 从env_code中提取年份（前4位）和地点（剩余部分）
pheno$year <- substr(pheno$env_code, 1, 4)
pheno$location <- substr(pheno$env_code, 5, nchar(pheno$env_code))
pheno$yloc <- pheno$env_code

# 转换为因子
pheno$year <- as.factor(pheno$year)
pheno$location <- as.factor(pheno$location)
pheno$yloc <- as.factor(pheno$yloc)
pheno$var <- as.factor(pheno$genotype)

# 移除TKW缺失值
pheno <- pheno[!is.na(pheno$TKW), ]
cat("最终分析观测值数量:", nrow(pheno), "\n")

# 第四部分：数据概览
# ==============================================================================
cat("\n=== 第四部分：数据概览 ===\n")
cat("品种数量:", length(levels(pheno$var)), "\n")
cat("年份数量:", length(levels(pheno$year)), "\n")
cat("地点数量:", length(levels(pheno$location)), "\n")
cat("年份-地点组合数量:", length(levels(pheno$yloc)), "\n")
cat("TKW平均值:", round(mean(pheno$TKW, na.rm = TRUE), 2), "\n")
cat("TKW标准差:", round(sd(pheno$TKW, na.rm = TRUE), 2), "\n")
cat("TKW范围:", round(min(pheno$TKW, na.rm = TRUE), 2), "-", 
    round(max(pheno$TKW, na.rm = TRUE), 2), "\n")

# 第五部分：设置贝叶斯分析参数
# ==============================================================================
cat("\n=== 第五部分：设置贝叶斯分析参数 ===\n")

nIter <- 12000
burnIn <- 2000
thin <- 10
set.seed(195021)

cat("总迭代次数:", nIter, "\n")
cat("预烧期:", burnIn, "\n")
cat("稀疏间隔:", thin, "\n")
cat("随机种子:", 195021, "\n")
cat("实际使用的后验样本数:", (nIter - burnIn) / thin, "\n")

# 第六部分：构建设计矩阵和拟合模型
# ==============================================================================
cat("\n=== 第六部分：构建设计矩阵和拟合模型 ===\n")

# 创建设计矩阵
ZY <- model.matrix(~ pheno$year - 1)
ZL <- model.matrix(~ pheno$location - 1)
ZYL <- model.matrix(~ pheno$yloc - 1)
Zv <- model.matrix(~ pheno$var - 1)

cat("\n设计矩阵维度:\n")
cat("ZY (年份):", dim(ZY), "\n")
cat("ZL (地点):", dim(ZL), "\n")
cat("ZYL (年份×地点):", dim(ZYL), "\n")
cat("Zv (品种):", dim(Zv), "\n")

# 构建ETA列表
Eta1 <- list(
  year = list(X = ZY, model = 'BRR', saveEffects = TRUE),
  loc = list(X = ZL, model = 'BRR', saveEffects = TRUE),
  yloc = list(X = ZYL, model = 'BRR', saveEffects = TRUE),
  v = list(X = Zv, model = 'BRR', saveEffects = TRUE)
)

cat("\n正在拟合贝叶斯模型，这可能需要一些时间...\n")
start_time <- Sys.time()

# 拟合BGLR模型
fm1 <- BGLR(y = pheno$TKW, ETA = Eta1, nIter = nIter, 
            burnIn = burnIn, thin = thin, saveAt = 'm1_', verbose = TRUE)

end_time <- Sys.time()
cat("模型拟合完成! 耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n")

# 第七部分：计算方差组分
# ==============================================================================
cat("\n=== 第七部分：计算方差组分 ===\n")

# Year-variance
cat("1. 计算年份方差...\n")
B <- readBinMat('m1_ETA_year_b.bin') # posterior samples
TMP <- tcrossprod(ZY, B)
vY <- apply(FUN = var, X = TMP, MARGIN = 2)
mean_vY <- mean(vY)
se_vY <- sd(vY) / sqrt(length(vY))
ci_vY <- quantile(vY, c(0.025, 0.975))

# Location variance
cat("2. 计算地点方差...\n")
B <- readBinMat('m1_ETA_loc_b.bin') # posterior samples 
TMP <- tcrossprod(ZL, B)
vL <- apply(FUN = var, X = TMP, MARGIN = 2)
mean_vL <- mean(vL)
se_vL <- sd(vL) / sqrt(length(vL))
ci_vL <- quantile(vL, c(0.025, 0.975))

# Year-Location variance
cat("3. 计算年份×地点互作方差...\n")
B <- readBinMat('m1_ETA_yloc_b.bin') # posterior samples
TMP <- tcrossprod(ZYL, B)
vYL <- apply(FUN = var, X = TMP, MARGIN = 2)
mean_vYL <- mean(vYL)
se_vYL <- sd(vYL) / sqrt(length(vYL))
ci_vYL <- quantile(vYL, c(0.025, 0.975))

# Cultivar variance
cat("4. 计算品种方差...\n")
B <- readBinMat('m1_ETA_v_b.bin') # posterior samples
TMP <- tcrossprod(Zv, B)
vG <- apply(FUN = var, X = TMP, MARGIN = 2)
mean_vG <- mean(vG)
se_vG <- sd(vG) / sqrt(length(vG))
ci_vG <- quantile(vG, c(0.025, 0.975))

# Residual variance
cat("5. 计算残差方差...\n")
varE <- scan("m1_varE.dat")
mean_varE <- mean(varE)
se_varE <- sd(varE) / sqrt(length(varE))
ci_varE <- quantile(varE, c(0.025, 0.975))

cat("方差组分计算完成!\n")

# 第八部分：检查向量长度并确保一致
# ==============================================================================
cat("\n=== 第八部分：检查向量长度 ===\n")

# 获取所有向量的长度
lengths <- c(
  vY = length(vY),
  vL = length(vL),
  vYL = length(vYL),
  vG = length(vG),
  varE = length(varE)
)

cat("各向量长度:\n")
print(lengths)

# 找到最小长度
min_length <- min(lengths)
cat("最小长度:", min_length, "\n")

# 如果长度不一致，截取所有向量到相同长度
if(length(unique(lengths)) > 1) {
  cat("警告: 向量长度不一致，将截取到最小长度", min_length, "\n")
  
  vY <- vY[1:min_length]
  vL <- vL[1:min_length]
  vYL <- vYL[1:min_length]
  vG <- vG[1:min_length]
  varE <- varE[1:min_length]
  
  # 重新计算均值和置信区间
  mean_vY <- mean(vY)
  se_vY <- sd(vY) / sqrt(length(vY))
  ci_vY <- quantile(vY, c(0.025, 0.975))
  
  mean_vL <- mean(vL)
  se_vL <- sd(vL) / sqrt(length(vL))
  ci_vL <- quantile(vL, c(0.025, 0.975))
  
  mean_vYL <- mean(vYL)
  se_vYL <- sd(vYL) / sqrt(length(vYL))
  ci_vYL <- quantile(vYL, c(0.025, 0.975))
  
  mean_vG <- mean(vG)
  se_vG <- sd(vG) / sqrt(length(vG))
  ci_vG <- quantile(vG, c(0.025, 0.975))
  
  mean_varE <- mean(varE)
  se_varE <- sd(varE) / sqrt(length(varE))
  ci_varE <- quantile(varE, c(0.025, 0.975))
}

# 第九部分：计算总环境方差
# ==============================================================================
cat("\n=== 第九部分：计算总环境方差 ===\n")

# 计算总环境方差（年份+地点+年份×地点互作）
total_YL_variance <- vY + vL + vYL
mean_total_YL <- mean(total_YL_variance)
se_total_YL <- sd(total_YL_variance) / sqrt(length(total_YL_variance))
ci_total_YL <- quantile(total_YL_variance, c(0.025, 0.975))

cat("总环境方差 (Y+L+YL):", round(mean_total_YL, 3), 
    "[", round(ci_total_YL[1], 3), ",", round(ci_total_YL[2], 3), "]\n")

# 第十部分：计算总方差和方差组分占比
# ==============================================================================
cat("\n=== 第十部分：计算方差组分占比 ===\n")

# 计算总表型方差（观测值的方差）
total_phenotypic_variance <- var(pheno$TKW)
cat("总表型方差:", round(total_phenotypic_variance, 3), "\n")

# 计算各组分在总表型方差中的占比
prop_year <- mean_vY / total_phenotypic_variance * 100
prop_loc <- mean_vL / total_phenotypic_variance * 100
prop_yloc <- mean_vYL / total_phenotypic_variance * 100
prop_var <- mean_vG / total_phenotypic_variance * 100
prop_residual <- mean_varE / total_phenotypic_variance * 100
prop_total_yl <- mean_total_YL / total_phenotypic_variance * 100

# 计算模型解释的总方差（1 - 残差占比）
total_explained <- 100 - prop_residual

cat("\n方差组分占比 (基于总表型方差):\n")
cat("年份效应:", round(prop_year, 2), "%\n")
cat("地点效应:", round(prop_loc, 2), "%\n")
cat("年份-地点互作效应:", round(prop_yloc, 2), "%\n")
cat("总环境效应:", round(prop_total_yl, 2), "%\n")
cat("品种效应:", round(prop_var, 2), "%\n")
cat("残差效应:", round(prop_residual, 2), "%\n")
cat("模型总解释率:", round(total_explained, 2), "%\n")

# 第十一部分：计算遗传力
# ==============================================================================
cat("\n=== 第十一部分：计算遗传参数 ===\n")

# 广义遗传力 (Broad-sense heritability)
H2 <- mean_vG / (mean_vG + mean_total_YL + mean_varE) * 100
cat("广义遗传力 (H², 基于模型方差组分):", round(H2, 2), "%\n")

# 另一种计算方式：基于总表型方差的遗传力
H2_alt <- prop_var
cat("广义遗传力 (H², 基于总表型方差):", round(H2_alt, 2), "%\n")

# 基因型×环境互作占比
GE_ratio <- prop_yloc / (prop_var + prop_yloc) * 100
cat("基因型×环境互作相对重要性:", round(GE_ratio, 2), "%\n")

# 计算环境方差与遗传方差的比例
env_gen_ratio <- mean_total_YL / mean_vG
cat("环境方差/遗传方差比例:", round(env_gen_ratio, 3), "\n")

# 第十二部分：创建输出目录
# ==============================================================================
cat("\n=== 第十二部分：创建输出目录 ===\n")

output_dir <- "baselineModel"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出目录:", output_dir, "\n")
} else {
  cat("输出目录已存在:", output_dir, "\n")
}

# 第十三部分：按照GW模型表格格式输出结果
# ==============================================================================
cat("\n=== 第十三部分：输出结果 ===\n")

# 创建与GW模型相同的表格格式
cat("Model                            Year                 Location             YxL                   Cultivar             Env. Cov.            SNPxEC              Error\n")
cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

# 主模型行 - 按照GW模型格式
tl_year_str <- sprintf("%.3f", mean_vY)
tl_location_str <- sprintf("%.3f", mean_vL)
tl_yxl_str <- sprintf("%.3f", mean_vYL)
tl_cultivar_str <- sprintf("%.3f", mean_vG)
tl_error_str <- sprintf("%.3f", mean_varE)

# 置信区间字符串
tl_year_ci <- sprintf("[%.3f,%.3f]", ci_vY[1], ci_vY[2])
tl_location_ci <- sprintf("[%.3f,%.3f]", ci_vL[1], ci_vL[2])
tl_yxl_ci <- sprintf("[%.3f,%.3f]", ci_vYL[1], ci_vYL[2])
tl_cultivar_ci <- sprintf("[%.3f,%.3f]", ci_vG[1], ci_vG[2])
tl_error_ci <- sprintf("[%.3f,%.3f]", ci_varE[1], ci_varE[2])

# 第一行：均值
cat(sprintf("%-32s%-21s%-21s%-22s%-21s%-21s%-21s%s\n", 
            "TL (baseline)", 
            tl_year_str,             # Year
            tl_location_str,         # Location  
            tl_yxl_str,              # YxL
            tl_cultivar_str,         # Cultivar
            "-",                     # Env. Cov.
            "-",                     # SNPxEC
            tl_error_str))           # Error

# 第二行：置信区间
cat(sprintf("%-32s%-21s%-21s%-22s%-21s%-21s%-21s%s\n", 
            "", 
            tl_year_ci,              # Year CI
            tl_location_ci,          # Location CI
            tl_yxl_ci,               # YxL CI
            tl_cultivar_ci,          # Cultivar CI
            "-",                     # Env. Cov. CI
            "-",                     # SNPxEC CI
            tl_error_ci))            # Error CI

# 第十四部分：保存所有结果到文件
# ==============================================================================
cat("\n=== 第十四部分：保存结果到文件 ===\n")

# 1. 创建与GW模型相同格式的结果表格
table_results <- data.frame(
  Model = c("TL (baseline)", ""),
  Year = c(tl_year_str, tl_year_ci),
  Location = c(tl_location_str, tl_location_ci),
  YxL = c(tl_yxl_str, tl_yxl_ci),
  Cultivar = c(tl_cultivar_str, tl_cultivar_ci),
  Env.Cov = c("-", "-"),
  SNPxEC = c("-", "-"),
  Error = c(tl_error_str, tl_error_ci)
)

write.csv(table_results, file.path(output_dir, "TL_variance_components_table.csv"), row.names = FALSE)
cat("  - 保存TL模型表格格式结果: TL_variance_components_table.csv\n")

# 2. 保存详细的数值结果
detailed_results <- data.frame(
  Component = c("Year", "Location", "Year_x_Location", "Cultivar", "Residual", "Total_Environmental", "Total_Phenotypic"),
  Mean = c(mean_vY, mean_vL, mean_vYL, mean_vG, mean_varE, mean_total_YL, total_phenotypic_variance),
  CI_Lower = c(ci_vY[1], ci_vL[1], ci_vYL[1], ci_vG[1], ci_varE[1], ci_total_YL[1], NA),
  CI_Upper = c(ci_vY[2], ci_vL[2], ci_vYL[2], ci_vG[2], ci_varE[2], ci_total_YL[2], NA),
  Proportion = c(prop_year, prop_loc, prop_yloc, prop_var, prop_residual, prop_total_yl, 100),
  SE = c(se_vY, se_vL, se_vYL, se_vG, se_varE, se_total_YL, NA)
)

write.csv(detailed_results, file.path(output_dir, "TL_detailed_variance_components.csv"), row.names = FALSE)
cat("  - 保存TL模型详细方差组分: TL_detailed_variance_components.csv\n")

# 3. 保存所有方差轨迹数据
variance_traces <- data.frame(
  Iteration = 1:length(vY),
  Year = vY,
  Location = vL,
  YearLocation = vYL,
  Genotype = vG,
  Residual = varE,
  Total_YL = total_YL_variance
)

write.csv(variance_traces, file.path(output_dir, "TL_variance_traces_data.csv"), row.names = FALSE)
cat("  - 保存TL模型方差轨迹数据: TL_variance_traces_data.csv\n")

# 4. 保存BGLR模型结果摘要
model_summary <- list(
  nIter = nIter,
  burnIn = burnIn,
  thin = thin,
  nPosterior = length(vY),
  data_info = list(
    n_observations = nrow(pheno),
    n_genotypes = length(levels(pheno$var)),
    n_years = length(levels(pheno$year)),
    n_locations = length(levels(pheno$location)),
    n_yloc = length(levels(pheno$yloc)),
    deleted_environments = env_to_remove
  ),
  variance_components = list(
    year = mean_vY,
    location = mean_vL,
    year_location = mean_vYL,
    cultivar = mean_vG,
    residual = mean_varE,
    total_environmental = mean_total_YL
  ),
  proportions = list(
    year = prop_year,
    location = prop_loc,
    year_location = prop_yloc,
    cultivar = prop_var,
    residual = prop_residual,
    total_environmental = prop_total_yl,
    total_explained = total_explained
  ),
  heritability = H2,
  GE_ratio = GE_ratio,
  env_gen_ratio = env_gen_ratio,
  computation_time = as.numeric(difftime(end_time, start_time, units = "mins"))
)

saveRDS(model_summary, file.path(output_dir, "TL_model_summary.rds"))
cat("  - 保存TL模型摘要: TL_model_summary.rds\n")

# 第十五部分：生成详细的分析报告
# ==============================================================================
cat("\n=== 第十五部分：生成详细分析报告 ===\n")

report <- paste0(
  "TL模型方差组分分析结果报告\n",
  "===========================\n\n",
  "分析日期: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  
  "1. 数据信息:\n",
  "   - 原始观测值数量: ", pheno_original, "\n",
  "   - 删除环境后观测值数量: ", pheno_after, "\n",
  "   - 最终分析观测值数量: ", nrow(pheno), "\n",
  "   - 品种数量: ", length(levels(pheno$var)), "\n",
  "   - 年份数量: ", length(levels(pheno$year)), "\n",
  "   - 地点数量: ", length(levels(pheno$location)), "\n",
  "   - 年份-地点组合数量: ", length(levels(pheno$yloc)), "\n",
  "   - TKW平均值: ", round(mean(pheno$TKW, na.rm = TRUE), 2), "\n",
  "   - TKW标准差: ", round(sd(pheno$TKW, na.rm = TRUE), 2), "\n",
  "   - TKW范围: ", round(min(pheno$TKW, na.rm = TRUE), 2), " - ", 
      round(max(pheno$TKW, na.rm = TRUE), 2), "\n",
  "   - 删除的环境: ", paste(env_to_remove, collapse = ", "), "\n\n",
  
  "2. 贝叶斯分析参数:\n",
  "   - 总迭代次数: ", nIter, "\n",
  "   - 预烧期: ", burnIn, "\n",
  "   - 稀疏间隔: ", thin, "\n",
  "   - 实际使用迭代次数: ", min_length, "\n",
  "   - 随机种子: 195021\n",
  "   - 模型拟合时间: ", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n\n",
  
  "3. 方差组分结果 (均值 [95%置信区间]):\n",
  "   - 年份方差: ", round(mean_vY, 4), " [", round(ci_vY[1], 4), ",", round(ci_vY[2], 4), "]\n",
  "   - 地点方差: ", round(mean_vL, 4), " [", round(ci_vL[1], 4), ",", round(ci_vL[2], 4), "]\n",
  "   - 年份×地点互作方差: ", round(mean_vYL, 4), " [", round(ci_vYL[1], 4), ",", round(ci_vYL[2], 4), "]\n",
  "   - 品种方差: ", round(mean_vG, 4), " [", round(ci_vG[1], 4), ",", round(ci_vG[2], 4), "]\n",
  "   - 残差方差: ", round(mean_varE, 4), " [", round(ci_varE[1], 4), ",", round(ci_varE[2], 4), "]\n",
  "   - 总环境方差: ", round(mean_total_YL, 4), " [", round(ci_total_YL[1], 4), ",", round(ci_total_YL[2], 4), "]\n",
  "   - 总表型方差: ", round(total_phenotypic_variance, 4), "\n\n",
  
  "4. 方差组分占比 (基于总表型方差):\n",
  "   - 年份效应: ", round(prop_year, 2), "%\n",
  "   - 地点效应: ", round(prop_loc, 2), "%\n",
  "   - 年份×地点互作: ", round(prop_yloc, 2), "%\n", 
  "   - 总环境效应: ", round(prop_total_yl, 2), "%\n",
  "   - 品种效应: ", round(prop_var, 2), "%\n",
  "   - 残差效应: ", round(prop_residual, 2), "%\n",
  "   - 模型总解释率: ", round(total_explained, 2), "%\n\n",
  
  "5. 遗传参数:\n",
  "   - 广义遗传力 (H², 基于模型方差组分): ", round(H2, 2), "%\n",
  "   - 广义遗传力 (H², 基于总表型方差): ", round(H2_alt, 2), "%\n",
  "   - 基因型×环境互作相对重要性: ", round(GE_ratio, 2), "%\n",
  "   - 环境方差/遗传方差比例: ", round(env_gen_ratio, 3), "\n\n",
  
  "6. 生成的输出文件:\n",
  "   - TL模型表格格式结果: TL_variance_components_table.csv\n",
  "   - TL模型详细方差组分: TL_detailed_variance_components.csv\n",
  "   - TL模型方差轨迹数据: TL_variance_traces_data.csv\n",
  "   - TL模型摘要: TL_model_summary.rds\n",
  "   - 本报告文件: TL_variance_analysis_report.txt\n"
)

writeLines(report, file.path(output_dir, "TL_variance_analysis_report.txt"))
cat("  - 保存TL模型详细报告: TL_variance_analysis_report.txt\n")

# 第十六部分：绘制图表
# ==============================================================================
cat("\n=== 第十六部分：绘制图表 ===\n")

# 1. 创建方差轨迹图
trace_plot_data <- data.frame(
  Iteration = rep(1:length(vY), 5),
  Variance = c(vY, vL, vYL, vG, varE),
  Component = factor(rep(c("Year", "Location", "Year×Location", "Genotype", "Residual"), 
                         each = length(vY)),
                     levels = c("Year", "Location", "Year×Location", "Genotype", "Residual"))
)

p1 <- ggplot(trace_plot_data, aes(x = Iteration, y = Variance, color = Component)) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  facet_wrap(~ Component, scales = "free_y", ncol = 1) +
  labs(title = "TL模型方差组分轨迹图",
       x = "迭代次数",
       y = "方差") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))

ggsave(file.path(output_dir, "TL_variance_trace_plots.png"), p1, width = 8, height = 10, dpi = 300)
cat("  - 保存TL模型方差轨迹图: TL_variance_trace_plots.png\n")

# 2. 创建方差组分占比饼图
prop_data <- data.frame(
  Component = c("Year", "Location", "Year×Location", "Genotype", "Residual"),
  Proportion = c(prop_year, prop_loc, prop_yloc, prop_var, prop_residual),
  Label = paste(c("Year", "Location", "Y×L", "Genotype", "Residual"), 
                sprintf("%.1f%%", c(prop_year, prop_loc, prop_yloc, prop_var, prop_residual)))
)

p2 <- ggplot(prop_data, aes(x = "", y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
            position = position_stack(vjust = 0.5), 
            size = 4) +
  labs(title = "TL模型方差组分占比",
       fill = "方差组分") +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "TL_variance_proportion_pie.png"), p2, width = 8, height = 6, dpi = 300)
cat("  - 保存TL模型方差组分占比饼图: TL_variance_proportion_pie.png\n")

# 3. 创建置信区间图
ci_data <- data.frame(
  Component = c("Year", "Location", "Year×Location", "Genotype", "Residual", "Total Env"),
  Mean = c(mean_vY, mean_vL, mean_vYL, mean_vG, mean_varE, mean_total_YL),
  Lower = c(ci_vY[1], ci_vL[1], ci_vYL[1], ci_vG[1], ci_varE[1], ci_total_YL[1]),
  Upper = c(ci_vY[2], ci_vL[2], ci_vYL[2], ci_vG[2], ci_varE[2], ci_total_YL[2])
)

p3 <- ggplot(ci_data, aes(x = Component, y = Mean)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "blue") +
  labs(title = "TL模型方差组分估计值及95%置信区间",
       x = "方差组分",
       y = "方差") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "TL_variance_CI_plots.png"), p3, width = 10, height = 6, dpi = 300)
cat("  - 保存TL模型方差置信区间图: TL_variance_CI_plots.png\n")

# 第十七部分：清理临时文件
# ==============================================================================
cat("\n=== 第十七部分：清理临时文件 ===\n")

temp_files <- c("m1_varE.dat", "m1_ETA_year_b.bin", "m1_ETA_loc_b.bin", 
                "m1_ETA_yloc_b.bin", "m1_ETA_v_b.bin",
                "m1_ETA_year_b.bin.txt", "m1_ETA_loc_b.bin.txt", 
                "m1_ETA_yloc_b.bin.txt", "m1_ETA_v_b.bin.txt")

temp_removed <- 0
for(file in temp_files) {
  if(file.exists(file)) {
    file.remove(file)
    cat("  - 删除临时文件:", file, "\n")
    temp_removed <- temp_removed + 1
  }
}

cat("临时文件清理完成! 共删除", temp_removed, "个临时文件。\n")

# 第十八部分：最终总结
# ==============================================================================
cat("\n=== 第十八部分：分析完成 ===\n")

cat("\n=== 分析流程总结 ===\n")
cat("1. 数据读取和预处理: 完成\n")
cat("2. 平均值计算: 完成\n")
cat("3. 数据预处理: 完成\n")
cat("4. 贝叶斯模型拟合: 完成\n")
cat("5. 方差组分计算: 完成\n")
cat("6. 结果分析和可视化: 完成\n")
cat("7. 结果保存: 完成\n\n")

cat("=== 主要结果 ===\n")
cat("- 年份方差: ", round(mean_vY, 3), " [", round(ci_vY[1], 3), ",", round(ci_vY[2], 3), "]\n")
cat("- 地点方差: ", round(mean_vL, 3), " [", round(ci_vL[1], 3), ",", round(ci_vL[2], 3), "]\n")
cat("- 年份×地点互作方差: ", round(mean_vYL, 3), " [", round(ci_vYL[1], 3), ",", round(ci_vYL[2], 3), "]\n")
cat("- 品种方差: ", round(mean_vG, 3), " [", round(ci_vG[1], 3), ",", round(ci_vG[2], 3), "]\n")
cat("- 残差方差: ", round(mean_varE, 3), " [", round(ci_varE[1], 3), ",", round(ci_varE[2], 3), "]\n")
cat("- 总环境方差: ", round(mean_total_YL, 3), " [", round(ci_total_YL[1], 3), ",", round(ci_total_YL[2], 3), "]\n")
cat("- 广义遗传力 (H²): ", round(H2, 2), "%\n")
cat("- 基因型×环境互作相对重要性: ", round(GE_ratio, 2), "%\n\n")

cat("=== 所有输出文件已保存到", output_dir, "目录 ===\n")
cat("1. TL_variance_components_table.csv - TL模型表格格式结果\n")
cat("2. TL_detailed_variance_components.csv - TL模型详细方差组分\n")
cat("3. TL_variance_traces_data.csv - TL模型方差轨迹数据\n")
cat("4. TL_model_summary.rds - TL模型摘要\n")
cat("5. TL_variance_analysis_report.txt - TL模型详细报告\n")
cat("6. TL_variance_trace_plots.png - TL模型方差轨迹图\n")
cat("7. TL_variance_proportion_pie.png - TL模型方差组分占比饼图\n")
cat("8. TL_variance_CI_plots.png - TL模型方差置信区间图\n\n")

cat("=== 数据处理统计 ===\n")
cat("- 原始观测值数量:", pheno_original, "\n")
cat("- 删除环境后观测值数量:", pheno_after, "\n")
cat("- 最终分析观测值数量:", nrow(pheno), "\n")
cat("- 品种数量:", length(levels(pheno$var)), "\n")
cat("- 年份数量:", length(levels(pheno$year)), "\n")
cat("- 地点数量:", length(levels(pheno$location)), "\n")
cat("- 年份-地点组合数量:", length(levels(pheno$yloc)), "\n")
cat("- 删除的环境:", paste(env_to_remove, collapse = ", "), "\n\n")

cat("=== 贝叶斯分析统计 ===\n")
cat("- 总迭代次数:", nIter, "\n")
cat("- 预烧期:", burnIn, "\n")
cat("- 稀疏间隔:", thin, "\n")
cat("- 实际使用迭代次数:", min_length, "\n")
cat("- 模型拟合时间:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n\n")

cat("所有分析已完成! 请检查", output_dir, "目录中的结果文件。\n")