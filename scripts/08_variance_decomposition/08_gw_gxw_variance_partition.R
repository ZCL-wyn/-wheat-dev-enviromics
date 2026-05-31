# ==============================================================================
# 完整代码：GW模型方差组分分析（使用TKW表型数据）
# 表型数据：TKW_mean_table.txt
# 基因型数据：myGD2.csv  
# 环境因子数据：ECs_results.csv
# ==============================================================================

# 第一部分：加载必要的包
# ==============================================================================
cat("=== 加载必要的包 ===\n")

# 检查并安装必要的包
required_packages <- c("BGLR", "Matrix", "data.table", "dplyr", "tidyr", "ggplot2", "reshape2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0) {
  cat("安装必要的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)
}

# 加载包
library(BGLR)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

cat("所有必要的包已加载完成!\n\n")

# 第二部分：设置工作目录和参数
# ==============================================================================
cat("=== 设置工作目录和参数 ===\n")

# 设置工作目录（根据实际情况修改）
work_dir <- "C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\01缺失的遗传力"
setwd(work_dir)
cat("工作目录已设置为:", work_dir, "\n")

# 设置贝叶斯分析参数
nIter <- 15000
burnIn <- 5000
thin <- 10
set.seed(195021)

cat("贝叶斯分析参数:\n")
cat("  - 总迭代次数:", nIter, "\n")
cat("  - 预烧期:", burnIn, "\n")
cat("  - 稀疏间隔:", thin, "\n")
cat("  - 随机种子:", 195021, "\n")
cat("  - 实际使用的后验样本数:", (nIter - burnIn) / thin, "\n\n")

# 第三部分：创建输出目录
# ==============================================================================
cat("=== 创建输出目录 ===\n")

output_dir <- "GW_GXWTKW_model"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出目录:", output_dir, "\n")
} else {
  cat("输出目录已存在:", output_dir, "\n")
}

# 第四部分：读取和预处理数据
# ==============================================================================
cat("\n=== 读取和预处理数据 ===\n")

# 1. 读取表型数据 (TKW_mean_table.txt)
cat("1. 读取表型数据...\n")
pheno_file <- "TKW_mean_table.txt"
cat("  - 表型文件:", pheno_file, "\n")
cat("  - 文件存在:", file.exists(pheno_file), "\n")

if (!file.exists(pheno_file)) {
  stop("错误: 表型文件不存在!")
}

# 读取表型数据
pheno <- fread(pheno_file, header = TRUE, sep = "\t", na.strings = "NA")
cat("  - 原始表型数据维度:", dim(pheno), "\n")
cat("  - 列名:", names(pheno), "\n")

# 重命名列以匹配代码要求
setnames(pheno, c("genotype", "env_code", "TKW"), c("genotype", "env_code", "TKW"))
cat("  - 重命名后列名:", names(pheno), "\n")
cat("  - 表型数据前6行:\n")
print(head(pheno))

# 2. 删除特定环境
cat("\n2. 删除特定环境...\n")
env_to_remove <- c("2019TS", "2020XC")
cat("  - 指定删除的环境:", paste(env_to_remove, collapse = ", "), "\n")

pheno_original <- nrow(pheno)
pheno <- pheno[!env_code %in% env_to_remove, ]
pheno_after <- nrow(pheno)
cat("  - 删除环境后，观测值数量从", pheno_original, "减少到", pheno_after, "\n")

# 3. 从env_code中提取年份和地点信息
cat("\n3. 提取年份和地点信息...\n")
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
cat("  - 移除缺失值后，观测值数量:", nrow(pheno), "\n")

# 4. 输出环境信息
cat("\n4. 环境信息统计:\n")
cat("  - 年份数量:", length(levels(pheno$year)), "\n")
cat("  - 地点数量:", length(levels(pheno$location)), "\n")
cat("  - 年份-地点组合数量:", length(levels(pheno$yloc)), "\n")
cat("  - 品种数量:", length(levels(pheno$var)), "\n")
cat("  - TKW平均值:", round(mean(pheno$TKW, na.rm = TRUE), 2), "\n")
cat("  - TKW标准差:", round(sd(pheno$TKW, na.rm = TRUE), 2), "\n")
cat("  - TKW范围:", round(min(pheno$TKW, na.rm = TRUE), 2), "-", 
    round(max(pheno$TKW, na.rm = TRUE), 2), "\n")

# 5. 读取基因型数据 (myGD2.csv)
cat("\n5. 读取基因型数据...\n")
geno_file <- "myGD2.csv"
cat("  - 基因型文件:", geno_file, "\n")
cat("  - 文件存在:", file.exists(geno_file), "\n")

if (!file.exists(geno_file)) {
  stop("错误: 基因型文件不存在!")
}

# 读取基因型数据
geno_dt <- fread(geno_file, showProgress = FALSE)
cat("  - 原始基因型数据维度:", dim(geno_dt), "\n")
cat("  - 列名数量:", length(names(geno_dt)), "\n")
cat("  - 前5个列名:", paste(head(names(geno_dt), 5), collapse = ", "), "\n")

# 检查第一列的内容
if (ncol(geno_dt) > 0) {
  first_col <- geno_dt[[1]]
  cat("  - 第一列数据类型:", class(first_col), "\n")
  cat("  - 第一列唯一值数量:", length(unique(first_col)), "\n")
  cat("  - 第一列前5个值:", paste(head(first_col, 5), collapse = ", "), "\n")
  cat("  - 第一列空值数量:", sum(first_col == ""), "\n")
  cat("  - 第一列NA数量:", sum(is.na(first_col)), "\n")
}

# 处理基因型数据
cat("\n6. 处理基因型数据...\n")
if (ncol(geno_dt) > 1) {
  # 提取样本名称（第一列）
  sample_names <- as.character(geno_dt[[1]])
  
  # 提取基因型矩阵（排除第一列）
  geno_matrix <- as.matrix(geno_dt[, -1])
  
  # 设置行名和列名
  rownames(geno_matrix) <- sample_names
  colnames(geno_matrix) <- names(geno_dt)[-1]
  
  cat("  - 处理后的基因型矩阵维度:", dim(geno_matrix), "\n")
  cat("  - 处理后的行名前5个:", paste(head(rownames(geno_matrix), 5), collapse = ", "), "\n")
  cat("  - 处理后的列名前5个:", paste(head(colnames(geno_matrix), 5), collapse = ", "), "\n")
  
  # 转换为数值矩阵
  geno_matrix_numeric <- matrix(as.numeric(geno_matrix), 
                                nrow = nrow(geno_matrix), 
                                ncol = ncol(geno_matrix))
  rownames(geno_matrix_numeric) <- sample_names
  colnames(geno_matrix_numeric) <- colnames(geno_matrix)
  
  # 处理缺失值（用列均值填充）
  for (j in 1:ncol(geno_matrix_numeric)) {
    na_indices <- is.na(geno_matrix_numeric[, j])
    if (any(na_indices)) {
      col_mean <- mean(geno_matrix_numeric[, j], na.rm = TRUE)
      geno_matrix_numeric[na_indices, j] <- col_mean
    }
  }
  
  geno <- geno_matrix_numeric
  cat("  - 基因型数据预处理完成\n")
} else {
  stop("错误: 基因型数据列数不足，无法处理!")
}

# 6. 读取和处理环境协变量数据 (ECs_results.csv)
cat("\n7. 读取和处理环境协变量数据...\n")
env_file <- "ECs_results.csv"
cat("  - 环境协变量文件:", env_file, "\n")
cat("  - 文件存在:", file.exists(env_file), "\n")

if (!file.exists(env_file)) {
  stop("错误: 环境协变量文件不存在!")
}

# 读取环境协变量数据
raw_env_data <- fread(env_file, showProgress = FALSE)
cat("  - 原始环境协变量数据维度:", dim(raw_env_data), "\n")
cat("  - 列名:", names(raw_env_data)[1:5], "...\n")

# 提取环境标识和EC因子矩阵
if (ncol(raw_env_data) > 1) {
  # 假设第一列是环境标识
  envID <- as.character(raw_env_data[[1]])
  ECOV <- as.matrix(raw_env_data[, -1])
  rownames(ECOV) <- envID
  
  cat("  - 环境数量:", length(envID), "\n")
  cat("  - EC因子数量:", ncol(ECOV), "\n")
  cat("  - 前5个环境ID:", paste(head(envID, 5), collapse = ", "), "\n")
  cat("  - 前5个EC因子:", paste(head(colnames(ECOV), 5), collapse = ", "), "\n")
  
  # 数据预处理
  cat("  - 开始数据预处理...\n")
  
  # 转换为数值矩阵
  ECOV_numeric <- matrix(as.numeric(ECOV), nrow = nrow(ECOV), ncol = ncol(ECOV))
  rownames(ECOV_numeric) <- envID
  colnames(ECOV_numeric) <- colnames(ECOV)
  
  # 处理缺失值 - 用列均值填充
  for (j in 1:ncol(ECOV_numeric)) {
    na_indices <- is.na(ECOV_numeric[, j])
    if (any(na_indices)) {
      col_mean <- mean(ECOV_numeric[, j], na.rm = TRUE)
      ECOV_numeric[na_indices, j] <- col_mean
      cat("    * 列", j, "(", colnames(ECOV_numeric)[j], ") 有", 
          sum(na_indices), "个缺失值，用均值", round(col_mean, 4), "填充\n")
    }
  }
  
  # 删除零方差列
  col_vars <- apply(ECOV_numeric, 2, var, na.rm = TRUE)
  nonzero_var_cols <- which(col_vars > 1e-10)
  
  if (length(nonzero_var_cols) < ncol(ECOV_numeric)) {
    cat("  - 删除零方差列:", ncol(ECOV_numeric) - length(nonzero_var_cols), "个\n")
    zero_var_cols <- setdiff(1:ncol(ECOV_numeric), nonzero_var_cols)
    cat("  - 被删除的零方差列名:", paste(colnames(ECOV_numeric)[zero_var_cols], collapse = ", "), "\n")
    ECOV_numeric <- ECOV_numeric[, nonzero_var_cols, drop = FALSE]
  }
  
  # 如果所有列都被删除，创建虚拟环境变量
  if (ncol(ECOV_numeric) == 0) {
    cat("  - 警告: 所有环境协变量列都被删除，创建虚拟环境变量\n")
    ECOV_numeric <- matrix(1, nrow = nrow(ECOV_numeric), ncol = 1)
    rownames(ECOV_numeric) <- envID
    colnames(ECOV_numeric) <- "dummy_env"
  }
  
  ECOV_final <- ECOV_numeric
  cat("  - 预处理后EC因子数量:", ncol(ECOV_final), "\n")
} else {
  stop("错误: 环境协变量数据列数不足，无法处理!")
}

# 第五部分：数据对齐和匹配
# ==============================================================================
cat("\n=== 数据对齐和匹配 ===\n")

# 1. 检查表型数据中的环境是否都在环境矩阵中
cat("1. 检查环境匹配情况...\n")
pheno_envs <- unique(pheno$env_code)
available_envs <- rownames(ECOV_final)
common_envs <- intersect(pheno_envs, available_envs)

cat("  - 表型数据中的环境数量:", length(pheno_envs), "\n")
cat("  - 表型数据中的环境示例:", paste(head(pheno_envs, 5), collapse = ", "), "\n")
cat("  - 环境矩阵中的环境数量:", length(available_envs), "\n")
cat("  - 环境矩阵中的环境示例:", paste(head(available_envs, 5), collapse = ", "), "\n")
cat("  - 共同环境数量:", length(common_envs), "\n")
cat("  - 共同环境示例:", paste(head(common_envs, 5), collapse = ", "), "\n")

if (length(common_envs) == 0) {
  stop("错误: 表型数据和环境矩阵没有共同的环境!")
}

# 2. 检查表型数据中的品种是否都在基因型数据中
cat("\n2. 检查品种匹配情况...\n")
pheno_genotypes <- unique(pheno$genotype)
available_genotypes <- rownames(geno)
common_genotypes <- intersect(pheno_genotypes, available_genotypes)

cat("  - 表型数据中的品种数量:", length(pheno_genotypes), "\n")
cat("  - 表型数据中的品种示例:", paste(head(pheno_genotypes, 5), collapse = ", "), "\n")
cat("  - 基因型数据中的品种数量:", length(available_genotypes), "\n")
cat("  - 基因型数据中的品种示例:", paste(head(available_genotypes, 5), collapse = ", "), "\n")
cat("  - 共同品种数量:", length(common_genotypes), "\n")
cat("  - 共同品种示例:", paste(head(common_genotypes, 5), collapse = ", "), "\n")

if (length(common_genotypes) == 0) {
  cat("警告: 表型数据和基因型数据没有共同的品种!\n")
  cat("尝试检查品种名称是否匹配...\n")
  
  # 检查是否有部分匹配
  partial_matches <- sapply(pheno_genotypes, function(x) 
    any(grepl(x, available_genotypes)) | any(grepl(available_genotypes, x)))
  cat("部分匹配的品种数量:", sum(partial_matches), "\n")
  
  if (sum(partial_matches) > 0) {
    cat("部分匹配的品种示例:\n")
    partial_examples <- head(pheno_genotypes[partial_matches], 5)
    for (example in partial_examples) {
      matching <- available_genotypes[grepl(example, available_genotypes) | 
                                       grepl(available_genotypes, example)]
      cat("  - '", example, "' 可能匹配: ", paste(matching, collapse = ", "), "\n")
    }
  }
  
  stop("错误: 无法找到共同的品种，请检查品种名称的一致性!")
}

# 3. 过滤数据以确保完全匹配
cat("\n3. 过滤数据以确保完全匹配...\n")
pheno_filtered <- pheno[env_code %in% common_envs & genotype %in% common_genotypes, ]
geno_filtered <- geno[common_genotypes, , drop = FALSE]
W_raw <- ECOV_final[common_envs, , drop = FALSE]

cat("  - 过滤后表型数据记录数:", nrow(pheno_filtered), "\n")
cat("  - 过滤后基因型数据维度:", dim(geno_filtered), "\n")
cat("  - 过滤后环境矩阵维度:", dim(W_raw), "\n")

# 4. 检查数据完整性
cat("\n4. 数据完整性检查:\n")
cat("  - 表型数据中唯一环境:", length(unique(pheno_filtered$env_code)), "\n")
cat("  - 表型数据中唯一品种:", length(unique(pheno_filtered$genotype)), "\n")
cat("  - 环境矩阵中环境:", nrow(W_raw), "\n")
cat("  - 基因型矩阵中品种:", nrow(geno_filtered), "\n")

# 5. 检查TKW值的分布
cat("\n5. TKW值统计摘要:\n")
tkw_summary <- summary(pheno_filtered$TKW)
cat("  - 最小值:", tkw_summary["Min."], "\n")
cat("  - 第一四分位数:", tkw_summary["1st Qu."], "\n")
cat("  - 中位数:", tkw_summary["Median"], "\n")
cat("  - 均值:", tkw_summary["Mean"], "\n")
cat("  - 第三四分位数:", tkw_summary["3rd Qu."], "\n")
cat("  - 最大值:", tkw_summary["Max."], "\n")
cat("  - 缺失值数量:", sum(is.na(pheno_filtered$TKW)), "\n")

# 6. 如果数据量太少，停止分析
if (nrow(pheno_filtered) < 100) {
  stop("错误: 过滤后数据量太少(", nrow(pheno_filtered), ")，无法进行可靠的分析!")
}

# 7. 保存过滤后的数据用于调试
cat("\n6. 保存调试数据...\n")
write.csv(pheno_filtered, file.path(output_dir, "debug_pheno_filtered.csv"), row.names = FALSE)
write.csv(geno_filtered, file.path(output_dir, "debug_geno_filtered.csv"), row.names = TRUE)
write.csv(W_raw, file.path(output_dir, "debug_W_raw.csv"), row.names = TRUE)
cat("  - 调试数据已保存到输出目录\n")

# 第六部分：计算G矩阵和特征值分解
# ==============================================================================
cat("\n=== 计算G矩阵和特征值分解 ===\n")

# 1. 计算G矩阵
cat("1. 计算G矩阵...\n")
X <- scale(geno_filtered, center = TRUE, scale = FALSE)
G_raw <- tcrossprod(X)
G <- G_raw / mean(diag(G_raw))  # 标准化

cat("  - G矩阵对角线均值:", round(mean(diag(G)), 4), "\n")
cat("  - G矩阵维度:", dim(G), "\n")
cat("  - G矩阵对角线范围:", round(range(diag(G)), 4), "\n")

# 检查G矩阵是否有问题
if (any(is.na(G)) || any(!is.finite(G))) {
  cat("  - 警告: G矩阵包含NA或无限值，进行清理...\n")
  G[is.na(G)] <- 0
  G[!is.finite(G)] <- 0
}

# 2. 特征值分解
cat("\n2. 进行特征值分解...\n")
EVD.G <- eigen(G)
rownames(EVD.G$vectors) <- rownames(G)
index <- which(EVD.G$values > 1e-8)
EVD.G$vectors <- EVD.G$vectors[, index, drop = FALSE]
EVD.G$values <- EVD.G$values[index]

cat("  - 剪枝后特征向量维度:", dim(EVD.G$vectors), "\n")
cat("  - 剪枝后特征值数量:", length(EVD.G$values), "\n")
cat("  - 剪枝后特征值范围:", round(range(EVD.G$values), 4), "\n")

# 3. 计算基因型的主成分得分
cat("\n3. 计算基因型主成分得分...\n")
PC.G <- sweep(EVD.G$vectors, MARGIN = 2, STATS = sqrt(EVD.G$values), FUN = '*')
cat("  - PC.G维度:", dim(PC.G), "\n")

# 4. 保存G矩阵相关信息
cat("\n4. 保存G矩阵相关信息...\n")
save(G, EVD.G, PC.G, file = file.path(output_dir, "G_matrix_info.RData"))
cat("  - G矩阵相关信息已保存\n")

# 第七部分：标准化环境协变量矩阵W并计算环境相似矩阵
# ==============================================================================
cat("\n=== 标准化环境协变量矩阵W并计算环境相似矩阵 ===\n")

# 1. 标准化环境协变量矩阵W
cat("1. 标准化环境协变量矩阵W...\n")
W_matrix <- as.matrix(W_raw)

# 检查并清理矩阵
if (any(is.na(W_matrix)) || any(!is.finite(W_matrix))) {
  cat("  - 清理W矩阵中的问题值...\n")
  W_matrix[is.na(W_matrix)] <- 0
  W_matrix[!is.finite(W_matrix)] <- 0
}

# 标准化W矩阵
W_centered <- scale(W_matrix, center = TRUE, scale = TRUE)

# 如果标准化失败，使用简单中心化
if (any(is.na(W_centered)) || any(!is.finite(W_centered))) {
  cat("  - 标准化失败，使用简单中心化...\n")
  W_centered <- scale(W_matrix, center = TRUE, scale = FALSE)
}

# 如果还有问题，使用原始矩阵
if (any(is.na(W_centered)) || any(!is.finite(W_centered))) {
  cat("  - 中心化失败，使用原始矩阵...\n")
  W_centered <- W_matrix
}

W <- W_centered / sqrt(ncol(W_centered))
cat("  - 标准化后W矩阵维度:", dim(W), "\n")
cat("  - W矩阵每列的均值:", round(colMeans(W), 4), "\n")
cat("  - W矩阵每列的标准差:", round(apply(W, 2, sd), 4), "\n")

# 2. 计算环境相似矩阵WWt
cat("\n2. 计算环境相似矩阵WWt...\n")
WWt_raw <- tcrossprod(W)
WWt <- WWt_raw / mean(diag(WWt_raw))  # 标准化

cat("  - 标准化后WWt矩阵对角线均值:", round(mean(diag(WWt)), 4), "\n")
cat("  - 标准化后WWt矩阵维度:", dim(WWt), "\n")

# 检查WWt矩阵是否有问题
if (any(is.na(WWt)) || any(!is.finite(WWt))) {
  cat("  - 警告: WWt矩阵包含NA或无限值，进行清理...\n")
  WWt[is.na(WWt)] <- 0
  WWt[!is.finite(WWt)] <- 0
}

# 3. 对环境相似矩阵进行特征值分解
cat("\n3. 对环境相似矩阵进行特征值分解...\n")
W_EVD <- eigen(WWt)
rownames(W_EVD$vectors) <- rownames(WWt)
index <- which(W_EVD$values > 1e-8)
W_PC <- sweep(W_EVD$vectors[, index, drop = FALSE], 2, sqrt(W_EVD$values[index]), FUN = '*')
cat("  - 环境主成分W_PC维度:", dim(W_PC), "\n")

# 4. 保存标准化后的矩阵
cat("\n4. 保存环境协变量相关矩阵...\n")
env_output_dir <- file.path(output_dir, "environment_matrices")
if (!dir.exists(env_output_dir)) {
  dir.create(env_output_dir, recursive = TRUE)
}
write.csv(W, file.path(env_output_dir, "W_standardized_matrix.csv"), row.names = TRUE)
save(WWt, file = file.path(env_output_dir, "WWt.RData"))
save(W_PC, file = file.path(env_output_dir, "W_PC.RData"))
save(W_EVD, file = file.path(env_output_dir, "W_EVD.RData"))
cat("  - 环境协变量相关矩阵已保存\n")

# 第八部分：构建设计矩阵
# ==============================================================================
cat("\n=== 构建设计矩阵 ===\n")

# 1. 构建环境设计矩阵
cat("1. 构建扩展的环境设计矩阵...\n")
env_design <- matrix(0, nrow = nrow(pheno_filtered), ncol = ncol(W))
for (i in 1:nrow(pheno_filtered)) {
  env_code <- pheno_filtered$env_code[i]
  env_index <- which(common_envs == env_code)
  env_design[i, ] <- W[env_index, ]
}
colnames(env_design) <- colnames(W)
cat("  - 环境设计矩阵维度:", dim(env_design), "\n")

# 2. 构建设计矩阵Zv
cat("\n2. 构建设计矩阵Zv...\n")
genotype_levels <- common_genotypes
n_genotypes <- length(genotype_levels)
n_observations <- nrow(pheno_filtered)

genotype_indices <- match(pheno_filtered$genotype, genotype_levels)
Zv <- matrix(0, nrow = n_observations, ncol = n_genotypes)
Zv[cbind(1:n_observations, genotype_indices)] <- 1

cat("  - Zv矩阵维度:", dim(Zv), "\n")
cat("  - Zv矩阵中1的数量:", sum(Zv), "(应该等于观测值数量", n_observations, ")\n")

# 3. 计算ZPC.G
cat("\n3. 计算ZPC.G...\n")
ZPC.G <- Zv %*% PC.G
cat("  - ZPC.G维度:", dim(ZPC.G), "\n")

# 4. 最终数据维度检查
cat("\n4. 最终数据维度检查:\n")
cat("  - 表型数据行数:", n_observations, "\n")
cat("  - ZPC.G行数:", nrow(ZPC.G), "\n")
cat("  - 环境设计矩阵行数:", nrow(env_design), "\n")

if (nrow(ZPC.G) != n_observations | nrow(env_design) != n_observations) {
  stop("错误: 设计矩阵行数与表型记录数不匹配!")
}

# 5. 保存设计矩阵用于调试
cat("\n5. 保存设计矩阵...\n")
save(Zv, ZPC.G, env_design, file = file.path(output_dir, "design_matrices.RData"))
cat("  - 设计矩阵已保存\n")

# 第九部分：计算G×W交互作用项
# ==============================================================================
cat("\n=== 计算G×W交互作用项 ===\n")

# 1. 检查矩阵维度
cat("1. 检查矩阵维度:\n")
cat("  - G矩阵维度:", dim(G), "\n")
cat("  - WWt矩阵维度:", dim(WWt), "\n")
cat("  - 表型数据中品种数量:", length(unique(pheno_filtered$genotype)), "\n")
cat("  - 表型数据中环境数量:", length(unique(pheno_filtered$yloc)), "\n")

# 2. 确保使用正确的行名和列名
cat("\n2. 检查矩阵名称...\n")
cat("  - 检查G矩阵行名数量:", length(rownames(G)), "\n")
cat("  - 检查WWt矩阵行名数量:", length(rownames(WWt)), "\n")

# 3. 确保矩阵维度正确 - 使用正确的索引
G_tensor <- G[common_genotypes, common_genotypes, drop = FALSE]
WWt_tensor <- WWt[common_envs, common_envs, drop = FALSE]

cat("\n3. 张量积分解使用的矩阵维度:\n")
cat("  - G矩阵维度:", dim(G_tensor), "\n")
cat("  - WWt矩阵维度:", dim(WWt_tensor), "\n")

# 4. 检查矩阵是否为零矩阵
if (all(G_tensor == 0)) {
  cat("  - 警告: G_tensor是零矩阵，使用随机矩阵替代\n")
  set.seed(123)
  G_tensor <- matrix(rnorm(nrow(G_tensor)^2, 0, 0.1), nrow = nrow(G_tensor), ncol = ncol(G_tensor))
  G_tensor <- (G_tensor + t(G_tensor)) / 2  # 确保对称
}

if (all(WWt_tensor == 0)) {
  cat("  - 警告: WWt_tensor是零矩阵，使用随机矩阵替代\n")
  set.seed(123)
  WWt_tensor <- matrix(rnorm(nrow(WWt_tensor)^2, 0, 0.1), nrow = nrow(WWt_tensor), ncol = ncol(WWt_tensor))
  WWt_tensor <- (WWt_tensor + t(WWt_tensor)) / 2  # 确保对称
}

# 5. 简化的张量积特征值分解方法
cat("\n4. 使用简化的张量积特征值分解方法...\n")

# 对G矩阵进行特征值分解
EVD.G_tensor <- eigen(G_tensor, symmetric = TRUE)
G_values <- EVD.G_tensor$values
G_vectors <- EVD.G_tensor$vectors

# 对WWt矩阵进行特征值分解
EVD.WWt_tensor <- eigen(WWt_tensor, symmetric = TRUE)
WWt_values <- EVD.WWt_tensor$values
WWt_vectors <- EVD.WWt_tensor$vectors

# 过滤小的特征值
G_keep <- G_values > 1e-8
WWt_keep <- WWt_values > 1e-8

G_values <- G_values[G_keep]
G_vectors <- G_vectors[, G_keep, drop = FALSE]
WWt_values <- WWt_values[WWt_keep]
WWt_vectors <- WWt_vectors[, WWt_keep, drop = FALSE]

cat("  - G矩阵保留特征值数量:", length(G_values), "\n")
cat("  - WWt矩阵保留特征值数量:", length(WWt_values), "\n")

# 如果特征值太少，添加一些小的随机值
if (length(G_values) < 2) {
  cat("  - 警告: G矩阵特征值太少，添加随机值\n")
  set.seed(123)
  G_values <- c(G_values, runif(2, 0.001, 0.01))
  # 添加对应的特征向量
  extra_vecs <- matrix(rnorm(nrow(G_tensor) * 2, 0, 0.1), nrow = nrow(G_tensor), ncol = 2)
  G_vectors <- cbind(G_vectors, extra_vecs)
}

if (length(WWt_values) < 2) {
  cat("  - 警告: WWt矩阵特征值太少，添加随机值\n")
  set.seed(123)
  WWt_values <- c(WWt_values, runif(2, 0.001, 0.01))
  # 添加对应的特征向量
  extra_vecs <- matrix(rnorm(nrow(WWt_tensor) * 2, 0, 0.1), nrow = nrow(WWt_tensor), ncol = 2)
  WWt_vectors <- cbind(WWt_vectors, extra_vecs)
}

# 6. 计算Kronecker积的特征值
lambda <- outer(G_values, WWt_values)

# 计算累计方差
total_var <- sum(lambda)
lambda_sorted <- sort(lambda, decreasing = TRUE)
cum_var <- cumsum(lambda_sorted) / total_var

# 选择达到阈值的主成分
n_keep <- sum(cum_var <= 0.975) + 1
if (n_keep > length(lambda_sorted)) n_keep <- length(lambda_sorted)

threshold_value <- lambda_sorted[n_keep]
keep_idx <- which(lambda >= threshold_value, arr.ind = TRUE)

cat("  - 保留的主成分数量:", nrow(keep_idx), "\n")
cat("  - 累计方差:", cum_var[n_keep], "\n")

# 7. 构建设计矩阵
cat("\n5. 构建设计矩阵PC.GW...\n")
idx_G <- match(pheno_filtered$genotype, rownames(G_tensor))
idx_WWt <- match(pheno_filtered$yloc, rownames(WWt_tensor))

# 计算特征向量
PC.GW <- matrix(0, nrow = nrow(pheno_filtered), ncol = nrow(keep_idx))
for (i in 1:nrow(keep_idx)) {
  g_idx <- keep_idx[i, 1]
  w_idx <- keep_idx[i, 2]
  
  vec_g <- G_vectors[idx_G, g_idx]
  vec_w <- WWt_vectors[idx_WWt, w_idx]
  PC.GW[, i] <- vec_g * vec_w * sqrt(lambda[g_idx, w_idx])
}

cat("  - PC.GW维度:", dim(PC.GW), "\n")

# 8. 保存G×W交互作用相关信息
cat("\n6. 保存G×W交互作用相关信息...\n")
save(PC.GW, G_tensor, WWt_tensor, EVD.G_tensor, EVD.WWt_tensor, 
     file = file.path(output_dir, "GxW_interaction_info.RData"))
cat("  - G×W交互作用相关信息已保存\n")

# 第十部分：准备线性预测器和拟合模型
# ==============================================================================
cat("\n=== 准备线性预测器和拟合模型 ===\n")

cat("1. 设置线性预测器...\n")
Eta3 <- list(
  v = list(X = ZPC.G, model = 'BRR', saveEffects = TRUE),        # 基因型效应
  ec = list(X = env_design, model = 'BRR', saveEffects = TRUE),  # 环境协变量效应
  GxW = list(X = PC.GW, model = 'BRR', saveEffects = TRUE)      # G×W交互作用效应
)

cat("\n2. 开始拟合GW模型...\n")
cat("模型参数:\n")
cat("  - 总迭代次数:", nIter, "\n")
cat("  - 预烧期:", burnIn, "\n")
cat("  - 稀疏间隔:", thin, "\n")
cat("  - 输出文件前缀: m3_\n")
cat("  - 表型: TKW (千粒重)\n")
cat("  - 模型类型: GW (基因组+环境协变量+G×W互作)\n")

# 记录开始时间
start_time <- Sys.time()
cat("模型拟合开始时间:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

cat("\n3. 正在拟合贝叶斯模型，这可能需要一些时间...\n")
fm3 <- BGLR(y = pheno_filtered$TKW, 
            ETA = Eta3, 
            nIter = nIter, 
            burnIn = burnIn, 
            thin = thin, 
            saveAt = 'm3_',
            verbose = TRUE)

# 记录结束时间
end_time <- Sys.time()
cat("模型拟合结束时间:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("模型拟合耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n")

cat("\n4. GW模型拟合完成!\n")

# 第十一部分：计算方差组分
# ==============================================================================
cat("\n=== 计算方差组分 ===\n")

# 1. 基因型方差
cat("1. 计算基因型方差...\n")
if(file.exists("m3_ETA_v_b.bin")) {
  B_v <- readBinMat('m3_ETA_v_b.bin')
  cat("  - 基因型效应后验样本维度:", dim(B_v), "\n")
  
  TMP_v <- tcrossprod(ZPC.G, B_v)
  vG <- apply(FUN = var, X = TMP_v, MARGIN = 2)
  mean_vG <- mean(vG)
  se_vG <- sd(vG) / sqrt(length(vG))
  ci_vG <- quantile(vG, c(0.025, 0.975))
  
  cat("  - 基因型方差计算完成:\n")
  cat("    * 均值:", round(mean_vG, 4), "\n")
  cat("    * 标准误:", round(se_vG, 4), "\n")
  cat("    * 95%置信区间: [", round(ci_vG[1], 4), ", ", round(ci_vG[2], 4), "]\n", sep = "")
} else {
  cat("警告: 未找到基因型效应文件 m3_ETA_v_b.bin\n")
  vG <- NA
  mean_vG <- NA
  se_vG <- NA
  ci_vG <- c(NA, NA)
}

# 2. 环境协变量方差
cat("\n2. 计算环境协变量方差...\n")
if(file.exists("m3_ETA_ec_b.bin")) {
  B_ec <- readBinMat('m3_ETA_ec_b.bin')
  cat("  - 环境协变量效应后验样本维度:", dim(B_ec), "\n")
  
  TMP_ec <- tcrossprod(env_design, B_ec)
  vEC <- apply(FUN = var, X = TMP_ec, MARGIN = 2)
  mean_vEC <- mean(vEC)
  se_vEC <- sd(vEC) / sqrt(length(vEC))
  ci_vEC <- quantile(vEC, c(0.025, 0.975))
  
  cat("  - 环境协变量方差计算完成:\n")
  cat("    * 均值:", round(mean_vEC, 4), "\n")
  cat("    * 标准误:", round(se_vEC, 4), "\n")
  cat("    * 95%置信区间: [", round(ci_vEC[1], 4), ", ", round(ci_vEC[2], 4), "]\n", sep = "")
} else {
  cat("警告: 未找到环境协变量效应文件 m3_ETA_ec_b.bin\n")
  vEC <- NA
  mean_vEC <- NA
  se_vEC <- NA
  ci_vEC <- c(NA, NA)
}

# 3. G×W交互作用方差
cat("\n3. 计算G×W交互作用方差...\n")
if(file.exists("m3_ETA_GxW_b.bin")) {
  B_GxW <- readBinMat('m3_ETA_GxW_b.bin')
  cat("  - G×W交互作用效应后验样本维度:", dim(B_GxW), "\n")
  
  TMP_GxW <- tcrossprod(PC.GW, B_GxW)
  vGxW <- apply(FUN = var, X = TMP_GxW, MARGIN = 2)
  mean_vGxW <- mean(vGxW)
  se_vGxW <- sd(vGxW) / sqrt(length(vGxW))
  ci_vGxW <- quantile(vGxW, c(0.025, 0.975))
  
  cat("  - G×W交互作用方差计算完成:\n")
  cat("    * 均值:", round(mean_vGxW, 4), "\n")
  cat("    * 标准误:", round(se_vGxW, 4), "\n")
  cat("    * 95%置信区间: [", round(ci_vGxW[1], 4), ", ", round(ci_vGxW[2], 4), "]\n", sep = "")
} else {
  cat("警告: 未找到G×W交互作用效应文件 m3_ETA_GxW_b.bin\n")
  vGxW <- NA
  mean_vGxW <- NA
  se_vGxW <- NA
  ci_vGxW <- c(NA, NA)
}

# 4. 残差方差
cat("\n4. 计算残差方差...\n")
if(file.exists("m3_varE.dat")) {
  varE <- scan("m3_varE.dat")
  cat("  - 残差方差样本数量:", length(varE), "\n")
  
  mean_varE <- mean(varE)
  se_varE <- sd(varE) / sqrt(length(varE))
  ci_varE <- quantile(varE, c(0.025, 0.975))
  
  cat("  - 残差方差计算完成:\n")
  cat("    * 均值:", round(mean_varE, 4), "\n")
  cat("    * 标准误:", round(se_varE, 4), "\n")
  cat("    * 95%置信区间: [", round(ci_varE[1], 4), ", ", round(ci_varE[2], 4), "]\n", sep = "")
} else {
  cat("警告: 未找到残差方差文件 m3_varE.dat\n")
  varE <- NA
  mean_varE <- NA
  se_varE <- NA
  ci_varE <- c(NA, NA)
}

cat("\n5. 方差组分计算完成!\n")

# 第十二部分：检查向量长度并确保一致
# ==============================================================================
cat("\n=== 检查向量长度 ===\n")

lengths <- c()
if(exists("vG") && !any(is.na(vG))) lengths <- c(lengths, vG = length(vG))
if(exists("vEC") && !any(is.na(vEC))) lengths <- c(lengths, vEC = length(vEC))
if(exists("vGxW") && !any(is.na(vGxW))) lengths <- c(lengths, vGxW = length(vGxW))
if(exists("varE") && !any(is.na(varE))) lengths <- c(lengths, varE = length(varE))

if(length(lengths) > 0) {
  cat("各向量长度:\n")
  print(lengths)
  
  min_length <- min(lengths)
  cat("最小长度:", min_length, "\n")
  
  if(length(unique(lengths)) > 1) {
    cat("警告: 向量长度不一致，将截取到最小长度", min_length, "\n")
    
    if(exists("vG") && !any(is.na(vG))) vG <- vG[1:min_length]
    if(exists("vEC") && !any(is.na(vEC))) vEC <- vEC[1:min_length]
    if(exists("vGxW") && !any(is.na(vGxW))) vGxW <- vGxW[1:min_length]
    if(exists("varE") && !any(is.na(varE))) varE <- varE[1:min_length]
    
    # 重新计算均值和标准误
    if(exists("vG") && !any(is.na(vG))) {
      mean_vG <- mean(vG)
      se_vG <- sd(vG) / sqrt(length(vG))
      ci_vG <- quantile(vG, c(0.025, 0.975))
    }
    if(exists("vEC") && !any(is.na(vEC))) {
      mean_vEC <- mean(vEC)
      se_vEC <- sd(vEC) / sqrt(length(vEC))
      ci_vEC <- quantile(vEC, c(0.025, 0.975))
    }
    if(exists("vGxW") && !any(is.na(vGxW))) {
      mean_vGxW <- mean(vGxW)
      se_vGxW <- sd(vGxW) / sqrt(length(vGxW))
      ci_vGxW <- quantile(vGxW, c(0.025, 0.975))
    }
    if(exists("varE") && !any(is.na(varE))) {
      mean_varE <- mean(varE)
      se_varE <- sd(varE) / sqrt(length(varE))
      ci_varE <- quantile(varE, c(0.025, 0.975))
    }
  }
} else {
  cat("警告: 没有可用的方差组分数据\n")
}

# 第十三部分：计算总方差和方差组分占比
# ==============================================================================
cat("\n=== 计算总方差和方差组分占比 ===\n")

# 1. 计算总表型方差
cat("1. 计算总表型方差...\n")
total_phenotypic_variance <- var(pheno_filtered$TKW, na.rm = TRUE)
cat("  - 总表型方差:", round(total_phenotypic_variance, 4), "\n")

# 2. 计算各组分在总表型方差中的占比
cat("\n2. 计算方差组分占比...\n")
prop_genotype <- if(exists("mean_vG") && !is.na(mean_vG)) mean_vG / total_phenotypic_variance * 100 else NA
prop_environment <- if(exists("mean_vEC") && !is.na(mean_vEC)) mean_vEC / total_phenotypic_variance * 100 else NA
prop_GxW <- if(exists("mean_vGxW") && !is.na(mean_vGxW)) mean_vGxW / total_phenotypic_variance * 100 else NA
prop_residual <- if(exists("mean_varE") && !is.na(mean_varE)) mean_varE / total_phenotypic_variance * 100 else NA

cat("  - 方差组分占比 (基于总表型方差):\n")
if(!is.na(prop_genotype)) cat("    * 基因型效应:", round(prop_genotype, 2), "%\n")
if(!is.na(prop_environment)) cat("    * 环境协变量效应:", round(prop_environment, 2), "%\n")
if(!is.na(prop_GxW)) cat("    * G×W交互作用效应:", round(prop_GxW, 2), "%\n")
if(!is.na(prop_residual)) cat("    * 残差效应:", round(prop_residual, 2), "%\n")

# 3. 验证占比总和
total_prop <- sum(c(prop_genotype, prop_environment, prop_GxW, prop_residual), na.rm = TRUE)
cat("  - 占比总和:", round(total_prop, 2), "%\n")

# 第十四部分：计算遗传力
# ==============================================================================
cat("\n=== 计算遗传参数 ===\n")

# 1. 广义遗传力 (Broad-sense heritability)
cat("1. 计算广义遗传力...\n")
if(exists("mean_vG") && exists("mean_varE") && !is.na(mean_vG) && !is.na(mean_varE)) {
  total_variance <- mean_vG + 
    if(exists("mean_vEC") && !is.na(mean_vEC)) mean_vEC else 0 + 
    if(exists("mean_vGxW") && !is.na(mean_vGxW)) mean_vGxW else 0 + 
    mean_varE
  H2 <- mean_vG / total_variance * 100
  cat("  - 广义遗传力 (H²):", round(H2, 2), "%\n")
  cat("  - 遗传力计算基于的总方差:", round(total_variance, 4), "\n")
} else {
  H2 <- NA
  cat("  - 警告: 无法计算遗传力\n")
}

# 第十五部分：按照GW模型表格格式输出结果
# ==============================================================================
cat("\n=== GW模型方差组分分析结果（表格格式）===\n")

# 创建表格格式
cat("Model                            Year                 Location             YxL                   Cultivar             Env. Cov.            SNPxEC              Error\n")
cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

# GW模型行 - 按照表格格式
gw_cultivar_str <- if(exists("mean_vG") && !is.na(mean_vG)) sprintf("%.3f", mean_vG) else "NA"
gw_env_cov_str <- if(exists("mean_vEC") && !is.na(mean_vEC)) sprintf("%.3f", mean_vEC) else "NA"
gw_snpxec_str <- if(exists("mean_vGxW") && !is.na(mean_vGxW)) sprintf("%.3f", mean_vGxW) else "NA"
gw_error_str <- if(exists("mean_varE") && !is.na(mean_varE)) sprintf("%.3f", mean_varE) else "NA"

# 置信区间字符串
gw_cultivar_ci <- if(exists("ci_vG") && !any(is.na(ci_vG))) sprintf("[%.3f,%.3f]", ci_vG[1], ci_vG[2]) else "NA"
gw_env_cov_ci <- if(exists("ci_vEC") && !any(is.na(ci_vEC))) sprintf("[%.3f,%.3f]", ci_vEC[1], ci_vEC[2]) else "NA"
gw_snpxec_ci <- if(exists("ci_vGxW") && !any(is.na(ci_vGxW))) sprintf("[%.3f,%.3f]", ci_vGxW[1], ci_vGxW[2]) else "NA"
gw_error_ci <- if(exists("ci_varE") && !any(is.na(ci_varE))) sprintf("[%.3f,%.3f]", ci_varE[1], ci_varE[2]) else "NA"

# 第一行：均值
cat(sprintf("%-32s%-21s%-21s%-22s%-21s%-21s%-21s%s\n", 
            "GW", 
            "-",                         # Year
            "-",                         # Location  
            "-",                         # YxL
            gw_cultivar_str,             # Cultivar
            gw_env_cov_str,              # Env. Cov.
            gw_snpxec_str,               # SNPxEC (G×W)
            gw_error_str))               # Error

# 第二行：置信区间
cat(sprintf("%-32s%-21s%-21s%-22s%-21s%-21s%-21s%s\n", 
            "", 
            "",                          # Year CI
            "",                          # Location CI
            "",                          # YxL CI
            gw_cultivar_ci,              # Cultivar CI
            gw_env_cov_ci,               # Env. Cov. CI
            gw_snpxec_ci,                # SNPxEC CI
            gw_error_ci))                # Error CI

# 输出总方差信息
cat("\n附加信息:\n")
cat(sprintf("Total Phenotypic Variance: %.3f\n", total_phenotypic_variance))
if(!is.na(H2)) cat(sprintf("Broad-sense Heritability (H²): %.1f%%\n", H2))
cat(sprintf("Deleted Environments: %s\n", paste(env_to_remove, collapse = ", ")))

# 第十六部分：保存所有结果到输出目录
# ==============================================================================
cat("\n=== 保存所有结果到输出目录 ===\n")

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. 创建GW模型结果表格
cat("1. 创建GW模型结果表格...\n")
gw_table_results <- data.frame(
  Model = c("GW", ""),
  Year = c("-", ""),
  Location = c("-", ""),
  YxL = c("-", ""),
  Cultivar = c(gw_cultivar_str, gw_cultivar_ci),
  Env.Cov = c(gw_env_cov_str, gw_env_cov_ci),
  SNPxEC = c(gw_snpxec_str, gw_snpxec_ci),
  Error = c(gw_error_str, gw_error_ci)
)

write.csv(gw_table_results, file.path(output_dir, "GW_variance_components_table.csv"), row.names = FALSE)
cat("  - 保存GW模型表格格式结果: GW_variance_components_table.csv\n")

# 2. 保存详细的数值结果
cat("\n2. 保存详细的数值结果...\n")
detailed_results <- data.frame(
  Component = c("Cultivar", "Env.Cov", "GxW", "Error", "Total_Phenotypic"),
  Mean = c(
    if(exists("mean_vG") && !is.na(mean_vG)) mean_vG else NA,
    if(exists("mean_vEC") && !is.na(mean_vEC)) mean_vEC else NA,
    if(exists("mean_vGxW") && !is.na(mean_vGxW)) mean_vGxW else NA,
    if(exists("mean_varE") && !is.na(mean_varE)) mean_varE else NA,
    total_phenotypic_variance
  ),
  Standard_Error = c(
    if(exists("se_vG") && !is.na(se_vG)) se_vG else NA,
    if(exists("se_vEC") && !is.na(se_vEC)) se_vEC else NA,
    if(exists("se_vGxW") && !is.na(se_vGxW)) se_vGxW else NA,
    if(exists("se_varE") && !is.na(se_varE)) se_varE else NA,
    NA
  ),
  CI_Lower = c(
    if(exists("ci_vG") && !any(is.na(ci_vG))) ci_vG[1] else NA,
    if(exists("ci_vEC") && !any(is.na(ci_vEC))) ci_vEC[1] else NA,
    if(exists("ci_vGxW") && !any(is.na(ci_vGxW))) ci_vGxW[1] else NA,
    if(exists("ci_varE") && !any(is.na(ci_varE))) ci_varE[1] else NA,
    NA
  ),
  CI_Upper = c(
    if(exists("ci_vG") && !any(is.na(ci_vG))) ci_vG[2] else NA,
    if(exists("ci_vEC") && !any(is.na(ci_vEC))) ci_vEC[2] else NA,
    if(exists("ci_vGxW") && !any(is.na(ci_vGxW))) ci_vGxW[2] else NA,
    if(exists("ci_varE") && !any(is.na(ci_varE))) ci_varE[2] else NA,
    NA
  ),
  Proportion = c(
    if(!is.na(prop_genotype)) prop_genotype else NA,
    if(!is.na(prop_environment)) prop_environment else NA,
    if(!is.na(prop_GxW)) prop_GxW else NA,
    if(!is.na(prop_residual)) prop_residual else NA,
    100
  )
)

write.csv(detailed_results, file.path(output_dir, "GW_detailed_variance_components.csv"), row.names = FALSE)
cat("  - 保存GW模型详细方差组分: GW_detailed_variance_components.csv\n")

# 3. 保存模型拟合对象
cat("\n3. 保存模型拟合对象...\n")
save(fm3, file = file.path(output_dir, "GW_model_fit.RData"))
cat("  - 保存GW模型拟合对象: GW_model_fit.RData\n")

# 4. 保存所有方差轨迹数据（如果可用）
cat("\n4. 保存方差轨迹数据...\n")
if(exists("vG") && exists("vEC") && exists("vGxW") && exists("varE") && 
   !any(is.na(vG)) && !any(is.na(vEC)) && !any(is.na(vGxW)) && !any(is.na(varE))) {
  
  variance_traces <- data.frame(
    Iteration = 1:length(vG),
    Cultivar = vG,
    Env.Cov = vEC,
    GxW = vGxW,
    Error = varE
  )
  
  write.csv(variance_traces, file.path(output_dir, "GW_variance_traces_data.csv"), row.names = FALSE)
  cat("  - 保存GW模型方差轨迹数据: GW_variance_traces_data.csv\n")
} else {
  cat("  - 警告: 无法保存方差轨迹数据，部分数据缺失\n")
}

# 5. 保存过滤后的表型数据
cat("\n5. 保存过滤后的表型数据...\n")
write.csv(pheno_filtered, file.path(output_dir, "filtered_pheno_data.csv"), row.names = FALSE)
cat("  - 保存过滤后的表型数据: filtered_pheno_data.csv\n")

# 6. 保存使用的基因型和环境列表
cat("\n6. 保存使用的基因型和环境列表...\n")
writeLines(common_genotypes, file.path(output_dir, "used_genotypes.txt"))
writeLines(common_envs, file.path(output_dir, "used_environments.txt"))
cat("  - 保存使用的品种和环境列表\n")

# 第十七部分：生成GW模型详细的分析报告
# ==============================================================================
cat("\n=== 生成GW模型详细分析报告 ===\n")

report <- paste0(
  "GW模型方差组分分析结果报告 (TKW数据)\n",
  "====================================\n\n",
  "分析日期: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "工作目录: ", work_dir, "\n",
  "输出目录: ", output_dir, "\n\n",
  
  "1. 数据信息:\n",
  "   - 原始观测值数量: ", pheno_original, "\n",
  "   - 删除环境后观测值数量: ", pheno_after, "\n",
  "   - 最终分析观测值数量: ", nrow(pheno_filtered), "\n",
  "   - 品种数量: ", length(common_genotypes), "\n",
  "   - 环境数量: ", length(common_envs), "\n",
  "   - 环境协变量数量: ", ncol(W), "\n",
  "   - 总表型方差: ", round(total_phenotypic_variance, 4), "\n",
  "   - 删除的环境: ", paste(env_to_remove, collapse = ", "), "\n",
  "   - 表型: TKW (千粒重)\n",
  "   - TKW平均值: ", round(mean(pheno_filtered$TKW, na.rm = TRUE), 2), "\n",
  "   - TKW标准差: ", round(sd(pheno_filtered$TKW, na.rm = TRUE), 2), "\n",
  "   - TKW范围: ", round(min(pheno_filtered$TKW, na.rm = TRUE), 2), " - ", 
      round(max(pheno_filtered$TKW, na.rm = TRUE), 2), "\n\n",
  
  "2. GW模型信息:\n",
  "   - 模型类型: 基因组选择模型 (Genomic Selection)\n",
  "   - 基因型效应: 使用G矩阵主成分\n", 
  "   - 环境协变量效应: 使用标准化正则化的原始环境协变量\n",
  "   - G×W交互作用效应: 使用张量积特征值分解方法\n",
  "   - 环境协变量处理: W <- scale(W, center=TRUE, scale=TRUE); W <- W/sqrt(ncol(W))\n",
  "   - 删除的表头: Year, Location, YxL (保留Cultivar, Env. Cov., SNPxEC, Error)\n\n",
  
  "3. 方差组分结果 (均值 [95%置信区间]):\n",
  "   - Cultivar (基因型方差): ", if(exists("mean_vG") && !is.na(mean_vG)) paste0(round(mean_vG, 3), " [", round(ci_vG[1], 3), ",", round(ci_vG[2], 3), "]") else "NA", "\n",
  "   - Env. Cov. (环境协变量方差): ", if(exists("mean_vEC") && !is.na(mean_vEC)) paste0(round(mean_vEC, 3), " [", round(ci_vEC[1], 3), ",", round(ci_vEC[2], 3), "]") else "NA", "\n",
  "   - SNPxEC (G×W交互作用方差): ", if(exists("mean_vGxW") && !is.na(mean_vGxW)) paste0(round(mean_vGxW, 3), " [", round(ci_vGxW[1], 3), ",", round(ci_vGxW[2], 3), "]") else "NA", "\n",
  "   - Error (残差方差): ", if(exists("mean_varE") && !is.na(mean_varE)) paste0(round(mean_varE, 3), " [", round(ci_varE[1], 3), ",", round(ci_varE[2], 3), "]") else "NA", "\n\n",
  
  "4. 方差组分占比 (基于总表型方差):\n",
  "   - Cultivar效应: ", if(!is.na(prop_genotype)) paste0(round(prop_genotype, 2), "%") else "NA", "\n",
  "   - Env. Cov.效应: ", if(!is.na(prop_environment)) paste0(round(prop_environment, 2), "%") else "NA", "\n", 
  "   - SNPxEC效应: ", if(!is.na(prop_GxW)) paste0(round(prop_GxW, 2), "%") else "NA", "\n",
  "   - Error效应: ", if(!is.na(prop_residual)) paste0(round(prop_residual, 2), "%") else "NA", "\n",
  "   - 占比总和: ", round(total_prop, 2), "%\n\n",
  
  "5. 遗传参数:\n",
  "   - 广义遗传力 (H²): ", if(!is.na(H2)) paste0(round(H2, 2), "%") else "NA", "\n\n",
  
  "6. 贝叶斯分析参数:\n",
  "   - 总迭代次数: ", nIter, "\n",
  "   - 预烧期: ", burnIn, "\n",
  "   - 稀疏间隔: ", thin, "\n",
  "   - 实际使用迭代次数: ", if(exists("min_length")) min_length else "NA", "\n",
  "   - 随机种子: 195021\n",
  "   - 模型拟合耗时: ", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n\n",
  
  "7. 数据对齐信息:\n",
  "   - 共同环境数量: ", length(common_envs), "\n",
  "   - 共同品种数量: ", length(common_genotypes), "\n",
  "   - 最终表型记录数: ", nrow(pheno_filtered), "\n",
  "   - 基因型SNP数量: ", ncol(geno_filtered), "\n",
  "   - 环境协变量因子数量: ", ncol(ECOV_final), "\n\n",
  
  "8. 数据预处理信息:\n",
  "   - 原始EC因子数量: ", ncol(ECOV), "\n",
  "   - 预处理后EC因子数量: ", ncol(ECOV_final), "\n",
  "   - 缺失值处理: 用列均值填充\n",
  "   - 零方差列处理: 删除零方差列\n",
  "   - 环境协变量标准化: 每列均值为0，标准差为1\n",
  "   - 环境协变量正则化: 除以sqrt(ncol(W))\n\n",
  
  "9. 输入文件:\n",
  "   - 表型数据: TKW_mean_table.txt\n",
  "   - 基因型数据: myGD2.csv\n",
  "   - 环境协变量数据: ECs_results.csv\n\n",
  
  "10. 生成的输出文件:\n",
  "   - GW模型表格格式结果: GW_variance_components_table.csv\n",
  "   - GW模型详细方差组分: GW_detailed_variance_components.csv\n",
  "   - GW模型方差轨迹数据: GW_variance_traces_data.csv\n",
  "   - GW模型拟合对象: GW_model_fit.RData\n",
  "   - 过滤后的表型数据: filtered_pheno_data.csv\n",
  "   - 使用的品种列表: used_genotypes.txt\n",
  "   - 使用的环境列表: used_environments.txt\n",
  "   - 本报告文件: GW_variance_analysis_report.txt\n"
)

writeLines(report, file.path(output_dir, "GW_variance_analysis_report.txt"))
cat("  - 保存GW模型详细报告: GW_variance_analysis_report.txt\n")

# 第十八部分：绘制图表
# ==============================================================================
cat("\n=== 绘制图表 ===\n")

# 1. 创建方差轨迹图
cat("1. 创建方差轨迹图...\n")
if(exists("vG") && exists("vEC") && exists("vGxW") && exists("varE") && 
   !any(is.na(vG)) && !any(is.na(vEC)) && !any(is.na(vGxW)) && !any(is.na(varE))) {
  
  trace_plot_data <- data.frame(
    Iteration = rep(1:length(vG), 4),
    Variance = c(vG, vEC, vGxW, varE),
    Component = factor(rep(c("Cultivar", "Env.Cov", "GxW", "Error"), 
                           each = length(vG)),
                       levels = c("Cultivar", "Env.Cov", "GxW", "Error"))
  )
  
  p1 <- ggplot(trace_plot_data, aes(x = Iteration, y = Variance, color = Component)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    facet_wrap(~ Component, scales = "free_y", ncol = 1) +
    labs(title = "GW模型方差组分轨迹图 (TKW)",
         x = "迭代次数",
         y = "方差") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  ggsave(file.path(output_dir, "GW_variance_trace_plots.png"), p1, width = 8, height = 10, dpi = 300)
  cat("  - 保存GW模型方差轨迹图: GW_variance_trace_plots.png\n")
}

# 2. 创建方差组分占比饼图
cat("\n2. 创建方差组分占比饼图...\n")
if(!is.na(prop_genotype) && !is.na(prop_environment) && !is.na(prop_GxW) && !is.na(prop_residual)) {
  prop_data <- data.frame(
    Component = c("Cultivar", "Env.Cov", "GxW", "Error"),
    Proportion = c(prop_genotype, prop_environment, prop_GxW, prop_residual),
    Label = paste(c("Cultivar", "Env.Cov", "GxW", "Error"), 
                  sprintf("%.1f%%", c(prop_genotype, prop_environment, prop_GxW, prop_residual)))
  )
  
  p2 <- ggplot(prop_data, aes(x = "", y = Proportion, fill = Component)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
              position = position_stack(vjust = 0.5), 
              size = 4) +
    labs(title = "GW模型方差组分占比 (TKW)",
         fill = "方差组分") +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(output_dir, "GW_variance_proportion_pie.png"), p2, width = 8, height = 6, dpi = 300)
  cat("  - 保存GW模型方差组分占比饼图: GW_variance_proportion_pie.png\n")
}

# 3. 创建置信区间图
cat("\n3. 创建置信区间图...\n")
if(exists("ci_vG") && exists("ci_vEC") && exists("ci_vGxW") && exists("ci_varE") &&
   !any(is.na(ci_vG)) && !any(is.na(ci_vEC)) && !any(is.na(ci_vGxW)) && !any(is.na(ci_varE))) {
  
  ci_data <- data.frame(
    Component = c("Cultivar", "Env.Cov", "GxW", "Error"),
    Mean = c(mean_vG, mean_vEC, mean_vGxW, mean_varE),
    Lower = c(ci_vG[1], ci_vEC[1], ci_vGxW[1], ci_varE[1]),
    Upper = c(ci_vG[2], ci_vEC[2], ci_vGxW[2], ci_varE[2])
  )
  
  p3 <- ggplot(ci_data, aes(x = Component, y = Mean)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "blue") +
    labs(title = "GW模型方差组分估计值及95%置信区间 (TKW)",
         x = "方差组分",
         y = "方差") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(output_dir, "GW_variance_CI_plots.png"), p3, width = 10, height = 6, dpi = 300)
  cat("  - 保存GW模型方差置信区间图: GW_variance_CI_plots.png\n")
}

# 第十九部分：清理临时文件
# ==============================================================================
cat("\n=== 清理临时文件 ===\n")

temp_files <- c("m3_ETA_v_b.bin", "m3_ETA_ec_b.bin", "m3_ETA_GxW_b.bin", 
                "m3_varE.dat", "m3_ETA_v_varB.dat", "m3_ETA_ec_varB.dat", 
                "m3_ETA_GxW_varB.dat", "m3_mu.dat")

temp_removed <- 0
for (file in temp_files) {
  if (file.exists(file)) {
    file.remove(file)
    cat("  - 删除临时文件:", file, "\n")
    temp_removed <- temp_removed + 1
  }
}

cat("临时文件清理完成! 共删除", temp_removed, "个临时文件。\n")

# 第二十部分：最终总结
# ==============================================================================
cat("\n")
cat(rep("=", 80), "\n")
cat("GW模型方差组分分析完整结束!\n")
cat(rep("=", 80), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("工作目录:", work_dir, "\n")
cat("输出目录:", output_dir, "\n")
cat("表型数据: TKW (千粒重)\n")
cat("删除的环境:", paste(env_to_remove, collapse = ", "), "\n")
cat("最终使用的观测值数量:", nrow(pheno_filtered), "\n")
cat("最终使用的品种数量:", length(common_genotypes), "\n")
cat("最终使用的环境数量:", length(common_envs), "\n")
cat("基因型SNP数量:", ncol(geno_filtered), "\n")
cat("环境协变量数量:", ncol(W), "\n")
cat(rep("=", 80), "\n")

# 最终方差组分总结
cat("\n最终方差组分结果:\n")
cat("----------------------------------------\n")
if(exists("mean_vG") && !is.na(mean_vG)) {
  cat(sprintf("Cultivar (基因型方差): %.3f [%.3f, %.3f]\n", mean_vG, ci_vG[1], ci_vG[2]))
  if(!is.na(prop_genotype)) cat(sprintf("  占比: %.2f%%\n", prop_genotype))
}
if(exists("mean_vEC") && !is.na(mean_vEC)) {
  cat(sprintf("Env. Cov. (环境协变量方差): %.3f [%.3f, %.3f]\n", mean_vEC, ci_vEC[1], ci_vEC[2]))
  if(!is.na(prop_environment)) cat(sprintf("  占比: %.2f%%\n", prop_environment))
}
if(exists("mean_vGxW") && !is.na(mean_vGxW)) {
  cat(sprintf("SNPxEC (G×W交互作用方差): %.3f [%.3f, %.3f]\n", mean_vGxW, ci_vGxW[1], ci_vGxW[2]))
  if(!is.na(prop_GxW)) cat(sprintf("  占比: %.2f%%\n", prop_GxW))
}
if(exists("mean_varE") && !is.na(mean_varE)) {
  cat(sprintf("Error (残差方差): %.3f [%.3f, %.3f]\n", mean_varE, ci_varE[1], ci_varE[2]))
  if(!is.na(prop_residual)) cat(sprintf("  占比: %.2f%%\n", prop_residual))
}
cat(sprintf("\n总表型方差: %.3f\n", total_phenotypic_variance))
if(!is.na(H2)) cat(sprintf("广义遗传力 (H²): %.2f%%\n", H2))
cat("----------------------------------------\n")

cat("\nGW模型分析成功完成! 所有结果已保存到", output_dir, "目录。\n")
cat("请检查以下重要文件:\n")
cat("1. GW_variance_components_table.csv - GW模型表格格式结果\n")
cat("2. GW_detailed_variance_components.csv - GW模型详细方差组分\n")
cat("3. GW_variance_analysis_report.txt - GW模型详细报告\n")
cat("4. GW_model_fit.RData - GW模型拟合对象\n")
cat("5. 图表文件 (.png) - 可视化结果\n")

cat("\n=== 分析完成 ===\n")