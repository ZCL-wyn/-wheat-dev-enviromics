#!/usr/bin/env Rscript

# ===========================================
# FW参数计算和GWAS分析完整代码 - 完整8参数版本
# 计算全部8个FW参数，但只对矫正截距和斜率进行GWAS
# 配置：
# 1. 工作目录：/mnt/7t_storage/zhangcl/TKW
# 2. 基因型文件：myGD2.csv
# 3. 图谱文件：myGM2.csv
# 4. 环境数据：EC8.csv
# 5. 表型数据：TKW_mean_table.txt
# 6. 表型表头：genotypeenv_codeTKW
# 7. 性状：TKW（千粒重）
# 8. 环境因子：DTR&mean&PostFlowering1_14
# ===========================================

# 设置工作目录
setwd("/mnt/7t_storage/zhangcl/TKW")

cat("=== FW参数计算和GWAS分析开始 ===\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ===========================================
# 1. 初始设置和包加载
# ===========================================

cat("\n=== 初始设置和包加载 ===\n")

num_cores <- parallel::detectCores()
cat("检测到CPU核心数:", num_cores, "\n")
used_cores <- max(1, floor(num_cores * 0.5))
cat("将使用", used_cores, "个CPU核心进行并行计算\n")

cat("加载必要的R包...\n")
required_packages <- c(
  "data.table", "BGLR", "glmnet", "caret", "ggplot2",
  "vcfR", "rrBLUP", "parallel", "doParallel", "foreach", 
  "Matrix", "dplyr", "tidyr", "reshape2", "qqman"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装包:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
    cat("成功加载包:", pkg, "\n")
  } else {
    cat("包已加载:", pkg, "\n")
  }
}

# 检查并安装BGLR包
if (!require("BGLR", character.only = TRUE, quietly = TRUE)) {
  cat("安装BGLR包...\n")
  install.packages("BGLR", dependencies = TRUE, repos = "https://cloud.r-project.org")
  library("BGLR")
  cat("成功加载BGLR包\n")
} else {
  cat("BGLR包已加载\n")
}

cat("注册并行后端...\n")
cl <- makeCluster(used_cores)
doParallel::registerDoParallel(cl)
cat("并行后端注册完成\n")

# ===========================================
# 2. 加载GAPIT函数
# ===========================================

cat("\n=== 加载GAPIT函数 ===\n")

cat("从网络加载GAPIT...\n")
tryCatch({
  source("http://zzlab.net/GAPIT/GAPIT.library.R")
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  cat("GAPIT函数加载成功\n")
}, error = function(e) {
  cat("网络GAPIT加载失败:", e$message, "\n")
  if (file.exists("GAPIT.library.R")) {
    source("GAPIT.library.R")
    cat("从本地加载GAPIT成功\n")
  } else {
    stop("无法加载GAPIT函数")
  }
})

# ===========================================
# 3. 数据准备
# ===========================================

cat("\n=== 数据准备 ===\n")

set.seed(195021)

# 创建输出目录
output_dir <- "FW_GWAS_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 读取表型数据
cat("读取表型数据...\n")
phenotype_data <- read.table("TKW_mean_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("原始表型数据维度:", dim(phenotype_data), "\n")

# 自动检测表型数据中的基因型ID列名
pheno_colnames <- colnames(phenotype_data)
cat("表型数据列名:", paste(pheno_colnames, collapse = ", "), "\n")

# 检查可能的基因型ID列名
geno_id_names <- c("line_code", "genotype", "taxa", "Line", "Accession", "ID", "Variety")
geno_id_col <- NULL

for (col_name in geno_id_names) {
  if (col_name %in% pheno_colnames) {
    geno_id_col <- col_name
    cat("检测到基因型ID列名:", geno_id_col, "\n")
    break
  }
}

if (is.null(geno_id_col)) {
  # 如果没有找到标准列名，假设第一列是基因型ID
  geno_id_col <- colnames(phenotype_data)[1]
  cat("未找到标准基因型ID列名，使用第一列:", geno_id_col, "\n")
}

# 重命名列以统一处理
colnames(phenotype_data)[colnames(phenotype_data) == geno_id_col] <- "genotype"
cat("统一后的表型数据列名:", colnames(phenotype_data), "\n")

# 找到TKW列
trait_col <- NULL
trait_names <- c("TKW", "tkw", "weight", "Value", "Phenotype")
for (col_name in trait_names) {
  if (col_name %in% pheno_colnames) {
    trait_col <- col_name
    cat("检测到性状列名:", trait_col, "\n")
    break
  }
}

if (is.null(trait_col)) {
  # 假设最后一列是性状值
  trait_col <- colnames(phenotype_data)[ncol(phenotype_data)]
  cat("未找到标准性状列名，使用最后一列:", trait_col, "\n")
}

# 确保性状列名正确
if (trait_col != "TKW") {
  colnames(phenotype_data)[colnames(phenotype_data) == trait_col] <- "TKW"
}

# 找到环境列
env_col <- NULL
env_names <- c("env_code", "Env", "Environment", "Location", "Site")
for (col_name in env_names) {
  if (col_name %in% pheno_colnames) {
    env_col <- col_name
    cat("检测到环境列名:", env_col, "\n")
    break
  }
}

if (is.null(env_col)) {
  # 假设第二列是环境
  env_col <- colnames(phenotype_data)[2]
  cat("未找到标准环境列名，使用第二列:", env_col, "\n")
}

# 确保环境列名正确
if (env_col != "env_code") {
  colnames(phenotype_data)[colnames(phenotype_data) == env_col] <- "env_code"
}

cat("最终表型数据结构:\n")
print(head(phenotype_data))

# 数据清洗
phenotype_data_clean <- phenotype_data[!is.na(phenotype_data$TKW), ]
cat("去除NA值后的表型数据维度:", dim(phenotype_data_clean), "\n")

# 读取环境参数
ecs_data <- read.csv("EC8.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(ecs_data)[1] <- "env_code"

# 环境过滤（根据需求调整）
env_to_remove <- c()  # 可根据需要添加要过滤的环境
if(length(env_to_remove) > 0){
  phenotype_data_filtered <- phenotype_data_clean[!phenotype_data_clean$env_code %in% env_to_remove, ]
  cat("环境过滤后的表型数据维度:", dim(phenotype_data_filtered), "\n")
} else {
  phenotype_data_filtered <- phenotype_data_clean
  cat("无环境过滤，使用全部表型数据\n")
}

# 读取基因型数据
cat("读取基因型数据...\n")
geno_dt <- fread("myGD2.csv", header = TRUE, data.table = FALSE)

# 自动检测基因型数据中的ID列名
geno_dt_colnames <- colnames(geno_dt)
cat("基因型数据列名（前5个）:", paste(head(geno_dt_colnames, 5), collapse = ", "), "\n")

# 基因型数据可能的ID列名
geno_data_id_names <- c("taxa", "line_code", "genotype", "ID", "Sample", "Accession")
geno_data_id_col <- NULL

# 检查第一列是否是基因型ID
first_col_values <- geno_dt[[1]]
first_col_name <- colnames(geno_dt)[1]

# 如果第一列看起来是基因型ID（不是SNP名字，且值类似品种名称）
if (first_col_name %in% geno_data_id_names) {
  geno_data_id_col <- first_col_name
  cat("检测到基因型ID列名（基于列名）:", geno_data_id_col, "\n")
} else {
  # 检查第一列的值是否看起来像品种名称（包含字母和数字的组合）
  sample_pattern <- sum(grepl("^[A-Za-z]+[0-9]*$", first_col_values[1:10])) > 5
  if (sample_pattern) {
    geno_data_id_col <- first_col_name
    cat("检测到基因型ID列名（基于值模式）:", geno_data_id_col, "\n")
  } else {
    # 遍历可能的列名
    for (col_name in geno_data_id_names) {
      if (col_name %in% geno_dt_colnames) {
        geno_data_id_col <- col_name
        cat("检测到基因型ID列名:", geno_data_id_col, "\n")
        break
      }
    }
  }
}

if (is.null(geno_data_id_col)) {
  # 如果都没找到，假设第一列是基因型ID
  geno_data_id_col <- colnames(geno_dt)[1]
  cat("未检测到标准基因型ID列名，使用第一列:", geno_data_id_col, "\n")
}

# 提取样本名称和基因型矩阵
sample_names <- geno_dt[[geno_data_id_col]]
geno_matrix <- as.matrix(geno_dt[, -which(colnames(geno_dt) == geno_data_id_col)])
rownames(geno_matrix) <- sample_names
cat("基因型矩阵维度:", dim(geno_matrix), "\n")
cat("基因型样本示例（前10个）:", head(sample_names, 10), "\n")

# 读取SNP图谱
cat("读取SNP图谱...\n")
map_data <- fread("myGM2.csv", header = TRUE, data.table = FALSE)

# 标准化SNP图谱列名
map_colnames <- colnames(map_data)
cat("SNP图谱原始列名:", paste(map_colnames, collapse = ", "), "\n")

# 可能的列名映射
map_name_mapping <- list(
  SNP = c("SNP", "Marker", "rs", "ID"),
  Chromosome = c("Chromosome", "Chr", "CHROM"),
  Position = c("Position", "Pos", "BP")
)

# 重命名列
for (target_name in names(map_name_mapping)) {
  possible_names <- map_name_mapping[[target_name]]
  for (possible_name in possible_names) {
    if (possible_name %in% map_colnames) {
      colnames(map_data)[colnames(map_data) == possible_name] <- target_name
      cat("将列", possible_name, "重命名为", target_name, "\n")
      break
    }
  }
}

# 如果缺少必需的列，尝试使用位置推断
if (!"SNP" %in% colnames(map_data)) {
  # 假设第一列是SNP名称
  colnames(map_data)[1] <- "SNP"
  cat("假设第一列为SNP名称\n")
}

if (!"Chromosome" %in% colnames(map_data)) {
  # 尝试找到染色体列
  for (i in 1:ncol(map_data)) {
    if (i != which(colnames(map_data) == "SNP")) {
      col_vals <- map_data[[i]]
      # 检查是否看起来像染色体编号
      if (is.character(col_vals) || is.numeric(col_vals)) {
        unique_vals <- unique(col_vals[1:100])
        if (length(unique_vals) <= 30 && all(grepl("^[0-9]+[A-Z]?$", as.character(unique_vals)))) {
          colnames(map_data)[i] <- "Chromosome"
          cat("假设列", i, "为染色体列\n")
          break
        }
      }
    }
  }
}

if (!"Position" %in% colnames(map_data)) {
  # 假设最后一列是位置
  colnames(map_data)[ncol(map_data)] <- "Position"
  cat("假设最后一列为位置列\n")
}

cat("SNP图谱数据维度:", dim(map_data), "\n")
cat("标准化后的SNP图谱列名:", colnames(map_data), "\n")

# 环境参数
kPara_Name <- "DTR&mean&PostFlowering1_14"
envMeanPara <- ecs_data[, c("env_code", kPara_Name)]
colnames(envMeanPara)[2] <- "kPara"
cat("选择的环境参数:", kPara_Name, "\n")

# 计算环境均值
env_mean_trait <- aggregate(TKW ~ env_code, data = phenotype_data_filtered, mean, na.rm = TRUE)
colnames(env_mean_trait)[2] <- "meanY"

# 清洗品种名称以确保匹配
clean_genotype_names <- function(names) {
  # 去除首尾空格
  names <- trimws(names)
  # 统一大小写（全部大写）
  names <- toupper(names)
  # 去除所有非字母数字字符，只保留字母和数字
  names <- gsub("[^A-Za-z0-9]", "", names)
  # 确保以字母开头
  names <- gsub("^([0-9]+)", "G\\1", names)  # 如果以数字开头，加G前缀
  return(names)
}

cat("\n=== 清洗和匹配品种名称 ===\n")

# 清洗表型数据中的品种名称
phenotype_data_filtered$genotype_clean <- clean_genotype_names(phenotype_data_filtered$genotype)
pheno_lines <- unique(phenotype_data_filtered$genotype_clean)
cat("表型数据清洗后品种数量:", length(pheno_lines), "\n")
cat("表型数据品种示例（前10个）:", head(pheno_lines, 10), "\n")

# 清洗基因型数据中的品种名称
geno_lines_clean <- clean_genotype_names(rownames(geno_matrix))
rownames(geno_matrix) <- geno_lines_clean
geno_lines <- geno_lines_clean
cat("基因型数据清洗后品种数量:", length(geno_lines), "\n")
cat("基因型数据品种示例（前10个）:", head(geno_lines, 10), "\n")

# 数据对齐
common_lines <- intersect(pheno_lines, geno_lines)
cat("共同品种数量:", length(common_lines), "\n")

if (length(common_lines) == 0) {
  cat("警告：没有共同品种！\n")
  
  # 尝试部分匹配
  cat("尝试部分匹配...\n")
  potential_matches <- list()
  
  for (pheno_line in pheno_lines[1:min(20, length(pheno_lines))]) {
    # 尝试多种匹配模式
    matches <- c()
    
    # 精确匹配（去除数字前的0）
    pheno_line_adj <- gsub("^G0+", "G", pheno_line)  # 去除G后面的0
    geno_line_adjs <- gsub("^G0+", "G", geno_lines)
    exact_match <- geno_lines[geno_line_adjs == pheno_line_adj]
    
    # 部分匹配（包含关系）
    partial_matches <- geno_lines[grep(pheno_line, geno_lines, ignore.case = TRUE)]
    
    # 反向部分匹配
    reverse_matches <- geno_lines[sapply(geno_lines, function(x) grepl(x, pheno_line, ignore.case = TRUE))]
    
    matches <- unique(c(exact_match, partial_matches, reverse_matches))
    
    if (length(matches) > 0) {
      potential_matches[[pheno_line]] <- matches
      cat("表型品种", pheno_line, "可能匹配基因型品种:", paste(matches, collapse = ", "), "\n")
    }
  }
  
  # 如果有潜在匹配，让用户选择
  if (length(potential_matches) > 0) {
    cat("\n找到", length(potential_matches), "个潜在匹配对\n")
    
    # 创建映射表（取第一个匹配）
    name_mapping <- data.frame(
      pheno_name = names(potential_matches),
      geno_name = sapply(potential_matches, function(x) x[1]),
      stringsAsFactors = FALSE
    )
    
    # 应用映射
    phenotype_data_filtered$genotype_clean <- 
      name_mapping$geno_name[match(phenotype_data_filtered$genotype_clean, name_mapping$pheno_name)]
    
    # 更新共同品种
    pheno_lines <- unique(phenotype_data_filtered$genotype_clean)
    common_lines <- intersect(pheno_lines, geno_lines)
    cat("应用映射后的共同品种数量:", length(common_lines), "\n")
    
    if (length(common_lines) == 0) {
      stop("错误：应用映射后仍然没有共同品种，请手动检查数据！")
    }
  } else {
    stop("错误：没有找到任何匹配的品种，请检查数据！")
  }
}

# 过滤数据
phenotype_data_final <- phenotype_data_filtered[phenotype_data_filtered$genotype_clean %in% common_lines, ]
geno_matrix_final <- geno_matrix[common_lines, ]
cat("最终表型数据维度:", dim(phenotype_data_final), "\n")
cat("最终基因型数据维度:", dim(geno_matrix_final), "\n")

# 检查SNP一致性
geno_snps <- colnames(geno_matrix_final)
map_snps <- map_data$SNP
common_snps <- intersect(geno_snps, map_snps)
cat("共同SNP数量:", length(common_snps), "\n")

if (length(common_snps) == 0) {
  # 如果SNP名称不匹配，尝试清洗SNP名称
  cat("SNP名称不匹配，尝试清洗SNP名称...\n")
  
  clean_snp_names <- function(names) {
    names <- gsub("[^[:alnum:]_]", "_", names)  # 替换特殊字符为下划线
    names <- gsub("_+", "_", names)  # 合并多个下划线
    names <- gsub("^_|_$", "", names)  # 去除首尾下划线
    return(names)
  }
  
  geno_snps_clean <- clean_snp_names(geno_snps)
  map_snps_clean <- clean_snp_names(map_snps)
  
  colnames(geno_matrix_final) <- geno_snps_clean
  map_data$SNP <- map_snps_clean
  
  common_snps <- intersect(geno_snps_clean, map_snps_clean)
  cat("清洗后的共同SNP数量:", length(common_snps), "\n")
}

# 保留共同SNP
geno_matrix_final <- geno_matrix_final[, common_snps, drop = FALSE]
map_data_final <- map_data[map_data$SNP %in% common_snps, ]
map_data_final <- map_data_final[match(colnames(geno_matrix_final), map_data_final$SNP), ]

# 质量控制
cat("\n=== 基因型数据质量控制 ===\n")

MAF.min <- 0.03
pNA.max <- 0.1

geno_matrix_qc <- geno_matrix_final

# MAF过滤
MAF_values <- apply(geno_matrix_qc, 2, function(x) mean(x, na.rm = TRUE) / 2)
drop_maf <- which(MAF_values < MAF.min | MAF_values > 1 - MAF.min)
if(length(drop_maf) > 0) {
    geno_matrix_qc <- geno_matrix_qc[, -drop_maf, drop = FALSE]
    map_data_final <- map_data_final[-drop_maf, ]
    cat("MAF过滤: 剔除了", length(drop_maf), "个位点\n")
}

# 缺失率过滤
missing_rates <- apply(geno_matrix_qc, 2, function(x) mean(is.na(x)))
drop_na <- which(missing_rates > pNA.max)
if(length(drop_na) > 0) {
    geno_matrix_qc <- geno_matrix_qc[, -drop_na, drop = FALSE]
    map_data_final <- map_data_final[-drop_na, ]
    cat("缺失率过滤: 剔除了", length(drop_na), "个位点\n")
}

cat("质量控制后基因型数据维度:", dim(geno_matrix_qc), "\n")

# 缺失值插补
cat("进行缺失值插补...\n")
for(j in 1:ncol(geno_matrix_qc)) {
    tmp <- geno_matrix_qc[, j]
    if(any(is.na(tmp))) {
        geno_matrix_qc[, j] <- ifelse(is.na(tmp), mean(tmp, na.rm = TRUE), tmp)
    }
}
cat("缺失值插补完成\n")

# ===========================================
# 4. 准备GAPIT格式数据
# ===========================================

cat("\n=== 准备GAPIT格式数据 ===\n")

# 标准化标记名称（进一步清洗）
standardize_snp_names <- function(names) {
  names <- gsub(";", "_", names)
  names <- gsub("-", "_", names)
  names <- gsub(":", "_", names)
  names <- gsub("/", "_", names)
  names <- gsub("\\|", "_", names)
  names <- gsub("\\(", "_", names)
  names <- gsub("\\)", "_", names)
  names <- gsub("\\[", "_", names)
  names <- gsub("\\]", "_", names)
  names <- gsub(" ", "_", names)
  names <- gsub("\\.", "_", names)
  names <- gsub("__+", "_", names)  # 合并多个下划线
  names <- gsub("^_|_$", "", names)  # 去除首尾下划线
  return(names)
}

new_geno_snp_names <- standardize_snp_names(colnames(geno_matrix_qc))
colnames(geno_matrix_qc) <- new_geno_snp_names
new_map_snp_names <- standardize_snp_names(map_data_final$SNP)
map_data_final$SNP <- new_map_snp_names

# 去除重复标记
keep_indices <- !duplicated(new_geno_snp_names)
geno_matrix_qc <- geno_matrix_qc[, keep_indices]
map_data_final <- map_data_final[keep_indices, ]

# 准备GD和GM
GD_all <- data.frame(taxa = rownames(geno_matrix_qc), geno_matrix_qc)
rownames(GD_all) <- GD_all$taxa
GM_all <- map_data_final
colnames(GM_all) <- c("SNP", "Chromosome", "Position")

# 将染色体转换为数值
GM_all$Chromosome <- as.numeric(as.character(GM_all$Chromosome))
GM_all$Position <- as.numeric(as.character(GM_all$Position))

cat("GAPIT GD数据维度:", dim(GD_all), "\n")
cat("GAPIT GM数据维度:", dim(GM_all), "\n")
cat("染色体分布:\n")
print(table(GM_all$Chromosome))

# ===========================================
# 5. 核心函数定义
# ===========================================

cat("\n=== 定义核心函数 ===\n")

# 5.1 FW参数计算函数 - 完整8参数版本
calculate_FW_parameters_full <- function(pheno_data, env_factors, env_means) {
  cat("计算完整的FW模型参数（8个参数）...\n")
  
  # 使用清洗后的品种名称
  colnames(pheno_data)[colnames(pheno_data) == "genotype_clean"] <- "line_code"
  colnames(pheno_data)[colnames(pheno_data) == "TKW"] <- "PH"
  
  line_codes <- unique(pheno_data$line_code)
  cat("品种数量:", length(line_codes), "\n")
  
  # 获取所有环境
  all_envs <- unique(pheno_data$env_code)
  cat("环境数量:", length(all_envs), "\n")
  
  # 计算环境参数平均值
  all_env_factors <- env_factors[env_factors$env_code %in% all_envs, ]
  mean_kPara <- mean(all_env_factors$kPara, na.rm = TRUE)
  cat("环境参数平均值:", mean_kPara, "\n")
  
  # 计算环境均值平均值
  all_env_means <- env_means[env_means$env_code %in% all_envs, ]
  mean_meanY <- mean(all_env_means$meanY, na.rm = TRUE)
  cat("环境均值平均值:", mean_meanY, "\n")

  results <- data.frame()
  
  for (line in line_codes) {
    line_data <- subset(pheno_data, line_code == line)
    
    if (nrow(line_data) >= 3) {
      # 合并当前品种的表型数据、环境参数数据和环境均值数据
      line_data_merged <- merge(line_data, env_factors, by = "env_code")
      line_data_merged <- merge(line_data_merged, env_means, by = "env_code")
      
      # 移除包含NA值的行
      line_data_clean <- line_data_merged[complete.cases(line_data_merged$PH, line_data_merged$meanY, line_data_merged$kPara), ]
      
      # 检查清理后的数据是否足够
      if (nrow(line_data_clean) >= 2) {
        tryCatch({
          # 1. 计算基于环境均值的线性回归模型
          lm_mean <- lm(PH ~ meanY, data = line_data_clean)
          
          # Intcp_mean = a_mean: 在环境均值平均值处的预测值
          a_mean <- as.numeric(predict(lm_mean, data.frame(meanY = mean_meanY)))
          
          # Slope_mean = b_mean: 回归斜率
          b_mean <- as.numeric(coef(lm_mean)[2])
          
          # R2_mean = r2_mean
          r2_mean <- summary(lm_mean)$r.squared
          
          # 2. 计算基于环境参数的线性回归模型
          lm_para <- lm(PH ~ kPara, data = line_data_clean)
          
          # Intcp_para_adj = a_para_adj: 在环境参数平均值处的预测值
          a_para_adj <- as.numeric(predict(lm_para, data.frame(kPara = mean_kPara)))
          
          # Intcp_para = a_para: 回归截距
          a_para <- as.numeric(coef(lm_para)[1])
          
          # Slope_para = b_para: 回归斜率
          b_para <- as.numeric(coef(lm_para)[2])
          
          # R2_para = r2_para
          r2_para <- summary(lm_para)$r.squared

          # 将结果添加到结果矩阵
          results <- rbind(results, data.frame(
            line_code = line,
            Intcp_mean = round(a_mean, 4),
            Slope_mean = round(b_mean, 4),
            Intcp_para_adj = round(a_para_adj, 4),
            Intcp_para = round(a_para, 4),
            Slope_para = round(b_para, 4),
            R2_mean = round(r2_mean, 4),
            R2_para = round(r2_para, 4),
            stringsAsFactors = FALSE
          ))
        }, error = function(e) {
          cat("  品种", line, "的FW参数计算失败:", e$message, "\n")
        })
      }
    }
  }

  cat("成功计算FW参数的品种数量:", nrow(results), "\n")
  
  return(list(
    FW_params = results,
    mean_kPara = mean_kPara,
    mean_meanY = mean_meanY
  ))
}

# 5.2 GAPIT GWAS函数
perform_gapit_gwas <- function(FW_params, GD_all, GM_all, trait, output_dir) {
  cat("对", trait, "进行GAPIT GWAS分析...\n")
  
  if (!is.data.frame(FW_params) || nrow(FW_params) == 0) {
    cat("错误：FW_params 无效\n")
    return(NULL)
  }
  
  # 准备表型数据
  pheno_df <- data.frame(
    Taxa = FW_params$line_code,
    Trait = FW_params[[trait]]
  )
  pheno_df <- pheno_df[complete.cases(pheno_df$Trait), ]
  
  cat(trait, "有效样本数量:", nrow(pheno_df), "\n")
  
  if (nrow(pheno_df) < 10) {
    cat("  有效样本不足，跳过", trait, "的GWAS\n")
    return(NULL)
  }
  
  # 共同基因型
  common_genotypes <- intersect(pheno_df$Taxa, GD_all$taxa)
  cat("共同基因型数量:", length(common_genotypes), "\n")
  
  if (length(common_genotypes) < 10) {
    cat("  共同基因型不足，跳过", trait, "的GWAS\n")
    return(NULL)
  }
  
  # 准备GAPIT输入
  pheno_gapit <- pheno_df[pheno_df$Taxa %in% common_genotypes, ]
  GD_gapit <- GD_all[GD_all$taxa %in% common_genotypes, ]
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  original_wd <- getwd()
  
  tryCatch({
    setwd(output_dir)
    
    cat("运行GAPIT Blink分析...\n")
    # 使用BLINK模型，PCA.total=2
    gapit_result <- GAPIT(
      Y = pheno_gapit,
      G = NULL,
      GD = GD_gapit,
      GM = GM_all,
      model = "Blink",
      PCA.total = 2,
      Multiple_analysis = FALSE,
      file.output = TRUE,
      Major.allele.zero = FALSE,
      SNP.MAF = 0.01,
      cutOff = 1,
      memo = paste("GWAS", trait, sep = "_")
    )
    
    setwd(original_wd)
    
    # 读取GWAS结果文件
    gwas_pattern <- paste0("^GAPIT\\.Association\\.GWAS_Results\\..*", trait, ".*\\.csv$")
    gwas_files <- list.files(output_dir, pattern = gwas_pattern, full.names = TRUE)
    
    if (length(gwas_files) == 0) {
      cat("警告：未找到GWAS结果文件\n")
      gwas_results_data <- NULL
    } else {
      cat("找到GWAS结果文件:", gwas_files[1], "\n")
      gwas_results_data <- read.csv(gwas_files[1], stringsAsFactors = FALSE)
      cat("GWAS结果文件读取成功，行数:", nrow(gwas_results_data), "\n")
    }
    
    return(list(
      gwas_results = gwas_results_data,
      GD = GD_gapit,
      GM = GM_all,
      phenotype_vector = pheno_gapit$Trait,
      taxa_names = pheno_gapit$Taxa,
      output_dir = output_dir
    ))
    
  }, error = function(e) {
    setwd(original_wd)
    cat("GAPIT分析失败:", e$message, "\n")
    return(NULL)
  })
}

# 5.3 显著标记选择函数
select_significant_markers <- function(gwas_results, GD_all, GM_all, 
                                       n_top_markers = 15, 
                                       p_threshold = NULL) {
  cat("从GWAS结果中选择显著标记...\n")
  
  if (is.null(gwas_results) || nrow(gwas_results) == 0) {
    cat("没有GWAS结果\n")
    return(character(0))
  }
  
  cat("GWAS结果行数:", nrow(gwas_results), "\n")
  
  # 确定P值列名
  pvalue_col <- NULL
  possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
  
  for (col in possible_pvalue_cols) {
    if (col %in% colnames(gwas_results)) {
      pvalue_col <- col
      break
    }
  }
  
  if (is.null(pvalue_col)) {
    cat("警告：未找到P值列，使用所有标记\n")
    if (nrow(gwas_results) <= n_top_markers) {
      return(gwas_results$SNP)
    } else {
      return(gwas_results$SNP[1:n_top_markers])
    }
  }
  
  # 使用0.05/10000作为Bonferroni校正阈值
  if (is.null(p_threshold)) {
    p_threshold <- 0.05 / 10000
    cat("Bonferroni校正p阈值(0.05/10000):", p_threshold, "\n")
  }
  
  # 筛选显著标记
  significant_markers <- gwas_results[gwas_results[[pvalue_col]] < p_threshold, ]
  cat("显著标记数量（p <", p_threshold, "）:", nrow(significant_markers), "\n")
  
  if (nrow(significant_markers) == 0) {
    cat("没有显著标记，使用所有标记中P值最小的前", n_top_markers, "个\n")
    gwas_results_sorted <- gwas_results[order(gwas_results[[pvalue_col]]), ]
    top_markers <- gwas_results_sorted$SNP[1:min(n_top_markers, nrow(gwas_results_sorted))]
    return(top_markers)
  }
  
  # 按P值升序排序
  significant_markers <- significant_markers[order(significant_markers[[pvalue_col]]), ]
  
  # 取前n_top_markers个标记
  n_to_select <- min(n_top_markers, nrow(significant_markers))
  top_markers <- significant_markers$SNP[1:n_to_select]
  
  cat("选择", n_to_select, "个显著标记\n")
  cat("前5个标记:", head(top_markers, 5), "\n")
  
  return(top_markers)
}

# ===========================================
# 6. 执行FW参数计算和GWAS分析
# ===========================================

cat("\n=== 执行FW参数计算和GWAS分析 ===\n")

# 6.1 计算完整的FW参数（8个参数）
cat("第一步：计算完整的FW参数（8个参数）...\n")
fw_result <- calculate_FW_parameters_full(phenotype_data_final, envMeanPara, env_mean_trait)

FW_params <- fw_result$FW_params
mean_kPara <- fw_result$mean_kPara
mean_meanY <- fw_result$mean_meanY

if (nrow(FW_params) < 5) {
  cat("FW参数不足，无法进行GWAS分析\n")
  stop("FW参数不足")
}

# 保存完整的FW参数
fw_file <- file.path(output_dir, "FW_parameters_full.csv")
write.csv(FW_params, fw_file, row.names = FALSE)
cat("完整FW参数已保存到:", fw_file, "\n")
cat("FW参数计算完成，品种数量:", nrow(FW_params), "\n")

# 输出参数统计摘要
cat("\nFW参数统计摘要:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat("环境参数平均值 (mean_kPara):", mean_kPara, "\n")
cat("环境均值平均值 (mean_meanY):", mean_meanY, "\n")
cat("\n各参数平均值:\n")
cat("Intcp_mean (基于环境均值的截距):", mean(FW_params$Intcp_mean, na.rm = TRUE), "\n")
cat("Slope_mean (基于环境均值的斜率):", mean(FW_params$Slope_mean, na.rm = TRUE), "\n")
cat("Intcp_para_adj (基于环境参数的矫正截距):", mean(FW_params$Intcp_para_adj, na.rm = TRUE), "\n")
cat("Intcp_para (基于环境参数的原始截距):", mean(FW_params$Intcp_para, na.rm = TRUE), "\n")
cat("Slope_para (基于环境参数的斜率):", mean(FW_params$Slope_para, na.rm = TRUE), "\n")
cat("R2_mean (基于环境均值回归的R2):", mean(FW_params$R2_mean, na.rm = TRUE), "\n")
cat("R2_para (基于环境参数回归的R2):", mean(FW_params$R2_para, na.rm = TRUE), "\n")

# 6.2 只对矫正截距和斜率进行GWAS
cat("\n第二步：只对矫正截距和斜率进行GWAS分析...\n")

# 只选择Intcp_para_adj和Slope_para进行GWAS
traits_for_gwas <- c("Intcp_para_adj", "Slope_para")
gwas_results <- list()
significant_markers_list <- list()

for (trait in traits_for_gwas) {
  cat("对", trait, "进行GWAS分析...\n")
  gwas_dir <- file.path(output_dir, paste0("GWAS_", trait))
  gwas_results[[trait]] <- perform_gapit_gwas(FW_params, GD_all, GM_all, trait, gwas_dir)
  
  if (!is.null(gwas_results[[trait]])) {
    # 从GWAS结果中选择显著标记
    cat("从", trait, "的GWAS结果中选择显著标记...\n")
    significant_markers <- select_significant_markers(
      gwas_results[[trait]]$gwas_results, 
      GD_all,
      GM_all,
      n_top_markers = 15
    )
    significant_markers_list[[trait]] <- significant_markers
    cat(trait, "选择的显著标记数量:", length(significant_markers), "\n")
    
    # 保存显著标记列表
    write.table(data.frame(SNP = significant_markers),
               file.path(gwas_dir, paste0("Significant_Markers_", trait, ".txt")),
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # 保存完整的GWAS结果
    if (!is.null(gwas_results[[trait]]$gwas_results)) {
      gwas_result_file <- file.path(gwas_dir, paste0("Complete_GWAS_Results_", trait, ".csv"))
      write.csv(gwas_results[[trait]]$gwas_results, gwas_result_file, row.names = FALSE)
      cat("完整GWAS结果已保存到:", gwas_result_file, "\n")
    }
  } else {
    cat(trait, "GWAS失败\n")
    significant_markers_list[[trait]] <- character(0)
  }
}

# 6.3 生成GWAS结果汇总
cat("\n第三步：生成GWAS结果汇总...\n")

# 创建汇总表格
summary_data <- data.frame()

for (trait in traits_for_gwas) {
  if (!is.null(gwas_results[[trait]]) && !is.null(gwas_results[[trait]]$gwas_results)) {
    gwas_data <- gwas_results[[trait]]$gwas_results
    
    # 确定P值列名
    pvalue_col <- NULL
    possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
    
    for (col in possible_pvalue_cols) {
      if (col %in% colnames(gwas_data)) {
        pvalue_col <- col
        break
      }
    }
    
    if (!is.null(pvalue_col)) {
      # 计算显著标记数量（不同阈值）
      p_threshold_1 <- 0.05 / nrow(gwas_data)  # Bonferroni校正
      p_threshold_2 <- 1e-5
      p_threshold_3 <- 1e-4
      
      sig_count_1 <- sum(gwas_data[[pvalue_col]] < p_threshold_1, na.rm = TRUE)
      sig_count_2 <- sum(gwas_data[[pvalue_col]] < p_threshold_2, na.rm = TRUE)
      sig_count_3 <- sum(gwas_data[[pvalue_col]] < p_threshold_3, na.rm = TRUE)
      
      summary_data <- rbind(summary_data, data.frame(
        Trait = trait,
        Total_SNPs = nrow(gwas_data),
        Significant_SNPs_Bonferroni = sig_count_1,
        Significant_SNPs_1e_5 = sig_count_2,
        Significant_SNPs_1e_4 = sig_count_3,
        Top_SNPs_Selected = length(significant_markers_list[[trait]]),
        Mean_P_value = mean(gwas_data[[pvalue_col]], na.rm = TRUE),
        Min_P_value = min(gwas_data[[pvalue_col]], na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(summary_data) > 0) {
  summary_file <- file.path(output_dir, "GWAS_Summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  cat("GWAS汇总已保存到:", summary_file, "\n")
  
  # 打印汇总
  cat("\nGWAS分析汇总:\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  print(summary_data)
}

# 6.4 创建曼哈顿图和QQ图
cat("\n第四步：创建可视化图表...\n")

tryCatch({
  library(ggplot2)
  library(qqman)
  
  for (trait in traits_for_gwas) {
    if (!is.null(gwas_results[[trait]]) && !is.null(gwas_results[[trait]]$gwas_results)) {
      gwas_data <- gwas_results[[trait]]$gwas_results
      
      # 确定P值列名
      pvalue_col <- NULL
      possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
      
      for (col in possible_pvalue_cols) {
        if (col %in% colnames(gwas_data)) {
          pvalue_col <- col
          break
        }
      }
      
      if (!is.null(pvalue_col) && "Chromosome" %in% colnames(gwas_data) && "Position" %in% colnames(gwas_data)) {
        # 准备曼哈顿图数据
        manhattan_data <- data.frame(
          SNP = gwas_data$SNP,
          CHR = as.numeric(gwas_data$Chromosome),
          BP = as.numeric(gwas_data$Position),
          P = as.numeric(gwas_data[[pvalue_col]])
        )
        
        manhattan_data <- manhattan_data[complete.cases(manhattan_data), ]
        
        if (nrow(manhattan_data) > 0) {
          # 曼哈顿图
          png_file <- file.path(output_dir, paste0("Manhattan_Plot_", trait, ".png"))
          png(png_file, width = 1200, height = 600)
          manhattan(manhattan_data, main = paste("Manhattan Plot -", trait))
          dev.off()
          cat("曼哈顿图已保存:", png_file, "\n")
          
          # QQ图
          qq_file <- file.path(output_dir, paste0("QQ_Plot_", trait, ".png"))
          png(qq_file, width = 600, height = 600)
          qq(manhattan_data$P, main = paste("QQ Plot -", trait))
          dev.off()
          cat("QQ图已保存:", qq_file, "\n")
        }
      }
    }
  }
}, error = function(e) {
  cat("可视化图表创建失败:", e$message, "\n")
})

# ===========================================
# 7. 结果保存和总结
# ===========================================

cat("\n=== 结果保存和总结 ===\n")

# 保存所有重要结果到R数据文件
save_data <- list(
  FW_params = FW_params,
  mean_kPara = mean_kPara,
  mean_meanY = mean_meanY,
  gwas_results = gwas_results,
  significant_markers = significant_markers_list,
  phenotype_data = phenotype_data_final,
  env_data = envMeanPara,
  env_means = env_mean_trait,
  GD_all = GD_all,
  GM_all = GM_all,
  analysis_time = Sys.time()
)

save_file <- file.path(output_dir, "FW_GWAS_Analysis_Results.RData")
save(save_data, file = save_file)
cat("所有分析结果已保存到R数据文件:", save_file, "\n")

# 生成分析报告
report_file <- file.path(output_dir, "Analysis_Report.txt")
sink(report_file)
cat("FW参数计算和GWAS分析报告\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("工作目录:", getwd(), "\n")
cat("输出目录:", output_dir, "\n")
cat("环境参数:", kPara_Name, "\n")
cat("\n")

cat("数据概览:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("表型数据品种数:", length(unique(phenotype_data_final$genotype_clean)), "\n")
cat("表型数据环境数:", length(unique(phenotype_data_final$env_code)), "\n")
cat("表型数据记录数:", nrow(phenotype_data_final), "\n")
cat("基因型数据品种数:", nrow(geno_matrix_qc), "\n")
cat("基因型数据SNP数:", ncol(geno_matrix_qc), "\n")
cat("环境参数平均值 (mean_kPara):", mean_kPara, "\n")
cat("环境均值平均值 (mean_meanY):", mean_meanY, "\n")
cat("\n")

cat("FW参数计算（8个参数）:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("成功计算FW参数的品种数:", nrow(FW_params), "\n")
cat("Intcp_mean (基于环境均值的截距)平均值:", mean(FW_params$Intcp_mean, na.rm = TRUE), "\n")
cat("Intcp_mean标准差:", sd(FW_params$Intcp_mean, na.rm = TRUE), "\n")
cat("Slope_mean (基于环境均值的斜率)平均值:", mean(FW_params$Slope_mean, na.rm = TRUE), "\n")
cat("Slope_mean标准差:", sd(FW_params$Slope_mean, na.rm = TRUE), "\n")
cat("Intcp_para_adj (基于环境参数的矫正截距)平均值:", mean(FW_params$Intcp_para_adj, na.rm = TRUE), "\n")
cat("Intcp_para_adj标准差:", sd(FW_params$Intcp_para_adj, na.rm = TRUE), "\n")
cat("Intcp_para (基于环境参数的原始截距)平均值:", mean(FW_params$Intcp_para, na.rm = TRUE), "\n")
cat("Intcp_para标准差:", sd(FW_params$Intcp_para, na.rm = TRUE), "\n")
cat("Slope_para (基于环境参数的斜率)平均值:", mean(FW_params$Slope_para, na.rm = TRUE), "\n")
cat("Slope_para标准差:", sd(FW_params$Slope_para, na.rm = TRUE), "\n")
cat("R2_mean (基于环境均值回归的R2)平均值:", mean(FW_params$R2_mean, na.rm = TRUE), "\n")
cat("R2_mean标准差:", sd(FW_params$R2_mean, na.rm = TRUE), "\n")
cat("R2_para (基于环境参数回归的R2)平均值:", mean(FW_params$R2_para, na.rm = TRUE), "\n")
cat("R2_para标准差:", sd(FW_params$R2_para, na.rm = TRUE), "\n")
cat("\n")

cat("GWAS分析结果（只对Intcp_para_adj和Slope_para）:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (trait in traits_for_gwas) {
  cat("\n性状:", trait, "\n")
  if (!is.null(gwas_results[[trait]]) && !is.null(gwas_results[[trait]]$gwas_results)) {
    cat("  GWAS分析SNP数:", nrow(gwas_results[[trait]]$gwas_results), "\n")
    cat("  选择的显著标记数:", length(significant_markers_list[[trait]]), "\n")
    if (length(significant_markers_list[[trait]]) > 0) {
      cat("  前5个显著标记:", paste(head(significant_markers_list[[trait]], 5), collapse = ", "), "\n")
    }
  } else {
    cat("  GWAS分析失败\n")
  }
}
sink()

cat("分析报告已保存到:", report_file, "\n")

# ===========================================
# 8. 清理和结束
# ===========================================

cat("\n=== 清理和结束 ===\n")

# 关闭并行
stopCluster(cl)
cat("并行集群已关闭\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("FW参数计算和GWAS分析完全结束!\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", output_dir, "\n")
cat("环境参数:", kPara_Name, "\n")
cat("\n")

cat("主要输出文件:\n")
cat("---------------\n")
cat("1. 完整FW参数 (8个参数): FW_parameters_full.csv\n")
cat("2. GWAS结果目录（只对Intcp_para_adj和Slope_para）: GWAS_Intcp_para_adj/ 和 GWAS_Slope_para/\n")
cat("3. 显著标记列表: Significant_Markers_*.txt\n")
cat("4. GWAS汇总: GWAS_Summary.csv\n")
cat("5. 曼哈顿图和QQ图: Manhattan_Plot_*.png 和 QQ_Plot_*.png\n")
cat("6. 完整R数据: FW_GWAS_Analysis_Results.RData\n")
cat("7. 分析报告: Analysis_Report.txt\n")

cat("\n最终总结:\n")
cat("1. FW参数计算成功，共计算了", nrow(FW_params), "个品种的8个FW参数\n")
cat("2. GWAS分析完成（只对矫正截距和斜率）:\n")
cat("   - Intcp_para_adj (矫正截距): 发现", length(significant_markers_list[["Intcp_para_adj"]]), "个显著标记\n")
cat("   - Slope_para (斜率): 发现", length(significant_markers_list[["Slope_para"]]), "个显著标记\n")
cat("3. 所有结果已保存到", output_dir, "目录\n")

cat("\n分析完成!\n")