#!/usr/bin/env Rscript

# ==================== 加载必需包 ====================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(vcfR)
  library(multcomp)
  library(car)
  library(RColorBrewer)
  library(openxlsx)  # 用于输出Excel文件
  library(kableExtra)  # 用于美化表格
  library(formattable)
  library(ggrepel)     # 用于标签
})

# ==================== 设置工作目录和路径 ====================
base_dir <- "/mnt/7t_storage/zhangcl/TKW"
setwd(base_dir)

# 输入文件路径
vcf_file <- file.path(base_dir, "extracted_G1_to_G339.vcf")
qtl_file <- file.path(base_dir, "QTL4.csv")
phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")
environment_file <- file.path(base_dir, "ECs_results.csv")
target_factor <- "PAR_TEMP&GS67"

# 输出目录
output_dir <- file.path(base_dir, "Haplotype_Analysis_Tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("========== 开始单倍型详细分析 ==========\n")
cat("工作目录:", base_dir, "\n")
cat("输出目录:", output_dir, "\n")

# ==================== 1. 读取和处理数据 ====================
cat("\n1. 读取和处理数据...\n")

# 读取QTL标记
qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)

# 自动识别SNP列
snp_col <- NULL
if ("SNP" %in% colnames(qtl_df)) {
  snp_col <- "SNP"
} else if ("Marker" %in% colnames(qtl_df)) {
  snp_col <- "Marker"
} else {
  snp_hit <- grep("SNP|snp|Marker|marker|ID|id", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)[1]
  if (!is.null(snp_hit)) {
    snp_col <- snp_hit
    message(sprintf("使用列名: '%s' 作为SNP列", snp_col))
  }
}

if (is.null(snp_col)) {
  stop("无法在QTL.csv中找到SNP/Marker列")
}

# 提取SNP列表并清理
qtl_snps <- unique(qtl_df[[snp_col]])
qtl_snps <- qtl_snps[!is.na(qtl_snps) & qtl_snps != ""]
cat(sprintf("从QTL文件中提取到 %d 个SNP标记\n", length(qtl_snps)))

# 读取VCF文件
cat("读取VCF文件...\n")
vcf_data <- read.vcfR(vcf_file, verbose = FALSE)
vcf_snp_ids <- vcf_data@fix[, "ID"]

# 匹配SNP
matched_indices <- which(vcf_snp_ids %in% qtl_snps)

# 如果通过ID匹配失败，尝试通过染色体+位置匹配
if (length(matched_indices) == 0) {
  cat("通过ID匹配失败，尝试通过染色体+位置匹配...\n")
  
  # 查找染色体和位置列
  chr_col <- grep("CHR|Chr|chr|chromosome|Chromosome", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)[1]
  pos_col <- grep("POS|Pos|pos|position|Position", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)[1]
  
  if (!is.null(chr_col) && !is.null(pos_col)) {
    qtl_pos_ids <- paste(qtl_df[[chr_col]], qtl_df[[pos_col]], sep = "_")
    vcf_pos_ids <- paste(vcf_data@fix[, "CHROM"], vcf_data@fix[, "POS"], sep = "_")
    matched_indices <- which(vcf_pos_ids %in% qtl_pos_ids)
    
    if (length(matched_indices) > 0) {
      cat(sprintf("通过染色体+位置匹配到 %d 个标记\n", length(matched_indices)))
    }
  }
}

if (length(matched_indices) == 0) {
  stop("在VCF文件中未找到匹配的SNP标记")
}

cat(sprintf("成功匹配 %d 个SNP标记\n", length(matched_indices)))

# 提取匹配的VCF数据
extracted_vcf <- vcf_data[matched_indices, ]

# ==================== 2. 提取基因型并构建单倍型 ====================
cat("\n2. 提取基因型并构建单倍型...\n")

# 提取基因型矩阵
gt_matrix <- extract.gt(extracted_vcf, element = "GT")

# 转换基因型为数值型 (0=纯合参考, 1=杂合, 2=纯合变异)
convert_gt_to_numeric <- function(gt_mat) {
  out <- matrix(NA_real_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  for (i in 1:nrow(gt_mat)) {
    for (j in 1:ncol(gt_mat)) {
      gt_str <- gt_mat[i, j]
      if (is.na(gt_str) || gt_str %in% c("./.", ".|.", ".")) {
        out[i, j] <- NA_real_
      } else {
        alleles <- strsplit(gt_str, "[|/]")[[1]]
        if (length(alleles) == 2) {
          a1 <- suppressWarnings(as.numeric(alleles[1]))
          a2 <- suppressWarnings(as.numeric(alleles[2]))
          if (is.finite(a1) && is.finite(a2)) {
            out[i, j] <- a1 + a2
          }
        }
      }
    }
  }
  rownames(out) <- rownames(gt_mat)
  colnames(out) <- colnames(gt_mat)
  out
}

gt_numeric <- convert_gt_to_numeric(gt_matrix)

# 构建单倍型数据框
haplotype_df <- as.data.frame(t(gt_numeric), stringsAsFactors = FALSE)
haplotype_df$Sample <- rownames(haplotype_df)

# 获取SNP列名
snp_cols <- setdiff(colnames(haplotype_df), "Sample")

# 过滤杂合子样本（只保留纯合子0或2）
haplotype_df$has_heterozygous <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
  any(x == 1, na.rm = TRUE)
})

cat(sprintf("原始样本数: %d\n", nrow(haplotype_df)))
haplotype_df <- haplotype_df[!haplotype_df$has_heterozygous, ]
cat(sprintf("过滤杂合子后样本数: %d\n", nrow(haplotype_df)))

# 构建单倍型字符串
haplotype_df$Haplotype <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
  paste(ifelse(is.na(x), "N", x), collapse = "-")
})

# 统计单倍型频率
haplotype_counts <- as.data.frame(table(haplotype_df$Haplotype), stringsAsFactors = FALSE)
colnames(haplotype_counts) <- c("Haplotype_String", "Frequency")
haplotype_counts <- haplotype_counts[order(haplotype_counts$Frequency, decreasing = TRUE), ]

# 为单倍型命名（使用罗马数字）
roman_nums <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")
n_haps <- nrow(haplotype_counts)

if (n_haps <= length(roman_nums)) {
  haplotype_names <- paste0("Hap", roman_nums[1:n_haps])
} else {
  haplotype_names <- paste0("Hap", 1:n_haps)
}

# 创建单倍型映射
haplotype_counts$Haplotype_Name <- haplotype_names
haplotype_counts$Proportion <- haplotype_counts$Frequency / sum(haplotype_counts$Frequency) * 100

# ==================== 3. 单倍型组成表 ====================
cat("\n3. 生成单倍型组成表...\n")

# 创建单倍型组成矩阵
haplotype_composition <- data.frame(
  Haplotype_Name = character(),
  Haplotype_String = character(),
  stringsAsFactors = FALSE
)

# 添加每个SNP的列
for (snp in snp_cols) {
  haplotype_composition[[snp]] <- character()
}

# 为每个单倍型填充组成信息
for (i in 1:nrow(haplotype_counts)) {
  hap_string <- haplotype_counts$Haplotype_String[i]
  hap_name <- haplotype_counts$Haplotype_Name[i]
  
  # 提取样本中的一行作为该单倍型的代表
  sample_row <- haplotype_df[haplotype_df$Haplotype == hap_string, ][1, ]
  
  # 创建新行
  new_row <- list(Haplotype_Name = hap_name, Haplotype_String = hap_string)
  
  # 添加每个SNP的等位基因
  for (snp in snp_cols) {
    gt_value <- sample_row[[snp]]
    if (!is.na(gt_value)) {
      if (gt_value == 0) {
        new_row[[snp]] <- "AA"
      } else if (gt_value == 2) {
        new_row[[snp]] <- "BB"
      } else {
        new_row[[snp]] <- "NA"
      }
    } else {
      new_row[[snp]] <- "NA"
    }
  }
  
  haplotype_composition <- rbind(haplotype_composition, new_row, stringsAsFactors = FALSE)
}

# ==================== 4. 读取表型和环境数据 ====================
cat("\n4. 读取表型和环境数据...\n")

# 读取表型数据
phenotype_data <- read.table(
  phenotype_file, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# 数据清理
phenotype_data$genotype <- as.character(phenotype_data$genotype)
phenotype_data$env_code <- as.character(phenotype_data$env_code)
phenotype_data$TKW <- as.numeric(phenotype_data$TKW)
colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
colnames(phenotype_data)[colnames(phenotype_data) == "TKW"] <- "PH"

# 过滤缺失值
phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]
cat(sprintf("表型数据行数: %d\n", nrow(phenotype_data)))

# 读取环境因子数据
ecs_data <- read.csv(environment_file, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(ecs_data)[1] <- "env_code"
ecs_data$env_code <- as.character(ecs_data$env_code)

# 检查目标因子列
if (!target_factor %in% colnames(ecs_data)) {
  cat("警告: 未找到目标因子列，尝试模糊匹配...\n")
  possible <- grep("PAR.*TEMP", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
  if (length(possible) > 0) {
    target_factor <- possible[1]
    cat(sprintf("使用: %s\n", target_factor))
  } else {
    stop("未找到匹配的环境因子列")
  }
}

# 提取环境因子数据
env_factor_data <- ecs_data[, c("env_code", target_factor)]
colnames(env_factor_data)[2] <- "env_factor_value"
env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)

# ==================== 5. 合并数据并拟合模型 ====================
cat("\n5. 合并数据并拟合模型...\n")

# 合并表型和单倍型数据
pheno_haplo <- merge(
  phenotype_data,
  haplotype_df[, c("Sample", "Haplotype")],
  by.x = "line_code", by.y = "Sample", all.x = FALSE, all.y = FALSE
)

# 合并环境因子
pheno_haplo_env <- merge(
  pheno_haplo,
  env_factor_data,
  by = "env_code", all.x = FALSE, all.y = FALSE
)

# 过滤缺失值
pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$env_factor_value), ]

# 添加单倍型名称映射
haplotype_mapping <- setNames(haplotype_counts$Haplotype_Name, haplotype_counts$Haplotype_String)
pheno_haplo_env$Haplotype_Name <- haplotype_mapping[pheno_haplo_env$Haplotype]

cat(sprintf("合并后数据行数: %d\n", nrow(pheno_haplo_env)))
cat(sprintf("株系数: %d\n", length(unique(pheno_haplo_env$line_code))))
cat(sprintf("单倍型数: %d\n", length(unique(pheno_haplo_env$Haplotype_Name))))

# ==================== 6. 拟合每个单倍型的线性模型 ====================
cat("\n6. 拟合每个单倍型的线性模型...\n")

# 为每个单倍型拟合模型
haplotype_models <- data.frame()
haplotype_stats <- data.frame()

for (hap_name in haplotype_names) {
  # 提取该单倍型的数据
  hap_data <- pheno_haplo_env[pheno_haplo_env$Haplotype_Name == hap_name, ]
  
  if (nrow(hap_data) >= 2) {
    # 拟合线性模型
    model <- lm(PH ~ env_factor_value, data = hap_data)
    model_summary <- summary(model)
    
    # 提取系数
    intercept <- coef(model)["(Intercept)"]
    slope <- coef(model)["env_factor_value"]
    
    # 提取统计信息
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    p_value <- model_summary$coefficients["env_factor_value", "Pr(>|t|)"]
    
    # 计算置信区间
    conf_int <- confint(model, level = 0.95)
    slope_ci_lower <- conf_int["env_factor_value", 1]
    slope_ci_upper <- conf_int["env_factor_value", 2]
    intercept_ci_lower <- conf_int["(Intercept)", 1]
    intercept_ci_upper <- conf_int["(Intercept)", 2]
    
    # 添加到模型结果表
    haplotype_models <- rbind(haplotype_models, data.frame(
      Haplotype_Name = hap_name,
      Intercept = round(intercept, 4),
      Intercept_CI_Lower = round(intercept_ci_lower, 4),
      Intercept_CI_Upper = round(intercept_ci_upper, 4),
      Slope = round(slope, 4),
      Slope_CI_Lower = round(slope_ci_lower, 4),
      Slope_CI_Upper = round(slope_ci_upper, 4),
      R_squared = round(r_squared, 4),
      Adj_R_squared = round(adj_r_squared, 4),
      P_value = format.pval(p_value, digits = 4),
      N_samples = nrow(hap_data),
      N_lines = length(unique(hap_data$line_code)),
      N_environments = length(unique(hap_data$env_code)),
      stringsAsFactors = FALSE
    ))
    
    # 添加到统计表
    haplotype_stats <- rbind(haplotype_stats, data.frame(
      Haplotype_Name = hap_name,
      N_lines = length(unique(hap_data$line_code)),
      N_samples = nrow(hap_data),
      Mean_PH = round(mean(hap_data$PH, na.rm = TRUE), 2),
      SD_PH = round(sd(hap_data$PH, na.rm = TRUE), 2),
      Min_PH = round(min(hap_data$PH, na.rm = TRUE), 2),
      Max_PH = round(max(hap_data$PH, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    ))
  } else {
    cat(sprintf("警告: 单倍型 %s 数据不足，跳过拟合\n", hap_name))
  }
}

# ==================== 7. 拟合每个株系的模型 ====================
cat("\n7. 拟合每个株系的模型...\n")

# 为每个株系拟合模型
line_models_detailed <- data.frame()
line_haplotype_info <- data.frame()

# 统计每个株系的环境数
line_env_counts <- pheno_haplo_env %>%
  group_by(line_code) %>%
  summarise(
    n_environments = n_distinct(env_code),
    n_observations = n(),
    stringsAsFactors = FALSE
  )

# 拟合每个株系的模型
for (i in 1:nrow(line_env_counts)) {
  line <- line_env_counts$line_code[i]
  n_env <- line_env_counts$n_environments[i]
  n_obs <- line_env_counts$n_observations[i]
  
  # 只拟合有足够数据的株系（至少2个环境）
  if (n_env >= 2) {
    line_data <- pheno_haplo_env[pheno_haplo_env$line_code == line, ]
    
    # 获取单倍型信息
    hap_name <- unique(line_data$Haplotype_Name)[1]
    hap_string <- unique(line_data$Haplotype)[1]
    
    # 拟合线性模型
    model <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
    
    if (!inherits(model, "try-error")) {
      model_summary <- summary(model)
      
      # 提取系数和统计信息
      intercept <- coef(model)["(Intercept)"]
      slope <- coef(model)["env_factor_value"]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      slope_p_value <- model_summary$coefficients["env_factor_value", "Pr(>|t|)"]
      
      # 计算置信区间
      conf_int <- try(confint(model, level = 0.95), silent = TRUE)
      if (!inherits(conf_int, "try-error")) {
        slope_ci_lower <- conf_int["env_factor_value", 1]
        slope_ci_upper <- conf_int["env_factor_value", 2]
        intercept_ci_lower <- conf_int["(Intercept)", 1]
        intercept_ci_upper <- conf_int["(Intercept)", 2]
      } else {
        slope_ci_lower <- NA
        slope_ci_upper <- NA
        intercept_ci_lower <- NA
        intercept_ci_upper <- NA
      }
      
      # 添加到株系模型详细表
      line_models_detailed <- rbind(line_models_detailed, data.frame(
        line_code = line,
        Haplotype_Name = hap_name,
        Haplotype_String = hap_string,
        n_environments = n_env,
        n_observations = n_obs,
        Intercept = round(intercept, 4),
        Intercept_CI_Lower = round(intercept_ci_lower, 4),
        Intercept_CI_Upper = round(intercept_ci_upper, 4),
        Slope = round(slope, 4),
        Slope_CI_Lower = round(slope_ci_lower, 4),
        Slope_CI_Upper = round(slope_ci_upper, 4),
        R_squared = round(r_squared, 4),
        Adj_R_squared = round(adj_r_squared, 4),
        Slope_P_value = format.pval(slope_p_value, digits = 4),
        stringsAsFactors = FALSE
      ))
      
      # 添加到株系-单倍型信息表
      line_haplotype_info <- rbind(line_haplotype_info, data.frame(
        line_code = line,
        Haplotype_Name = hap_name,
        Haplotype_String = hap_string,
        n_environments = n_env,
        n_observations = n_obs,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 排序
line_models_detailed <- line_models_detailed[order(line_models_detailed$Haplotype_Name, line_models_detailed$line_code), ]
line_haplotype_info <- line_haplotype_info[order(line_haplotype_info$Haplotype_Name, line_haplotype_info$line_code), ]

cat(sprintf("成功拟合 %d 个株系的模型\n", nrow(line_models_detailed)))

# ==================== 8. 差异分析（斜率和截距） ====================
cat("\n8. 进行斜率和截距的差异分析...\n")

# 准备用于ANOVA的数据
line_models_anova <- line_models_detailed[, c("line_code", "Haplotype_Name", "Slope", "Intercept")]

# 执行ANOVA和Tukey HSD检验（斜率和截距）
anova_results <- list()
tukey_results <- list()

if (nrow(line_models_anova) >= 5 && length(unique(line_models_anova$Haplotype_Name)) >= 2) {
  # 将单倍型转换为因子
  line_models_anova$Haplotype_Name <- as.factor(line_models_anova$Haplotype_Name)
  
  # ==================== 斜率差异分析 ====================
  cat("  - 斜率差异分析...\n")
  slope_anova <- aov(Slope ~ Haplotype_Name, data = line_models_anova)
  slope_anova_summary <- summary(slope_anova)
  
  # 整理斜率ANOVA结果
  anova_results$Slope <- data.frame(
    Source = c("Haplotype", "Residuals"),
    Df = c(slope_anova_summary[[1]]$Df[1], slope_anova_summary[[1]]$Df[2]),
    Sum_Sq = c(slope_anova_summary[[1]]$`Sum Sq`[1], slope_anova_summary[[1]]$`Sum Sq`[2]),
    Mean_Sq = c(slope_anova_summary[[1]]$`Mean Sq`[1], slope_anova_summary[[1]]$`Mean Sq`[2]),
    F_value = c(slope_anova_summary[[1]]$`F value`[1], NA),
    P_value = c(slope_anova_summary[[1]]$`Pr(>F)`[1], NA),
    stringsAsFactors = FALSE
  )
  
  # 斜率Tukey HSD检验
  if (length(unique(line_models_anova$Haplotype_Name)) >= 2) {
    slope_tukey <- glht(slope_anova, linfct = mcp(Haplotype_Name = "Tukey"))
    slope_tukey_summary <- summary(slope_tukey)
    
    slope_tukey_df <- data.frame(
      Comparison = names(slope_tukey_summary$test$coefficients),
      Difference = as.numeric(slope_tukey_summary$test$coefficients),
      Std_Error = as.numeric(slope_tukey_summary$test$sigma),
      t_value = as.numeric(slope_tukey_summary$test$tstat),
      P_value = as.numeric(slope_tukey_summary$test$pvalues),
      Significant = ifelse(slope_tukey_summary$test$pvalues < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    tukey_results$Slope <- slope_tukey_df
  }
  
  # ==================== 截距差异分析 ====================
  cat("  - 截距差异分析...\n")
  intercept_anova <- aov(Intercept ~ Haplotype_Name, data = line_models_anova)
  intercept_anova_summary <- summary(intercept_anova)
  
  # 整理截距ANOVA结果
  anova_results$Intercept <- data.frame(
    Source = c("Haplotype", "Residuals"),
    Df = c(intercept_anova_summary[[1]]$Df[1], intercept_anova_summary[[1]]$Df[2]),
    Sum_Sq = c(intercept_anova_summary[[1]]$`Sum Sq`[1], intercept_anova_summary[[1]]$`Sum Sq`[2]),
    Mean_Sq = c(intercept_anova_summary[[1]]$`Mean Sq`[1], intercept_anova_summary[[1]]$`Mean Sq`[2]),
    F_value = c(intercept_anova_summary[[1]]$`F value`[1], NA),
    P_value = c(intercept_anova_summary[[1]]$`Pr(>F)`[1], NA),
    stringsAsFactors = FALSE
  )
  
  # 截距Tukey HSD检验
  if (length(unique(line_models_anova$Haplotype_Name)) >= 2) {
    intercept_tukey <- glht(intercept_anova, linfct = mcp(Haplotype_Name = "Tukey"))
    intercept_tukey_summary <- summary(intercept_tukey)
    
    intercept_tukey_df <- data.frame(
      Comparison = names(intercept_tukey_summary$test$coefficients),
      Difference = as.numeric(intercept_tukey_summary$test$coefficients),
      Std_Error = as.numeric(intercept_tukey_summary$test$sigma),
      t_value = as.numeric(intercept_tukey_summary$test$tstat),
      P_value = as.numeric(intercept_tukey_summary$test$pvalues),
      Significant = ifelse(intercept_tukey_summary$test$pvalues < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    tukey_results$Intercept <- intercept_tukey_df
  }
  
  cat("斜率和截距差异分析完成\n")
} else {
  cat("警告: 数据不足，无法进行ANOVA/Tukey HSD检验\n")
  anova_results$Slope <- data.frame()
  anova_results$Intercept <- data.frame()
  tukey_results$Slope <- data.frame()
  tukey_results$Intercept <- data.frame()
}

# ==================== 9. 合并所有统计结果 ====================
cat("\n9. 合并所有统计结果...\n")

# 合并单倍型频率、统计和模型结果
haplotype_summary <- merge(
  haplotype_counts[, c("Haplotype_Name", "Frequency", "Proportion")],
  haplotype_stats,
  by = "Haplotype_Name",
  all = TRUE
)

haplotype_summary <- merge(
  haplotype_summary,
  haplotype_models[, c("Haplotype_Name", "Intercept", "Slope", "R_squared", "P_value")],
  by = "Haplotype_Name",
  all = TRUE
)

# 重新排序列
haplotype_summary <- haplotype_summary %>%
  arrange(desc(Frequency)) %>%
  select(
    Haplotype_Name, Frequency, Proportion, N_lines, N_samples,
    Mean_PH, SD_PH, Min_PH, Max_PH,
    Intercept, Slope, R_squared, P_value
  )

# 格式化百分比
haplotype_summary$Proportion <- sprintf("%.2f%%", haplotype_summary$Proportion)

# ==================== 10. 生成株系模型统计表 ====================
cat("\n10. 生成株系模型统计表...\n")

# 计算每个株系的性能排名
line_models_detailed$Slope_Rank <- rank(-line_models_detailed$Slope, ties.method = "min")  # 斜率越大排名越高
line_models_detailed$Intercept_Rank <- rank(-line_models_detailed$Intercept, ties.method = "min")  # 截距越大排名越高
line_models_detailed$R_squared_Rank <- rank(-line_models_detailed$R_squared, ties.method = "min")  # R²越大排名越高

# 添加综合排名（平均值排名）
line_models_detailed$Composite_Rank <- round(
  (line_models_detailed$Slope_Rank + line_models_detailed$Intercept_Rank + line_models_detailed$R_squared_Rank) / 3,
  1
)

# 重新排序
line_models_detailed <- line_models_detailed[order(line_models_detailed$Composite_Rank), ]

# ==================== 11. 生成单倍型热图（可视化） ====================
cat("\n11. 生成单倍型热图...\n")

# 准备热图数据
heatmap_data <- haplotype_composition[, -c(1:2)]  # 移除前两列（名称和字符串）
heatmap_data <- as.matrix(heatmap_data)

# 将AA/BB转换为数值
heatmap_numeric <- matrix(0, nrow = nrow(heatmap_data), ncol = ncol(heatmap_data))
rownames(heatmap_numeric) <- haplotype_composition$Haplotype_Name
colnames(heatmap_numeric) <- colnames(heatmap_data)

for (i in 1:nrow(heatmap_data)) {
  for (j in 1:ncol(heatmap_data)) {
    if (heatmap_data[i, j] == "AA") {
      heatmap_numeric[i, j] <- 0
    } else if (heatmap_data[i, j] == "BB") {
      heatmap_numeric[i, j] <- 1
    } else {
      heatmap_numeric[i, j] <- NA
    }
  }
}

# 绘制热图
pdf(file.path(output_dir, "Haplotype_Heatmap.pdf"), width = 12, height = 8)

# 创建热图数据框
heatmap_df <- as.data.frame(heatmap_numeric)
heatmap_df$Haplotype <- rownames(heatmap_df)

# 转换为长格式
heatmap_long <- heatmap_df %>%
  pivot_longer(cols = -Haplotype, names_to = "SNP", values_to = "Genotype") %>%
  mutate(
    Haplotype = factor(Haplotype, levels = rev(haplotype_names)),
    SNP = factor(SNP, levels = colnames(heatmap_numeric)),
    Genotype = factor(Genotype, levels = c(0, 1), labels = c("AA", "BB"))
  )

# 绘制热图
p_heatmap <- ggplot(heatmap_long, aes(x = SNP, y = Haplotype, fill = Genotype)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("AA" = "#4E79A7", "BB" = "#F28E2B"), na.value = "gray90") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(
    title = "单倍型组成热图",
    fill = "等位基因"
  ) +
  coord_equal()

print(p_heatmap)
dev.off()

# ==================== 12. 生成斜率-截距散点图 ====================
cat("\n12. 生成斜率-截距散点图...\n")

pdf(file.path(output_dir, "Slope_Intercept_Scatter.pdf"), width = 10, height = 8)

# 使用株系数据创建散点图
p_scatter <- ggplot(line_models_detailed, aes(x = Intercept, y = Slope)) +
  geom_point(aes(color = Haplotype_Name, size = R_squared), alpha = 0.7) +
  scale_size_continuous(range = c(2, 6), name = "R²") +
  geom_text_repel(
    data = line_models_detailed %>% filter(Composite_Rank <= 10),  # 只标注排名前10的株系
    aes(label = line_code),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray50",
    max.overlaps = 20
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  ) +
  labs(
    title = "株系斜率-截距关系图",
    x = "截距 (Intercept)",
    y = "斜率 (Slope)",
    color = "单倍型"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = mean(line_models_detailed$Intercept, na.rm = TRUE), 
             linetype = "dashed", color = "gray50")

print(p_scatter)
dev.off()

# ==================== 13. 生成株系性能排名图 ====================
cat("\n13. 生成株系性能排名图...\n")

pdf(file.path(output_dir, "Line_Performance_Ranking.pdf"), width = 14, height = 10)

# 选择排名前30的株系进行可视化
top_lines <- head(line_models_detailed[order(line_models_detailed$Composite_Rank), ], 30)

# 转换为长格式以便绘图
ranking_long <- top_lines %>%
  select(line_code, Haplotype_Name, Slope, Intercept, R_squared) %>%
  pivot_longer(cols = c(Slope, Intercept, R_squared), 
               names_to = "Parameter", 
               values_to = "Value")

# 创建排名图
p_ranking <- ggplot(ranking_long, aes(x = line_code, y = Value, fill = Parameter)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Parameter, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(
    title = "株系性能参数排名（前30名）",
    x = "株系",
    y = "参数值"
  )

print(p_ranking)
dev.off()

# ==================== 14. 保存所有结果为Excel文件 ====================
cat("\n14. 保存所有结果为Excel文件...\n")

# 创建Excel工作簿
wb <- createWorkbook()

# 1. 单倍型统计总表
addWorksheet(wb, "Haplotype_Summary")
writeData(wb, "Haplotype_Summary", haplotype_summary, rowNames = FALSE)

# 2. 单倍型组成表
addWorksheet(wb, "Haplotype_Composition")
writeData(wb, "Haplotype_Composition", haplotype_composition, rowNames = FALSE)

# 3. 单倍型模型参数
addWorksheet(wb, "Haplotype_Models")
writeData(wb, "Haplotype_Models", haplotype_models, rowNames = FALSE)

# 4. 株系模型详细结果（包含斜率和截距）
addWorksheet(wb, "Line_Models_Detailed")
writeData(wb, "Line_Models_Detailed", line_models_detailed, rowNames = FALSE)

# 5. 株系-单倍型对应表
addWorksheet(wb, "Line_Haplotype_Mapping")
writeData(wb, "Line_Haplotype_Mapping", line_haplotype_info, rowNames = FALSE)

# 6. 斜率ANOVA结果
if (nrow(anova_results$Slope) > 0) {
  addWorksheet(wb, "ANOVA_Slope")
  writeData(wb, "ANOVA_Slope", anova_results$Slope, rowNames = FALSE)
}

# 7. 截距ANOVA结果
if (nrow(anova_results$Intercept) > 0) {
  addWorksheet(wb, "ANOVA_Intercept")
  writeData(wb, "ANOVA_Intercept", anova_results$Intercept, rowNames = FALSE)
}

# 8. 斜率Tukey HSD结果
if (nrow(tukey_results$Slope) > 0) {
  addWorksheet(wb, "Tukey_Slope")
  writeData(wb, "Tukey_Slope", tukey_results$Slope, rowNames = FALSE)
}

# 9. 截距Tukey HSD结果
if (nrow(tukey_results$Intercept) > 0) {
  addWorksheet(wb, "Tukey_Intercept")
  writeData(wb, "Tukey_Intercept", tukey_results$Intercept, rowNames = FALSE)
}

# 10. 原始数据
addWorksheet(wb, "Raw_Data")
writeData(wb, "Raw_Data", pheno_haplo_env, rowNames = FALSE)

# 11. SNP信息
snp_info <- data.frame(
  SNP = snp_cols,
  Chromosome = vcf_data@fix[match(snp_cols, rownames(vcf_data@fix)), "CHROM"],
  Position = vcf_data@fix[match(snp_cols, rownames(vcf_data@fix)), "POS"],
  REF = vcf_data@fix[match(snp_cols, rownames(vcf_data@fix)), "REF"],
  ALT = vcf_data@fix[match(snp_cols, rownames(vcf_data@fix)), "ALT"],
  stringsAsFactors = FALSE
)
addWorksheet(wb, "SNP_Information")
writeData(wb, "SNP_Information", snp_info, rowNames = FALSE)

# 保存Excel文件
excel_file <- file.path(output_dir, "Haplotype_Analysis_Results.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)

# ==================== 15. 保存为CSV文件 ====================
cat("\n15. 保存为CSV文件...\n")

write.csv(haplotype_summary, 
          file.path(output_dir, "Haplotype_Summary_Table.csv"), 
          row.names = FALSE, quote = FALSE)

write.csv(haplotype_composition, 
          file.path(output_dir, "Haplotype_Composition_Table.csv"), 
          row.names = FALSE, quote = FALSE)

write.csv(haplotype_models, 
          file.path(output_dir, "Haplotype_Model_Parameters.csv"), 
          row.names = FALSE, quote = FALSE)

write.csv(line_models_detailed, 
          file.path(output_dir, "Line_Models_Detailed.csv"), 
          row.names = FALSE, quote = FALSE)

write.csv(line_haplotype_info, 
          file.path(output_dir, "Line_Haplotype_Mapping.csv"), 
          row.names = FALSE, quote = FALSE)

if (nrow(anova_results$Slope) > 0) {
  write.csv(anova_results$Slope, 
            file.path(output_dir, "ANOVA_Slope_Results.csv"), 
            row.names = FALSE, quote = FALSE)
}

if (nrow(anova_results$Intercept) > 0) {
  write.csv(anova_results$Intercept, 
            file.path(output_dir, "ANOVA_Intercept_Results.csv"), 
            row.names = FALSE, quote = FALSE)
}

if (nrow(tukey_results$Slope) > 0) {
  write.csv(tukey_results$Slope, 
            file.path(output_dir, "Tukey_HSD_Slope_Results.csv"), 
            row.names = FALSE, quote = FALSE)
}

if (nrow(tukey_results$Intercept) > 0) {
  write.csv(tukey_results$Intercept, 
            file.path(output_dir, "Tukey_HSD_Intercept_Results.csv"), 
            row.names = FALSE, quote = FALSE)
}

# ==================== 16. 生成HTML格式的报告 ====================
cat("\n16. 生成HTML格式的报告...\n")

html_file <- file.path(output_dir, "Haplotype_Analysis_Report.html")

# 创建HTML内容
html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>单倍型分析报告</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #34495e; margin-top: 30px; }
    h3 { color: #2c3e50; margin-top: 20px; }
    table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
    th { background-color: #3498db; color: white; padding: 10px; text-align: left; }
    td { padding: 8px; border: 1px solid #ddd; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
    .highlight { background-color: #fff3cd; padding: 5px; border-radius: 3px; }
    .highlight-red { background-color: #f8d7da; color: #721c24; padding: 5px; border-radius: 3px; }
    .highlight-green { background-color: #d4edda; color: #155724; padding: 5px; border-radius = 3px; }
    .figures { display: flex; flex-wrap: wrap; gap: 20px; margin-top: 20px; }
    .figure { flex: 1 1 300px; }
    .figure img { max-width: 100%; height: auto; border: 1px solid #ddd; }
    .section { margin-bottom: 30px; }
    .note { font-style: italic; color: #6c757d; margin-top: 10px; }
    .significance-yes { color: #28a745; font-weight: bold; }
    .significance-no { color: #dc3545; font-weight: bold; }
    .top-line { background-color: #e8f4fd; }
  </style>
</head>
<body>
  <h1>单倍型分析报告</h1>
  
  <div class="summary">
    <h2>分析摘要</h2>
    <p><strong>分析日期:</strong> ', format(Sys.Date(), "%Y年%m月%d日"), '</p>
    <p><strong>总单倍型数:</strong> ', n_haps, '</p>
    <p><strong>总样本数:</strong> ', sum(haplotype_counts$Frequency), '</p>
    <p><strong>总株系数:</strong> ', length(unique(pheno_haplo_env$line_code)), '</p>
    <p><strong>成功拟合模型的株系数:</strong> ', nrow(line_models_detailed), '</p>
    <p><strong>总环境数:</strong> ', length(unique(pheno_haplo_env$env_code)), '</p>
    <p><strong>目标环境因子:</strong> ', target_factor, '</p>
  </div>
  
  <div class="section">
    <h2>单倍型统计总表</h2>
    <p>下表显示了每个单倍型的样本数、比例、表型统计和模型参数：</p>
')

# 添加单倍型统计表
html_content <- paste0(html_content, 
                       kable(head(haplotype_summary, 10), format = "html", row.names = FALSE) %>%
                         kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                       full_width = TRUE) %>%
                         as.character())

html_content <- paste0(html_content, '
    <p class="note">注：完整表格包含 ', nrow(haplotype_summary), ' 行数据，可在Excel文件中查看。</p>
  </div>
  
  <div class="section">
    <h2>株系-单倍型对应表（前20行）</h2>
    <p>下表显示了每个株系所属的单倍型：</p>')

# 添加株系-单倍型对应表
if (nrow(line_haplotype_info) > 0) {
  html_content <- paste0(html_content, 
                         kable(head(line_haplotype_info, 20), format = "html", row.names = FALSE) %>%
                           kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                         full_width = TRUE) %>%
                           as.character())
  
  html_content <- paste0(html_content, '
    <p class="note">注：完整表格包含 ', nrow(line_haplotype_info), ' 个株系，可在Excel文件中查看。</p>
  </div>
  
  <div class="section">
    <h2>株系模型参数排名（前10名）</h2>
    <p>下表显示了综合排名前10的株系的模型参数：</p>')
  
  # 添加株系模型参数排名表
  if (nrow(line_models_detailed) > 0) {
    top_10_lines <- head(line_models_detailed[, c("line_code", "Haplotype_Name", "Intercept", "Slope", "R_squared", "Composite_Rank")], 10)
    html_content <- paste0(html_content, 
                           kable(top_10_lines, format = "html", row.names = FALSE) %>%
                             kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                           full_width = TRUE) %>%
                             as.character())
    
    html_content <- paste0(html_content, '
    <p class="note">注：完整表格包含 ', nrow(line_models_detailed), ' 个株系的详细模型参数，可在Excel文件中查看。</p>
  </div>')
  }
}

# ==================== 添加差异分析结果 ====================
if (nrow(anova_results$Slope) > 0 || nrow(anova_results$Intercept) > 0) {
  html_content <- paste0(html_content, '
  <div class="section">
    <h2>差异分析结果</h2>')
  
  # ==================== 斜率ANOVA分析 ====================
  if (nrow(anova_results$Slope) > 0) {
    slope_p_value <- anova_results$Slope$P_value[1]
    slope_significant <- ifelse(slope_p_value < 0.05, "显著", "不显著")
    slope_class <- ifelse(slope_p_value < 0.05, "highlight-green", "highlight-red")
    
    html_content <- paste0(html_content, '
    <h3>斜率ANOVA分析</h3>
    <p>斜率在不同单倍型间差异<span class="', slope_class, '">', slope_significant, 
    ifelse(slope_p_value < 0.05, ' (P < 0.05)', ' (P ≥ 0.05)'), '</span></p>')
    
    html_content <- paste0(html_content, 
                           kable(anova_results$Slope, format = "html", row.names = FALSE) %>%
                             kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                           full_width = TRUE) %>%
                             as.character())
  }
  
  # ==================== 截距ANOVA分析 ====================
  if (nrow(anova_results$Intercept) > 0) {
    intercept_p_value <- anova_results$Intercept$P_value[1]
    intercept_significant <- ifelse(intercept_p_value < 0.05, "显著", "不显著")
    intercept_class <- ifelse(intercept_p_value < 0.05, "highlight-green", "highlight-red")
    
    html_content <- paste0(html_content, '
    <h3>截距ANOVA分析</h3>
    <p>截距在不同单倍型间差异<span class="', intercept_class, '">', intercept_significant, 
    ifelse(intercept_p_value < 0.05, ' (P < 0.05)', ' (P ≥ 0.05)'), '</span></p>')
    
    html_content <- paste0(html_content, 
                           kable(anova_results$Intercept, format = "html", row.names = FALSE) %>%
                             kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                           full_width = TRUE) %>%
                             as.character())
  }
  
  html_content <- paste0(html_content, '</div>')
}

# ==================== 添加图形展示部分 ====================
html_content <- paste0(html_content, '
  <div class="section">
    <h2>图形展示</h2>
    <div class="figures">
      <div class="figure">
        <h3>单倍型组成热图</h3>
        <p>AA: 参考等位基因纯合; BB: 变异等位基因纯合</p>
        <img src="Haplotype_Heatmap.pdf" alt="单倍型热图">
      </div>
      <div class="figure">
        <h3>株系斜率-截距关系图</h3>
        <p>点的大小表示R²，颜色表示单倍型</p>
        <img src="Slope_Intercept_Scatter.pdf" alt="斜率-截距散点图">
      </div>
      <div class="figure">
        <h3>株系性能排名图</h3>
        <p>显示排名前30株系的斜率、截距和R²</p>
        <img src="Line_Performance_Ranking.pdf" alt="株系性能排名图">
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>数据文件</h2>
    <p>详细数据可在以下文件中找到：</p>
    <ul>
      <li><a href="Haplotype_Analysis_Results.xlsx">Haplotype_Analysis_Results.xlsx</a> - 完整的Excel分析结果</li>
      <li><a href="Haplotype_Summary_Table.csv">Haplotype_Summary_Table.csv</a> - 单倍型统计表</li>
      <li><a href="Haplotype_Composition_Table.csv">Haplotype_Composition_Table.csv</a> - 单倍型组成表</li>
      <li><a href="Haplotype_Model_Parameters.csv">Haplotype_Model_Parameters.csv</a> - 单倍型模型参数表</li>
      <li><a href="Line_Models_Detailed.csv">Line_Models_Detailed.csv</a> - 株系模型详细结果（包含斜率和截距）</li>
      <li><a href="Line_Haplotype_Mapping.csv">Line_Haplotype_Mapping.csv</a> - 株系-单倍型对应表</li>
')

# 添加差异分析文件链接
if (nrow(anova_results$Slope) > 0) {
  html_content <- paste0(html_content, '
      <li><a href="ANOVA_Slope_Results.csv">ANOVA_Slope_Results.csv</a> - 斜率ANOVA结果</li>')
}
if (nrow(anova_results$Intercept) > 0) {
  html_content <- paste0(html_content, '
      <li><a href="ANOVA_Intercept_Results.csv">ANOVA_Intercept_Results.csv</a> - 截距ANOVA结果</li>')
}

html_content <- paste0(html_content, '
    </ul>
  </div>
  
  <hr>
  <p><em>分析完成于: ', format(Sys.time(), "%Y年%m月%d日 %H:%M:%S"), '</em></p>
</body>
</html>
')

# 写入HTML文件
writeLines(html_content, html_file)

# ==================== 17. 输出关键结果到控制台 ====================
cat("\n17. 输出关键结果到控制台...\n")

cat("\n========== 单倍型统计摘要 ==========\n")
cat(sprintf("检测到单倍型总数: %d\n", n_haps))
cat(sprintf("总样本数: %d\n", sum(haplotype_counts$Frequency)))
cat(sprintf("样本最多的单倍型: %s (%d 个样本, %.1f%%)\n", 
            haplotype_counts$Haplotype_Name[1],
            haplotype_counts$Frequency[1],
            haplotype_counts$Proportion[1]))

cat("\n========== 株系模型拟合摘要 ==========\n")
cat(sprintf("成功拟合模型的株系数: %d\n", nrow(line_models_detailed)))
cat(sprintf("平均每个株系的环境数: %.1f\n", mean(line_models_detailed$n_environments, na.rm = TRUE)))
cat(sprintf("平均R²: %.3f\n", mean(line_models_detailed$R_squared, na.rm = TRUE)))

cat("\n========== 株系性能排名前5 ==========\n")
top_5 <- head(line_models_detailed[, c("line_code", "Haplotype_Name", "Intercept", "Slope", "R_squared", "Composite_Rank")], 5)
print(top_5)

cat("\n========== 单倍型-株系分布 ==========\n")
hap_line_distribution <- line_haplotype_info %>%
  group_by(Haplotype_Name) %>%
  summarise(
    n_lines = n(),
    proportion = sprintf("%.1f%%", n() / nrow(line_haplotype_info) * 100),
    stringsAsFactors = FALSE
  ) %>%
  arrange(desc(n_lines))

print(hap_line_distribution)

if (nrow(anova_results$Slope) > 0) {
  cat("\n========== 差异分析摘要 ==========\n")
  
  # 斜率结果
  slope_p_value <- anova_results$Slope$P_value[1]
  cat(sprintf("斜率ANOVA: F = %.3f, P = %s\n", 
              anova_results$Slope$F_value[1],
              format.pval(slope_p_value, digits = 3)))
  
  if (slope_p_value < 0.05) {
    cat("斜率在不同单倍型间存在显著差异 (P < 0.05)\n")
  } else {
    cat("斜率在不同单倍型间无显著差异\n")
  }
  
  # 截距结果
  if (nrow(anova_results$Intercept) > 0) {
    intercept_p_value <- anova_results$Intercept$P_value[1]
    cat(sprintf("\n截距ANOVA: F = %.3f, P = %s\n", 
                anova_results$Intercept$F_value[1],
                format.pval(intercept_p_value, digits = 3)))
    
    if (intercept_p_value < 0.05) {
      cat("截距在不同单倍型间存在显著差异 (P < 0.05)\n")
    } else {
      cat("截距在不同单倍型间无显著差异\n")
    }
  }
}

# ==================== 18. 完成信息 ====================
cat("\n========== 分析完成 ==========\n")
cat(sprintf("输出文件保存在: %s\n", output_dir))
cat("\n生成的主要文件:\n")
cat(sprintf("1. Excel汇总文件: %s\n", excel_file))
cat(sprintf("2. HTML报告: %s\n", html_file))
cat(sprintf("3. 株系模型详细结果: %s\n", file.path(output_dir, "Line_Models_Detailed.csv")))
cat(sprintf("4. 株系-单倍型对应表: %s\n", file.path(output_dir, "Line_Haplotype_Mapping.csv")))
cat(sprintf("5. 单倍型组成热图: %s\n", file.path(output_dir, "Haplotype_Heatmap.pdf")))
cat(sprintf("6. 株系斜率-截距图: %s\n", file.path(output_dir, "Slope_Intercept_Scatter.pdf")))
cat(sprintf("7. 株系性能排名图: %s\n", file.path(output_dir, "Line_Performance_Ranking.pdf")))
cat(sprintf("8. CSV数据表: %s 等\n", file.path(output_dir, "Haplotype_Summary_Table.csv")))

cat("\n📊 重要发现:\n")
if (nrow(hap_line_distribution) > 0) {
  cat(sprintf("1. 最常见的单倍型: %s (%s 个株系)\n", 
              hap_line_distribution$Haplotype_Name[1],
              hap_line_distribution$n_lines[1]))
}
if (nrow(top_5) > 0) {
  cat(sprintf("2. 综合排名第一的株系: %s (单倍型: %s, R²=%.3f)\n",
              top_5$line_code[1], top_5$Haplotype_Name[1], top_5$R_squared[1]))
}
if (nrow(anova_results$Slope) > 0 && anova_results$Slope$P_value[1] < 0.05) {
  cat("3. 斜率在不同单倍型间存在显著差异，表明环境响应性受遗传控制\n")
}
if (nrow(anova_results$Intercept) > 0 && anova_results$Intercept$P_value[1] < 0.05) {
  cat("4. 截距在不同单倍型间存在显著差异，表明基础产量潜力受遗传控制\n")
}

cat("\n✅ 单倍型详细分析完成！\n")