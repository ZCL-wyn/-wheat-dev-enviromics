#!/usr/bin/env Rscript

# ==================== 0. 初始设置 ====================
cat("========== Thousand Kernel Weight QTL Epistasis Interaction Plot (QTKW-2D × TaGW2L-7D) ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# 设置基础目录
base_dir <- "/mnt/7t_storage/zhangcl/TKW"

# ==================== 1. 加载必需的R包 ====================
cat("\n1. Loading required R packages...\n")

required_packages <- c(
  "vcfR",          # 用于读取VCF文件
  "ggplot2",       # 用于绘图
  "RColorBrewer",  # 用于颜色调色板
  "dplyr",         # 用于数据处理
  "tidyr",         # 用于数据整理
  "ggpubr",        # 用于统计标注和图形组合
  "car",           # 用于方差分析
  "emmeans",       # 用于事后检验
  "patchwork",     # 用于图形组合
  "scales",        # 用于更好的图形比例
  "cowplot",       # 用于图形组合
  "reshape2"       # 用于数据重塑
)

# 安装缺失的包
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

# 加载所有包
cat("Loading packages:\n")
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
  cat(sprintf("  ✓ %s (version: %s)\n", pkg, packageVersion(pkg)))
}

# ==================== 2. 定义辅助函数 ====================
cat("\n2. Defining helper functions...\n")

# 函数: 将基因型转换为数值 (0/2) - 纯合子编码
convert_gt_to_numeric <- function(gt_mat) {
  cat("  Converting genotypes to numeric (0/1/2)...\n")
  out <- matrix(NA_real_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  
  for (i in 1:nrow(gt_mat)) {
    for (j in 1:ncol(gt_mat)) {
      gt_str <- gt_mat[i, j]
      
      # 处理缺失值
      if (is.na(gt_str) || gt_str %in% c("./.", ".|.", ".", "NA", "")) {
        out[i, j] <- NA_real_
      } else {
        # 分割等位基因
        alleles <- strsplit(gt_str, "[|/]")[[1]]
        
        if (length(alleles) == 2) {
          # 尝试转换为数值
          a1 <- suppressWarnings(as.numeric(alleles[1]))
          a2 <- suppressWarnings(as.numeric(alleles[2]))
          
          if (is.finite(a1) && is.finite(a2)) {
            out[i, j] <- a1 + a2
          } else {
            # 如果是字符型等位基因（如A/T/C/G），转换为数值
            if (alleles[1] == alleles[2]) {
              # 纯合子
              out[i, j] <- 0  # 默认设为0，代表参考等位基因
            } else {
              # 杂合子
              out[i, j] <- 1
            }
          }
        } else {
          out[i, j] <- NA_real_
        }
      }
    }
  }
  
  rownames(out) <- rownames(gt_mat)
  colnames(out) <- colnames(gt_mat)
  
  cat(sprintf("  Conversion completed: %d rows × %d columns\n", nrow(out), ncol(out)))
  return(out)
}

# 函数: 拟合线性模型并提取截距和斜率
extract_line_effects <- function(line_data) {
  if (nrow(line_data) < 2) {
    return(NULL)
  }
  
  # 拟合线性模型: PH ~ 环境因子
  model <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
  
  if (inherits(model, "try-error")) {
    return(NULL)
  }
  
  # 提取系数
  coefs <- coef(model)
  if (length(coefs) < 2) {
    return(NULL)
  }
  
  # 提取统计信息
  fit_summary <- summary(model)
  
  return(list(
    Intercept = coefs[1],
    Slope = coefs[2],
    Intercept_SE = fit_summary$coefficients[1, 2],
    Slope_SE = fit_summary$coefficients[2, 2],
    Intercept_p = fit_summary$coefficients[1, 4],
    Slope_p = fit_summary$coefficients[2, 4],
    R_squared = fit_summary$r.squared,
    n_points = nrow(line_data)
  ))
}

# 函数: 执行两因素方差分析并返回结果
perform_two_way_anova <- function(data, y_var, factor1, factor2) {
  # 构建公式
  formula <- as.formula(paste(y_var, "~", factor1, "*", factor2))
  
  # 执行ANOVA
  anova_result <- aov(formula, data = data)
  anova_summary <- summary(anova_result)
  
  # 提取效应大小
  ss <- anova_summary[[1]]$`Sum Sq`
  df <- anova_summary[[1]]$Df
  f_val <- anova_summary[[1]]$`F value`
  p_val <- anova_summary[[1]]$`Pr(>F)`
  
  # 计算偏η²
  total_ss <- sum(ss, na.rm = TRUE)
  partial_eta2 <- ss / total_ss
  
  # 创建结果数据框
  result_df <- data.frame(
    Source = c(factor1, factor2, paste(factor1, ":", factor2, sep = ""), "Residuals"),
    SS = ss,
    DF = df,
    F_value = f_val,
    P_value = p_val,
    Partial_eta2 = partial_eta2,
    stringsAsFactors = FALSE
  )
  
  # 移除Residuals行的F值和P值
  result_df$F_value[nrow(result_df)] <- NA
  result_df$P_value[nrow(result_df)] <- NA
  
  return(list(
    anova = anova_result,
    summary = anova_summary,
    table = result_df[1:(nrow(result_df)-1), ],  # 移除Residuals
    formula = formula
  ))
}

# ==================== 3. 主分析函数 ====================
cat("\n3. Defining main analysis function for QTKW-2D × TaGW2L-7D interaction...\n")

analyze_qtl_interaction <- function() {
  cat("\n========== Starting QTKW-2D × TaGW2L-7D Interaction Analysis ==========\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # ==================== 1. 路径设置 ====================
  cat("\n[Step 1] Setting up paths...\n")
  output_dir <- file.path(base_dir, "QTKW_TaGW2L_Interaction_Analysis")
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    cat(sprintf("  Creating output directory: %s\n", output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
  } else {
    cat(sprintf("  Output directory already exists: %s\n", output_dir))
  }
  
  # 定义输入文件路径
  vcf_file <- file.path(base_dir, "extracted_G1_to_G339.vcf")
  qtl_file <- file.path(base_dir, "QTL4.csv")
  phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")
  environment_file <- file.path(base_dir, "ECs_results.csv")
  
  cat(sprintf("  VCF file: %s\n", vcf_file))
  cat(sprintf("  QTL file: %s\n", qtl_file))
  cat(sprintf("  Phenotype file: %s\n", phenotype_file))
  cat(sprintf("  Environment factor file: %s\n", environment_file))
  
  # ==================== 2. 读取QTL标记信息（只关注两个目标QTL） ====================
  cat("\n[Step 2] Reading QTL marker information for QTKW-2D and TaGW2L-7D...\n")
  
  if (!file.exists(qtl_file)) {
    stop(sprintf("Error: QTL file does not exist: %s", qtl_file))
  }
  
  # 读取QTL数据
  qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  Read QTL data: %d rows × %d columns\n", nrow(qtl_df), ncol(qtl_df)))
  
  # 提取目标QTL：QTKW-2D和TaGW2L-7D
  target_qtls <- c("QTKW-2D", "TaGW2L-7D")
  
  # 检查目标QTL是否存在
  found_qtls <- target_qtls[target_qtls %in% qtl_df$QTLname]
  if (length(found_qtls) < 2) {
    missing_qtls <- setdiff(target_qtls, found_qtls)
    stop(sprintf("Error: Some target QTLs not found in QTL4.csv: %s", paste(missing_qtls, collapse = ", ")))
  }
  
  # 筛选目标QTL
  target_qtl_df <- qtl_df[qtl_df$QTLname %in% target_qtls, ]
  
  cat("\n  Target QTL Information:\n")
  for (i in 1:nrow(target_qtl_df)) {
    cat(sprintf("    %s: %s (Chr%s: %s, P=%.2e, effect=%.3f, PVE=%.2f%%)\n", 
                target_qtl_df$QTLname[i], target_qtl_df$SNP[i], 
                target_qtl_df$Chromosome[i], target_qtl_df$Position[i],
                target_qtl_df$P.value[i], target_qtl_df$effect[i],
                target_qtl_df$`Phenotype_Variance_Explained(%)`[i]))
  }
  
  # 提取SNP IDs和QTL名称
  target_snps <- target_qtl_df$SNP
  target_qtl_names <- target_qtl_df$QTLname
  
  # 创建有效的R变量名（替换连字符为下划线）
  target_qtl_names_safe <- make.names(target_qtl_names)
  
  # 创建SNP到QTL名称的映射
  snp_to_qtl <- setNames(target_qtl_names, target_snps)
  snp_to_qtl_safe <- setNames(target_qtl_names_safe, target_snps)
  
  n_target_qtls <- length(target_snps)
  cat(sprintf("\n  Found %d target QTL markers\n", n_target_qtls))
  
  # ==================== 3. 读取VCF并提取目标SNPs ====================
  cat("\n[Step 3] Reading VCF and extracting target SNPs...\n")
  
  if (!file.exists(vcf_file)) {
    stop(sprintf("Error: VCF file does not exist: %s", vcf_file))
  }
  
  # 读取VCF文件
  vcf_data <- read.vcfR(vcf_file, verbose = FALSE)
  cat(sprintf("  Read VCF: %d variants × %d samples\n", 
              nrow(vcf_data@fix), ncol(vcf_data@gt)))
  
  # 匹配SNP IDs
  vcf_snp_ids <- vcf_data@fix[, "ID"]
  matched_indices <- which(vcf_snp_ids %in% target_snps)
  
  cat(sprintf("  Matched %d markers by ID\n", length(matched_indices)))
  
  if (length(matched_indices) != n_target_qtls) {
    cat("  Warning: Not all target QTL SNPs found in VCF\n")
    cat("  Found:", paste(vcf_snp_ids[matched_indices], collapse = ", "), "\n")
  }
  
  if (length(matched_indices) == 0) {
    stop("Error: No matching SNP markers found in VCF")
  }
  
  # 提取目标SNPs
  extracted_vcf <- vcf_data[matched_indices, ]
  cat(sprintf("  Extracted %d matching SNP markers\n", length(matched_indices)))
  
  # ==================== 4. 提取基因型并转换为0/1/2 ====================
  cat("\n[Step 4] Extracting and converting genotypes...\n")
  
  # 提取基因型
  gt_matrix <- extract.gt(extracted_vcf, element = "GT")
  cat(sprintf("  Genotype matrix: %d rows × %d columns\n", nrow(gt_matrix), ncol(gt_matrix)))
  
  # 转换为数值型（0/1/2）
  gt_numeric <- convert_gt_to_numeric(gt_matrix)
  
  # 显示基因型统计
  genotype_counts <- table(gt_numeric, useNA = "always")
  cat("  Genotype distribution:\n")
  for (val in names(genotype_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", 
                val, genotype_counts[val], 
                genotype_counts[val]/sum(!is.na(gt_numeric))*100))
  }
  
  # ==================== 5. 过滤杂合子并构建QTL基因型数据框 ====================
  cat("\n[Step 5] Filtering heterozygotes and constructing QTL genotype data frame...\n")
  
  # 转换为数据框并过滤杂合子
  geno_df <- as.data.frame(t(gt_numeric), stringsAsFactors = FALSE)
  geno_df$Sample <- rownames(geno_df)
  
  # 过滤杂合子（基因型为1的样本）
  snp_cols <- setdiff(colnames(geno_df), "Sample")
  geno_df$has_heterozygous <- apply(geno_df[, snp_cols, drop = FALSE], 1, function(x) {
    any(x == 1, na.rm = TRUE)
  })
  
  n_heterozygous <- sum(geno_df$has_heterozygous)
  cat(sprintf("  Found %d heterozygous samples (%.1f%%)\n", 
              n_heterozygous, n_heterozygous/nrow(geno_df)*100))
  
  # 移除杂合子
  geno_df <- geno_df[!geno_df$has_heterozygous, ]
  cat(sprintf("  Number of samples after filtering heterozygotes: %d\n", nrow(geno_df)))
  
  # 将基因型转换为0/1编码（0=参考等位基因，1=替代等位基因）
  # 注意：在gt_numeric中，0=参考等位基因纯合，2=替代等位基因纯合
  # 我们需要将2转换为1
  for (col in snp_cols) {
    geno_df[[col]] <- ifelse(geno_df[[col]] == 2, 1, geno_df[[col]])
  }
  
  # 使用安全的QTL名称重命名列
  colnames(geno_df)[1:length(snp_cols)] <- target_qtl_names_safe
  
  # 将QTL基因型转换为因子
  qtl_cols_safe <- target_qtl_names_safe
  for (col in qtl_cols_safe) {
    geno_df[[col]] <- as.factor(geno_df[[col]])
  }
  
  # 检查每个QTL的等位基因频率
  cat("\n  QTL allele frequencies after filtering heterozygotes:\n")
  for (i in 1:length(qtl_cols_safe)) {
    col <- qtl_cols_safe[i]
    original_name <- target_qtl_names[i]
    freq <- table(geno_df[[col]])
    cat(sprintf("    %s (%s): 0=%d (%.1f%%), 1=%d (%.1f%%)\n", 
                original_name, col,
                ifelse("0" %in% names(freq), freq["0"], 0),
                ifelse("0" %in% names(freq), freq["0"]/nrow(geno_df)*100, 0),
                ifelse("1" %in% names(freq), freq["1"], 0),
                ifelse("1" %in% names(freq), freq["1"]/nrow(geno_df)*100, 0)))
  }
  
  # ==================== 6. 读取表型和环境因子数据 ====================
  cat("\n[Step 6] Reading phenotype and environment factor data...\n")
  
  # 检查文件是否存在
  if (!file.exists(phenotype_file)) {
    stop(sprintf("Error: Phenotype file does not exist: %s", phenotype_file))
  }
  
  if (!file.exists(environment_file)) {
    stop(sprintf("Error: Environment factor file does not exist: %s", environment_file))
  }
  
  # 读取表型数据
  cat("  Reading phenotype data...\n")
  phenotype_data <- read.table(
    phenotype_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, na.strings = c("", "NA")
  )
  
  cat(sprintf("  Phenotype data: %d rows × %d columns\n", nrow(phenotype_data), ncol(phenotype_data)))
  
  # 数据清理和重命名
  phenotype_data$genotype <- as.character(phenotype_data$genotype)
  phenotype_data$env_code <- as.character(phenotype_data$env_code)
  phenotype_data$TKW <- as.numeric(phenotype_data$TKW)
  
  # 重命名列
  colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
  colnames(phenotype_data)[colnames(phenotype_data) == "TKW"] <- "PH"
  
  # 移除缺失值
  n_before <- nrow(phenotype_data)
  phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]
  n_after <- nrow(phenotype_data)
  
  cat(sprintf("  Removing missing values: %d → %d rows (removed %d rows)\n", 
              n_before, n_after, n_before - n_after))
  
  # 读取环境因子数据
  cat("  Reading environment factor data...\n")
  ecs_data <- read.csv(environment_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  
  # 重命名第一列
  colnames(ecs_data)[1] <- "env_code"
  ecs_data$env_code <- as.character(ecs_data$env_code)
  
  cat(sprintf("  Environment factor data: %d rows × %d columns\n", nrow(ecs_data), ncol(ecs_data)))
  
  # 使用PAR_TEMP&GS67作为目标环境因子
  target_factor <- "PAR_TEMP&GS67"
  
  # 确认目标环境因子存在
  if (!target_factor %in% colnames(ecs_data)) {
    cat("  Warning: Target factor not found, trying fuzzy matching...\n")
    possible <- grep("PAR.*TEMP", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
    if (length(possible) > 0) {
      target_factor <- possible[1]
      cat(sprintf("  Using: %s\n", target_factor))
    } else {
      stop("Error: Cannot find matching environment factor")
    }
  }
  
  cat(sprintf("  Target environment factor: %s\n", target_factor))
  
  # 提取环境因子数据
  env_factor_data <- ecs_data[, c("env_code", target_factor)]
  colnames(env_factor_data)[2] <- "env_factor_value"
  env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)
  
  # 移除缺失值
  n_before_env <- nrow(env_factor_data)
  env_factor_data <- env_factor_data[!is.na(env_factor_data$env_factor_value), ]
  n_after_env <- nrow(env_factor_data)
  
  cat(sprintf("  Environment factor data: %d → %d rows (removed %d missing values)\n", 
              n_before_env, n_after_env, n_before_env - n_after_env))
  
  # ==================== 7. 合并数据并计算株系响应 ====================
  cat("\n[Step 7] Merging data and calculating line responses...\n")
  
  # 合并表型和环境因子
  pheno_env_data <- merge(phenotype_data, env_factor_data, by = "env_code", all.x = TRUE)
  pheno_env_data <- pheno_env_data[!is.na(pheno_env_data$env_factor_value), ]
  
  cat(sprintf("  Merged phenotype and environment data: %d rows\n", nrow(pheno_env_data)))
  
  # 计算每个株系的截距和斜率
  cat("  Calculating line-specific slopes and intercepts...\n")
  
  line_effects <- data.frame()
  
  for (line in unique(pheno_env_data$line_code)) {
    line_data <- pheno_env_data[pheno_env_data$line_code == line, ]
    
    if (nrow(line_data) >= 2) {  # 至少需要2个点来拟合直线
      effects <- extract_line_effects(line_data)
      
      if (!is.null(effects)) {
        line_effects <- rbind(line_effects, data.frame(
          line_code = line,
          Intercept = effects$Intercept,
          Slope = effects$Slope,
          Intercept_SE = effects$Intercept_SE,
          Slope_SE = effects$Slope_SE,
          Intercept_p = effects$Intercept_p,
          Slope_p = effects$Slope_p,
          R_squared = effects$R_squared,
          n_points = effects$n_points,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  cat(sprintf("  Calculated effects for %d lines\n", nrow(line_effects)))
  
  if (nrow(line_effects) == 0) {
    stop("Error: No line effects calculated. Check if lines have sufficient data points.")
  }
  
  # 合并QTL基因型到株系效应数据
  final_data <- merge(line_effects, geno_df, by.x = "line_code", by.y = "Sample", all.x = TRUE)
  
  # 移除缺失值
  final_data <- final_data[complete.cases(final_data[, qtl_cols_safe]), ]
  
  cat(sprintf("  Final analysis data: %d lines with complete QTL genotype information\n", nrow(final_data)))
  
  # 添加原始QTL名称列以便于解释
  for (i in 1:length(target_qtl_names)) {
    final_data[[target_qtl_names[i]]] <- final_data[[target_qtl_names_safe[i]]]
  }
  
  # ==================== 8. 检查数据分布 ====================
  cat("\n[Step 8] Checking data distribution...\n")
  
  # 基本统计
  cat("  Basic statistics for Intercept:\n")
  cat(sprintf("    Mean: %.2f ± %.2f (SD)\n", mean(final_data$Intercept), sd(final_data$Intercept)))
  cat(sprintf("    Range: %.2f to %.2f\n", min(final_data$Intercept), max(final_data$Intercept)))
  
  cat("  Basic statistics for Slope:\n")
  cat(sprintf("    Mean: %.4f ± %.4f (SD)\n", mean(final_data$Slope), sd(final_data$Slope)))
  cat(sprintf("    Range: %.4f to %.4f\n", min(final_data$Slope), max(final_data$Slope)))
  
  # 检查每个QTL基因型组合的样本数
  cat("\n  Sample size for each QTL genotype combination:\n")
  genotype_counts <- final_data %>%
    group_by(!!sym(target_qtl_names_safe[1]), !!sym(target_qtl_names_safe[2])) %>%
    summarise(n = n(), .groups = "drop")
  
  print(as.data.frame(genotype_counts))
  
  # ==================== 9. 执行两因素方差分析 ====================
  cat("\n[Step 9] Performing two-way ANOVA for intercept...\n")
  
  # 获取QTL名称
  qtl1_safe <- target_qtl_names_safe[1]
  qtl2_safe <- target_qtl_names_safe[2]
  qtl1_original <- target_qtl_names[1]
  qtl2_original <- target_qtl_names[2]
  
  cat(sprintf("  Testing interaction between %s and %s on intercept\n", 
              qtl1_original, qtl2_original))
  
  # 执行ANOVA
  anova_result <- perform_two_way_anova(final_data, "Intercept", qtl1_safe, qtl2_safe)
  
  # 打印ANOVA结果
  cat("\n  Two-way ANOVA results for Intercept:\n")
  print(anova_result$table)
  
  # 提取互作项的p值
  interaction_row <- anova_result$table[grepl(":", anova_result$table$Source), ]
  interaction_p <- interaction_row$P_value
  
  cat(sprintf("\n  Interaction p-value: %.4f\n", interaction_p))
  
  # ==================== 10. 执行事后检验（如果互作显著） ====================
  if (interaction_p < 0.05) {
    cat("\n  Significant interaction detected! Performing post-hoc tests...\n")
    
    # 使用emmeans进行事后检验
    model <- lm(Intercept ~ get(qtl1_safe) * get(qtl2_safe), data = final_data)
    emm_interaction <- emmeans(model, specs = pairwise ~ get(qtl1_safe) | get(qtl2_safe))
    emm_simple <- emmeans(model, specs = pairwise ~ get(qtl1_safe) : get(qtl2_safe))
    
    # 保存事后检验结果
    posthoc_results <- list(
      interaction_contrasts = summary(emm_interaction$contrasts),
      simple_contrasts = summary(emm_simple$contrasts)
    )
    
    cat("\n  Post-hoc contrasts (QTKW-2D within each TaGW2L-7D genotype):\n")
    print(posthoc_results$interaction_contrasts)
    
  } else {
    cat("\n  No significant interaction detected at p < 0.05 level\n")
    posthoc_results <- NULL
  }
  
  # ==================== 11. 绘制互作图 ====================
  cat("\n[Step 10] Creating interaction plots...\n")
  
  # 11.1 创建箱型图
  cat("  Creating boxplot for interaction...\n")
  
  # 为基因型创建更好的标签
  final_data$TaGW2L_7D_label <- factor(
    final_data[[qtl2_safe]],
    levels = c("0", "1"),
    labels = c(paste0(qtl2_original, ": Ref"), paste0(qtl2_original, ": Alt"))
  )
  
  final_data$QTKW_2D_label <- factor(
    final_data[[qtl1_safe]],
    levels = c("0", "1"),
    labels = c(paste0(qtl1_original, ": Ref"), paste0(qtl1_original, ": Alt"))
  )
  
  # 创建组合基因型标签
  final_data$Genotype_combination <- interaction(
    final_data$TaGW2L_7D_label,
    final_data$QTKW_2D_label,
    sep = " × "
  )
  
  # 计算每个基因型组合的平均值和标准误
  genotype_summary <- final_data %>%
    group_by(Genotype_combination, !!sym(qtl1_safe), !!sym(qtl2_safe)) %>%
    summarise(
      n = n(),
      mean_intercept = mean(Intercept, na.rm = TRUE),
      sd_intercept = sd(Intercept, na.rm = TRUE),
      se_intercept = sd_intercept / sqrt(n),
      ci_lower = mean_intercept - 1.96 * se_intercept,
      ci_upper = mean_intercept + 1.96 * se_intercept,
      .groups = "drop"
    )
  
  # 创建箱型图
  p_boxplot <- ggplot(final_data, aes(x = TaGW2L_7D_label, y = Intercept, fill = QTKW_2D_label)) +
    geom_boxplot(
      alpha = 0.7,
      outlier.shape = 16,
      outlier.size = 2,
      outlier.alpha = 0.5,
      width = 0.7
    ) +
    geom_point(
      data = genotype_summary,
      aes(x = TaGW2L_7D_label, y = mean_intercept, group = QTKW_2D_label),
      position = position_dodge(width = 0.7),
      shape = 23,
      size = 4,
      fill = "white",
      color = "black"
    ) +
    geom_errorbar(
      data = genotype_summary,
      aes(x = TaGW2L_7D_label, y = mean_intercept, 
          ymin = mean_intercept - se_intercept, 
          ymax = mean_intercept + se_intercept, group = QTKW_2D_label),
      position = position_dodge(width = 0.7),
      width = 0.2,
      color = "black",
      size = 0.8
    ) +
    scale_fill_brewer(palette = "Set2", name = qtl1_original) +
    labs(
      title = sprintf("Interaction between %s and %s on TKW Intercept", 
                      qtl1_original, qtl2_original),
      subtitle = sprintf("Two-way ANOVA: Interaction p = %.4f", interaction_p),
      x = qtl2_original,
      y = "Intercept (Baseline TKW, g)",
      caption = sprintf("n = %d lines | Environment factor: %s", 
                       nrow(final_data), target_factor)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
  
  # 11.2 创建小提琴图
  cat("  Creating violin plot...\n")
  
  p_violin <- ggplot(final_data, aes(x = TaGW2L_7D_label, y = Intercept, fill = QTKW_2D_label)) +
    geom_violin(
      alpha = 0.6,
      trim = TRUE,
      scale = "width",
      width = 0.8
    ) +
    geom_boxplot(
      width = 0.15,
      alpha = 0.8,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_point(
      data = genotype_summary,
      aes(x = TaGW2L_7D_label, y = mean_intercept, group = QTKW_2D_label),
      position = position_dodge(width = 0.8),
      shape = 23,
      size = 4,
      fill = "white",
      color = "black"
    ) +
    geom_errorbar(
      data = genotype_summary,
      aes(x = TaGW2L_7D_label, y = mean_intercept, 
          ymin = mean_intercept - se_intercept, 
          ymax = mean_intercept + se_intercept, group = QTKW_2D_label),
      position = position_dodge(width = 0.8),
      width = 0.15,
      color = "black",
      size = 0.8
    ) +
    scale_fill_brewer(palette = "Set1", name = qtl1_original) +
    labs(
      title = sprintf("Distribution of TKW Intercept by %s and %s", 
                      qtl2_original, qtl1_original),
      x = qtl2_original,
      y = "Intercept (Baseline TKW, g)",
      caption = "Violin plots show data distribution, boxplots show quartiles"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # 11.3 创建点图（带误差线）
  cat("  Creating point plot with error bars...\n")
  
  p_point <- ggplot(genotype_summary, 
                    aes(x = TaGW2L_7D_label, y = mean_intercept, 
                        color = QTKW_2D_label, group = QTKW_2D_label)) +
    geom_line(position = position_dodge(width = 0.3), size = 1, alpha = 0.6) +
    geom_point(position = position_dodge(width = 0.3), size = 4) +
    geom_errorbar(
      aes(ymin = mean_intercept - se_intercept, 
          ymax = mean_intercept + se_intercept),
      position = position_dodge(width = 0.3),
      width = 0.2,
      size = 1
    ) +
    geom_text(
      aes(label = sprintf("n=%d", n), y = mean_intercept + se_intercept + 0.5),
      position = position_dodge(width = 0.3),
      size = 3.5,
      vjust = 0
    ) +
    scale_color_brewer(palette = "Dark2", name = qtl1_original) +
    labs(
      title = sprintf("Interaction Plot: %s × %s", qtl1_original, qtl2_original),
      subtitle = sprintf("Interaction p-value: %.4f", interaction_p),
      x = qtl2_original,
      y = "Mean Intercept (Baseline TKW, g) ± SE",
      caption = "Error bars represent standard error"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  # 11.4 创建柱状图
  cat("  Creating bar plot...\n")
  
  p_bar <- ggplot(genotype_summary, 
                  aes(x = TaGW2L_7D_label, y = mean_intercept, 
                      fill = QTKW_2D_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             width = 0.7, alpha = 0.8) +
    geom_errorbar(
      aes(ymin = mean_intercept - se_intercept, 
          ymax = mean_intercept + se_intercept),
      position = position_dodge(width = 0.8),
      width = 0.25,
      color = "black",
      size = 0.7
    ) +
    geom_text(
      aes(label = sprintf("%.1f\n(n=%d)", mean_intercept, n), 
          y = mean_intercept + se_intercept + 0.5),
      position = position_dodge(width = 0.8),
      size = 3.5,
      vjust = 0
    ) +
    scale_fill_brewer(palette = "Set3", name = qtl1_original) +
    labs(
      title = sprintf("Mean TKW Intercept by %s and %s", 
                      qtl2_original, qtl1_original),
      x = qtl2_original,
      y = "Mean Intercept (Baseline TKW, g) ± SE",
      caption = "Error bars represent standard error"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # ==================== 12. 创建综合图形 ====================
  cat("\n[Step 11] Creating comprehensive summary figure...\n")
  
  # 组合所有图形
  combined_plot <- plot_grid(
    p_boxplot + theme(legend.position = "bottom"),
    p_point + theme(legend.position = "none"),
    p_violin + theme(legend.position = "none"),
    p_bar + theme(legend.position = "none"),
    ncol = 2,
    nrow = 2,
    labels = c("A", "B", "C", "D"),
    align = "hv"
  )
  
  # 添加总标题
  title_plot <- ggdraw() +
    draw_label(
      sprintf("QTL Interaction Analysis: %s × %s on Thousand Kernel Weight", 
              qtl1_original, qtl2_original),
      fontface = "bold",
      size = 16,
      x = 0.5,
      hjust = 0.5
    ) +
    draw_label(
      sprintf("Environment factor: %s | Number of lines: %d", 
              target_factor, nrow(final_data)),
      size = 12,
      x = 0.5,
      y = 0.7,
      hjust = 0.5
    )
  
  final_plot <- plot_grid(
    title_plot,
    combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 0.9)
  )
  
  # ==================== 13. 保存所有结果 ====================
  cat("\n[Step 12] Saving all results...\n")
  
  # 保存数据文件
  write.csv(final_data, file.path(output_dir, "Line_Effects_QTL_Genotypes.csv"), row.names = FALSE)
  write.csv(genotype_summary, file.path(output_dir, "Genotype_Summary_Statistics.csv"), row.names = FALSE)
  
  # 保存ANOVA结果
  sink(file.path(output_dir, "ANOVA_Results.txt"))
  cat("========== Two-Way ANOVA Results ==========\n\n")
  cat(sprintf("Dependent variable: Intercept (Baseline TKW)\n"))
  cat(sprintf("Independent variables: %s and %s\n\n", qtl1_original, qtl2_original))
  cat("ANOVA Table:\n")
  print(anova_result$table)
  cat("\n")
  
  if (!is.null(posthoc_results)) {
    cat("Post-hoc Contrasts (QTKW-2D within each TaGW2L-7D genotype):\n")
    print(posthoc_results$interaction_contrasts)
    cat("\n")
    
    cat("Pairwise Comparisons (All genotype combinations):\n")
    print(posthoc_results$simple_contrasts)
  }
  
  cat("\n=============================================\n")
  cat("Interpretation:\n")
  cat(sprintf("- Main effect of %s: p = %.4f\n", 
              qtl1_original, 
              anova_result$table$P_value[1]))
  cat(sprintf("- Main effect of %s: p = %.4f\n", 
              qtl2_original, 
              anova_result$table$P_value[2]))
  cat(sprintf("- Interaction effect: p = %.4f\n", interaction_p))
  
  if (interaction_p < 0.05) {
    cat("- Significant interaction detected: The effect of QTKW-2D on intercept depends on TaGW2L-7D genotype.\n")
  } else {
    cat("- No significant interaction: Effects of QTKW-2D and TaGW2L-7D are additive.\n")
  }
  
  sink()
  
  # 保存图形
  cat("  Saving individual plots...\n")
  ggsave(file.path(output_dir, "Interaction_Boxplot.png"), 
         p_boxplot, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "Interaction_Violin.png"), 
         p_violin, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "Interaction_Point.png"), 
         p_point, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "Interaction_Bar.png"), 
         p_bar, width = 10, height = 7, dpi = 300, bg = "white")
  
  cat("  Saving comprehensive figure...\n")
  ggsave(file.path(output_dir, "Comprehensive_Interaction_Plot.png"), 
         final_plot, width = 14, height = 12, dpi = 300, bg = "white")
  
  # ==================== 14. 创建总结报告 ====================
  cat("\n[Step 13] Creating summary report...\n")
  
  summary_text <- c(
    "=========================================================",
    "  QTL INTERACTION ANALYSIS: QTKW-2D × TaGW2L-7D",
    "=========================================================",
    sprintf("Analysis time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Output directory: %s", output_dir),
    "",
    "1. ANALYSIS PARAMETERS",
    sprintf("  Environment factor: %s", target_factor),
    sprintf("  Target QTLs: %s and %s", qtl1_original, qtl2_original),
    sprintf("  Number of lines analyzed: %d", nrow(final_data)),
    "",
    "2. QTL INFORMATION",
    ""
  )
  
  # 添加QTL信息
  for (i in 1:nrow(target_qtl_df)) {
    summary_text <- c(summary_text,
                     sprintf("  %s:", target_qtl_df$QTLname[i]),
                     sprintf("    SNP: %s", target_qtl_df$SNP[i]),
                     sprintf("    Chromosome: %s, Position: %s", 
                             target_qtl_df$Chromosome[i], target_qtl_df$Position[i]),
                     sprintf("    P-value: %.2e, Effect: %.3f", 
                             target_qtl_df$P.value[i], target_qtl_df$effect[i]),
                     sprintf("    PVE: %.2f%%", 
                             target_qtl_df$`Phenotype_Variance_Explained(%)`[i]),
                     "")
  )
  
  summary_text <- c(summary_text,
                    "3. SAMPLE SIZES BY GENOTYPE COMBINATION",
                    ""
  )
  
  # 添加基因型组合样本数
  for (i in 1:nrow(genotype_summary)) {
    summary_text <- c(summary_text,
                     sprintf("  %s: n = %d lines, mean intercept = %.2f ± %.2f", 
                             genotype_summary$Genotype_combination[i],
                             genotype_summary$n[i],
                             genotype_summary$mean_intercept[i],
                             genotype_summary$se_intercept[i]))
  }
  
  summary_text <- c(summary_text,
                    "",
                    "4. STATISTICAL RESULTS",
                    sprintf("  Two-way ANOVA for Intercept:"),
                    sprintf("    Main effect of %s: F = %.2f, p = %.4f", 
                            qtl1_original, 
                            anova_result$table$F_value[1],
                            anova_result$table$P_value[1]),
                    sprintf("    Main effect of %s: F = %.2f, p = %.4f", 
                            qtl2_original, 
                            anova_result$table$F_value[2],
                            anova_result$table$P_value[2]),
                    sprintf("    Interaction effect: F = %.2f, p = %.4f", 
                            anova_result$table$F_value[3],
                            interaction_p),
                    ""
  )
  
  # 添加互作解释
  if (interaction_p < 0.05) {
    summary_text <- c(summary_text,
                     "  INTERPRETATION:",
                     "    ✓ Significant interaction detected!",
                     "    ✓ The effect of QTKW-2D on TKW intercept depends on TaGW2L-7D genotype.",
                     "    ✓ The genetic effects are non-additive (epistasis present).",
                     "")
    
    # 添加简单效应描述
    summary_text <- c(summary_text,
                     "  SIMPLE EFFECTS ANALYSIS:",
                     "    The effect of QTKW-2D genotype on intercept:",
                     sprintf("    - When TaGW2L-7D = Ref: QTKW-2D effect = %.2f", 
                             genotype_summary$mean_intercept[genotype_summary[[qtl1_safe]] == "1" & 
                                                              genotype_summary[[qtl2_safe]] == "0"] - 
                             genotype_summary$mean_intercept[genotype_summary[[qtl1_safe]] == "0" & 
                                                              genotype_summary[[qtl2_safe]] == "0"]),
                     sprintf("    - When TaGW2L-7D = Alt: QTKW-2D effect = %.2f", 
                             genotype_summary$mean_intercept[genotype_summary[[qtl1_safe]] == "1" & 
                                                              genotype_summary[[qtl2_safe]] == "1"] - 
                             genotype_summary$mean_intercept[genotype_summary[[qtl1_safe]] == "0" & 
                                                              genotype_summary[[qtl2_safe]] == "1"]),
                     "")
  } else {
    summary_text <- c(summary_text,
                     "  INTERPRETATION:",
                     "    ✗ No significant interaction detected.",
                     "    ✓ The effects of QTKW-2D and TaGW2L-7D are additive.",
                     "    ✓ Each QTL contributes independently to TKW intercept.",
                     "")
  }
  
  summary_text <- c(summary_text,
                    "5. KEY FINDINGS",
                    "",
                    sprintf("  Best genotype combination for high TKW intercept: %s", 
                            genotype_summary$Genotype_combination[which.max(genotype_summary$mean_intercept)]),
                    sprintf("    Mean intercept: %.2f g", max(genotype_summary$mean_intercept)),
                    "",
                    sprintf("  Worst genotype combination for TKW intercept: %s", 
                            genotype_summary$Genotype_combination[which.min(genotype_summary$mean_intercept)]),
                    sprintf("    Mean intercept: %.2f g", min(genotype_summary$mean_intercept)),
                    "",
                    sprintf("  Difference between best and worst: %.2f g", 
                            max(genotype_summary$mean_intercept) - min(genotype_summary$mean_intercept)),
                    "",
                    "6. GENERATED FILES",
                    "  1. Line_Effects_QTL_Genotypes.csv - Complete dataset",
                    "  2. Genotype_Summary_Statistics.csv - Summary statistics",
                    "  3. ANOVA_Results.txt - Full statistical output",
                    "  4. Comprehensive_Interaction_Plot.png - All plots combined",
                    "  5. Interaction_Boxplot.png - Boxplot with interaction",
                    "  6. Interaction_Violin.png - Violin plot",
                    "  7. Interaction_Point.png - Point plot with error bars",
                    "  8. Interaction_Bar.png - Bar plot",
                    "",
                    "=========================================================",
                    "ANALYSIS COMPLETE",
                    "=========================================================")
  
  writeLines(summary_text, file.path(output_dir, "Final_Summary_Report.txt"))
  
  # ==================== 15. 完成 ====================
  cat("\n========== Analysis Completed ==========\n")
  cat(sprintf("✅ All results saved to: %s\n", output_dir))
  
  # 列出生成的文件
  generated_files <- list.files(output_dir, pattern = "\\.(csv|txt|png)$", full.names = TRUE)
  cat(sprintf("\nGenerated %d files:\n", length(generated_files)))
  for (file in generated_files) {
    file_size <- file.info(file)$size / 1024  # KB
    cat(sprintf("  - %s (%.1f KB)\n", basename(file), file_size))
  }
  
  cat("\nKey findings:\n")
  if (interaction_p < 0.05) {
    cat(sprintf("  ✓ SIGNIFICANT INTERACTION detected between %s and %s (p = %.4f)\n", 
                qtl1_original, qtl2_original, interaction_p))
    cat("  ✓ Epistasis present: Effect of QTKW-2D depends on TaGW2L-7D genotype\n")
  } else {
    cat(sprintf("  ✗ No significant interaction (p = %.4f)\n", interaction_p))
    cat("  ✓ Additive genetic effects only\n")
  }
  
  cat(sprintf("\nBest genotype combination: %s (%.2f g)\n", 
              genotype_summary$Genotype_combination[which.max(genotype_summary$mean_intercept)],
              max(genotype_summary$mean_intercept)))
  
  cat(sprintf("Worst genotype combination: %s (%.2f g)\n", 
              genotype_summary$Genotype_combination[which.min(genotype_summary$mean_intercept)],
              min(genotype_summary$mean_intercept)))
  
  cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("Total runtime:", format(Sys.time() - start_time, digits = 2), "\n")
  cat("========================================\n")
  
  return(list(
    final_data = final_data,
    genotype_summary = genotype_summary,
    anova_result = anova_result,
    interaction_p = interaction_p,
    plots = list(
      boxplot = p_boxplot,
      violin = p_violin,
      point = p_point,
      bar = p_bar,
      combined = final_plot
    ),
    qtl_names = list(
      qtl1_original = qtl1_original,
      qtl2_original = qtl2_original,
      qtl1_safe = qtl1_safe,
      qtl2_safe = qtl2_safe
    )
  ))
}

# ==================== 4. 主程序 ====================
cat("\n4. Starting main QTL interaction analysis program...\n")

# 记录开始时间
start_time <- Sys.time()

# 执行主函数
tryCatch({
  analysis_results <- analyze_qtl_interaction()
  cat("\n✅ QTL interaction analysis completed successfully!\n")
  
  # 显示关键结果
  cat("\n=== KEY RESULTS SUMMARY ===\n")
  
  cat(sprintf("Interaction between %s and %s on TKW intercept:\n",
              analysis_results$qtl_names$qtl1_original,
              analysis_results$qtl_names$qtl2_original))
  
  cat(sprintf("  Two-way ANOVA interaction p-value: %.4f\n", analysis_results$interaction_p))
  
  if (analysis_results$interaction_p < 0.05) {
    cat("  ➤ SIGNIFICANT INTERACTION DETECTED!\n")
    cat("  ➤ Epistasis present: Genetic effects are non-additive\n")
  } else {
    cat("  ➤ No significant interaction\n")
    cat("  ➤ Genetic effects are additive\n")
  }
  
  cat("\nGenotype combination effects:\n")
  for (i in 1:nrow(analysis_results$genotype_summary)) {
    cat(sprintf("  %s: %.2f ± %.2f g (n=%d)\n",
                analysis_results$genotype_summary$Genotype_combination[i],
                analysis_results$genotype_summary$mean_intercept[i],
                analysis_results$genotype_summary$se_intercept[i],
                analysis_results$genotype_summary$n[i]))
  }
  
  cat("\nAnalysis completed successfully!\n")
  
}, error = function(e) {
  cat(sprintf("\n❌ Program execution failed: %s\n", e$message))
  cat("\nError details:\n")
  print(e)
  cat("\nPlease check if the input file paths and formats are correct.\n")
  quit(status = 1)
})

# ==================== 5. 程序结束 ====================
cat("\n5. Program finished\n")
cat("Total runtime:", format(Sys.time() - start_time, digits = 2), "\n")

cat("\n========== QTL Interaction Analysis Completed ==========\n")