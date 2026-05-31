# ==================== 专门计算图d拟合统计量 ====================
cat("\n=== 开始计算图d (QTL Effect by PAR_TEMP&GS67) 的拟合统计量 ===\n")

# 设置路径
base_dir <- "/mnt/7t_storage/zhangcl/TKW"
slope_qtl_path <- file.path(base_dir, "QTL4.csv")
base_env_path  <- file.path(base_dir, "GAPIT_BLINK_Results/Single_Environment")
ecs_path       <- file.path(base_dir, "ECs_results.csv")

# 1. 读取数据并检查结构
cat("1. 读取和检查数据...\n")

# 1.1 读取QTL数据
if (file.exists(slope_qtl_path)) {
  slope_qtl <- read.csv(slope_qtl_path, stringsAsFactors = FALSE)
  cat(sprintf("  读取QTL数据: %d 行, %d 列\n", nrow(slope_qtl), ncol(slope_qtl)))
  cat("  列名:", paste(colnames(slope_qtl), collapse=", "), "\n")
} else {
  stop(paste("QTL文件不存在:", slope_qtl_path))
}

# 1.2 读取环境协变量数据
if (file.exists(ecs_path)) {
  ecs_data <- read.csv(ecs_path, stringsAsFactors = FALSE)
  cat(sprintf("  读取环境协变量数据: %d 行, %d 列\n", nrow(ecs_data), ncol(ecs_data)))
  cat("  列名:", paste(colnames(ecs_data), collapse=", "), "\n")
} else {
  stop(paste("环境协变量文件不存在:", ecs_path))
}

# 2. 识别环境代码列
cat("\n2. 识别环境代码列...\n")
env_code_candidates <- c("env", "env_code", "environment", "Environment", "Env")
env_code_col <- NULL
for (candidate in env_code_candidates) {
  if (candidate %in% colnames(ecs_data)) {
    env_code_col <- candidate
    break
  }
}
if (is.null(env_code_col)) {
  # 尝试其他可能列名
  env_pattern <- grep("env", colnames(ecs_data), ignore.case = TRUE, value = TRUE)
  if (length(env_pattern) > 0) {
    env_code_col <- env_pattern[1]
  } else {
    # 如果没有找到环境代码列，尝试使用第一列
    env_code_col <- colnames(ecs_data)[1]
    cat(sprintf("  警告: 未找到环境代码列，使用第一列 '%s' 作为环境代码\n", env_code_col))
  }
}
cat(sprintf("  使用环境代码列: %s\n", env_code_col))
colnames(ecs_data)[colnames(ecs_data) == env_code_col] <- "env_code"

# 3. 识别PAR_TEMP&GS67列
cat("\n3. 识别PAR_TEMP&GS67列...\n")
par_temp_candidates <- c("PAR_TEMP&GS67", "PAR_TEMP_GS67", "PAR_TEMP", "TEMP", "PAR_TEMP.GS67")
par_temp_col <- NULL
for (candidate in par_temp_candidates) {
  if (candidate %in% colnames(ecs_data)) {
    par_temp_col <- candidate
    break
  }
}
if (is.null(par_temp_col)) {
  # 尝试模糊匹配
  par_pattern <- grep("PAR|TEMP|GS", colnames(ecs_data), ignore.case = TRUE, value = TRUE)
  if (length(par_pattern) > 0) {
    par_temp_col <- par_pattern[1]
    cat(sprintf("  警告: 未找到精确的PAR_TEMP&GS67列，使用 '%s' 代替\n", par_temp_col))
  } else {
    stop("错误: 找不到环境协变量列")
  }
}
cat(sprintf("  使用环境协变量列: %s\n", par_temp_col))
colnames(ecs_data)[colnames(ecs_data) == par_temp_col] <- "Xvar"

# 4. 提取QTL效应值
cat("\n4. 提取各环境QTL效应值...\n")

# 检查QTL数据中的必要列
required_cols <- c("QTLname", "SNP")
missing_cols <- setdiff(required_cols, colnames(slope_qtl))
if (length(missing_cols) > 0) {
  cat("  警告: QTL数据缺少列:", paste(missing_cols, collapse=", "), "\n")
  # 尝试查找SNP列的其他可能名称
  snp_candidates <- c("SNP", "snp", "Marker", "marker", "ID", "id")
  for (candidate in snp_candidates) {
    if (candidate %in% colnames(slope_qtl) && !"SNP" %in% colnames(slope_qtl)) {
      colnames(slope_qtl)[colnames(slope_qtl) == candidate] <- "SNP"
      cat(sprintf("  将列 '%s' 重命名为 'SNP'\n", candidate))
    }
  }
}

# 检查QTLname列
if (!"QTLname" %in% colnames(slope_qtl)) {
  qtlname_candidates <- c("QTLname", "QTL", "qtl", "Name", "name")
  for (candidate in qtlname_candidates) {
    if (candidate %in% colnames(slope_qtl)) {
      colnames(slope_qtl)[colnames(slope_qtl) == candidate] <- "QTLname"
      cat(sprintf("  将列 '%s' 重命名为 'QTLname'\n", candidate))
      break
    }
  }
}

# 提取唯一的QTL-SNP组合
all_effects <- unique(slope_qtl[, c("QTLname", "SNP")])

# 添加Chromosome和Position列（如果存在）
if ("Chromosome" %in% colnames(slope_qtl)) {
  all_effects$Chromosome <- slope_qtl$Chromosome[match(all_effects$SNP, slope_qtl$SNP)]
}
if ("Position" %in% colnames(slope_qtl)) {
  all_effects$Position <- slope_qtl$Position[match(all_effects$SNP, slope_qtl$SNP)]
}

# 5. 从各环境文件中提取效应值
cat("\n5. 从各环境文件中提取效应值...\n")
env_dirs <- c("2024NY", "2024PY", "2024YL", "2024ZMD",
              "2025YLW", "2025YLZ", "2025ZMDW", "2025ZMDZ")

for (env in env_dirs) {
  file_path <- file.path(base_env_path, env,
                         "GAPIT.Association.GWAS_Results.TKW.BLINK.Trait(NYC).csv")
  if (file.exists(file_path)) {
    env_data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # 检查Effect列
    if ("Effect" %in% colnames(env_data)) {
      # 检查SNP列名
      if (!"SNP" %in% colnames(env_data)) {
        snp_candidates <- c("SNP", "snp", "Marker", "marker", "ID", "id", "rs")
        for (candidate in snp_candidates) {
          if (candidate %in% colnames(env_data)) {
            colnames(env_data)[colnames(env_data) == candidate] <- "SNP"
            break
          }
        }
      }
      
      if ("SNP" %in% colnames(env_data)) {
        env_effects <- env_data[, c("SNP", "Effect")]
        colnames(env_effects)[2] <- env
        all_effects <- merge(all_effects, env_effects, by = "SNP", all.x = TRUE)
        cat(sprintf("  环境 %s: 读取 %d 个SNP效应值\n", env, sum(!is.na(env_effects$Effect))))
      } else {
        cat(sprintf("  环境 %s: 文件中没有SNP列\n", env))
        all_effects[[env]] <- NA
      }
    } else {
      cat(sprintf("  环境 %s: 文件中没有Effect列\n", env))
      all_effects[[env]] <- NA
    }
  } else {
    cat(sprintf("  环境 %s: 文件不存在\n", env))
    all_effects[[env]] <- NA
  }
}

# 6. 转换成长格式并合并环境协变量
cat("\n6. 数据转换和合并...\n")

long_format <- data.frame()
for (i in 1:nrow(all_effects)) {
  row_data <- all_effects[i, ]
  qtl_name <- row_data$QTLname
  
  for (env in env_dirs) {
    effect_val <- row_data[[env]]
    if (!is.na(effect_val)) {
      long_format <- rbind(long_format, data.frame(
        QTLname = qtl_name,
        env_code = env,
        Effect = effect_val,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("  转换后的数据: %d 行\n", nrow(long_format)))

# 合并环境协变量
env_xvar <- ecs_data[, c("env_code", "Xvar")]
long_format <- merge(long_format, env_xvar, by = "env_code", all.x = TRUE)
cat(sprintf("  合并环境协变量后: %d 行\n", nrow(long_format)))

# 7. 按QTL和环境分组计算平均效应
cat("\n7. 按QTL和环境分组计算平均效应...\n")

qtl_summary <- data.frame()
if (nrow(long_format) > 0) {
  for (qtl in unique(long_format$QTLname)) {
    qtl_data <- long_format[long_format$QTLname == qtl, ]
    for (env in unique(qtl_data$env_code)) {
      env_data <- qtl_data[qtl_data$env_code == env, ]
      if (nrow(env_data) > 0) {
        mean_effect <- mean(env_data$Effect, na.rm = TRUE)
        sd_effect <- sd(env_data$Effect, na.rm = TRUE)
        n_snps <- nrow(env_data)
        xvar <- unique(env_data$Xvar)
        
        if (length(xvar) > 0 && !is.na(xvar[1])) {
          qtl_summary <- rbind(qtl_summary, data.frame(
            QTLname = qtl,
            env_code = env,
            Xvar = xvar[1],
            mean_effect = mean_effect,
            sd_effect = sd_effect,
            n_snps = n_snps,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

cat(sprintf("  汇总数据: %d 行, %d 个QTL\n", nrow(qtl_summary), length(unique(qtl_summary$QTLname))))

# 8. 对每个QTL进行线性回归分析
cat("\n8. 进行线性回归分析...\n")

regression_results <- data.frame()

if (nrow(qtl_summary) > 0) {
  qtl_list <- unique(qtl_summary$QTLname)
  for (qtl in qtl_list) {
    qtl_data <- qtl_summary[qtl_summary$QTLname == qtl, ]
    
    if (nrow(qtl_data) >= 3) {  # 至少需要3个点才能拟合
      # 检查Xvar是否有变异
      if (sd(qtl_data$Xvar, na.rm = TRUE) > 0) {
        # 进行线性回归
        lm_model <- lm(mean_effect ~ Xvar, data = qtl_data)
        lm_summary <- summary(lm_model)
        
        # 提取系数
        coefficients <- coef(lm_summary)
        
        # 提取回归统计量
        r_squared <- lm_summary$r.squared
        adj_r_squared <- lm_summary$adj.r.squared
        
        if (!is.null(lm_summary$fstatistic)) {
          f_statistic <- lm_summary$fstatistic[1]
          df_model <- lm_summary$fstatistic[2]
          df_residual <- lm_summary$fstatistic[3]
          p_value <- pf(f_statistic, df_model, df_residual, lower.tail = FALSE)
        } else {
          f_statistic <- NA
          df_model <- NA
          df_residual <- NA
          p_value <- NA
        }
        
        # 提取斜率和截距的详细信息
        if ("(Intercept)" %in% rownames(coefficients)) {
          intercept <- coefficients["(Intercept)", "Estimate"]
          intercept_se <- coefficients["(Intercept)", "Std. Error"]
          intercept_t <- coefficients["(Intercept)", "t value"]
          intercept_p <- coefficients["(Intercept)", "Pr(>|t|)"]
        } else {
          intercept <- NA
          intercept_se <- NA
          intercept_t <- NA
          intercept_p <- NA
        }
        
        if ("Xvar" %in% rownames(coefficients)) {
          slope <- coefficients["Xvar", "Estimate"]
          slope_se <- coefficients["Xvar", "Std. Error"]
          slope_t <- coefficients["Xvar", "t value"]
          slope_p <- coefficients["Xvar", "Pr(>|t|)"]
        } else {
          slope <- NA
          slope_se <- NA
          slope_t <- NA
          slope_p <- NA
        }
        
        # 计算置信区间（95%）
        conf_int <- tryCatch({
          confint(lm_model, level = 0.95)
        }, error = function(e) {
          matrix(NA, nrow = 2, ncol = 2, 
                 dimnames = list(c("(Intercept)", "Xvar"), c("2.5 %", "97.5 %")))
        })
        
        # 存储结果
        regression_results <- rbind(regression_results, data.frame(
          QTLname = qtl,
          n_points = nrow(qtl_data),
          # 截距信息
          Intercept = intercept,
          Intercept_SE = intercept_se,
          Intercept_CI_lower = if ("(Intercept)" %in% rownames(conf_int)) conf_int["(Intercept)", 1] else NA,
          Intercept_CI_upper = if ("(Intercept)" %in% rownames(conf_int)) conf_int["(Intercept)", 2] else NA,
          Intercept_t = intercept_t,
          Intercept_p = intercept_p,
          # 斜率信息
          Slope = slope,
          Slope_SE = slope_se,
          Slope_CI_lower = if ("Xvar" %in% rownames(conf_int)) conf_int["Xvar", 1] else NA,
          Slope_CI_upper = if ("Xvar" %in% rownames(conf_int)) conf_int["Xvar", 2] else NA,
          Slope_t = slope_t,
          Slope_p = slope_p,
          # 模型整体信息
          R_squared = r_squared,
          Adjusted_R_squared = adj_r_squared,
          F_statistic = f_statistic,
          DF_model = df_model,
          DF_residual = df_residual,
          Model_p_value = p_value,
          stringsAsFactors = FALSE
        ))
        
        cat(sprintf("  QTL: %-15s | 斜率: %8.4f | R²: %6.4f | p值: %8.4f\n", 
                    qtl, slope, r_squared, p_value))
      } else {
        cat(sprintf("  QTL: %-15s | Xvar无变异，跳过回归分析\n", qtl))
      }
    } else {
      cat(sprintf("  QTL: %-15s | 数据点不足 (n=%d)，跳过回归分析\n", qtl, nrow(qtl_data)))
    }
  }
} else {
  cat("  警告: 没有汇总数据可用于回归分析\n")
}

# 9. 计算残差分析
cat("\n9. 计算残差分析...\n")

residual_analysis <- data.frame()
if (nrow(regression_results) > 0) {
  for (i in 1:nrow(regression_results)) {
    qtl <- regression_results$QTLname[i]
    qtl_data <- qtl_summary[qtl_summary$QTLname == qtl, ]
    
    if (nrow(qtl_data) >= 3) {
      lm_model <- lm(mean_effect ~ Xvar, data = qtl_data)
      residuals <- resid(lm_model)
      
      # 计算残差统计量
      residual_stats <- data.frame(
        QTLname = qtl,
        Residual_mean = mean(residuals),
        Residual_sd = sd(residuals),
        Residual_min = min(residuals),
        Residual_max = max(residuals),
        stringsAsFactors = FALSE
      )
      
      # 尝试计算Durbin-Watson检验（自相关）
      if (requireNamespace("car", quietly = TRUE)) {
        tryCatch({
          dw_test <- car::durbinWatsonTest(lm_model)
          residual_stats$DW_statistic <- dw_test$dw
          residual_stats$DW_p_value <- dw_test$p
        }, error = function(e) {
          residual_stats$DW_statistic <- NA
          residual_stats$DW_p_value <- NA
        })
      } else {
        residual_stats$DW_statistic <- NA
        residual_stats$DW_p_value <- NA
      }
      
      # Shapiro-Wilk正态性检验
      if (length(residuals) >= 3 && length(residuals) <= 5000) {
        tryCatch({
          shapiro_result <- shapiro.test(residuals)
          residual_stats$Shapiro_W <- shapiro_result$statistic
          residual_stats$Shapiro_p <- shapiro_result$p.value
        }, error = function(e) {
          residual_stats$Shapiro_W <- NA
          residual_stats$Shapiro_p <- NA
        })
      } else {
        residual_stats$Shapiro_W <- NA
        residual_stats$Shapiro_p <- NA
      }
      
      residual_analysis <- rbind(residual_analysis, residual_stats)
    }
  }
}

# 10. 计算每个QTL的效应值范围和环境响应
cat("\n10. 计算QTL效应范围和环境响应...\n")

qtl_range_analysis <- data.frame()
if (nrow(qtl_summary) > 0) {
  for (qtl in unique(qtl_summary$QTLname)) {
    qtl_data <- qtl_summary[qtl_summary$QTLname == qtl, ]
    
    if (nrow(qtl_data) > 0) {
      # 计算效应值范围
      effect_range <- max(qtl_data$mean_effect) - min(qtl_data$mean_effect)
      effect_sd <- sd(qtl_data$mean_effect)
      effect_mean <- mean(qtl_data$mean_effect)
      effect_cv <- ifelse(effect_mean != 0, effect_sd / abs(effect_mean) * 100, NA)
      
      # 计算环境协变量范围
      xvar_range <- max(qtl_data$Xvar) - min(qtl_data$Xvar)
      
      # 计算单位环境变化引起的效应变化（如果有斜率）
      slope_info <- regression_results[regression_results$QTLname == qtl, ]
      if (nrow(slope_info) > 0) {
        slope_value <- slope_info$Slope[1]
        effect_change_per_unit <- slope_value * xvar_range
      } else {
        slope_value <- NA
        effect_change_per_unit <- NA
      }
      
      qtl_range_analysis <- rbind(qtl_range_analysis, data.frame(
        QTLname = qtl,
        n_environments = nrow(qtl_data),
        mean_effect_overall = effect_mean,
        min_effect = min(qtl_data$mean_effect),
        max_effect = max(qtl_data$mean_effect),
        effect_range = effect_range,
        effect_sd = effect_sd,
        effect_cv_percent = effect_cv,
        Xvar_range = xvar_range,
        slope_if_available = slope_value,
        effect_change_over_Xrange = effect_change_per_unit,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 11. 整合所有结果并保存
cat("\n11. 整合和保存结果...\n")

if (nrow(regression_results) > 0) {
  # 合并所有结果
  if (nrow(residual_analysis) > 0) {
    final_results <- merge(regression_results, residual_analysis, by = "QTLname", all = TRUE)
  } else {
    final_results <- regression_results
  }
  
  if (nrow(qtl_range_analysis) > 0) {
    final_results <- merge(final_results, qtl_range_analysis, by = "QTLname", all = TRUE)
  }
  
  # 添加显著性标记
  final_results$Significance <- ifelse(final_results$Model_p_value < 0.001, "***",
                                       ifelse(final_results$Model_p_value < 0.01, "**",
                                              ifelse(final_results$Model_p_value < 0.05, "*", "NS")))
  
  final_results$Slope_significance <- ifelse(final_results$Slope_p < 0.001, "***",
                                             ifelse(final_results$Slope_p < 0.01, "**",
                                                    ifelse(final_results$Slope_p < 0.05, "*", "NS")))
  
  # 排序结果（按R²从高到低）
  final_results <- final_results[order(-final_results$R_squared, na.last = TRUE), ]
  
  # 保存所有结果到CSV文件
  output_file <- file.path(base_dir, "QTL_Effect_Regression_Analysis.csv")
  write.csv(final_results, output_file, row.names = FALSE, na = "")
  
  # 保存简化版本（主要结果）
  keep_cols <- c("QTLname", "n_points", "Slope", "Slope_SE", "Slope_p", "Slope_significance",
                 "Intercept", "R_squared", "Adjusted_R_squared", "Model_p_value", "Significance",
                 "effect_range", "effect_cv_percent")
  
  keep_cols <- keep_cols[keep_cols %in% colnames(final_results)]
  simplified_results <- final_results[, keep_cols, drop = FALSE]
  
  simplified_file <- file.path(base_dir, "QTL_Effect_Regression_Summary.csv")
  write.csv(simplified_results, simplified_file, row.names = FALSE, na = "")
  
  # 12. 输出关键统计摘要
  cat("\n=== 图d拟合统计量计算完成 ===\n")
  cat(sprintf("分析完成！共分析了 %d 个QTL\n", nrow(final_results)))
  cat(sprintf("结果已保存至:\n"))
  cat(sprintf("  完整结果: %s\n", output_file))
  cat(sprintf("  简化摘要: %s\n", simplified_file))
  
  cat("\n=== 主要发现 ===\n")
  
  if ("Significance" %in% colnames(final_results)) {
    cat("1. 回归模型显著性:\n")
    sig_counts <- table(final_results$Significance, useNA = "ifany")
    for (sig in names(sig_counts)) {
      cat(sprintf("   %s: %d 个QTL\n", sig, sig_counts[sig]))
    }
  }
  
  if ("Slope_significance" %in% colnames(final_results)) {
    cat("\n2. 斜率显著性:\n")
    slope_sig_counts <- table(final_results$Slope_significance, useNA = "ifany")
    for (sig in names(slope_sig_counts)) {
      cat(sprintf("   %s: %d 个QTL\n", sig, slope_sig_counts[sig]))
    }
  }
  
  if ("R_squared" %in% colnames(final_results)) {
    cat("\n3. 决定系数分布:\n")
    r2_vals <- final_results$R_squared[!is.na(final_results$R_squared)]
    if (length(r2_vals) > 0) {
      r2_summary <- summary(r2_vals)
      cat(sprintf("   最小值: %.4f\n", r2_summary[1]))
      cat(sprintf("   中位数: %.4f\n", r2_summary[3]))
      cat(sprintf("   平均值: %.4f\n", r2_summary[4]))
      cat(sprintf("   最大值: %.4f\n", r2_summary[6]))
    } else {
      cat("   无有效的R²值\n")
    }
  }
  
  if ("Slope" %in% colnames(final_results)) {
    cat("\n4. 斜率范围:\n")
    slope_vals <- final_results$Slope[!is.na(final_results$Slope)]
    if (length(slope_vals) > 0) {
      slope_summary <- summary(slope_vals)
      cat(sprintf("   最小值: %.4f\n", slope_summary[1]))
      cat(sprintf("   中位数: %.4f\n", slope_summary[3]))
      cat(sprintf("   平均值: %.4f\n", slope_summary[4]))
      cat(sprintf("   最大值: %.4f\n", slope_summary[6]))
    } else {
      cat("   无有效的斜率值\n")
    }
  }
  
  # 找出最显著和最不显著的QTL
  if (nrow(final_results) > 0 && "Model_p_value" %in% colnames(final_results)) {
    valid_rows <- !is.na(final_results$Model_p_value)
    if (sum(valid_rows) > 0) {
      final_results_valid <- final_results[valid_rows, ]
      
      if (nrow(final_results_valid) > 0) {
        most_sig <- final_results_valid[which.min(final_results_valid$Model_p_value), ]
        
        cat("\n5. 关键QTL:\n")
        cat(sprintf("   最显著QTL: %s (p=%.2e", most_sig$QTLname, most_sig$Model_p_value))
        if ("R_squared" %in% colnames(most_sig)) {
          cat(sprintf(", R²=%.3f)\n", most_sig$R_squared))
        } else {
          cat(")\n")
        }
      }
    }
  }
} else {
  cat("警告：没有回归结果可保存\n")
}

cat("\n=== 分析完成 ===\n")

# 输出数据概览
cat("\n=== 数据概览 ===\n")
cat(sprintf("QTL数据行数: %d\n", nrow(slope_qtl)))
cat(sprintf("环境协变量数据行数: %d\n", nrow(ecs_data)))
cat(sprintf("合并后的长格式数据行数: %d\n", nrow(long_format)))
cat(sprintf("汇总数据行数: %d\n", nrow(qtl_summary)))
cat(sprintf("回归分析结果数: %d\n", nrow(regression_results)))