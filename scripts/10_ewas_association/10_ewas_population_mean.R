# 完整修复版 - 修复select错误和图形警告
library(dplyr)
library(ggplot2)
library(tidyr)

# 基础设置
data_dir <- "C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\01缺失的遗传力\\"
pheno_file <- "TKW_mean_table.txt"
env_file <- "ECs_results.csv"
output_dir <- paste0(data_dir, "EWAS_robust_analysis\\")

# 要删除的环境列表
envs_to_remove <- c()

# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. 加载数据并删除指定环境
load_and_filter_data <- function(data_dir, pheno_file, env_file, envs_to_remove) {
  cat("📊 加载数据并删除指定环境...\n")
  
  # 构建完整文件路径
  pheno_full_path <- file.path(data_dir, pheno_file)
  env_full_path <- file.path(data_dir, env_file)
  
  cat("   表型文件路径:", pheno_full_path, "\n")
  cat("   环境因子文件路径:", env_full_path, "\n")
  
  # 检查文件是否存在
  if (!file.exists(pheno_full_path)) {
    stop("❌ 表型文件不存在: ", pheno_full_path)
  }
  if (!file.exists(env_full_path)) {
    stop("❌ 环境因子文件不存在: ", env_full_path)
  }
  
  # 加载表型数据
  pheno <- read.table(pheno_full_path, header = TRUE, sep = "\t", 
                     stringsAsFactors = FALSE, na.strings = c("", "NA"))
  
  # 显示表型数据的前几行
  cat("   表型数据维度:", dim(pheno), "\n")
  cat("   表型数据列名:", paste(colnames(pheno), collapse = ", "), "\n")
  
  # 重命名genotype列为line_code
  pheno <- pheno %>% rename(line_code = genotype)
  
  # 确认表型列名
  pheno_col_name <- "TKW"
  
  # 检查表型列是否存在
  if (!pheno_col_name %in% colnames(pheno)) {
    cat("   警告: 列名", pheno_col_name, "不存在。可用列名:", paste(colnames(pheno), collapse = ", "), "\n")
    possible_pheno_cols <- setdiff(colnames(pheno), c("line_code", "env_code", "ID", "id"))
    if (length(possible_pheno_cols) > 0) {
      pheno_col_name <- possible_pheno_cols[1]
      cat("   自动选择表型列:", pheno_col_name, "\n")
    } else {
      stop("❌ 无法识别表型列")
    }
  }
  
  # 有效株系筛选（至少2个环境有数据）
  valid_lines <- pheno %>%
    group_by(line_code) %>%
    summarise(n_env = sum(!is.na(.data[[pheno_col_name]]))) %>%
    filter(n_env >= 2) %>%
    pull(line_code)
  
  pheno_filtered <- pheno %>% filter(line_code %in% valid_lines)
  
  # 计算环境表型均值（应用环境删除）
  env_means <- pheno_filtered %>%
    filter(!env_code %in% envs_to_remove) %>%
    group_by(env_code) %>%
    summarise(pheno_mean = mean(.data[[pheno_col_name]], na.rm = TRUE),
              n_lines = n(),
              pheno_sd = sd(.data[[pheno_col_name]], na.rm = TRUE)) %>%
    filter(!is.na(pheno_mean)) %>%
    arrange(env_code)
  
  # 加载环境因子数据
  env_factors <- read.csv(env_full_path, 
                         header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  cat("   环境因子数据维度:", dim(env_factors), "\n")
  cat("   环境因子数据列数:", ncol(env_factors), "\n")
  
  # 保持原始列名
  original_colnames <- colnames(env_factors)
  env_col <- original_colnames[1]
  env_factors <- env_factors %>% rename(env_code = !!env_col)
  
  # 环境因子数据也要删除相同的环境
  env_factors_filtered <- env_factors %>%
    filter(!env_code %in% envs_to_remove) %>%
    arrange(env_code)
  
  # 环境匹配
  common_envs <- intersect(env_means$env_code, env_factors_filtered$env_code)
  
  # 最终筛选和排序
  env_means_final <- env_means %>% 
    filter(env_code %in% common_envs) %>%
    arrange(env_code)
  
  env_factors_final <- env_factors_filtered %>% 
    filter(env_code %in% common_envs) %>%
    arrange(env_code)
  
  cat("   初始环境数:", length(unique(pheno_filtered$env_code)), "\n")
  cat("   删除环境数:", length(envs_to_remove), "\n")
  cat("   剩余环境数:", length(common_envs), "\n")
  cat("   剩余环境:", paste(common_envs, collapse = ", "), "\n")
  cat("   环境因子数:", ncol(env_factors_final) - 1, "\n")
  
  return(list(
    env_means = env_means_final,
    env_factors = env_factors_final,
    env_factors_original = env_factors,
    common_envs = common_envs,
    pheno_col_name = pheno_col_name
  ))
}

# 2. 预筛选环境因子
prefilter_environmental_factors <- function(env_factors, env_means, output_dir) {
  cat("\n🔍 预筛选环境因子（过滤缺失值和零方差）...\n")
  
  # 合并数据以检查
  combined_data <- merge(env_means, env_factors, by = "env_code")
  
  # 提取环境因子列名
  factor_cols <- setdiff(colnames(env_factors), "env_code")
  
  # 记录筛选结果
  filter_results <- data.frame(
    factor = character(),
    total_envs = integer(),
    missing_count = integer(),
    missing_rate = numeric(),
    zero_variance = logical(),
    keep = logical(),
    stringsAsFactors = FALSE
  )
  
  valid_factors <- c()
  
  cat("   检查", length(factor_cols), "个环境因子...\n")
  pb <- txtProgressBar(min = 0, max = length(factor_cols), style = 3)
  
  for (i in 1:length(factor_cols)) {
    factor <- factor_cols[i]
    
    tryCatch({
      # 获取该因子的值
      factor_values <- combined_data[[factor]]
      
      # 转换为数值型
      factor_values_numeric <- suppressWarnings(as.numeric(as.character(factor_values)))
      
      # 计算缺失值情况
      missing_count <- sum(is.na(factor_values_numeric))
      total_envs <- length(factor_values_numeric)
      
      # 防止除以零
      if (total_envs > 0) {
        missing_rate <- missing_count / total_envs
      } else {
        missing_rate <- 1
      }
      
      # 检查方差
      non_missing_values <- factor_values_numeric[!is.na(factor_values_numeric)]
      zero_var <- FALSE
      
      if (length(non_missing_values) >= 2) {
        variance <- var(non_missing_values, na.rm = TRUE)
        zero_var <- is.na(variance) || variance == 0
      } else {
        zero_var <- TRUE
      }
      
      # 筛选条件
      condition1 <- !is.na(missing_rate) && missing_rate <= 0.5
      condition2 <- !is.na(length(non_missing_values)) && length(non_missing_values) >= 3
      condition3 <- !is.na(zero_var) && !zero_var
      
      keep_factor <- condition1 && condition2 && condition3
      
      # 记录结果
      result_row <- data.frame(
        factor = factor,
        total_envs = total_envs,
        missing_count = missing_count,
        missing_rate = round(missing_rate, 3),
        zero_variance = zero_var,
        keep = keep_factor,
        stringsAsFactors = FALSE
      )
      
      filter_results <- rbind(filter_results, result_row)
      
      if (!is.na(keep_factor) && keep_factor) {
        valid_factors <- c(valid_factors, factor)
      }
      
    }, error = function(e) {
      result_row <- data.frame(
        factor = factor,
        total_envs = length(combined_data[[factor]]),
        missing_count = NA_integer_,
        missing_rate = NA_real_,
        zero_variance = NA,
        keep = FALSE,
        stringsAsFactors = FALSE
      )
      
      filter_results <- rbind(filter_results, result_row)
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 保存筛选结果
  write.csv(filter_results, 
            file.path(output_dir, "environment_factor_prefilter_results.csv"),
            row.names = FALSE)
  
  # 统计信息
  cat("\n   预筛选统计:\n")
  cat("   总环境因子数:", length(factor_cols), "\n")
  cat("   保留的环境因子数:", length(valid_factors), "\n")
  cat("   过滤的环境因子数:", length(factor_cols) - length(valid_factors), "\n")
  
  filtered_by_missing <- sum(filter_results$missing_rate > 0.5, na.rm = TRUE)
  filtered_by_samples <- sum((filter_results$total_envs - filter_results$missing_count) < 3, na.rm = TRUE)
  filtered_by_variance <- sum(filter_results$zero_variance, na.rm = TRUE)
  
  cat("   过滤原因:\n")
  cat("     - 缺失率>50%:", filtered_by_missing, "\n")
  cat("     - 非缺失值<3个:", filtered_by_samples, "\n")
  cat("     - 零方差:", filtered_by_variance, "\n")
  
  return(valid_factors)
}

# 3. 筛选有效环境因子
filter_valid_factors <- function(env_factors, valid_factors) {
  cols_to_keep <- c("env_code", valid_factors)
  env_factors_valid <- env_factors[, cols_to_keep, drop = FALSE]
  return(env_factors_valid)
}

# 4. 批量计算相关性
calculate_correlations_batch <- function(env_means, env_factors, output_dir) {
  cat("\n📈 批量计算环境因子与表型的相关性...\n")
  
  # 合并环境均值与环境因子数据
  combined_data <- merge(env_means, env_factors, by = "env_code")
  
  # 提取环境因子列名
  factor_cols <- setdiff(colnames(env_factors), "env_code")
  
  cat("   将分析", length(factor_cols), "个环境因子...\n")
  
  # 初始化结果数据框
  results <- data.frame(
    factor = character(0),
    n_env = integer(0),
    pearson_cor = numeric(0),
    pearson_p = numeric(0),
    spearman_cor = numeric(0),
    spearman_p = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # 创建显著结果目录
  significant_dir <- file.path(output_dir, "significant_correlations")
  if (!dir.exists(significant_dir)) dir.create(significant_dir)
  
  # 进度条
  pb <- txtProgressBar(min = 0, max = length(factor_cols), style = 3)
  
  for (i in 1:length(factor_cols)) {
    factor <- factor_cols[i]
    
    tryCatch({
      # 检查是否有缺失值
      complete_cases <- complete.cases(combined_data$pheno_mean, combined_data[[factor]])
      n_complete <- sum(complete_cases)
      
      if (n_complete >= 3) {
        # 转换为数值型
        factor_values <- suppressWarnings(as.numeric(as.character(combined_data[[factor]])))
        
        # 检查数据有效性
        non_missing_vals <- factor_values[complete_cases]
        if (length(non_missing_vals) >= 3 && var(non_missing_vals, na.rm = TRUE) > 0) {
          
          # 计算Pearson相关性
          pearson_result <- tryCatch({
            cor.test(combined_data$pheno_mean[complete_cases], 
                     factor_values[complete_cases], 
                     method = "pearson")
          }, error = function(e) NULL)
          
          # 计算Spearman相关性
          spearman_result <- tryCatch({
            cor.test(combined_data$pheno_mean[complete_cases], 
                     factor_values[complete_cases], 
                     method = "spearman", exact = FALSE)
          }, error = function(e) NULL)
          
          # 保存结果
          if (!is.null(pearson_result) && !is.null(spearman_result)) {
            result_row <- data.frame(
              factor = factor,
              n_env = n_complete,
              pearson_cor = pearson_result$estimate,
              pearson_p = pearson_result$p.value,
              spearman_cor = spearman_result$estimate,
              spearman_p = spearman_result$p.value,
              stringsAsFactors = FALSE
            )
            
            results <- rbind(results, result_row)
            
            # 如果p值<0.05，绘制散点图
            if (pearson_result$p.value < 0.05 || spearman_result$p.value < 0.05) {
              plot_data <- data.frame(
                x = factor_values[complete_cases],
                y = combined_data$pheno_mean[complete_cases]
              )
              
              plot_data <- plot_data[is.finite(plot_data$x) & is.finite(plot_data$y), ]
              
              if (nrow(plot_data) >= 3) {
                clean_factor_name <- gsub("[^[:alnum:]]", "_", factor)
                plot_file <- file.path(significant_dir, paste0("correlation_", clean_factor_name, ".png"))
                
                p <- ggplot(plot_data, aes(x = x, y = y)) +
                  geom_point(size = 3, alpha = 0.7) +
                  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
                  labs(title = paste0("Correlation between ", factor, " and Phenotype Mean"),
                       x = factor,
                       y = "Phenotype Mean",
                       subtitle = sprintf("Pearson r = %.3f (p = %.3e), Spearman ρ = %.3f (p = %.3e)",
                                         pearson_result$estimate, pearson_result$p.value,
                                         spearman_result$estimate, spearman_result$p.value)) +
                  theme_minimal() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))
                
                ggsave(plot_file, p, width = 8, height = 6, dpi = 300)
              }
            }
          }
        }
      }
    }, error = function(e) {
      # 静默处理错误
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 按Pearson相关性绝对值排序并重置行名
  if (nrow(results) > 0) {
    results <- results %>% 
      arrange(desc(abs(pearson_cor))) %>%
      as.data.frame()  # 确保是普通数据框
    
    # 重置行名
    rownames(results) <- NULL
  }
  
  cat("   完成相关性分析，共成功分析", nrow(results), "个环境因子\n")
  
  return(results)
}

# 5. 结果筛选和多重比较校正 - 修复版本
filter_and_report_results <- function(results) {
  cat("\n🎯 筛选显著相关性结果并应用多重比较校正...\n")
  
  if (nrow(results) == 0) {
    cat("   没有分析结果可用于筛选\n")
    return(data.frame())
  }
  
  # 重置行名以避免问题
  rownames(results) <- NULL
  
  # 应用多重比较校正（FDR校正）
  results$pearson_p_adj <- p.adjust(results$pearson_p, method = "BH")
  results$spearman_p_adj <- p.adjust(results$spearman_p, method = "BH")
  
  # 定义显著标准（FDR < 0.05）
  significant_results <- results %>%
    filter(pearson_p_adj < 0.05 | spearman_p_adj < 0.05)
  
  # 重置行名
  rownames(significant_results) <- NULL
  
  # 打印前20个最显著的结果
  if (nrow(significant_results) > 0) {
    cat("   发现", nrow(significant_results), "个显著相关环境因子（FDR < 0.05）:\n")
    
    # 创建要显示的数据框
    display_df <- significant_results
    
    # 选择要显示的列
    cols_to_display <- c("factor", "n_env", "pearson_cor", "pearson_p", "pearson_p_adj", 
                         "spearman_cor", "spearman_p", "spearman_p_adj")
    
    # 只保留实际存在的列
    cols_to_display <- cols_to_display[cols_to_display %in% colnames(display_df)]
    
    if (length(cols_to_display) > 0) {
      display_df <- display_df[, cols_to_display, drop = FALSE]
      
      # 只显示前20行
      if (nrow(display_df) > 20) {
        display_df <- display_df[1:20, , drop = FALSE]
      }
      
      # 打印结果
      print(display_df)
    } else {
      cat("   警告：没有找到预期的列名\n")
      print(colnames(significant_results))
    }
  } else {
    # 如果没有FDR显著的，显示名义上显著的结果（p < 0.05）
    nominally_significant <- results %>%
      filter(pearson_p < 0.05 | spearman_p < 0.05)
    
    rownames(nominally_significant) <- NULL
    
    if (nrow(nominally_significant) > 0) {
      cat("   未发现FDR显著相关，但发现", nrow(nominally_significant), "个名义上显著的环境因子（p < 0.05）:\n")
      
      # 创建要显示的数据框
      display_df <- nominally_significant
      
      # 选择要显示的列
      cols_to_display <- c("factor", "n_env", "pearson_cor", "pearson_p", "pearson_p_adj", 
                           "spearman_cor", "spearman_p", "spearman_p_adj")
      
      # 只保留实际存在的列
      cols_to_display <- cols_to_display[cols_to_display %in% colnames(display_df)]
      
      if (length(cols_to_display) > 0) {
        display_df <- display_df[, cols_to_display, drop = FALSE]
        
        # 只显示前20行
        if (nrow(display_df) > 20) {
          display_df <- display_df[1:20, , drop = FALSE]
        }
        
        # 打印结果
        print(display_df)
      } else {
        cat("   警告：没有找到预期的列名\n")
        print(colnames(nominally_significant))
      }
      
      # 返回名义上显著的结果
      return(nominally_significant)
    } else {
      cat("   未发现任何显著相关环境因子（即使名义上也不显著）\n")
      return(data.frame())
    }
  }
  
  return(significant_results)
}

# 6. 保存分析结果
save_analysis_results <- function(results, significant_results, data_info, output_dir) {
  cat("\n💾 保存分析结果...\n")
  
  # 重置行名
  if (nrow(results) > 0) {
    rownames(results) <- NULL
  }
  
  if (nrow(significant_results) > 0) {
    rownames(significant_results) <- NULL
  }
  
  # 保存完整结果
  if (nrow(results) > 0) {
    write.csv(results, 
              file.path(output_dir, "complete_correlation_results.csv"),
              row.names = FALSE, quote = FALSE)
    cat("   完整结果已保存: complete_correlation_results.csv\n")
    
    # 保存显著结果
    if (nrow(significant_results) > 0) {
      write.csv(significant_results, 
                file.path(output_dir, "significant_correlation_results.csv"),
                row.names = FALSE, quote = FALSE)
      cat("   显著结果已保存: significant_correlation_results.csv\n")
    }
    
    # 保存数据汇总信息
    summary_info <- data.frame(
      item = c("初始环境数", "删除环境数", "最终环境数", 
               "总环境因子数", "预筛选后因子数", "成功分析因子数", 
               "名义显著因子数(p<0.05)", "FDR显著因子数(FDR<0.05)"),
      value = c(data_info$initial_env_count,
                length(data_info$envs_removed),
                length(data_info$common_envs),
                ncol(data_info$env_factors_original) - 1,
                length(data_info$valid_factors),
                nrow(results),
                sum(results$pearson_p < 0.05 | results$spearman_p < 0.05),
                nrow(significant_results))
    )
    
    write.csv(summary_info, 
              file.path(output_dir, "analysis_summary.csv"),
              row.names = FALSE)
    cat("   分析汇总已保存: analysis_summary.csv\n")
    
    # 保存前100个最相关的结果用于快速查看
    if (nrow(results) > 0) {
      top_100 <- results
      if (nrow(top_100) > 100) {
        top_100 <- top_100[1:100, , drop = FALSE]
      }
      
      write.csv(top_100,
                file.path(output_dir, "top_100_correlations.csv"),
                row.names = FALSE, quote = FALSE)
      cat("   前100个最相关结果已保存: top_100_correlations.csv\n")
    }
    
    # 保存R数据对象
    analysis_data <- list(
      results = results,
      significant_results = significant_results,
      summary = summary_info,
      env_info = list(
        removed_envs = data_info$envs_removed,
        common_envs = data_info$common_envs,
        valid_factors = data_info$valid_factors
      )
    )
    
    saveRDS(analysis_data, file.path(output_dir, "analysis_data.rds"))
    cat("   完整分析数据已保存: analysis_data.rds\n")
  } else {
    cat("   没有分析结果可保存\n")
  }
}

# 7. 生成相关性分布图 - 修复版本
generate_correlation_plots <- function(results, output_dir) {
  cat("\n📊 生成相关性分析图...\n")
  
  if (nrow(results) == 0) {
    cat("   没有数据生成图形\n")
    return()
  }
  
  # 1. Pearson相关性分布图
  # 移除非有限值
  pearson_data <- results[is.finite(results$pearson_cor), ]
  
  if (nrow(pearson_data) > 0) {
    p1 <- ggplot(pearson_data, aes(x = pearson_cor)) +
      geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7, na.rm = TRUE) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Distribution of Pearson Correlation Coefficients",
           x = "Pearson Correlation Coefficient",
           y = "Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path(output_dir, "pearson_correlation_distribution.png"), 
           p1, width = 8, height = 6, dpi = 300)
  }
  
  # 2. Spearman相关性分布图
  spearman_data <- results[is.finite(results$spearman_cor), ]
  
  if (nrow(spearman_data) > 0) {
    p2 <- ggplot(spearman_data, aes(x = spearman_cor)) +
      geom_histogram(bins = 50, fill = "lightgreen", color = "black", alpha = 0.7, na.rm = TRUE) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Distribution of Spearman Correlation Coefficients",
           x = "Spearman Correlation Coefficient",
           y = "Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path(output_dir, "spearman_correlation_distribution.png"), 
           p2, width = 8, height = 6, dpi = 300)
  }
  
  # 3. P值分布图
  valid_p_data <- results[is.finite(results$pearson_p) & is.finite(results$spearman_p), ]
  
  if (nrow(valid_p_data) > 0) {
    results_long <- valid_p_data %>%
      dplyr::select(factor, pearson_p, spearman_p) %>%
      pivot_longer(cols = c(pearson_p, spearman_p), 
                   names_to = "method", 
                   values_to = "p_value")
    
    # 移除无限值
    results_long <- results_long[is.finite(-log10(results_long$p_value)), ]
    
    if (nrow(results_long) > 0) {
      p3 <- ggplot(results_long, aes(x = -log10(p_value), fill = method)) +
        geom_histogram(alpha = 0.6, position = "identity", bins = 50, na.rm = TRUE) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
        labs(title = "Distribution of -log10(p-values)",
             x = "-log10(p-value)",
             y = "Frequency",
             fill = "Correlation Method") +
        scale_fill_manual(values = c("pearson_p" = "skyblue", "spearman_p" = "lightgreen")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom")
      
      ggsave(file.path(output_dir, "pvalue_distribution.png"), 
             p3, width = 8, height = 6, dpi = 300)
    }
  }
  
  # 4. 相关性散点图（Pearson vs Spearman）
  valid_corr_data <- results[is.finite(results$pearson_cor) & is.finite(results$spearman_cor), ]
  
  if (nrow(valid_corr_data) > 1) {
    p4 <- ggplot(valid_corr_data, aes(x = pearson_cor, y = spearman_cor)) +
      geom_point(alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Pearson vs Spearman Correlation Coefficients",
           x = "Pearson Correlation Coefficient",
           y = "Spearman Correlation Coefficient") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path(output_dir, "pearson_vs_spearman.png"), 
           p4, width = 8, height = 6, dpi = 300)
  }
  
  cat("   相关性分析图已保存\n")
}

# 8. 生成详细报告
generate_detailed_report <- function(results, significant_results, data_info, output_dir) {
  cat("\n📄 生成详细分析报告...\n")
  
  # 创建报告文件
  report_file <- file.path(output_dir, "detailed_analysis_report.txt")
  
  sink(report_file)
  
  cat("===================================================\n")
  cat("环境因子相关性分析详细报告\n")
  cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("===================================================\n\n")
  
  # 1. 数据概览
  cat("1. 数据概览\n")
  cat("   ==================================\n")
  cat("   初始环境数:", data_info$initial_env_count, "\n")
  cat("   删除环境数:", length(data_info$envs_removed), "\n")
  if (length(data_info$envs_removed) > 0) {
    cat("   删除的环境:", paste(data_info$envs_removed, collapse = ", "), "\n")
  }
  cat("   最终环境数:", length(data_info$common_envs), "\n")
  cat("   分析的环境:", paste(data_info$common_envs, collapse = ", "), "\n")
  cat("   总环境因子数:", ncol(data_info$env_factors_original) - 1, "\n")
  cat("   预筛选后因子数:", length(data_info$valid_factors), "\n")
  cat("   成功分析因子数:", nrow(results), "\n")
  cat("   分析成功率:", round(nrow(results)/(ncol(data_info$env_factors_original)-1)*100, 1), "%\n\n")
  
  # 2. 相关性分析结果概览
  cat("2. 相关性分析结果概览\n")
  cat("   ==================================\n")
  if (nrow(results) > 0) {
    valid_pearson <- results$pearson_cor[is.finite(results$pearson_cor)]
    valid_spearman <- results$spearman_cor[is.finite(results$spearman_cor)]
    
    if (length(valid_pearson) > 0) {
      cat("   Pearson相关性范围:", round(min(valid_pearson, na.rm = TRUE), 3), 
          "到", round(max(valid_pearson, na.rm = TRUE), 3), "\n")
      cat("   Pearson相关性中位数:", round(median(valid_pearson, na.rm = TRUE), 3), "\n")
    }
    
    if (length(valid_spearman) > 0) {
      cat("   Spearman相关性范围:", round(min(valid_spearman, na.rm = TRUE), 3), 
          "到", round(max(valid_spearman, na.rm = TRUE), 3), "\n")
      cat("   Spearman相关性中位数:", round(median(valid_spearman, na.rm = TRUE), 3), "\n")
    }
    
    # 计算正负相关比例
    pos_pearson <- sum(results$pearson_cor > 0, na.rm = TRUE)
    neg_pearson <- sum(results$pearson_cor < 0, na.rm = TRUE)
    total_pearson <- pos_pearson + neg_pearson
    
    if (total_pearson > 0) {
      cat("   正Pearson相关因子数:", pos_pearson, "(", round(pos_pearson/total_pearson*100, 1), "%)\n")
      cat("   负Pearson相关因子数:", neg_pearson, "(", round(neg_pearson/total_pearson*100, 1), "%)\n\n")
    }
    
    # 显著性统计
    cat("   名义显著因子数(p < 0.05):", sum(results$pearson_p < 0.05 | results$spearman_p < 0.05, na.rm = TRUE), "\n")
    cat("   其中Pearson显著:", sum(results$pearson_p < 0.05, na.rm = TRUE), "\n")
    cat("   其中Spearman显著:", sum(results$spearman_p < 0.05, na.rm = TRUE), "\n")
    cat("   FDR显著因子数(FDR < 0.05):", nrow(significant_results), "\n\n")
  } else {
    cat("   没有成功的相关性分析结果\n\n")
  }
  
  sink()
  
  cat("   详细报告已保存: detailed_analysis_report.txt\n")
}

# 9. 主函数
main <- function() {
  cat("==========================================\n")
  cat("🌍 环境因子相关性分析 - 完整修复版本\n")
  cat("   修复select错误和图形警告\n")
  cat("==========================================\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    # 加载数据并删除指定环境
    data <- load_and_filter_data(data_dir, pheno_file, env_file, envs_to_remove)
    
    # 预筛选环境因子
    valid_factors <- prefilter_environmental_factors(
      data$env_factors, 
      data$env_means, 
      output_dir
    )
    
    if (length(valid_factors) == 0) {
      cat("❌ 没有通过预筛选的环境因子，停止分析\n")
      return()
    }
    
    # 筛选有效环境因子
    env_factors_valid <- filter_valid_factors(data$env_factors, valid_factors)
    
    # 计算相关性
    cat("\n⚠️ 注意：接下来将分析", length(valid_factors), "个预筛选后的环境因子\n")
    cat("   这可能需要一些时间，请耐心等待...\n")
    
    results <- calculate_correlations_batch(data$env_means, env_factors_valid, output_dir)
    
    # 结果筛选和多重比较校正
    significant_results <- filter_and_report_results(results)
    
    # 生成图形
    generate_correlation_plots(results, output_dir)
    
    # 保存结果
    data_info <- list(
      initial_env_count = length(unique(data$env_factors_original$env_code)),
      envs_removed = envs_to_remove,
      common_envs = data$common_envs,
      env_factors_original = data$env_factors_original,
      valid_factors = valid_factors
    )
    
    save_analysis_results(results, significant_results, data_info, output_dir)
    
    # 生成详细报告
    generate_detailed_report(results, significant_results, data_info, output_dir)
    
    # 最终报告
    duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    cat("\n==========================================\n")
    cat("✅ 分析完成!\n")
    cat("   耗时:", duration, "分钟\n")
    cat("   初始环境:", data_info$initial_env_count, "个\n")
    cat("   删除环境:", length(envs_to_remove), "个\n")
    cat("   最终环境:", length(data$common_envs), "个\n")
    cat("   总环境因子:", ncol(data$env_factors_original) - 1, "个\n")
    cat("   预筛选后因子数:", length(valid_factors), "个\n")
    cat("   成功分析:", nrow(results), "个\n")
    
    if (nrow(results) > 0) {
      valid_pearson <- results$pearson_cor[is.finite(results$pearson_cor)]
      valid_spearman <- results$spearman_cor[is.finite(results$spearman_cor)]
      
      if (length(valid_pearson) > 0) {
        cat("   Pearson相关范围:", round(min(valid_pearson, na.rm = TRUE), 3), 
            "到", round(max(valid_pearson, na.rm = TRUE), 3), "\n")
      }
      
      if (length(valid_spearman) > 0) {
        cat("   Spearman相关范围:", round(min(valid_spearman, na.rm = TRUE), 3), 
            "到", round(max(valid_spearman, na.rm = TRUE), 3), "\n")
      }
      
      cat("   名义显著(p<0.05):", sum(results$pearson_p < 0.05 | results$spearman_p < 0.05, na.rm = TRUE), "个\n")
    }
    
    if (nrow(significant_results) > 0) {
      cat("   FDR显著(FDR<0.05):", nrow(significant_results), "个\n")
      cat("   最显著的因子:", significant_results$factor[1], "\n")
      cat("   Pearson r =", round(significant_results$pearson_cor[1], 3), 
          " (p.adj =", format(significant_results$pearson_p_adj[1], scientific = TRUE, digits = 3), ")\n")
    } else if (sum(results$pearson_p < 0.05 | results$spearman_p < 0.05, na.rm = TRUE) > 0) {
      cat("   ⚠️ 有", sum(results$pearson_p < 0.05 | results$spearman_p < 0.05, na.rm = TRUE), "个名义显著的因子(p<0.05)\n")
    }
    
    cat("   结果保存目录:", output_dir, "\n")
    cat("   生成的文件包括:\n")
    cat("     - 环境因子预筛选结果: environment_factor_prefilter_results.csv\n")
    cat("     - 完整相关性结果: complete_correlation_results.csv\n")
    cat("     - 显著相关性结果: significant_correlation_results.csv\n")
    cat("     - 前100个最相关结果: top_100_correlations.csv\n")
    cat("     - 分析汇总: analysis_summary.csv\n")
    cat("     - 详细分析报告: detailed_analysis_report.txt\n")
    cat("     - 各种统计图形\n")
    cat("==========================================\n")
    
  }, error = function(e) {
    cat("\n❌ 分析过程中出现错误:\n")
    cat("   错误信息:", e$message, "\n")
    cat("   调用堆栈:\n")
    print(e$call)
  })
}

# 执行分析
main()