# 表型分析代码（BLUE作为估计值与各环境比较）- 环境匹配修复版
# 修复环境过滤后的数据匹配问题
# 修改为：小麦千粒重分析 - 单一环境因子分析版
# 只分析目标环境因子：PAR_TEMP&GS67

# 1. 加载所需库
required_pkgs <- c("ggplot2", "dplyr", "lme4", "emmeans", "gridExtra", "grid", 
                   "colorspace", "tidyr", "purrr", "Cairo", "svglite")
lapply(required_pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
})

# 2. 设置路径参数及环境删除配置
data_dir <- "/mnt/7t_storage/zhangcl/TKW/"
output_dir <- paste0(data_dir, "loocv_single_factor\\")

# -------------------------- 环境删除配置 --------------------------
envs_to_remove <- c()
# ------------------------------------------------------------------

# 3. 创建输出目录
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("已创建输出目录：", output_dir, "\n")
}

# 4. 数据加载与环境匹配修复
cat("开始数据加载与环境匹配...\n")

# 4.1 加载表型数据
phenotype_data <- read.table(
    file = paste0(data_dir, "TKW_mean_table.txt"), 
    header = TRUE, sep = "\t", 
    stringsAsFactors = FALSE, na.strings = c("", "NA")
) %>% 
    rename(
        line_code = genotype,  # 将genotype列重命名为line_code以匹配后续代码
        PH = TKW  # 将TKW列重命名为PH以匹配后续代码（为了代码复用）
    ) %>%
    mutate(
        line_code = as.character(line_code), 
        env_code = as.character(env_code),
        PH = as.numeric(PH)
    ) %>%
    filter(!is.na(PH))

# 4.2 加载环境因子数据
ecs_data <- read.csv(
    file = paste0(data_dir, "ECs_results.csv"), 
    header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
)

# 重命名第一列为env_code
first_col_name <- colnames(ecs_data)[1]
ecs_data <- ecs_data %>%
    rename(env_code = 1) %>%
    mutate(env_code = as.character(env_code))

# 4.3 环境匹配与过滤
cat("执行环境匹配与过滤...\n")

# 获取两个数据集的环境列表
pheno_envs <- unique(phenotype_data$env_code)
ecs_envs <- unique(ecs_data$env_code)

cat("表型数据环境数:", length(pheno_envs), "\n")
cat("环境因子数据环境数:", length(ecs_envs), "\n")

# 找出共同的环境
common_envs <- intersect(pheno_envs, ecs_envs)
cat("共同环境数:", length(common_envs), "\n")

# 找出需要删除的环境（在共同环境中的）
envs_to_remove_final <- envs_to_remove[envs_to_remove %in% common_envs]
cat("实际要删除的环境:", paste(envs_to_remove_final, collapse = ", "), "\n")

# 最终保留的环境
final_envs <- setdiff(common_envs, envs_to_remove_final)
cat("最终保留的环境数:", length(final_envs), "\n")

# 4.4 过滤两个数据集，只保留共同且不需要删除的环境
phenotype_data_filtered <- phenotype_data %>%
    filter(env_code %in% final_envs)

ecs_data_filtered <- ecs_data %>%
    filter(env_code %in% final_envs)

cat("过滤后表型数据环境数:", length(unique(phenotype_data_filtered$env_code)), "\n")
cat("过滤后环境因子数据环境数:", length(unique(ecs_data_filtered$env_code)), "\n")

# 5. 设置目标环境因子
# 只分析 PAR_TEMP&GS67
target_factor <- "PAR_TEMP&GS67"

# 检查目标因子是否存在
if (!target_factor %in% colnames(ecs_data_filtered)) {
    cat("\n错误：目标环境因子", target_factor, "不存在于数据中！\n")
    cat("可用的环境因子：\n")
    print(colnames(ecs_data_filtered)[-1])
    stop("请检查环境因子名称是否正确")
}

cat("\n目标环境因子：", target_factor, "\n")

# 6. 预处理（计算BLUE和有效品种）
cat("\n开始数据预处理...\n")

# 筛选有效品种
valid_lines <- phenotype_data_filtered %>%
    group_by(line_code) %>%
    summarise(n_obs = sum(!is.na(PH))) %>%
    filter(n_obs >= 4) %>%
    pull(line_code)

cat("有效品种筛选完成 - 数量：", length(valid_lines), "\n")

env_codes_all <- unique(phenotype_data_filtered$env_code)
n_env <- length(env_codes_all)
cat("总环境数量：", n_env, "\n")

# 环境-颜色映射
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
color_mapping <- setNames(gg_color_hue(length(env_codes_all)), env_codes_all)

# 计算环境均值（只需一次）
all_env_means <- phenotype_data_filtered %>%
    group_by(env_code) %>%
    summarise(meanY = mean(PH, na.rm = TRUE)) %>%
    ungroup()

# 7. 计算BLUE值的函数（用于LOOCV）
calculate_blue_for_env <- function(train_envs, phenotype_data, valid_lines) {
    # 使用训练环境计算BLUE
    blue_data <- phenotype_data %>%
        filter(line_code %in% valid_lines, 
               env_code %in% train_envs, 
               !is.na(PH)) %>%
        mutate(
            line_code = as.factor(line_code), 
            env_code = as.factor(env_code),
            PH = as.numeric(PH)
        )
    
    if (nrow(blue_data) == 0 || 
        length(unique(blue_data$env_code)) < 2 || 
        length(unique(blue_data$line_code)) < 2) {
        cat("警告：训练数据不足，无法计算BLUE\n")
        return(NULL)
    }
    
    tryCatch({
        blue_model <- lmer(PH ~ (1 | line_code) + (1 | env_code), data = blue_data)
        
        blue_values <- ranef(blue_model)$line_code %>%
            as.data.frame() %>%
            mutate(
                line_code = rownames(.),
                blue_estimate = `(Intercept)` + fixef(blue_model)[1]
            ) %>%
            select(line_code, blue_estimate)
        
        return(blue_values)
        
    }, error = function(e) {
        cat("lmer模型失败，尝试使用lm模型...\n")
        tryCatch({
            blue_model <- lm(PH ~ line_code + env_code, data = blue_data)
            
            coefs <- coef(blue_model)
            line_coefs <- coefs[grepl("line_code", names(coefs))]
            
            blue_values <- data.frame(
                line_code = gsub("line_code", "", names(line_coefs)),
                blue_estimate = coefs[1] + line_coefs,
                stringsAsFactors = FALSE
            )
            
            base_line <- setdiff(levels(blue_data$line_code), 
                                gsub("line_code", "", names(line_coefs)))[1]
            blue_values <- rbind(
                data.frame(line_code = base_line, blue_estimate = coefs[1]),
                blue_values
            )
            return(blue_values)
        }, error = function(e2) {
            cat("lm模型也失败，使用品种均值...\n")
            line_means <- blue_data %>%
                group_by(line_code) %>%
                summarise(blue_estimate = mean(PH, na.rm = TRUE))
            return(line_means)
        })
    })
}

# 8. 交叉验证函数（修复版）- 增加BLUE的LOOCV
run_loocv <- function(phenotype_data, env_data, env_factor_data, all_env_means, 
                      method = c("env_factor", "env_mean", "blue"), 
                      min_train = 2, valid_lines = NULL, target_factor = NULL) {  
    
    method <- match.arg(method)
    env_codes <- unique(phenotype_data$env_code)
    results <- data.frame(
        env_code = character(), line_code = character(),
        observed = numeric(), predicted = numeric(),
        factor_value = numeric(), stringsAsFactors = FALSE
    )
    
    if (is.null(valid_lines)) valid_lines <- unique(phenotype_data$line_code)
    
    cat("开始", method, "LOOCV，测试环境数量：", length(env_codes), "\n")
    
    for (test_env in env_codes) {
        cat("  处理测试环境:", test_env, "...")
        
        if (method != "blue" && !(test_env %in% env_factor_data$env_code)) {
            cat("警告：环境", test_env, "无", target_factor, "数据，跳过\n")
            next
        }
        
        if (method != "blue") {
            current_factor_value <- env_factor_data[env_factor_data$env_code == test_env, "factor_value"]
        } else {
            current_factor_value <- NA  # BLUE方法不需要环境因子值
        }
        
        train_envs <- setdiff(env_codes, test_env)
        test_data <- filter(phenotype_data, env_code == test_env, line_code %in% valid_lines, !is.na(PH))
        
        if (nrow(test_data) == 0) {
            cat("无有效测试数据，跳过\n")
            next
        }
        
        if (method == "env_factor") {
            # 使用基础R方法创建环境参数数据
            env_params <- data.frame(
                env_code = env_data$env_code,
                kPara = as.numeric(env_data[[target_factor]])
            )
            
            test_data <- merge(test_data, env_params, by = "env_code")
            
            for (i in 1:nrow(test_data)) {
                line <- test_data$line_code[i]
                train_data <- filter(phenotype_data, line_code == line, env_code %in% train_envs, !is.na(PH)) %>%
                    merge(env_params, by = "env_code")
                
                if (nrow(train_data) >= min_train && nrow(train_data) > 1) {
                    tryCatch({
                        model <- lm(PH ~ kPara, data = train_data)
                        pred <- predict(model, data.frame(kPara = test_data$kPara[i]))
                        results <- rbind(results, data.frame(
                            env_code = test_env, line_code = line,
                            observed = test_data$PH[i], predicted = pred,
                            factor_value = current_factor_value
                        ))
                    }, error = function(e) {
                        pred <- mean(train_data$PH, na.rm = TRUE)
                        results <- rbind(results, data.frame(
                            env_code = test_env, line_code = line,
                            observed = test_data$PH[i], predicted = pred,
                            factor_value = current_factor_value
                        ))
                    })
                }
            }
        }
        
        if (method == "env_mean") {
            env_means <- filter(all_env_means, env_code %in% train_envs)
            test_mean <- all_env_means[all_env_means$env_code == test_env, "meanY"]
            
            if (length(test_mean) == 0 || is.na(test_mean)) {
                cat("警告：测试环境", test_env, "无均值数据，跳过\n")
                next
            }
            
            for (i in 1:nrow(test_data)) {
                line <- test_data$line_code[i]
                train_data <- filter(phenotype_data, line_code == line, env_code %in% train_envs, !is.na(PH)) %>%
                    inner_join(env_means, by = "env_code")
                
                if (nrow(train_data) >= min_train && nrow(train_data) > 1) {
                    tryCatch({
                        model <- lm(PH ~ meanY, data = train_data)
                        pred <- predict(model, data.frame(meanY = test_mean))
                        results <- rbind(results, data.frame(
                            env_code = test_env, line_code = line,
                            observed = test_data$PH[i], predicted = pred,
                            factor_value = current_factor_value
                        ))
                    }, error = function(e) {
                        pred <- mean(train_data$PH, na.rm = TRUE)
                        results <- rbind(results, data.frame(
                            env_code = test_env, line_code = line,
                            observed = test_data$PH[i], predicted = pred,
                            factor_value = current_factor_value
                        ))
                    })
                }
            }
        }
        
        if (method == "blue") {
            # BLUE的LOOCV：使用训练环境计算BLUE，预测测试环境
            blue_values <- calculate_blue_for_env(train_envs, phenotype_data, valid_lines)
            
            if (!is.null(blue_values)) {
                for (i in 1:nrow(test_data)) {
                    line <- test_data$line_code[i]
                    if (line %in% blue_values$line_code) {
                        pred <- blue_values$blue_estimate[blue_values$line_code == line]
                    } else {
                        # 如果训练环境中没有该品种，使用训练环境中所有品种的平均表现
                        train_data_all <- filter(phenotype_data, env_code %in% train_envs, !is.na(PH))
                        pred <- mean(train_data_all$PH, na.rm = TRUE)
                    }
                    
                    results <- rbind(results, data.frame(
                        env_code = test_env, line_code = line,
                        observed = test_data$PH[i], predicted = pred,
                        factor_value = current_factor_value
                    ))
                }
            }
        }
        
        cat("完成 - 获得", nrow(results), "个预测\n")
    }
    
    cat(method, "LOOCV完成 - 总预测样本：", nrow(results), "\n")
    return(results %>% filter(!is.na(predicted) & !is.na(observed)))
}

# 9. 生成标签数据函数
generate_uniform_label_data <- function(df, color_map, ref_obs_range, ref_pred_range) {
    if (nrow(df) == 0) {
        return(data.frame())
    }
    
    obs_extended_min <- ref_obs_range[1]
    obs_extended_max <- ref_obs_range[2]
    pred_extended_min <- ref_pred_range[1]
    pred_extended_max <- ref_pred_range[2]
    
    n_labels <- n_distinct(df$env_code)
    
    # 动态计算间距，适应更多环境
    spacing <- (pred_extended_max - pred_extended_min) / max(n_labels * 1.1, 1)
    
    # 计算每个环境的相关系数和RMSE
    env_stats <- df %>%
        group_by(env_code) %>%
        summarise(
            cor_val = round(cor(observed, predicted, use = "complete.obs"), 2),
            rmse_val = round(sqrt(mean((observed - predicted)^2, na.rm = TRUE)), 2)
        )
    
    env_stats %>%
        mutate(
            label = paste0("\u25A0 ", env_code, " (r=", cor_val, ", RMSE=", rmse_val, ")"),
            point_color = color_map[env_code],
            x_pos = obs_extended_min + (obs_extended_max - obs_extended_min) * 0.005,
            y_pos = pred_extended_max - (row_number() - 1) * spacing - (pred_extended_max - pred_extended_min) * 0.005
        )
}

# 10. 保存多种格式图像的函数
save_multiple_formats <- function(plot, file_path, width = 20, height = 20, dpi = 600) {
    # 确保目录存在
    dir_path <- dirname(file_path)
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    
    # 提取基础文件名
    base_name <- tools::file_path_sans_ext(file_path)
    
    cat("保存图片到:", base_name, "\n")
    
    # 1. 保存为TIFF格式 (高分辨率，适合出版)
    tiff_file <- paste0(base_name, ".tiff")
    tryCatch({
        ggsave(tiff_file, plot, width = width, height = height, dpi = dpi, 
               device = "tiff", compression = "lzw", bg = "white")
        cat("  TIFF格式已保存:", tiff_file, "\n")
    }, error = function(e) {
        cat("  TIFF保存失败:", e$message, "\n")
    })
    
    # 2. 保存为PDF格式 (矢量图形，适合出版)
    pdf_file <- paste0(base_name, ".pdf")
    tryCatch({
        ggsave(pdf_file, plot, width = width, height = height, device = "pdf", 
               bg = "white")
        cat("  PDF格式已保存:", pdf_file, "\n")
    }, error = function(e) {
        cat("  PDF保存失败:", e$message, "\n")
    })
    
    # 3. 保存为JPG格式 (600 PPI)
    jpg_file <- paste0(base_name, ".jpg")
    tryCatch({
        ggsave(jpg_file, plot, width = width, height = height, dpi = dpi,
               device = "jpeg", quality = 100, bg = "white")
        cat("  JPG格式已保存:", jpg_file, " (600 PPI)\n")
    }, error = function(e) {
        cat("  JPG保存失败:", e$message, "\n")
    })
    
    # 4. 保存为PNG格式 (高分辨率)
    png_file <- paste0(base_name, ".png")
    tryCatch({
        ggsave(png_file, plot, width = width, height = height, dpi = dpi,
               device = "png", bg = "white")
        cat("  PNG格式已保存:", png_file, "\n")
    }, error = function(e) {
        cat("  PNG保存失败:", e$message, "\n")
    })
    
    # 5. 保存为SVG格式 (矢量图形)
    svg_file <- paste0(base_name, ".svg")
    tryCatch({
        ggsave(svg_file, plot, width = width, height = height, device = svglite)
        cat("  SVG格式已保存:", svg_file, "\n")
    }, error = function(e) {
        cat("  SVG保存失败:", e$message, "\n")
    })
    
    # 6. 保存为EPS格式 (适合Latex)
    eps_file <- paste0(base_name, ".eps")
    tryCatch({
        ggsave(eps_file, plot, width = width, height = height, device = "eps", 
               bg = "white")
        cat("  EPS格式已保存:", eps_file, "\n")
    }, error = function(e) {
        cat("  EPS保存失败:", e$message, "\n")
    })
}

# 11. 主要分析函数（针对单一环境因子）
analyze_single_factor <- function(target_factor, phenotype_data_filtered, ecs_data_filtered, 
                                  all_env_means, valid_lines, color_mapping, 
                                  output_dir, env_codes_all, n_env) {
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
    cat("开始分析环境因子:", target_factor, "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
    
    # 为当前环境因子创建子目录
    factor_output_dir <- paste0(output_dir, gsub("[^[:alnum:]]", "_", target_factor), "\\")
    if (!dir.exists(factor_output_dir)) {
        dir.create(factor_output_dir, recursive = TRUE)
        cat("已创建环境因子输出目录：", factor_output_dir, "\n")
    }
    
    # 创建环境因子数据
    env_factor_data <- data.frame(
        env_code = ecs_data_filtered$env_code,
        factor_value = ecs_data_filtered[[target_factor]]
    )
    
    # 清理数据：移除NA值并转换为数值型
    env_factor_data <- env_factor_data %>%
        filter(!is.na(factor_value)) %>%
        mutate(factor_value = as.numeric(factor_value))
    
    cat("目标环境因子", target_factor, "数据加载完成 - 有效环境数：", nrow(env_factor_data), "\n")
    
    # 运行LOOCV - 三种方法
    cat("\n开始运行环境因子LOOCV...\n")
    env_factor_cv <- tryCatch({
        run_loocv(
            phenotype_data = phenotype_data_filtered, 
            env_data = ecs_data_filtered,
            env_factor_data = env_factor_data, 
            all_env_means = all_env_means,
            method = "env_factor", 
            min_train = 2, 
            valid_lines = valid_lines,
            target_factor = target_factor
        )
    }, error = function(e) {
        cat("环境因子LOOCV失败：", e$message, "\n")
        data.frame(env_code = character(), line_code = character(),
                   observed = numeric(), predicted = numeric(),
                   factor_value = numeric())
    })
    
    cat("环境因子LOOCV完成 - 有效样本：", nrow(env_factor_cv), "\n")
    
    cat("\n开始运行环境均值LOOCV...\n")
    env_mean_cv <- tryCatch({
        run_loocv(
            phenotype_data = phenotype_data_filtered,
            env_data = ecs_data_filtered,
            env_factor_data = env_factor_data, 
            all_env_means = all_env_means,
            method = "env_mean", 
            min_train = 2, 
            valid_lines = valid_lines,
            target_factor = target_factor
        )
    }, error = function(e) {
        cat("环境均值LOOCV失败：", e$message, "\n")
        data.frame(env_code = character(), line_code = character(),
                   observed = numeric(), predicted = numeric(),
                   factor_value = numeric())
    })
    
    cat("环境均值LOOCV完成 - 有效样本：", nrow(env_mean_cv), "\n")
    
    cat("\n开始运行BLUE LOOCV...\n")
    blue_cv <- tryCatch({
        run_loocv(
            phenotype_data = phenotype_data_filtered,
            env_data = ecs_data_filtered,
            env_factor_data = env_factor_data, 
            all_env_means = all_env_means,
            method = "blue", 
            min_train = 2, 
            valid_lines = valid_lines,
            target_factor = target_factor
        )
    }, error = function(e) {
        cat("BLUE LOOCV失败：", e$message, "\n")
        data.frame(env_code = character(), line_code = character(),
                   observed = numeric(), predicted = numeric(),
                   factor_value = numeric())
    })
    
    cat("BLUE LOOCV完成 - 有效样本：", nrow(blue_cv), "\n")
    
    # 检查是否有足够数据继续
    if (nrow(env_factor_cv) == 0 && nrow(env_mean_cv) == 0 && nrow(blue_cv) == 0) {
        cat("警告：所有方法都未能生成有效结果，跳过此环境因子\n")
        return(NULL)
    }
    
    # 保存结果
    write.csv(env_factor_cv, 
              file = paste0(factor_output_dir, "loocv_env_factor.csv"), row.names = FALSE)
    write.csv(env_mean_cv, 
              file = paste0(factor_output_dir, "loocv_env_mean.csv"), row.names = FALSE)
    write.csv(blue_cv, 
              file = paste0(factor_output_dir, "loocv_blue.csv"), row.names = FALSE)
    
    cv_summary <- bind_rows(
        env_factor_cv %>% mutate(method = "env_factor"),
        env_mean_cv %>% mutate(method = "env_mean"),
        blue_cv %>% mutate(method = "blue")
    )
    write.csv(cv_summary, 
              file = paste0(factor_output_dir, "all_summary.csv"), row.names = FALSE)
    cat("分析结果已保存至：", factor_output_dir, "\n")
    
    # 计算参考范围
    if (nrow(env_factor_cv) > 0) {
        ref_obs_range <- range(env_factor_cv$observed, na.rm = TRUE)
        ref_pred_range <- range(env_factor_cv$predicted, na.rm = TRUE)
    } else if (nrow(env_mean_cv) > 0) {
        ref_obs_range <- range(env_mean_cv$observed, na.rm = TRUE)
        ref_pred_range <- range(env_mean_cv$predicted, na.rm = TRUE)
    } else {
        ref_obs_range <- range(blue_cv$observed, na.rm = TRUE)
        ref_pred_range <- range(blue_cv$predicted, na.rm = TRUE)
    }
    
    ref_obs_extended_range <- c(
        ref_obs_range[1] - diff(ref_obs_range) * 0.25,
        ref_obs_range[2] + diff(ref_obs_range) * 0.25
    )
    ref_pred_extended_range <- c(
        ref_pred_range[1] - diff(ref_pred_range) * 0.25,
        ref_pred_range[2] + diff(ref_pred_range) * 0.25
    )
    
    # 生成标签数据
    label_data_env_factor <- generate_uniform_label_data(env_factor_cv, color_mapping, 
                                                         ref_obs_extended_range, ref_pred_extended_range)
    label_data_env_mean <- generate_uniform_label_data(env_mean_cv, color_mapping, 
                                                       ref_obs_extended_range, ref_pred_extended_range)
    label_data_blue <- generate_uniform_label_data(blue_cv, color_mapping, 
                                                   ref_obs_extended_range, ref_pred_extended_range)
    
    # 准备第四张图数据
    env_cor_data <- all_env_means %>%
        inner_join(env_factor_data, by = "env_code") %>%
        dplyr::select(env_code, meanY, factor_value) %>%
        rename(env_tkw_mean = meanY) %>%
        filter(!is.na(env_tkw_mean), !is.na(factor_value))
    
    env_cor <- if (nrow(env_cor_data) > 1) {
        round(cor(env_cor_data$factor_value, env_cor_data$env_tkw_mean, use = "complete.obs"), 2)
    } else {
        0
    }
    
    # 计算整体相关性和RMSE
    get_cor <- function(df) {
        if (nrow(df) > 1) {
            round(cor(df$observed, df$predicted, use = "complete.obs"), 2)
        } else {
            0
        }
    }
    
    get_rmse <- function(df) {
        if (nrow(df) > 0) {
            round(sqrt(mean((df$observed - df$predicted)^2, na.rm = TRUE)), 2)
        } else {
            NA
        }
    }
    
    correlations <- list(
        env_factor = get_cor(env_factor_cv),
        env_mean = get_cor(env_mean_cv),
        blue = get_cor(blue_cv)
    )
    
    rmses <- list(
        env_factor = get_rmse(env_factor_cv),
        env_mean = get_rmse(env_mean_cv),
        blue = get_rmse(blue_cv)
    )
    
    # 结果可视化 - 统一主题（添加abcd标签，去掉了主标题）
    uniform_theme <- theme_bw() + theme(
        plot.title = element_text(size = 16, hjust = 0, face = "bold", margin = margin(0,0,10,0)),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold", margin = margin(5,5,0,0)),
        panel.grid = element_blank(),
        plot.margin = unit(c(5, 5, 5, 5), "mm"),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1)
    )
    
    # 创建图形（添加abcd标签）
    plots <- list()
    
    # 子图1：环境因子预测 (a)
    if (nrow(env_factor_cv) > 0) {
        p1 <- ggplot(env_factor_cv, aes(x = observed, y = predicted, color = env_code)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
            geom_point(size = 2.5, alpha = 0.8) +
            scale_color_manual(values = color_mapping) +
            xlim(ref_obs_extended_range) + ylim(ref_pred_extended_range) +
            coord_fixed(ratio = 1) +
            labs(x = "Observed TKW (g)", y = "Predicted TKW (g)", 
                 title = "a) Env Factor (LOOCV)") +
            uniform_theme
        
        if (nrow(label_data_env_factor) > 0) {
            p1 <- p1 + geom_text(
                data = label_data_env_factor,
                aes(x = x_pos, y = y_pos, label = label),
                color = "black",
                hjust = 0, vjust = 1, size = 3.5, fontface = "bold"
            )
        }
        
        p1 <- p1 + annotate("text", 
                            x = ref_obs_extended_range[1] + diff(ref_obs_extended_range) * 0.05,
                            y = ref_pred_extended_range[2] - diff(ref_pred_extended_range) * 0.05,
                            label = sprintf("Overall: r = %.2f, RMSE = %.2f", 
                                          correlations$env_factor, rmses$env_factor), 
                            size = 5, fontface = "bold", color = "black", hjust = 0)
        
        plots$p1 <- p1
    }
    
    # 子图2：环境均值预测 (b)
    if (nrow(env_mean_cv) > 0) {
        p2 <- ggplot(env_mean_cv, aes(x = observed, y = predicted, color = env_code)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
            geom_point(size = 2.5, alpha = 0.8) +
            scale_color_manual(values = color_mapping) +
            xlim(ref_obs_extended_range) + ylim(ref_pred_extended_range) +
            coord_fixed(ratio = 1) +
            labs(x = "Observed TKW (g)", y = "Predicted TKW (g)", 
                 title = "b) Env Mean (LOOCV)") +
            uniform_theme
        
        if (nrow(label_data_env_mean) > 0) {
            p2 <- p2 + geom_text(
                data = label_data_env_mean,
                aes(x = x_pos, y = y_pos, label = label),
                color = "black",
                hjust = 0, vjust = 1, size = 3.5, fontface = "bold"
            )
        }
        
        p2 <- p2 + annotate("text", 
                            x = ref_obs_extended_range[1] + diff(ref_obs_extended_range) * 0.05,
                            y = ref_pred_extended_range[2] - diff(ref_pred_extended_range) * 0.05,
                            label = sprintf("Overall: r = %.2f, RMSE = %.2f", 
                                          correlations$env_mean, rmses$env_mean), 
                            size = 5, fontface = "bold", color = "black", hjust = 0)
        
        plots$p2 <- p2
    }
    
    # 子图3：BLUE LOOCV预测 (c)
    if (nrow(blue_cv) > 0) {
        p3 <- ggplot(blue_cv, aes(x = observed, y = predicted, color = env_code)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
            geom_point(size = 2.5, alpha = 0.8) +
            scale_color_manual(values = color_mapping) +
            xlim(ref_obs_extended_range) + ylim(ref_pred_extended_range) +
            coord_fixed(ratio = 1) +
            labs(x = "Observed TKW (g)", y = "BLUE Predicted (g)", 
                 title = "c) BLUE (LOOCV)") +
            uniform_theme
        
        if (nrow(label_data_blue) > 0) {
            p3 <- p3 + geom_text(
                data = label_data_blue,
                aes(x = x_pos, y = y_pos, label = label),
                color = "black",
                hjust = 0, vjust = 1, size = 3.5, fontface = "bold"
            )
        }
        
        p3 <- p3 + annotate("text", 
                            x = ref_obs_extended_range[1] + diff(ref_obs_extended_range) * 0.05,
                            y = ref_pred_extended_range[2] - diff(ref_pred_extended_range) * 0.05,
                            label = sprintf("Overall: r = %.2f, RMSE = %.2f", 
                                          correlations$blue, rmses$blue), 
                            size = 5, fontface = "bold", color = "black", hjust = 0)
        
        plots$p3 <- p3
    }
    
    # 子图4：环境因子与均值相关性 (d)
    if (nrow(env_cor_data) > 0) {
        x_range <- range(env_cor_data$factor_value, na.rm = TRUE)
        y_range <- range(env_cor_data$env_tkw_mean, na.rm = TRUE)
        x_extended_range <- c(x_range[1] - diff(x_range) * 0.2, x_range[2] + diff(x_range) * 0.2)
        y_extended_range <- c(y_range[1] - diff(y_range) * 0.2, y_range[2] + diff(y_range) * 0.2)
        
        p4 <- ggplot(env_cor_data, aes(x = factor_value, y = env_tkw_mean, color = env_code)) +
            geom_point(size = 4.5, alpha = 0.8) +
            geom_text(
                aes(label = env_code),
                size = 3.2, fontface = "bold", vjust = -1.2, color = "black"
            ) +
            geom_smooth(method = "lm", se = TRUE, fill = "lightblue", color = "darkblue", 
                       linewidth = 1.5, alpha = 0.2) +
            scale_color_manual(values = color_mapping) +
            xlim(x_extended_range) + ylim(y_extended_range) +
            coord_fixed(ratio = diff(x_extended_range)/diff(y_extended_range)) +
            annotate("text", 
                     x = x_extended_range[1] + diff(x_extended_range) * 0.05,
                     y = y_extended_range[2] - diff(y_extended_range) * 0.05,
                     label = sprintf("Correlation: r = %.2f", env_cor), 
                     size = 5, fontface = "bold", color = "black", hjust = 0) +
            labs(x = target_factor, y = "Environment TKW Mean (g)", 
                 title = "d) Environment Mean vs Factor") +
            uniform_theme
        
        plots$p4 <- p4
    }
    
    # 组合图形并保存多种格式（去掉主标题和脚注）
    if (length(plots) > 0) {
        n_plots <- length(plots)
        
        if (n_plots == 4) {
            # 使用arrangeGrob但去掉顶部和底部的文本
            combined_plot <- arrangeGrob(
                plots$p1, plots$p2, plots$p3, plots$p4,
                ncol = 2,
                padding = unit(0.5, "cm")  # 减少内边距
            )
        } else {
            combined_plot <- arrangeGrob(
                grobs = plots,
                ncol = 2,
                padding = unit(0.5, "cm")
            )
        }
        
        # 保存多种格式
        base_output_plot <- paste0(factor_output_dir, "TKW_comparison_all_methods_LOOCV")
        cat("\n开始保存多种格式图像...\n")
        
        # 保存组合图
        save_multiple_formats(combined_plot, base_output_plot, width = 20, height = 20, dpi = 600)
        
        # 同时保存每个子图为单独的高分辨率图像
        cat("\n开始保存单独的子图图像...\n")
        
        if (!is.null(plots$p1)) {
            save_multiple_formats(plots$p1, paste0(factor_output_dir, "a_env_factor_LOOCV"), 
                                 width = 10, height = 10, dpi = 600)
        }
        if (!is.null(plots$p2)) {
            save_multiple_formats(plots$p2, paste0(factor_output_dir, "b_env_mean_LOOCV"), 
                                 width = 10, height = 10, dpi = 600)
        }
        if (!is.null(plots$p3)) {
            save_multiple_formats(plots$p3, paste0(factor_output_dir, "c_blue_LOOCV"), 
                                 width = 10, height = 10, dpi = 600)
        }
        if (!is.null(plots$p4)) {
            save_multiple_formats(plots$p4, paste0(factor_output_dir, "d_env_correlation"), 
                                 width = 10, height = 10, dpi = 600)
        }
        
        cat("所有图像已保存！\n")
    } else {
        cat("警告：没有足够的有效数据生成图形\n")
    }
    
    # 输出性能指标
    performance <- data.frame(
        Method = character(),
        样本量 = integer(),
        相关系数r = numeric(),
        R平方 = numeric(),
        RMSE = numeric(),
        stringsAsFactors = FALSE
    )
    
    if (nrow(env_factor_cv) > 0) {
        performance <- rbind(performance, data.frame(
            Method = "环境因子预测 (LOOCV)",
            样本量 = nrow(env_factor_cv),
            相关系数r = correlations$env_factor,
            R平方 = correlations$env_factor^2,
            RMSE = rmses$env_factor
        ))
    }
    
    if (nrow(env_mean_cv) > 0) {
        performance <- rbind(performance, data.frame(
            Method = "环境均值预测 (LOOCV)",
            样本量 = nrow(env_mean_cv),
            相关系数r = correlations$env_mean,
            R平方 = correlations$env_mean^2,
            RMSE = rmses$env_mean
        ))
    }
    
    if (nrow(blue_cv) > 0) {
        performance <- rbind(performance, data.frame(
            Method = "BLUE预测 (LOOCV)",
            样本量 = nrow(blue_cv),
            相关系数r = correlations$blue,
            R平方 = correlations$blue^2,
            RMSE = rmses$blue
        ))
    }
    
    # 保存性能指标
    if (nrow(performance) > 0) {
        # 保存为CSV
        write.csv(performance, 
                  file = paste0(factor_output_dir, "performance_summary_LOOCV.csv"), row.names = FALSE)
        
        # 保存为TXT（更易读）
        perf_txt_file <- paste0(factor_output_dir, "performance_summary_LOOCV.txt")
        sink(perf_txt_file)
        cat("======================================================\n")
        cat("小麦千粒重(TKW)预测性能汇总\n")
        cat("======================================================\n")
        cat("分析日期:", format(Sys.Date(), "%Y年%m月%d日"), "\n")
        cat("目标环境因子:", target_factor, "\n")
        cat("环境数量:", n_env, "\n")
        cat("有效品种数量:", length(valid_lines), "\n")
        cat("环境因子与TKW均值相关性:", round(env_cor, 3), "\n")
        cat("\n--- LOOCV预测性能汇总 ---\n")
        cat("\n")
        print(performance, row.names = FALSE)
        cat("\n")
        cat("最优方法（基于相关系数）:", performance$Method[which.max(performance$相关系数r)], "\n")
        cat("最优方法（基于RMSE）:", performance$Method[which.min(performance$RMSE)], "\n")
        sink()
        
        best_method <- performance$Method[which.max(performance$相关系数r)]
        best_method_rmse <- performance$Method[which.min(performance$RMSE)]
        
        cat("\n=== LOOCV性能汇总 ===\n")
        print(performance, row.names = FALSE)
        cat("\n最优方法（基于相关系数）：", best_method, "\n")
        cat("最优方法（基于RMSE）：", best_method_rmse, "\n")
        cat("环境因子与TKW均值相关性：", env_cor, "\n")
        
        return(list(
            factor = target_factor,
            performance = performance,
            best_by_cor = best_method,
            best_by_rmse = best_method_rmse,
            env_factor_r = ifelse(nrow(env_factor_cv) > 0, correlations$env_factor, NA),
            env_mean_r = ifelse(nrow(env_mean_cv) > 0, correlations$env_mean, NA),
            blue_r = ifelse(nrow(blue_cv) > 0, correlations$blue, NA),
            env_cor = env_cor,
            n_env = n_env,
            n_lines = length(valid_lines)
        ))
    } else {
        cat("没有可用的性能数据\n")
        return(NULL)
    }
}

# 12. 运行单一环境因子分析
cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("开始分析目标环境因子：", target_factor, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

tryCatch({
    result <- analyze_single_factor(
        target_factor = target_factor,
        phenotype_data_filtered = phenotype_data_filtered,
        ecs_data_filtered = ecs_data_filtered,
        all_env_means = all_env_means,
        valid_lines = valid_lines,
        color_mapping = color_mapping,
        output_dir = output_dir,
        env_codes_all = env_codes_all,
        n_env = n_env
    )
    
    if (!is.null(result)) {
        cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
        cat("分析完成！结果已保存至：", output_dir, "\n")
        cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")
        
        # 生成详细报告
        report_file <- paste0(output_dir, "analysis_report_", 
                            format(Sys.Date(), "%Y%m%d"), ".txt")
        sink(report_file)
        cat("======================================================\n")
        cat("小麦千粒重(TKW)环境因子分析详细报告\n")
        cat("======================================================\n")
        cat("分析日期：", format(Sys.Date(), "%Y年%m月%d日"), "\n")
        cat("分析时间：", format(Sys.time(), "%H:%M:%S"), "\n")
        cat("目标环境因子：", result$factor, "\n")
        cat("环境数量：", result$n_env, "\n")
        cat("有效品种数量：", result$n_lines, "\n")
        cat("环境因子与TKW均值相关性：", round(result$env_cor, 4), "\n")
        cat("\n--- LOOCV预测性能详细汇总 ---\n")
        cat("\n")
        
        # 格式化输出性能表
        perf_table <- result$performance
        colnames(perf_table) <- c("Method", "N", "r", "R²", "RMSE")
        print(perf_table, row.names = FALSE)
        
        cat("\n--- 关键结果摘要 ---\n")
        cat("1. 最优预测方法（基于相关系数）：", result$best_by_cor, "\n")
        cat("   相关系数r =", round(perf_table$r[which.max(perf_table$r)], 4), "\n")
        cat("   R² =", round(perf_table$`R²`[which.max(perf_table$r)], 4), "\n")
        cat("   RMSE =", round(perf_table$RMSE[which.max(perf_table$r)], 4), "\n")
        
        cat("\n2. 最优预测方法（基于RMSE）：", result$best_by_rmse, "\n")
        cat("   RMSE =", round(min(perf_table$RMSE, na.rm = TRUE), 4), "\n")
        cat("   相关系数r =", round(perf_table$r[which.min(perf_table$RMSE)], 4), "\n")
        
        cat("\n3. 环境因子解释能力：\n")
        if (abs(result$env_cor) > 0.7) {
            cat("   环境因子", result$factor, "与TKW均值有强相关关系\n")
        } else if (abs(result$env_cor) > 0.5) {
            cat("   环境因子", result$factor, "与TKW均值有中等相关关系\n")
        } else if (abs(result$env_cor) > 0.3) {
            cat("   环境因子", result$factor, "与TKW均值有弱相关关系\n")
        } else {
            cat("   环境因子", result$factor, "与TKW均值相关性较弱\n")
        }
        
        cat("\n4. 图像文件保存情况：\n")
        cat("   组合图：TKW_comparison_all_methods_LOOCV.[tiff/pdf/jpg/png/svg/eps]\n")
        cat("   单独子图：a_env_factor_LOOCV, b_env_mean_LOOCV, c_blue_LOOCV, d_env_correlation\n")
        
        cat("\n--- 技术细节 ---\n")
        cat("分析方法：留一环境交叉验证（LOOCV）\n")
        cat("预测方法：环境因子回归、环境均值回归、BLUE预测\n")
        cat("有效品种筛选：至少4个环境有观测数据\n")
        cat("BLUE计算：线性混合模型（lmer）或线性模型（lm）\n")
        cat("图像格式：TIFF、PDF、JPG（600 PPI）、PNG、SVG、EPS\n")
        
        cat("\n======================================================\n")
        cat("报告结束\n")
        cat("======================================================\n")
        sink()
        
        cat("详细分析报告已保存至：", report_file, "\n")
        
        # 显示简要结果
        cat("\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
        cat("分析完成！关键结果：\n")
        cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
        cat("环境因子：", result$factor, "\n")
        cat("最佳预测方法：", result$best_by_cor, "\n")
        cat("最佳相关系数：", round(max(result$performance$相关系数r, na.rm = TRUE), 4), "\n")
        cat("最佳R²：", round(max(result$performance$R平方, na.rm = TRUE), 4), "\n")
        cat("最低RMSE：", round(min(result$performance$RMSE, na.rm = TRUE), 4), "\n")
        cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
        
    } else {
        cat("分析失败：没有生成有效结果\n")
    }
}, error = function(e) {
    cat("\n分析环境因子", target_factor, "时出错:", e$message, "\n")
    cat("错误追踪：\n")
    print(traceback())
})

# 13. 最终总结
cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("所有分析完成！\n")
cat("结果保存在：", output_dir, "\n")
cat("生成的图像格式：TIFF、PDF、JPG（600 PPI）、PNG、SVG、EPS\n")
cat("子图标签：a)、b)、c)、d)\n")
cat("已去除主标题和脚注\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# 14. 清理工作空间（可选）
cat("\n清理工作空间...\n")
rm(list = ls(pattern = "^temp_"))
gc()
cat("清理完成。\n")