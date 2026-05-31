# 清除环境残留与控制台输出
rm(list = ls(all = TRUE))
cat("\014")
message("✅ 前置环境清理完成")

# Set working directory
setwd("C:/Users/Lenovo/Desktop/小麦千粒重文章/001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割")

# 1. 加载必需库
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)

# 2. 核心路径配置
growthstage_dir <- "C:/Users/Lenovo/Desktop/小麦千粒重文章/001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割/03FutureAPSIMgrowthstage"

# 3. 检查文件并读取数据
message("=== 检查并读取数据文件 ===")

# 3.1 检查修复后的数据文件
fixed_data_path <- file.path(growthstage_dir, "weather_phenology_fixed.csv")
validation_report_dir <- file.path(growthstage_dir, "validation_reports")

if (file.exists(fixed_data_path)) {
  message("读取修复后的数据文件...")
  phenology_data <- read_csv(fixed_data_path, show_col_types = FALSE)
} else {
  # 如果修复后的数据不存在，读取原始数据
  original_data_path <- file.path(growthstage_dir, "weather_phenology_final.csv")
  if (file.exists(original_data_path)) {
    message("读取原始数据文件...")
    phenology_data <- read_csv(original_data_path, show_col_types = FALSE)
  } else {
    stop("❌ 未找到任何数据文件！")
  }
}

# 标准化数据格式
phenology_data <- phenology_data %>%
  mutate(
    env = as.character(env),
    date = as.Date(date),
    Wheat.Phenology.CurrentStageName = as.character(Wheat.Phenology.CurrentStageName),
    Wheat.Phenology.Zadok.Stage = as.numeric(Wheat.Phenology.Zadok.Stage)
  )

message("数据维度：", paste(dim(phenology_data), collapse = " x "))
message("环境数量：", length(unique(phenology_data$env)))
message("日期范围：", min(phenology_data$date), " 到 ", max(phenology_data$date))

# 4. 定义理论物候阶段顺序
message("\n=== 定义理论物候阶段顺序 ===")

theory_stage_order <- tibble(
  Stage = c(
    "Sowing",              # 播种 (Zadok 0)
    "Germination",         # 发芽 (Zadok 5-7)
    "Emergence",           # 出苗 (Zadok 10-11)
    "DoubleRidge",         # 双棱期 (Zadok 11-14)
    "MaximumSpikeletPrimordia", # 最大小穗原基 (Zadok 15-17)
    "FlagLeafAppearance",  # 旗叶出现 (Zadok 39)
    "HeadEmergence",       # 抽穗 (Zadok 55-57)
    "Anthesis",            # 开花 (Zadok 65-66)
    "MaximumGrainLength",  # 最大籽粒长度 (Zadok 71-72)
    "EndGrainFill",        # 灌浆结束 (Zadok 87-88)
    "HarvestRipe"          # 成熟 (Zadok 90)
  ),
  Order = 1:11,
  TheoreticalZadokMin = c(0, 5, 10, 11, 15, 39, 55, 65, 71, 87, 90),
  TheoreticalZadokMax = c(1, 7, 11, 14, 17, 39, 57, 66, 72, 88, 90),
  TheoreticalZadokTypical = c(0.5, 6, 10.5, 12.5, 16, 39, 56, 65.5, 71.5, 87.5, 90)
)

message("理论阶段顺序：")
print(theory_stage_order)

# 5. 分析缺失阶段问题
message("\n=== 分析物候阶段完整性 ===")

# 5.1 计算每个环境的物候阶段数量
stage_counts <- phenology_data %>%
  filter(!is.na(Wheat.Phenology.CurrentStageName) & 
           Wheat.Phenology.CurrentStageName != "") %>%
  group_by(env) %>%
  summarise(
    stage_count = n_distinct(Wheat.Phenology.CurrentStageName),
    stages_present = list(unique(Wheat.Phenology.CurrentStageName)),
    .groups = "drop"
  )

# 5.2 识别缺失的阶段
identify_missing_stages <- function(stages_present, theory_stages) {
  missing <- setdiff(theory_stages, stages_present)
  return(list(missing = missing, count = length(missing)))
}

stage_counts <- stage_counts %>%
  rowwise() %>%
  mutate(
    missing_info = list(identify_missing_stages(stages_present, theory_stage_order$Stage)),
    missing_stages = list(missing_info$missing),
    missing_count = missing_info$count
  ) %>%
  select(-missing_info) %>%
  ungroup()

# 5.3 统计信息
completeness_stats <- stage_counts %>%
  summarise(
    total_envs = n(),
    avg_stages = mean(stage_count, na.rm = TRUE),
    median_stages = median(stage_count, na.rm = TRUE),
    min_stages = min(stage_count, na.rm = TRUE),
    max_stages = max(stage_count, na.rm = TRUE),
    envs_with_11_stages = sum(stage_count == 11),
    envs_missing_stages = sum(missing_count > 0),
    total_missing_stages = sum(missing_count)
  )

message("📊 物候阶段完整性统计：")
print(completeness_stats)

# 5.4 显示缺失阶段的分布
missing_stage_distribution <- stage_counts %>%
  filter(missing_count > 0) %>%
  unnest(missing_stages) %>%
  group_by(missing_stages) %>%
  summarise(
    env_count = n(),
    percentage = round(n() / nrow(stage_counts) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(env_count))

message("\n📋 缺失阶段分布：")
print(missing_stage_distribution)

# 6. 检查验证问题
message("\n=== 检查验证问题 ===")

if (dir.exists(validation_report_dir)) {
  issues_path <- file.path(validation_report_dir, "validation_issues.csv")
  
  if (file.exists(issues_path)) {
    validation_issues <- read_csv(issues_path, show_col_types = FALSE)
    
    if (nrow(validation_issues) > 0) {
      message("找到验证问题文件，包含 ", nrow(validation_issues), " 个问题")
      
      # 按问题类型分组
      issue_summary <- validation_issues %>%
        group_by(IssueType) %>%
        summarise(
          count = n(),
          envs_affected = n_distinct(env),
          .groups = "drop"
        ) %>%
        arrange(desc(count))
      
      message("\n📋 验证问题统计：")
      print(issue_summary)
      
      # 获取有问题的环境
      problem_envs <- unique(validation_issues$env)
      message("\n⚠️ 有问题的环境数量：", length(problem_envs))
      message("前5个有问题的环境：")
      print(head(problem_envs, 5))
      
      # 显示具体问题示例
      if (nrow(validation_issues) > 0) {
        message("\n📄 问题示例：")
        print(head(validation_issues, 5))
      }
    } else {
      message("验证问题文件为空")
    }
  } else {
    message("未找到验证问题文件")
  }
} else {
  message("验证报告目录不存在")
}

# 7. 智能填补缺失阶段
message("\n=== 开始智能填补缺失阶段 ===")

# 7.1 智能填补函数
impute_missing_stages_intelligently <- function(env_data, env_id, theory_order) {
  # 提取现有阶段
  existing_stages <- env_data %>%
    filter(!is.na(Wheat.Phenology.CurrentStageName) & 
             Wheat.Phenology.CurrentStageName != "") %>%
    arrange(date) %>%
    distinct(Wheat.Phenology.CurrentStageName, .keep_all = TRUE)
  
  if (nrow(existing_stages) == 0) {
    message("环境 ", env_id, ": 没有现有阶段，无法填补")
    # 确保返回的数据有填补标记列
    if (!"imputation_flag" %in% colnames(env_data)) {
      env_data$imputation_flag <- FALSE
      env_data$imputation_method <- NA
      env_data$imputation_notes <- NA
    }
    return(env_data)
  }
  
  # 识别缺失阶段
  existing_stage_names <- existing_stages$Wheat.Phenology.CurrentStageName
  missing_stages <- setdiff(theory_order$Stage, existing_stage_names)
  
  if (length(missing_stages) == 0) {
    # 如果没有缺失阶段，确保返回的数据有填补标记列
    if (!"imputation_flag" %in% colnames(env_data)) {
      env_data$imputation_flag <- FALSE
      env_data$imputation_method <- NA
      env_data$imputation_notes <- NA
    }
    return(env_data)
  }
  
  message("环境 ", env_id, ": 缺失 ", length(missing_stages), " 个阶段")
  
  # 创建现有阶段的理论顺序映射
  existing_with_order <- existing_stages %>%
    left_join(theory_order %>% select(Stage, Order), 
              by = c("Wheat.Phenology.CurrentStageName" = "Stage")) %>%
    arrange(Order, date)
  
  # 为每个缺失阶段进行填补
  new_rows <- list()
  
  for (missing_stage in missing_stages) {
    # 获取缺失阶段的理论信息
    stage_info <- theory_order %>% filter(Stage == missing_stage)
    stage_order <- stage_info$Order
    theory_zadok <- stage_info$TheoreticalZadokTypical
    
    # 找到前后阶段
    prev_stages <- existing_with_order %>% filter(Order < stage_order)
    next_stages <- existing_with_order %>% filter(Order > stage_order)
    
    if (nrow(prev_stages) > 0 && nrow(next_stages) > 0) {
      # 前后都有阶段，使用线性插值
      prev_stage <- prev_stages[nrow(prev_stages), ]  # 最近的先前阶段
      next_stage <- next_stages[1, ]  # 最近的后续阶段
      
      # 计算插值权重
      order_range <- next_stage$Order - prev_stage$Order
      if (order_range > 0) {
        weight1 <- (stage_order - prev_stage$Order) / order_range
        weight2 <- (next_stage$Order - stage_order) / order_range
        
        # 插值日期
        days_between <- as.numeric(difftime(next_stage$date, prev_stage$date, units = "days"))
        imputed_date <- prev_stage$date + days(round(weight1 * days_between))
        
        # 插值Zadok值
        imputed_zadok <- prev_stage$Wheat.Phenology.Zadok.Stage * weight2 + 
                         next_stage$Wheat.Phenology.Zadok.Stage * weight1
        
        imputation_method <- "线性插值"
        
      } else {
        # 如果顺序相同（不应该发生），使用中间值
        imputed_date <- prev_stage$date + days(1)
        imputed_zadok <- theory_zadok
        imputation_method <- "理论值填补"
      }
      
    } else if (nrow(prev_stages) > 0) {
      # 只有先前阶段，使用外推
      prev_stage <- prev_stages[nrow(prev_stages), ]
      
      # 计算平均阶段间隔
      if (nrow(prev_stages) >= 2) {
        intervals <- numeric(nrow(prev_stages) - 1)
        for (i in 2:nrow(prev_stages)) {
          intervals[i-1] <- prev_stages$Order[i] - prev_stages$Order[i-1]
        }
        avg_interval <- mean(intervals, na.rm = TRUE)
      } else {
        avg_interval <- 1  # 默认间隔
      }
      
      # 外推日期和Zadok值
      order_diff <- stage_order - prev_stage$Order
      imputed_date <- prev_stage$date + days(round(7 * order_diff))  # 假设每个阶段间隔约7天
      imputed_zadok <- prev_stage$Wheat.Phenology.Zadok.Stage + 
                       (theory_zadok - prev_stage$TheoreticalZadokTypical)
      imputation_method <- "向前外推"
      
    } else if (nrow(next_stages) > 0) {
      # 只有后续阶段，使用反向外推
      next_stage <- next_stages[1, ]
      
      # 计算平均阶段间隔
      if (nrow(next_stages) >= 2) {
        intervals <- numeric(nrow(next_stages) - 1)
        for (i in 2:nrow(next_stages)) {
          intervals[i-1] <- next_stages$Order[i] - next_stages$Order[i-1]
        }
        avg_interval <- mean(intervals, na.rm = TRUE)
      } else {
        avg_interval <- 1
      }
      
      # 反向外推
      order_diff <- next_stage$Order - stage_order
      imputed_date <- next_stage$date - days(round(7 * order_diff))
      imputed_zadok <- next_stage$Wheat.Phenology.Zadok.Stage - 
                       (next_stage$TheoreticalZadokTypical - theory_zadok)
      imputation_method <- "向后外推"
      
    } else {
      # 没有参考阶段，使用理论值
      # 找到该环境的播种日期（如果存在）
      sowing_data <- env_data %>% 
        filter(Wheat.Phenology.CurrentStageName == "Sowing") %>%
        slice(1)
      
      if (nrow(sowing_data) > 0) {
        # 基于播种日期估算
        base_date <- sowing_data$date
        days_to_add <- (stage_order - 1) * 7  # 假设每个阶段间隔7天
        imputed_date <- base_date + days(days_to_add)
      } else {
        # 使用环境的最小日期
        base_date <- min(env_data$date, na.rm = TRUE)
        days_to_add <- (stage_order - 1) * 7
        imputed_date <- base_date + days(days_to_add)
      }
      
      imputed_zadok <- theory_zadok
      imputation_method <- "理论估算"
    }
    
    # 确保日期在合理范围内
    env_min_date <- min(env_data$date, na.rm = TRUE)
    env_max_date <- max(env_data$date, na.rm = TRUE)
    imputed_date <- max(env_min_date, min(imputed_date, env_max_date))
    
    # 确保Zadok值在合理范围内
    imputed_zadok <- max(0, min(100, imputed_zadok))
    
    # 创建新行
    new_row <- env_data[1, ]  # 使用第一行作为模板
    
    # 填充新行的物候信息
    new_row$date <- imputed_date
    new_row$Wheat.Phenology.CurrentStageName <- missing_stage
    new_row$Wheat.Phenology.Zadok.Stage <- round(imputed_zadok, 1)
    
    # 添加填补标记
    new_row$imputation_flag <- TRUE
    new_row$imputation_method <- imputation_method
    new_row$imputation_notes <- paste("填补缺失阶段:", missing_stage)
    
    new_rows[[length(new_rows) + 1]] <- new_row
    
    message(sprintf("  ✓ 填补阶段 %s: %s (Zadok=%.1f, 方法=%s)", 
                    missing_stage, as.character(imputed_date), 
                    imputed_zadok, imputation_method))
  }
  
  # 合并新行和原始数据
  if (length(new_rows) > 0) {
    new_rows_df <- bind_rows(new_rows)
    
    # 确保原始数据有填补标记列
    if (!"imputation_flag" %in% colnames(env_data)) {
      env_data$imputation_flag <- FALSE
      env_data$imputation_method <- NA
      env_data$imputation_notes <- NA
    }
    
    # 合并数据
    combined_data <- bind_rows(env_data, new_rows_df) %>%
      arrange(date)
    
    message(sprintf("  填补完成，新增 %d 行", length(new_rows)))
    
    return(combined_data)
  }
  
  # 确保返回的数据有填补标记列
  if (!"imputation_flag" %in% colnames(env_data)) {
    env_data$imputation_flag <- FALSE
    env_data$imputation_method <- NA
    env_data$imputation_notes <- NA
  }
  
  return(env_data)
}

# 7.2 执行智能填补
message("\n开始对所有环境进行智能填补...")

all_envs <- unique(phenology_data$env)
total_envs <- length(all_envs)

# 存储填补后的数据
filled_data_list <- list()
imputation_log <- tibble(
  env = character(),
  stage = character(),
  imputed_date = as.Date(character()),
  imputed_zadok = numeric(),
  theory_zadok = numeric(),
  method = character(),
  notes = character()
)

for (i in 1:total_envs) {
  env_id <- all_envs[i]
  
  if (i %% 50 == 0) {
    message(sprintf("正在处理第 %d/%d 个环境: %s", i, total_envs, env_id))
  }
  
  # 获取环境数据
  env_data <- phenology_data %>% filter(env == env_id)
  
  # 执行智能填补
  filled_data <- impute_missing_stages_intelligently(env_data, env_id, theory_stage_order)
  
  # 记录填补信息 - 修复版本：检查列是否存在
  if ("imputation_flag" %in% colnames(filled_data)) {
    # 筛选填补的行
    imputed_rows <- filled_data %>%
      filter(imputation_flag == TRUE) %>%
      select(env, Wheat.Phenology.CurrentStageName, date, 
             Wheat.Phenology.Zadok.Stage, imputation_method)
    
    if (nrow(imputed_rows) > 0) {
      for (j in 1:nrow(imputed_rows)) {
        row <- imputed_rows[j, ]
        stage_info <- theory_stage_order %>% filter(Stage == row$Wheat.Phenology.CurrentStageName)
        
        imputation_log <- imputation_log %>%
          add_row(
            env = env_id,
            stage = row$Wheat.Phenology.CurrentStageName,
            imputed_date = row$date,
            imputed_zadok = row$Wheat.Phenology.Zadok.Stage,
            theory_zadok = ifelse(nrow(stage_info) > 0, stage_info$TheoreticalZadokTypical, NA),
            method = row$imputation_method,
            notes = paste("智能填补:", row$Wheat.Phenology.CurrentStageName)
          )
      }
    }
  }
  
  filled_data_list[[env_id]] <- filled_data
}

# 合并所有填补后的数据
phenology_filled <- bind_rows(filled_data_list) %>%
  arrange(env, date)

message("\n✅ 智能填补完成！")
message("填补后的数据维度：", paste(dim(phenology_filled), collapse = " x "))

# 8. 验证填补结果
message("\n=== 验证填补结果 ===")

# 8.1 重新计算阶段数量
stage_counts_filled <- phenology_filled %>%
  filter(!is.na(Wheat.Phenology.CurrentStageName) & 
           Wheat.Phenology.CurrentStageName != "") %>%
  group_by(env) %>%
  summarise(
    stage_count = n_distinct(Wheat.Phenology.CurrentStageName),
    .groups = "drop"
  )

completeness_stats_filled <- stage_counts_filled %>%
  summarise(
    total_envs = n(),
    avg_stages = mean(stage_count, na.rm = TRUE),
    median_stages = median(stage_count, na.rm = TRUE),
    min_stages = min(stage_count, na.rm = TRUE),
    max_stages = max(stage_count, na.rm = TRUE),
    envs_with_11_stages = sum(stage_count == 11),
    envs_with_10_stages = sum(stage_count == 10),
    envs_with_9_stages = sum(stage_count <= 9)
  )

message("📊 填补后物候阶段完整性统计：")
print(completeness_stats_filled)

# 8.2 检查填补的质量
if (nrow(imputation_log) > 0) {
  message("\n📊 填补统计：")
  
  imputation_summary <- imputation_log %>%
    summarise(
      total_imputations = n(),
      unique_envs = n_distinct(env),
      unique_stages = n_distinct(stage),
      avg_zadok_diff = mean(abs(imputed_zadok - theory_zadok), na.rm = TRUE),
      .groups = "drop"
    )
  
  print(imputation_summary)
  
  # 按填补方法统计
  method_stats <- imputation_log %>%
    group_by(method) %>%
    summarise(
      count = n(),
      percentage = round(n() / nrow(imputation_log) * 100, 1),
      avg_zadok_diff = mean(abs(imputed_zadok - theory_zadok), na.rm = TRUE),
      .groups = "drop"
    )
  
  message("\n📋 填补方法统计：")
  print(method_stats)
}

# 9. 验证物候阶段顺序
message("\n=== 验证物候阶段顺序 ===")

# 9.1 顺序验证函数
validate_phenology_sequence <- function(env_data, theory_order) {
  env_id <- unique(env_data$env)
  
  # 提取物候阶段
  stage_data <- env_data %>%
    filter(!is.na(Wheat.Phenology.CurrentStageName) & 
             Wheat.Phenology.CurrentStageName != "") %>%
    arrange(date) %>%
    distinct(Wheat.Phenology.CurrentStageName, .keep_all = TRUE)
  
  if (nrow(stage_data) == 0) {
    return(list(
      env = env_id,
      valid = FALSE,
      issues = "没有物候阶段数据",
      stage_count = 0
    ))
  }
  
  # 检查阶段顺序
  issues <- character()
  stage_names <- stage_data$Wheat.Phenology.CurrentStageName
  stage_dates <- stage_data$date
  
  # 检查日期顺序
  for (i in 2:length(stage_dates)) {
    if (stage_dates[i] <= stage_dates[i-1]) {
      issues <- c(issues, 
                 sprintf("阶段顺序错误: %s (%s) 不晚于 %s (%s)",
                        stage_names[i], stage_dates[i],
                        stage_names[i-1], stage_dates[i-1]))
    }
  }
  
  # 检查Zadok值顺序
  zadok_values <- stage_data$Wheat.Phenology.Zadok.Stage
  for (i in 2:length(zadok_values)) {
    if (!is.na(zadok_values[i]) && !is.na(zadok_values[i-1])) {
      if (zadok_values[i] < zadok_values[i-1]) {
        issues <- c(issues,
                   sprintf("Zadok值递减: %s (%.1f) -> %s (%.1f)",
                          stage_names[i-1], zadok_values[i-1],
                          stage_names[i], zadok_values[i]))
      }
    }
  }
  
  valid <- (length(issues) == 0)
  
  return(list(
    env = env_id,
    valid = valid,
    issues = if (length(issues) > 0) paste(issues, collapse = "; ") else "无",
    stage_count = nrow(stage_data)
  ))
}

# 9.2 执行顺序验证
validation_results <- list()

for (env_id in unique(phenology_filled$env)) {
  env_data <- phenology_filled %>% filter(env == env_id)
  validation <- validate_phenology_sequence(env_data, theory_stage_order)
  validation_results[[env_id]] <- validation
}

# 9.3 汇总验证结果
validation_summary <- bind_rows(lapply(validation_results, as.data.frame))

validation_stats <- validation_summary %>%
  summarise(
    total_envs = n(),
    valid_envs = sum(valid),
    invalid_envs = sum(!valid),
    valid_percentage = round(valid_envs / total_envs * 100, 1),
    avg_stage_count = mean(stage_count, na.rm = TRUE),
    .groups = "drop"
  )

message("📊 物候阶段顺序验证统计：")
print(validation_stats)

# 9.4 显示有问题的环境
problem_envs <- validation_summary %>% filter(!valid) %>% pull(env)

if (length(problem_envs) > 0) {
  message("\n⚠️ 有问题的环境数量：", length(problem_envs))
  message("前5个有问题的环境：")
  
  for (i in 1:min(5, length(problem_envs))) {
    env_id <- problem_envs[i]
    result <- validation_results[[env_id]]
    message(sprintf("  %s: %s", env_id, result$issues))
  }
} else {
  message("✅ 所有环境物候阶段顺序正确！")
}

# 10. 保存最终结果
message("\n=== 保存最终结果 ===")

# 10.1 创建最终结果目录
final_results_dir <- file.path(growthstage_dir, "final_results")
if (!dir.exists(final_results_dir)) dir.create(final_results_dir, recursive = TRUE)

# 10.2 保存填补后的数据
final_data_path <- file.path(final_results_dir, "phenology_complete_filled.csv")
write_csv(phenology_filled, final_data_path)
message("✅ 已保存填补后的完整数据：", final_data_path)

# 10.3 保存填补日志
if (nrow(imputation_log) > 0) {
  imputation_log_path <- file.path(final_results_dir, "imputation_log_detailed.csv")
  write_csv(imputation_log, imputation_log_path)
  message("✅ 已保存详细填补日志：", imputation_log_path)
}

# 10.4 保存验证结果
validation_summary_path <- file.path(final_results_dir, "final_validation_summary.csv")
write_csv(validation_summary, validation_summary_path)
message("✅ 已保存验证结果：", validation_summary_path)

# 10.5 保存物候阶段摘要
# 辅助函数：获取第一个非NA值
first_non_na <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA)
  }
  return(x[!is.na(x)][1])
}

phenology_summary <- phenology_filled %>%
  filter(!is.na(Wheat.Phenology.CurrentStageName) & 
           Wheat.Phenology.CurrentStageName != "") %>%
  group_by(env) %>%
  summarise(
    site = first_non_na(site),
    scenario = first_non_na(scenario),
    stage_count = n_distinct(Wheat.Phenology.CurrentStageName),
    date_range = paste(min(date, na.rm = TRUE), "到", max(date, na.rm = TRUE)),
    duration_days = as.numeric(difftime(max(date, na.rm = TRUE), min(date, na.rm = TRUE), units = "days")) + 1,
    first_stage = first(Wheat.Phenology.CurrentStageName[order(date)]),
    last_stage = last(Wheat.Phenology.CurrentStageName[order(date)]),
    has_sowing = "Sowing" %in% Wheat.Phenology.CurrentStageName,
    has_harvestripe = "HarvestRipe" %in% Wheat.Phenology.CurrentStageName,
    .groups = "drop"
  )

phenology_summary_path <- file.path(final_results_dir, "phenology_summary_final.csv")
write_csv(phenology_summary, phenology_summary_path)
message("✅ 已保存物候阶段摘要：", phenology_summary_path)

# 11. 生成最终报告
message("\n=== 生成最终报告 ===")

final_report <- tibble(
  项目 = c(
    "处理前的环境数量",
    "处理前的平均物候阶段数",
    "填补的阶段总数",
    "填补后的环境数量", 
    "填补后的平均物候阶段数",
    "拥有完整11个阶段的环境数",
    "物候阶段顺序正确的环境数",
    "物候阶段顺序正确的比例"
  ),
  值 = c(
    total_envs,
    round(completeness_stats$avg_stages, 1),
    ifelse(nrow(imputation_log) > 0, nrow(imputation_log), 0),
    n_distinct(phenology_filled$env),
    round(completeness_stats_filled$avg_stages, 1),
    completeness_stats_filled$envs_with_11_stages,
    validation_stats$valid_envs,
    paste0(validation_stats$valid_percentage, "%")
  )
)

message("📋 最终处理报告：")
print(final_report)

# 保存最终报告
final_report_path <- file.path(final_results_dir, "final_processing_report.csv")
write_csv(final_report, final_report_path)
message("✅ 已保存最终处理报告：", final_report_path)

# 12. 可视化最终结果
message("\n=== 生成可视化图表 ===")

plots_final_dir <- file.path(final_results_dir, "final_plots")
if (!dir.exists(plots_final_dir)) dir.create(plots_final_dir, recursive = TRUE)

# 12.1 物候阶段数量分布
stage_count_plot <- stage_counts_filled %>%
  ggplot(aes(x = stage_count)) +
  geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 11, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 11.5, y = max(table(stage_counts_filled$stage_count))/2, 
           label = "理论值: 11", color = "red", hjust = 0, size = 4) +
  labs(
    title = "填补后各环境物候阶段数量分布",
    subtitle = paste("平均阶段数:", round(completeness_stats_filled$avg_stages, 1)),
    x = "物候阶段数量",
    y = "环境数量"
  ) +
  theme_minimal()

ggsave(file.path(plots_final_dir, "stage_count_distribution_final.png"), 
       stage_count_plot, width = 10, height = 6, dpi = 300)

# 12.2 验证状态分布
validation_status_plot <- validation_summary %>%
  mutate(Status = ifelse(valid, "有效", "无效")) %>%
  ggplot(aes(x = Status, fill = Status)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("有效" = "green", "无效" = "red")) +
  labs(
    title = "物候阶段顺序验证状态",
    subtitle = paste("有效环境比例:", validation_stats$valid_percentage, "%"),
    x = "验证状态",
    y = "环境数量"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plots_final_dir, "validation_status_final.png"), 
       validation_status_plot, width = 8, height = 6, dpi = 300)

# 12.3 填补方法分布（如果有填补）
if (nrow(imputation_log) > 0) {
  method_plot <- imputation_log %>%
    group_by(method) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ggplot(aes(x = reorder(method, -count), y = count, fill = method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(count, "\n(", round(percentage, 1), "%)")), 
              vjust = -0.3, size = 3) +
    labs(
      title = "填补方法分布",
      x = "填补方法",
      y = "填补次数"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(plots_final_dir, "imputation_methods_distribution.png"), 
         method_plot, width = 10, height = 6, dpi = 300)
}

message("✅ 已保存最终可视化图表到：", plots_final_dir)

message("\n🎉 物候阶段完整性填补与验证完成！")
message("📊 主要成果：")
message("  1. 填补前平均阶段数：", round(completeness_stats$avg_stages, 1))
message("  2. 填补后平均阶段数：", round(completeness_stats_filled$avg_stages, 1))
message("  3. 完整环境数（11个阶段）：", completeness_stats_filled$envs_with_11_stages)
message("  4. 顺序正确环境比例：", validation_stats$valid_percentage, "%")
message("📁 所有结果保存在：", final_results_dir)