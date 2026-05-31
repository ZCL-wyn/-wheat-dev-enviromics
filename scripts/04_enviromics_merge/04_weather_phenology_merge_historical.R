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
library(stringr)

# 2. 核心路径配置
future_growthstage_dir <- "C:/Users/Lenovo/Desktop/小麦千粒重文章/001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割/03FutureAPSIMgrowthstage"
historical_dir <- "C:/Users/Lenovo/Desktop/小麦千粒重文章/001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割/03APSIMgrowthstage"

# 历史数据文件路径
historical_file_path <- file.path(historical_dir, "final_weather_cut_by_phenology.csv")

# 未来数据文件路径
future_file_path <- file.path(future_growthstage_dir, "final_weather_cut_by_phenology.csv")

# 输出文件路径
merged_output_path <- file.path(future_growthstage_dir, "merged_historical_future_growthstage_data.csv")

# 3. 读取数据
message("=== 开始读取历史数据和未来数据 ===")

# 检查文件是否存在
if (!file.exists(historical_file_path)) {
  stop(paste("❌ 历史数据文件不存在：", historical_file_path))
}

if (!file.exists(future_file_path)) {
  stop(paste("❌ 未来数据文件不存在：", future_file_path))
}

# 读取历史数据
message(paste("读取历史数据：", historical_file_path))
historical_data <- read_csv(historical_file_path, show_col_types = FALSE)

# 读取未来数据
message(paste("读取未来数据：", future_file_path))
future_data <- read_csv(future_file_path, show_col_types = FALSE)

# 4. 查看数据结构
message("\n📊 数据结构：")
message("历史数据维度：", paste(dim(historical_data), collapse = " x "))
message("未来数据维度：", paste(dim(future_data), collapse = " x "))

message("\n历史数据列名：")
print(colnames(historical_data))

message("\n未来数据列名：")
print(colnames(future_data))

# 5. 标准化列名和数据结构
message("\n=== 标准化数据结构 ===")

# 找出两个数据集共有的列
common_cols <- intersect(colnames(historical_data), colnames(future_data))
message("共有的列数：", length(common_cols))
message("共有列：")
print(common_cols)

# 找出只在历史数据中存在的列
historical_only <- setdiff(colnames(historical_data), colnames(future_data))
message("\n只在历史数据中存在的列（", length(historical_only), "个）：")
print(historical_only)

# 找出只在未来数据中存在的列
future_only <- setdiff(colnames(future_data), colnames(historical_data))
message("\n只在未来数据中存在的列（", length(future_only), "个）：")
print(future_only)

# 6. 检查关键列是否存在（移除Wheat.Phenology.ThermalTime）
required_cols <- c("env", "Wheat.Phenology.Zadok.Stage", "Wheat.Phenology.CurrentStageName", 
                   "YYYYMMDD", "TMEAN", "TMAX", "TMIN")

missing_in_historical <- setdiff(required_cols, colnames(historical_data))
missing_in_future <- setdiff(required_cols, colnames(future_data))

if (length(missing_in_historical) > 0) {
  message("⚠️ 历史数据中缺失的列：", paste(missing_in_historical, collapse = ", "))
}

if (length(missing_in_future) > 0) {
  message("⚠️ 未来数据中缺失的列：", paste(missing_in_future, collapse = ", "))
}

# 7. 移除GDD_r重命名逻辑（因为无Wheat.Phenology.ThermalTime列）
message("\n=== 跳过GDD_r重命名（无Wheat.Phenology.ThermalTime列）===")
message("✅ 已跳过不必要的列重命名步骤")

# 8. 处理YYYYMMDD列
message("\n=== 处理日期列 ===")

# 检查历史数据中是否有YYYYMMDD列，如果没有则尝试其他日期列
if (!"YYYYMMDD" %in% colnames(historical_data)) {
  # 查找可能的日期列
  date_patterns <- c("Date", "DATE", "date", "Time", "TIME", "time", "day", "Day", "DAY")
  date_cols <- colnames(historical_data)[sapply(colnames(historical_data), function(x) any(grepl(paste(date_patterns, collapse = "|"), x)))]
  
  if (length(date_cols) > 0) {
    # 使用第一个找到的日期列
    date_col <- date_cols[1]
    historical_data <- historical_data %>%
      rename(YYYYMMDD = !!date_col)
    message(paste("✅ 历史数据中重命名", date_col, "为YYYYMMDD"))
  } else {
    message("⚠️ 历史数据中没有找到日期列")
  }
}

# 确保YYYYMMDD是日期格式
if ("YYYYMMDD" %in% colnames(historical_data)) {
  historical_data <- historical_data %>%
    mutate(YYYYMMDD = as.Date(YYYYMMDD))
  message("✅ 历史数据YYYYMMDD已转换为日期格式")
}

if ("YYYYMMDD" %in% colnames(future_data)) {
  future_data <- future_data %>%
    mutate(YYYYMMDD = as.Date(YYYYMMDD))
  message("✅ 未来数据YYYYMMDD已转换为日期格式")
}

# 9. 确保两个数据集有相同的列
message("\n=== 统一两个数据集的列 ===")

# 获取两个数据集的所有列
all_historical_cols <- colnames(historical_data)
all_future_cols <- colnames(future_data)
all_cols <- unique(c(all_historical_cols, all_future_cols))

message("总列数：", length(all_cols))

# 为历史数据添加缺失的列
missing_in_historical <- setdiff(all_cols, all_historical_cols)
if (length(missing_in_historical) > 0) {
  message("为历史数据添加缺失的列（", length(missing_in_historical), "个）...")
  for (col in missing_in_historical) {
    historical_data[[col]] <- NA
  }
  message("✅ 历史数据缺失列已添加")
}

# 为未来数据添加缺失的列
missing_in_future <- setdiff(all_cols, all_future_cols)
if (length(missing_in_future) > 0) {
  message("为未来数据添加缺失的列（", length(missing_in_future), "个）...")
  for (col in missing_in_future) {
    future_data[[col]] <- NA
  }
  message("✅ 未来数据缺失列已添加")
}

# 确保两个数据集的列顺序一致
historical_data <- historical_data[, all_cols]
future_data <- future_data[, all_cols]

message("✅ 列已统一")

# 10. 合并数据
message("\n=== 合并历史数据和未来数据 ===")

# 添加数据源标识
historical_data$Data_Source <- "Historical"
future_data$Data_Source <- "Future"

# 使用bind_rows合并
merged_data <- bind_rows(historical_data, future_data)

message("合并后数据维度：", paste(dim(merged_data), collapse = " x "))
message("历史数据行数：", nrow(historical_data))
message("未来数据行数：", nrow(future_data))
message("合并后总行数：", nrow(merged_data))

# 11. 检查合并后数据质量（移除GDD_r相关统计）
message("\n=== 检查合并后数据质量 ===")

# 检查数据源分布
source_dist <- merged_data %>%
  group_by(Data_Source) %>%
  summarise(
    Rows = n(),
    Unique_Envs = n_distinct(env, na.rm = TRUE),
    Date_Range = ifelse(
      all(is.na(YYYYMMDD)), 
      "No dates", 
      paste(min(YYYYMMDD, na.rm = TRUE), "to", max(YYYYMMDD, na.rm = TRUE))
    ),
    .groups = "drop"
  )

message("数据源分布：")
print(source_dist)

# 检查关键列的缺失值（移除GDD_r）
if ("env" %in% colnames(merged_data)) {
  missing_stats <- merged_data %>%
    summarise(
      env_NA = sum(is.na(env)),
      YYYYMMDD_NA = sum(is.na(YYYYMMDD)),
      Zadok_NA = sum(is.na(Wheat.Phenology.Zadok.Stage)),
      StageName_NA = sum(is.na(Wheat.Phenology.CurrentStageName)),
      TMEAN_NA = ifelse("TMEAN" %in% colnames(merged_data), sum(is.na(TMEAN)), NA),
      .groups = "drop"
    )
  
  message("\n关键列缺失值统计：")
  print(missing_stats)
}

# 12. 重新排列列顺序（移除GDD_r）
message("\n=== 重新排列列顺序 ===")

# 定义优先显示的列（移除GDD_r）
priority_cols <- c("Data_Source", "env", "YYYYMMDD", 
                   "Wheat.Phenology.Zadok.Stage", "Wheat.Phenology.CurrentStageName")

# 气象数据列（按您要求的顺序）
weather_cols <- c("TMEAN", "TMAX", "TMIN", "T2MDEW", "ASSD", "APAR", "DL", 
                  "PR", "QV2M", "RH", "WS2M", "VPD", "GDD_base0", "GDD_simple", 
                  "FRUE", "GDD_cardinal", "DTR")

# 只保留实际存在的列
priority_cols <- priority_cols[priority_cols %in% colnames(merged_data)]
weather_cols <- weather_cols[weather_cols %in% colnames(merged_data)]

# 其他列
other_cols <- setdiff(colnames(merged_data), c(priority_cols, weather_cols))

# 重新排列列顺序
merged_data <- merged_data %>%
  select(all_of(priority_cols), all_of(weather_cols), all_of(other_cols))

message("列顺序已重新排列")
message("合并后数据列顺序：")
print(colnames(merged_data))

# 13. 保存合并后的数据
message("\n=== 保存合并后的数据 ===")

write.csv(merged_data, merged_output_path, row.names = FALSE)
message(paste("✅ 已保存合并后的数据：", merged_output_path))

# 14. 生成数据摘要报告（移除Mean_GDD_r）
message("\n=== 生成数据摘要报告 ===")

# 按数据源和环境统计
env_summary <- merged_data %>%
  group_by(Data_Source, env) %>%
  summarise(
    Start_Date = min(YYYYMMDD, na.rm = TRUE),
    End_Date = max(YYYYMMDD, na.rm = TRUE),
    Days = as.numeric(difftime(End_Date, Start_Date, units = "days")) + 1,
    Min_Zadok = min(Wheat.Phenology.Zadok.Stage, na.rm = TRUE),
    Max_Zadok = max(Wheat.Phenology.Zadok.Stage, na.rm = TRUE),
    Mean_Temp = ifelse("TMEAN" %in% colnames(merged_data), mean(TMEAN, na.rm = TRUE), NA),
    Total_Precip = ifelse("PR" %in% colnames(merged_data), sum(PR, na.rm = TRUE), NA),
    .groups = "drop"
  )

# 保存环境摘要
env_summary_path <- file.path(future_growthstage_dir, "merged_data_environment_summary.csv")
write.csv(env_summary, env_summary_path, row.names = FALSE)
message(paste("✅ 已保存环境摘要：", env_summary_path))

# 15. 物候阶段统计（移除Avg_GDD_r）
message("\n=== 物候阶段统计 ===")

phenology_stats <- merged_data %>%
  filter(!is.na(Wheat.Phenology.CurrentStageName)) %>%
  group_by(Data_Source, Wheat.Phenology.CurrentStageName) %>%
  summarise(
    Count = n(),
    Avg_Zadok = mean(Wheat.Phenology.Zadok.Stage, na.rm = TRUE),
    Min_Date = min(YYYYMMDD, na.rm = TRUE),
    Max_Date = max(YYYYMMDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Data_Source, Avg_Zadok)

# 保存物候统计
phenology_stats_path <- file.path(future_growthstage_dir, "merged_data_phenology_statistics.csv")
write.csv(phenology_stats, phenology_stats_path, row.names = FALSE)
message(paste("✅ 已保存物候统计：", phenology_stats_path))

# 16. 生成最终报告
message("\n📊 ========== 合并数据最终报告 ==========")
message("📁 输出目录：", future_growthstage_dir)
message("📄 生成的文件：")
message("  1. merged_historical_future_growthstage_data.csv - 合并后的完整数据")
message("  2. merged_data_environment_summary.csv - 环境摘要")
message("  3. merged_data_phenology_statistics.csv - 物候统计")

message("\n📈 数据统计：")
message(paste("  - 历史数据行数：", nrow(historical_data)))
message(paste("  - 未来数据行数：", nrow(future_data)))
message(paste("  - 合并后总行数：", nrow(merged_data)))
message(paste("  - 历史环境数：", source_dist$Unique_Envs[source_dist$Data_Source == "Historical"]))
message(paste("  - 未来环境数：", source_dist$Unique_Envs[source_dist$Data_Source == "Future"]))
message(paste("  - 总环境数：", n_distinct(merged_data$env)))

message("\n🔧 注意事项：")
message("  - 因无Wheat.Phenology.ThermalTime列，已移除所有GDD_r相关逻辑")
message("  - 所有统计均基于现有数据列完成")

message("\n🎉 数据合并完成！")

# 17. 显示前几行数据示例
message("\n=== 合并后数据前10行示例 ===")
print(head(merged_data, 10))