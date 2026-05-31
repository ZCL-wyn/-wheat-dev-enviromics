########################## 1. 前置检查：自动安装并加载必要包 ##########################
if (!require("lubridate", quietly = TRUE)) {
  install.packages("lubridate", dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(lubridate)
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(dplyr)
}
if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr", dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(tidyr)
}
if (!require("zoo", quietly = TRUE)) {
  install.packages("zoo", dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(zoo)
}

# 清除环境残留与控制台输出
rm(list = ls(all = TRUE))
cat("\014")
message("✅ 前置环境清理完成")


########################## 2. 配置工作目录与数据读取 ##########################
work_dir <- "C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割"
if (!dir.exists(work_dir)) {
  stop(paste("❌ 工作目录不存在：", work_dir))
}
setwd(work_dir)
message(paste("✅ 工作目录：", work_dir))

# 读取站点信息
loc_data <- data.frame(
  site = c("NY", "PY", "YL", "ZMD"),
  latitude = c(32.9126, 35.7627, 34.2879, 33.014),
  longitude = c(112.4626, 115.0292, 108.0001, 114.0219),
  stringsAsFactors = FALSE
)
message("✅ 站点数据：")
print(loc_data)

# 读取原始气象数据
weather_file <- "all_sites_ssp126_ssp585_2001_2050.csv"
if (!file.exists(weather_file)) {
  stop(paste("❌ 气象数据缺失：", weather_file))
}

original_weather <- read.csv(weather_file, stringsAsFactors = FALSE)
message(paste("✅ 读取原始气象数据：", nrow(original_weather), "行"))
message("原始数据集列名：")
print(names(original_weather))

# 确保时间列格式正确
if ("time" %in% names(original_weather)) {
  original_weather$time <- as.Date(original_weather$time)
} else if ("date" %in% names(original_weather)) {
  original_weather$time <- as.Date(original_weather$date)
} else {
  stop("❌ 数据集中没有找到时间列（time或date）")
}

message(paste("📅 时间范围：", min(original_weather$time), "至", max(original_weather$time)))
message(paste("📊 站点：", paste(unique(original_weather$site), collapse = ", ")))
message(paste("📊 情景：", paste(unique(original_weather$scenario), collapse = ", ")))


########################## 3. 创建数据矫正和插补函数 ##########################
correct_and_interpolate_weather <- function(weather_df, site_name, scenario_name) {
  # 筛选特定站点和情景的数据
  site_scenario_data <- weather_df %>%
    filter(site == site_name, scenario == scenario_name) %>%
    arrange(time)
  
  if (nrow(site_scenario_data) == 0) {
    warning(paste("❌ 没有数据：站点", site_name, "情景", scenario_name))
    return(NULL)
  }
  
  message(paste("  处理站点", site_name, "情景", scenario_name, 
                "（", nrow(site_scenario_data), "行数据）"))
  
  # 1. 创建完整的日期序列（从第一天到最后一天）
  full_date_seq <- seq.Date(min(site_scenario_data$time), 
                           max(site_scenario_data$time), 
                           by = "day")
  
  # 2. 检查是否有缺失日期
  missing_dates <- setdiff(full_date_seq, site_scenario_data$time)
  if (length(missing_dates) > 0) {
    message(paste("  发现", length(missing_dates), "个缺失日期，将进行插补"))
    
    # 创建包含所有日期的完整数据框
    full_data <- data.frame(
      time = full_date_seq,
      site = site_name,
      scenario = scenario_name,
      stringsAsFactors = FALSE
    )
    
    # 合并数据
    full_data <- left_join(full_data, site_scenario_data, 
                          by = c("time", "site", "scenario"))
  } else {
    full_data <- site_scenario_data
  }
  
  # 3. 单位转换和矫正
  
  # 3.1 温度数据矫正
  # 检查是否有温度数据
  if ("tasmax" %in% names(full_data) && "tasmin" %in% names(full_data)) {
    message("  矫正温度数据...")
    
    # 检查单位：如果平均值>100，可能是开尔文(K)
    if (mean(full_data$tasmax, na.rm = TRUE) > 100) {
      full_data$tasmax_corrected <- full_data$tasmax - 273.15
      full_data$tasmin_corrected <- full_data$tasmin - 273.15
    } else {
      full_data$tasmax_corrected <- full_data$tasmax
      full_data$tasmin_corrected <- full_data$tasmin
    }
    
    # 确保最低温度不高于最高温度
    for (i in 1:nrow(full_data)) {
      if (!is.na(full_data$tasmax_corrected[i]) && !is.na(full_data$tasmin_corrected[i])) {
        if (full_data$tasmin_corrected[i] > full_data$tasmax_corrected[i]) {
          # 交换值
          temp <- full_data$tasmin_corrected[i]
          full_data$tasmin_corrected[i] <- full_data$tasmax_corrected[i]
          full_data$tasmax_corrected[i] <- temp
        }
      }
    }
    
    # 计算日平均温度
    full_data$tas_corrected <- (full_data$tasmax_corrected + full_data$tasmin_corrected) / 2
  } else if ("tas" %in% names(full_data)) {
    message("  使用平均温度计算最高最低温度...")
    
    # 检查单位
    if (mean(full_data$tas, na.rm = TRUE) > 100) {
      full_data$tas_corrected <- full_data$tas - 273.15
    } else {
      full_data$tas_corrected <- full_data$tas
    }
    
    # 使用日较差（10°C）计算最高最低温度
    full_data$tasmax_corrected <- full_data$tas_corrected + 5
    full_data$tasmin_corrected <- full_data$tas_corrected - 5
  } else {
    warning("  缺少温度数据，使用默认值")
    full_data$tas_corrected <- 15
    full_data$tasmax_corrected <- 20
    full_data$tasmin_corrected <- 10
  }
  
  # 3.2 降水数据矫正
  if ("pr" %in% names(full_data)) {
    message("  矫正降水数据...")
    
    # 检查单位：如果平均值很小，可能是kg/m²/s
    avg_pr <- mean(full_data$pr, na.rm = TRUE)
    if (avg_pr < 0.01) {
      # kg/m²/s 到 mm/day: * 86400
      full_data$pr_corrected <- full_data$pr * 86400
      message(paste("    降水单位转换：kg/m²/s → mm/day (平均", round(avg_pr*86400, 2), "mm/day)"))
    } else {
      full_data$pr_corrected <- full_data$pr
      message(paste("    降水单位：mm/day (平均", round(avg_pr, 2), "mm/day)"))
    }
    
    # 确保降水不为负
    full_data$pr_corrected <- pmax(0, full_data$pr_corrected)
  } else {
    warning("  缺少降水数据，使用默认值")
    full_data$pr_corrected <- 0
  }
  
  # 3.3 太阳辐射数据矫正
  if ("rsds" %in% names(full_data)) {
    message("  矫正太阳辐射数据...")
    
    # W/m² 到 MJ/m²/day: * 0.0864
    full_data$rsds_corrected <- full_data$rsds * 0.0864
    
    avg_rsds <- mean(full_data$rsds, na.rm = TRUE)
    message(paste("    太阳辐射单位转换：W/m² → MJ/m²/day (平均", 
                  round(avg_rsds*0.0864, 2), "MJ/m²/day)"))
    
    # 确保辐射不为负
    full_data$rsds_corrected <- pmax(0, full_data$rsds_corrected)
  } else {
    warning("  缺少太阳辐射数据，使用默认值")
    full_data$rsds_corrected <- 15
  }
  
  # 3.4 相对湿度数据矫正
  if ("hurs" %in% names(full_data)) {
    message("  矫正相对湿度数据...")
    full_data$hurs_corrected <- full_data$hurs
    
    # 确保相对湿度在合理范围内
    full_data$hurs_corrected <- pmin(100, pmax(0, full_data$hurs_corrected))
  } else {
    message("  估算相对湿度...")
    # 基于温度估算相对湿度
    tavg <- full_data$tas_corrected
    full_data$hurs_corrected <- 80 - (tavg - 15) * 2  # 15°C时80%，每升高1度减少2%
    full_data$hurs_corrected <- pmax(30, pmin(100, full_data$hurs_corrected))
  }
  
  # 3.5 风速数据矫正
  if ("sfcwind" %in% names(full_data)) {
    message("  矫正风速数据...")
    full_data$sfcwind_corrected <- full_data$sfcwind
    
    # 确保风速不为负
    full_data$sfcwind_corrected <- pmax(0, full_data$sfcwind_corrected)
  } else {
    warning("  缺少风速数据，使用默认值")
    full_data$sfcwind_corrected <- 2.0
  }
  
  # 3.6 气压数据矫正
  if ("ps" %in% names(full_data)) {
    message("  矫正气压数据...")
    
    # 帕(Pa) 到 千帕(kPa): / 1000
    full_data$ps_corrected <- full_data$ps / 1000
    
    avg_ps <- mean(full_data$ps, na.rm = TRUE)
    message(paste("    气压单位转换：Pa → kPa (平均", round(avg_ps/1000, 2), "kPa)"))
    
    # 确保气压在合理范围内
    full_data$ps_corrected <- pmax(80, pmin(110, full_data$ps_corrected))
  } else {
    message("  估算气压...")
    # 基于海拔估算气压（标准海平面气压101.3 kPa，每升高100米降低1.2 kPa）
    # 假设平均海拔50米
    full_data$ps_corrected <- 101.3 - (50 / 100 * 1.2)
  }
  
  # 4. 插补缺失数据
  message("  插补缺失数据...")
  
  # 需要插补的变量
  vars_to_interpolate <- c(
    "tasmax_corrected", "tasmin_corrected", "tas_corrected",
    "pr_corrected", "rsds_corrected", "hurs_corrected",
    "sfcwind_corrected", "ps_corrected"
  )
  
  # 检查哪些变量存在
  existing_vars <- intersect(vars_to_interpolate, names(full_data))
  
  for (var in existing_vars) {
    # 计算缺失比例
    missing_ratio <- sum(is.na(full_data[[var]])) / nrow(full_data)
    
    if (missing_ratio > 0) {
      message(paste("    插补变量", var, "（缺失", 
                   round(missing_ratio * 100, 1), "%）"))
      
      # 使用线性插值（前后各7天窗口）
      full_data[[var]] <- na.approx(full_data[[var]], 
                                    na.rm = FALSE, 
                                    maxgap = 7)
      
      # 如果还有缺失，使用前后值的平均值
      if (any(is.na(full_data[[var]]))) {
        message(paste("      使用前后值平均填补剩余缺失"))
        full_data[[var]] <- na.fill(full_data[[var]], "extend")
      }
    }
  }
  
  # 5. 计算水汽压（vapr）
  message("  计算水汽压...")
  
  # 使用修正后的温度计算饱和水汽压
  es_tmax <- 0.6108 * exp((17.27 * full_data$tasmax_corrected) / 
                           (full_data$tasmax_corrected + 237.3))
  es_tmin <- 0.6108 * exp((17.27 * full_data$tasmin_corrected) / 
                           (full_data$tasmin_corrected + 237.3))
  es <- (es_tmax + es_tmin) / 2
  
  # 计算实际水汽压
  full_data$vapr_corrected <- es * (full_data$hurs_corrected / 100)
  
  # 确保水汽压在合理范围内
  full_data$vapr_corrected <- pmax(0.1, pmin(10, full_data$vapr_corrected))
  
  message(paste("✅ 站点", site_name, "情景", scenario_name, "数据矫正完成"))
  
  return(full_data)
}


########################## 4. 批量矫正所有站点和情景的数据 ##########################
message("\n🔧 开始矫正气象数据...")

# 获取所有站点和情景组合
sites <- unique(loc_data$site)
scenarios <- unique(original_weather$scenario)

corrected_data_list <- list()
counter <- 0
total_combinations <- length(sites) * length(scenarios)

for (site in sites) {
  for (scenario in scenarios) {
    counter <- counter + 1
    message(paste("\n[", counter, "/", total_combinations, "] 矫正", site, "-", scenario))
    
    corrected_data <- correct_and_interpolate_weather(
      weather_df = original_weather,
      site_name = site,
      scenario_name = scenario
    )
    
    if (!is.null(corrected_data)) {
      corrected_data_list[[paste(site, scenario, sep = "_")]] <- corrected_data
    }
  }
}

# 合并所有矫正后的数据
if (length(corrected_data_list) > 0) {
  corrected_weather <- bind_rows(corrected_data_list)
  
  # 保存矫正后的数据
  corrected_file <- "all_sites_ssp126_ssp585_2001_2050_corrected.csv"
  write.csv(corrected_weather, corrected_file, row.names = FALSE)
  message(paste("\n✅ 矫正后的数据保存：", corrected_file))
  message(paste("📊 数据行数：", nrow(corrected_weather)))
} else {
  stop("❌ 没有成功矫正任何数据")
}


########################## 5. 数据质量检查与统计 ##########################
message("\n📊 矫正数据质量检查...")

# 基本统计信息
cat("\n=== 矫正数据基本统计 ===\n")
cat("时间范围:", as.character(min(corrected_weather$time)), "至", 
    as.character(max(corrected_weather$time)), "\n")
cat("站点数量:", length(unique(corrected_weather$site)), "\n")
cat("气候情景:", paste(unique(corrected_weather$scenario), collapse = ", "), "\n")
cat("总数据行数:", nrow(corrected_weather), "\n\n")

# 各变量统计
if ("tasmax_corrected" %in% names(corrected_weather)) {
  cat("最高温度(°C):\n")
  cat("  最小值:", round(min(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), "\n")
  cat("  平均值:", round(mean(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), "\n")
  cat("  最大值:", round(max(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), "\n")
  cat("  缺失值:", sum(is.na(corrected_weather$tasmax_corrected)), "\n\n")
}

if ("tasmin_corrected" %in% names(corrected_weather)) {
  cat("最低温度(°C):\n")
  cat("  最小值:", round(min(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), "\n")
  cat("  平均值:", round(mean(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), "\n")
  cat("  最大值:", round(max(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), "\n")
  cat("  缺失值:", sum(is.na(corrected_weather$tasmin_corrected)), "\n\n")
}

if ("pr_corrected" %in% names(corrected_weather)) {
  cat("降水(mm/day):\n")
  cat("  最小值:", round(min(corrected_weather$pr_corrected, na.rm = TRUE), 1), "\n")
  cat("  平均值:", round(mean(corrected_weather$pr_corrected, na.rm = TRUE), 2), "\n")
  cat("  最大值:", round(max(corrected_weather$pr_corrected, na.rm = TRUE), 1), "\n")
  cat("  总降水量:", round(sum(corrected_weather$pr_corrected, na.rm = TRUE), 0), "mm\n")
  cat("  缺失值:", sum(is.na(corrected_weather$pr_corrected)), "\n\n")
}

if ("rsds_corrected" %in% names(corrected_weather)) {
  cat("太阳辐射(MJ/m²/day):\n")
  cat("  最小值:", round(min(corrected_weather$rsds_corrected, na.rm = TRUE), 1), "\n")
  cat("  平均值:", round(mean(corrected_weather$rsds_corrected, na.rm = TRUE), 1), "\n")
  cat("  最大值:", round(max(corrected_weather$rsds_corrected, na.rm = TRUE), 1), "\n")
  cat("  缺失值:", sum(is.na(corrected_weather$rsds_corrected)), "\n\n")
}

if ("hurs_corrected" %in% names(corrected_weather)) {
  cat("相对湿度(%):\n")
  cat("  最小值:", round(min(corrected_weather$hurs_corrected, na.rm = TRUE), 1), "\n")
  cat("  平均值:", round(mean(corrected_weather$hurs_corrected, na.rm = TRUE), 1), "\n")
  cat("  最大值:", round(max(corrected_weather$hurs_corrected, na.rm = TRUE), 1), "\n")
  cat("  缺失值:", sum(is.na(corrected_weather$hurs_corrected)), "\n\n")
}

# 检查数据逻辑错误
message("\n🔍 数据逻辑检查...")

# 检查最高温度是否低于最低温度
if (all(c("tasmax_corrected", "tasmin_corrected") %in% names(corrected_weather))) {
  logic_errors <- sum(corrected_weather$tasmax_corrected < corrected_weather$tasmin_corrected, na.rm = TRUE)
  if (logic_errors > 0) {
    warning(paste("⚠️  发现", logic_errors, "行数据最高温度低于最低温度"))
  } else {
    message("✅ 最高/最低温度逻辑正确")
  }
}

# 检查负值
negative_vars <- c()
for (var in c("pr_corrected", "rsds_corrected", "hurs_corrected", "sfcwind_corrected")) {
  if (var %in% names(corrected_weather)) {
    neg_count <- sum(corrected_weather[[var]] < 0, na.rm = TRUE)
    if (neg_count > 0) {
      negative_vars <- c(negative_vars, paste(var, "(", neg_count, "个负值)"))
    }
  }
}

if (length(negative_vars) > 0) {
  warning(paste("⚠️  发现负值：", paste(negative_vars, collapse = ", ")))
} else {
  message("✅ 所有变量无负值")
}

# 检查缺失值
missing_summary <- data.frame(
  Variable = character(),
  Missing_Count = integer(),
  Missing_Percent = numeric(),
  stringsAsFactors = FALSE
)

for (var in c("tasmax_corrected", "tasmin_corrected", "pr_corrected", 
              "rsds_corrected", "hurs_corrected", "sfcwind_corrected")) {
  if (var %in% names(corrected_weather)) {
    missing_count <- sum(is.na(corrected_weather[[var]]))
    missing_percent <- round(missing_count / nrow(corrected_weather) * 100, 2)
    missing_summary <- rbind(missing_summary, 
                            data.frame(Variable = var, 
                                       Missing_Count = missing_count,
                                       Missing_Percent = missing_percent))
  }
}

if (nrow(missing_summary) > 0) {
  cat("\n缺失值统计:\n")
  print(missing_summary)
  
  # 检查是否有未插补的缺失值
  unprocessed_missing <- missing_summary[missing_summary$Missing_Count > 0, ]
  if (nrow(unprocessed_missing) > 0) {
    warning("⚠️  仍有未处理的缺失值，建议检查数据插补过程")
  } else {
    message("✅ 所有缺失值已成功插补")
  }
}

# 保存数据质量报告
quality_report <- file.path(getwd(), "data_quality_report.txt")
sink(quality_report)

cat("气象数据矫正质量报告\n")
cat("==================\n\n")
cat("报告生成时间:", as.character(Sys.time()), "\n")
cat("原始数据文件:", weather_file, "\n")
cat("矫正数据文件:", corrected_file, "\n\n")

cat("1. 数据概况\n")
cat("   时间范围:", as.character(min(corrected_weather$time)), "至", 
    as.character(max(corrected_weather$time)), "\n")
cat("   站点数量:", length(unique(corrected_weather$site)), "\n")
cat("   气候情景:", paste(unique(corrected_weather$scenario), collapse = ", "), "\n")
cat("   总数据行数:", nrow(corrected_weather), "\n\n")

cat("2. 单位转换说明\n")
cat("   - 温度: 如原始数据>100°C则视为开尔文(K)，已转换为摄氏温度(°C)\n")
cat("   - 降水: 如平均值<0.01则视为kg/m²/s，已转换为mm/day\n")
cat("   - 太阳辐射: W/m² 转换为 MJ/m²/day (×0.0864)\n")
cat("   - 气压: Pa 转换为 kPa (÷1000)\n")
cat("   - 水汽压: 基于温度和相对湿度计算\n\n")

cat("3. 数据逻辑检查结果\n")
if (exists("logic_errors")) {
  cat("   最高温度低于最低温度:", logic_errors, "行\n")
}
cat("   负值检查:", ifelse(length(negative_vars) > 0, 
                         paste("发现", length(negative_vars), "个变量有负值"), 
                         "无负值"), "\n")
cat("   缺失值插补: 使用线性插值(7天窗口)和前后值平均填补\n\n")

cat("4. 各变量统计摘要\n")
if ("tasmax_corrected" %in% names(corrected_weather)) {
  cat("   最高温度: ", round(mean(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), 
      "°C (", round(min(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), 
      "~", round(max(corrected_weather$tasmax_corrected, na.rm = TRUE), 1), ")\n")
}
if ("tasmin_corrected" %in% names(corrected_weather)) {
  cat("   最低温度: ", round(mean(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), 
      "°C (", round(min(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), 
      "~", round(max(corrected_weather$tasmin_corrected, na.rm = TRUE), 1), ")\n")
}
if ("pr_corrected" %in% names(corrected_weather)) {
  cat("   降水: ", round(mean(corrected_weather$pr_corrected, na.rm = TRUE), 2), 
      "mm/day (总量", round(sum(corrected_weather$pr_corrected, na.rm = TRUE), 0), "mm)\n")
}
if ("rsds_corrected" %in% names(corrected_weather)) {
  cat("   太阳辐射: ", round(mean(corrected_weather$rsds_corrected, na.rm = TRUE), 1), 
      "MJ/m²/day\n")
}
if ("hurs_corrected" %in% names(corrected_weather)) {
  cat("   相对湿度: ", round(mean(corrected_weather$hurs_corrected, na.rm = TRUE), 1), 
      "%\n")
}

cat("\n5. 处理步骤\n")
cat("   1) 缺失日期检测与填补\n")
cat("   2) 单位转换与标准化\n")
cat("   3) 数据逻辑检查与修正\n")
cat("   4) 缺失数据线性插补\n")
cat("   5) 水汽压计算\n")
cat("   6) 数据质量验证\n")

sink()

message(paste("\n✅ 数据质量报告保存：", quality_report))

# 输出最终提示
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("✅ 气象数据矫正完成！\n")
cat("📊 主要成果：\n")
cat("1. ✅ 矫正后的气象数据：", corrected_file, "\n")
cat("2. ✅ 数据质量报告：", quality_report, "\n")
cat("\n📋 矫正内容包括：\n")
cat("• 单位转换（温度、降水、辐射、气压）\n")
cat("• 数据逻辑检查与修正\n")
cat("• 缺失数据插补\n")
cat("• 衍生变量计算（水汽压）\n")
cat("• 数据质量验证\n")
cat("\n⚠️  请检查数据质量报告，确保数据符合APSIM输入要求。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")