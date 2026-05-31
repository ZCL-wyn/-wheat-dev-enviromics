###############################################################################
########################## NASA数据多套日积温计算脚本 ##########################
###############################################################################
# 作者：AI助手
# 版本：v2.0
# 日期：2024年
# 描述：读取NASA数据，计算多套日积温方案
# 注意：不做met文件，只计算日积温，不计算累计积温
###############################################################################

# 清除环境变量
rm(list = ls(all = TRUE))
cat("\014")

# 设置工作目录
work_dir <- "C:/Users/Lenovo/Desktop/小麦千粒重文章/001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割"
setwd(work_dir)

message("==================================================")
message("📁 工作目录：", work_dir)
message("📅 当前时间：", Sys.time())
message("==================================================")

###############################################################################
########################## 1. 加载必要包 ######################################
###############################################################################
message("📦 加载必要的R包...")

# 检查并安装缺失的包
required_packages <- c("dplyr", "lubridate", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0) {
  message("安装缺失的包：", paste(new_packages, collapse = ", "))
  install.packages(new_packages, dependencies = TRUE)
}

# 加载包
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
})

message("✅ 包加载完成")

###############################################################################
########################## 2. 积温计算函数 ####################################
###############################################################################

# 1) 四参数温度响应 FRUE：0~1
calc_FRUE <- function(T, Tbase1 = 0, Topt1 = 15, Topt2 = 25, Tbase2 = 36) {
  # 保证数值
  T <- as.numeric(T)
  
  # 初始化结果向量
  fr <- numeric(length(T))
  
  # 上升段：Tbase1 -> Topt1
  idx1 <- !is.na(T) & T > Tbase1 & T < Topt1
  fr[idx1] <- (T[idx1] - Tbase1) / (Topt1 - Tbase1)
  
  # 平台段：Topt1 -> Topt2
  idx2 <- !is.na(T) & T >= Topt1 & T <= Topt2
  fr[idx2] <- 1
  
  # 下降段：Topt2 -> Tbase2
  idx3 <- !is.na(T) & T > Topt2 & T < Tbase2
  fr[idx3] <- (Tbase2 - T[idx3]) / (Tbase2 - Topt2)
  
  # 其他情况保持0
  pmin(1, pmax(0, fr))
}

# 2) 简单积温（基温+封顶）
calc_GDD_simple <- function(Tmean, Tbase1 = 0, Tbase2 = 36) {
  Tmean <- as.numeric(Tmean)
  Tcap <- pmin(Tmean, Tbase2)  # 封顶温度
  pmax(0, Tcap - Tbase1)       # 减去基温，确保非负
}

# 3) 三基点积温（加权积温）
calc_GDD_cardinal <- function(Tmean, Tbase1 = 0, Topt1 = 15, Topt2 = 25, Tbase2 = 36) {
  Tmean <- as.numeric(Tmean)
  
  # 计算简单积温
  gdd_simple <- calc_GDD_simple(Tmean, Tbase1, Tbase2)
  
  # 计算温度响应
  frue <- calc_FRUE(Tmean, Tbase1, Topt1, Topt2, Tbase2)
  
  # 加权积温
  gdd_simple * frue
}

###############################################################################
########################## 3. 主处理函数 ######################################
###############################################################################
process_nasa_data <- function(input_file) {
  
  message("\n🔄 开始处理NASA数据...")
  message("📂 输入文件：", input_file)
  
  # 检查输入文件是否存在
  if (!file.exists(input_file)) {
    stop(paste("❌ 输入文件不存在：", input_file))
  }
  
  # 读取数据
  message("📊 读取数据...")
  nasa_data <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # 检查必要列
  required_cols <- c("YYYYMMDD", "TMAX", "TMIN", "LAT")
  missing_cols <- setdiff(required_cols, names(nasa_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("❌ 输入数据缺失必要列：", paste(missing_cols, collapse = ", ")))
  }
  
  # 转换日期格式
  if(is.numeric(nasa_data$YYYYMMDD)) {
    nasa_data$date <- as.Date(as.character(nasa_data$YYYYMMDD), format = "%Y%m%d")
  } else {
    nasa_data$date <- as.Date(nasa_data$YYYYMMDD)
  }
  
  message("✅ 数据读取完成")
  message(paste("  行数：", nrow(nasa_data)))
  message(paste("  列数：", ncol(nasa_data)))
  message(paste("  时间范围：", min(nasa_data$date), "至", max(nasa_data$date)))
  message(paste("  纬度范围：", min(nasa_data$LAT, na.rm = TRUE), "至", max(nasa_data$LAT, na.rm = TRUE)))
  
  #############################################################################
  # 日积温计算
  #############################################################################
  message("\n🌡️ 计算日积温相关变量...")
  
  # 计算日平均温度（用于积温计算）- 使用(TMAX+TMIN)/2
  nasa_data$Tmean_for_gdd <- (nasa_data$TMAX + nasa_data$TMIN) / 2
  
  # 覆盖原有的GDD和FRUE（如果需要）
  # 计算新的GDD（按您的要求，使用简单积温）
  nasa_data$GDD <- round(calc_GDD_simple(nasa_data$Tmean_for_gdd, Tbase1 = 0, Tbase2 = 36), 2)
  
  # 计算新的FRUE
  nasa_data$FRUE <- round(calc_FRUE(nasa_data$Tmean_for_gdd, Tbase1 = 0, Topt1 = 15, 
                                    Topt2 = 25, Tbase2 = 36), 3)
  
  # 计算新增的三列积温
  # 1. GDD_base0: 0℃基准，无封顶
  nasa_data$GDD_base0 <- round(pmax(0, nasa_data$Tmean_for_gdd - 0), 2)
  
  # 2. GDD_simple: 简单积温（0℃基准，36℃封顶）
  nasa_data$GDD_simple <- round(calc_GDD_simple(nasa_data$Tmean_for_gdd, Tbase1 = 0, Tbase2 = 36), 2)
  
  # 3. GDD_cardinal: 三基点积温（加权积温）
  nasa_data$GDD_cardinal <- round(
    calc_GDD_cardinal(nasa_data$Tmean_for_gdd, Tbase1 = 0, Topt1 = 15, 
                      Topt2 = 25, Tbase2 = 36), 
    2
  )
  
  # 移除中间变量
  nasa_data$Tmean_for_gdd <- NULL
  nasa_data$date <- NULL
  
  #############################################################################
  # 整理输出列顺序
  #############################################################################
  message("\n📋 整理输出数据...")
  
  # 定义输出列顺序（按照您的要求）
  final_columns <- c(
    "YYYYMMDD", "LON", "LAT", "env", "DOY.x", "daysFromStart",
    "TMEAN", "TMAX", "TMIN", "T2MDEW", "ASLD", "ASSD", "APAR",
    "PR", "EVPTRNS", "QV2M", "RH", "GW", "WS2M", "P_ETP",
    "VPD", "DL", "RTA", "PAR_TEMP", "ETP", "PETP",
    "GDD", "FRUE", "DTR", 
    "GDD_base0", "GDD_simple", "GDD_cardinal"
  )
  
  # 只保留实际存在于数据中的列
  existing_columns <- final_columns[final_columns %in% names(nasa_data)]
  missing_in_output <- setdiff(final_columns, existing_columns)
  extra_columns <- setdiff(names(nasa_data), final_columns)
  
  if(length(missing_in_output) > 0) {
    message("⚠️  以下列在输出顺序中指定但数据中不存在：")
    message(paste("  ", paste(missing_in_output, collapse = ", ")))
  }
  
  if(length(extra_columns) > 0) {
    message("⚠️  以下列存在于数据中但不在输出顺序中，将被放在最后：")
    message(paste("  ", paste(extra_columns, collapse = ", ")))
    existing_columns <- c(existing_columns, extra_columns)
  }
  
  # 重新排列列
  final_data <- nasa_data %>%
    select(all_of(existing_columns))
  
  message("✅ 数据处理完成")
  
  return(final_data)
}

###############################################################################
########################## 4. 执行数据处理 ####################################
###############################################################################

# 设置输入文件名
input_filename <- "A.csv"

# 检查文件是否存在
if (!file.exists(input_filename)) {
  # 列出当前目录下的CSV文件
  csv_files <- list.files(pattern = "\\.csv$")
  if(length(csv_files) > 0) {
    message("当前目录下的CSV文件：")
    for(i in seq_along(csv_files)) {
      message(sprintf("  %d. %s", i, csv_files[i]))
    }
    input_filename <- readline(prompt = "请输入要处理的文件编号或文件名: ")
    if(grepl("^[0-9]+$", input_filename)) {
      idx <- as.integer(input_filename)
      if(idx >= 1 && idx <= length(csv_files)) {
        input_filename <- csv_files[idx]
      }
    }
  } else {
    stop("❌ 当前目录下没有CSV文件")
  }
}

message(paste("📂 将处理文件：", input_filename))

# 执行数据处理
processed_data <- process_nasa_data(input_file = input_filename)

###############################################################################
########################## 5. 保存结果 #########################################
###############################################################################
message("\n💾 保存处理结果...")

# 生成输出文件名（基于输入文件名）
output_filename <- gsub("\\.csv$", "_with_thermal.csv", input_filename)
if (output_filename == input_filename) {
  output_filename <- paste0(tools::file_path_sans_ext(input_filename), "_with_thermal.csv")
}

# 保存为CSV文件
write.csv(processed_data, output_filename, row.names = FALSE)

message(paste("✅ 结果保存至：", output_filename))
message(paste("📊 输出数据维度：", nrow(processed_data), "行 ×", ncol(processed_data), "列"))

###############################################################################
########################## 6. 数据质量检查 #####################################
###############################################################################
message("\n🔍 数据质量检查...")

# 1. 基本统计信息
message("📈 基本统计信息：")
cat(paste("  总记录数：", nrow(processed_data), "\n"))
cat(paste("  总变量数：", ncol(processed_data), "\n"))

# 检查是否有日期列并显示时间范围
if("YYYYMMDD" %in% names(processed_data)) {
  if(is.numeric(processed_data$YYYYMMDD)) {
    dates <- as.Date(as.character(processed_data$YYYYMMDD), format = "%Y%m%d")
  } else {
    dates <- as.Date(processed_data$YYYYMMDD)
  }
  cat(paste("  时间范围：", min(dates, na.rm = TRUE), "至", max(dates, na.rm = TRUE), "\n"))
}

# 2. 检查TMIN和TMAX
if(all(c("TMIN", "TMAX") %in% names(processed_data))) {
  tmin_tmax_check <- processed_data %>%
    summarise(
      TMIN_gt_TMAX = sum(TMIN > TMAX, na.rm = TRUE),
      TMIN_min = round(min(TMIN, na.rm = TRUE), 1),
      TMAX_max = round(max(TMAX, na.rm = TRUE), 1),
      TMEAN_min = ifelse("TMEAN" %in% names(.), round(min(TMEAN, na.rm = TRUE), 1), NA),
      TMEAN_max = ifelse("TMEAN" %in% names(.), round(max(TMEAN, na.rm = TRUE), 1), NA)
    )
  
  message("🌡️ 温度检查：")
  cat(paste("  TMIN > TMAX的记录数：", tmin_tmax_check$TMIN_gt_TMAX, "\n"))
  cat(paste("  TMIN最小值：", tmin_tmax_check$TMIN_min, "°C\n"))
  cat(paste("  TMAX最大值：", tmin_tmax_check$TMAX_max, "°C\n"))
  if(!is.na(tmin_tmax_check$TMEAN_min)) {
    cat(paste("  TMEAN范围：", tmin_tmax_check$TMEAN_min, "-", tmin_tmax_check$TMEAN_max, "°C\n"))
  }
}

# 3. 积温相关统计
thermal_stats <- processed_data %>%
  summarise(
    GDD_mean = ifelse("GDD" %in% names(.), round(mean(GDD, na.rm = TRUE), 2), NA),
    GDD_min = ifelse("GDD" %in% names(.), round(min(GDD, na.rm = TRUE), 2), NA),
    GDD_max = ifelse("GDD" %in% names(.), round(max(GDD, na.rm = TRUE), 2), NA),
    FRUE_mean = ifelse("FRUE" %in% names(.), round(mean(FRUE, na.rm = TRUE), 3), NA),
    FRUE_min = ifelse("FRUE" %in% names(.), round(min(FRUE, na.rm = TRUE), 3), NA),
    FRUE_max = ifelse("FRUE" %in% names(.), round(max(FRUE, na.rm = TRUE), 3), NA),
    GDD_simple_mean = ifelse("GDD_simple" %in% names(.), round(mean(GDD_simple, na.rm = TRUE), 2), NA),
    GDD_cardinal_mean = ifelse("GDD_cardinal" %in% names(.), round(mean(GDD_cardinal, na.rm = TRUE), 2), NA)
  )

message("🌡️ 日积温统计：")
if(!is.na(thermal_stats$GDD_mean)) {
  cat(paste("  GDD平均值：", thermal_stats$GDD_mean, "°C-day\n"))
  cat(paste("  GDD范围：", thermal_stats$GDD_min, "-", thermal_stats$GDD_max, "°C-day\n"))
}
if(!is.na(thermal_stats$FRUE_mean)) {
  cat(paste("  FRUE平均值：", thermal_stats$FRUE_mean, "\n"))
  cat(paste("  FRUE范围：", thermal_stats$FRUE_min, "-", thermal_stats$FRUE_max, "\n"))
}
if(!is.na(thermal_stats$GDD_simple_mean)) {
  cat(paste("  GDD_simple平均值：", thermal_stats$GDD_simple_mean, "°C-day\n"))
}
if(!is.na(thermal_stats$GDD_cardinal_mean)) {
  cat(paste("  GDD_cardinal平均值：", thermal_stats$GDD_cardinal_mean, "°C-day\n"))
}

# 4. 数据预览
message("\n👀 数据预览（前6行）：")
print(head(processed_data))

# 显示新增的列
new_thermal_cols <- c("GDD_base0", "GDD_simple", "GDD_cardinal")
existing_new_cols <- new_thermal_cols[new_thermal_cols %in% names(processed_data)]

if(length(existing_new_cols) > 0) {
  message("\n✨ 新增的积温列：")
  for(col in existing_new_cols) {
    cat(paste("  • ", col, "\n"))
  }
}

###############################################################################
########################## 7. 生成处理日志 #####################################
###############################################################################
message("\n📝 生成处理日志...")

log_filename <- gsub("\\.csv$", "_thermal_log.txt", output_filename)
if (log_filename == output_filename) {
  log_filename <- paste0(tools::file_path_sans_ext(output_filename), "_log.txt")
}

sink(log_filename)

cat("NASA数据日积温计算处理日志\n")
cat("==================================================\n\n")
cat("处理时间：", as.character(Sys.time()), "\n")
cat("工作目录：", work_dir, "\n")
cat("输入文件：", input_filename, "\n")
cat("输出文件：", output_filename, "\n\n")

cat("数据概览：\n")
cat("--------------------------------------------\n")
cat(paste("行数：", nrow(processed_data), "\n"))
cat(paste("列数：", ncol(processed_data), "\n"))
cat(paste("时间范围：", if(exists("dates")) paste(min(dates, na.rm = TRUE), "至", max(dates, na.rm = TRUE)) else "未知", "\n\n"))

cat("温度参数设置：\n")
cat("--------------------------------------------\n")
cat("使用小麦默认温度参数：\n")
cat(paste("Tbase1（基温）：", 0, "°C\n"))
cat(paste("Topt1（最适温度下限）：", 15, "°C\n"))
cat(paste("Topt2（最适温度上限）：", 25, "°C\n"))
cat(paste("Tbase2（上限温度）：", 36, "°C\n\n"))

cat("积温计算说明：\n")
cat("--------------------------------------------\n")
cat("1. 日均温计算：使用 (TMAX + TMIN) / 2\n")
cat("2. GDD：简单积温，max(0, min(Tmean, 36) - 0)\n")
cat("3. FRUE：四参数温度响应函数（0-1）\n")
cat("4. GDD_base0：0℃基准无封顶积温，max(0, Tmean - 0)\n")
cat("5. GDD_simple：与GDD相同，简单积温\n")
cat("6. GDD_cardinal：三基点积温，GDD_simple × FRUE\n\n")

cat("新增的积温变量：\n")
cat("--------------------------------------------\n")
for(col in existing_new_cols) {
  cat(paste(col, "\n"))
}
cat("\n")

cat("数据质量检查结果：\n")
cat("--------------------------------------------\n")
if(all(c("TMIN", "TMAX") %in% names(processed_data))) {
  cat(paste("TMIN > TMAX的记录数：", tmin_tmax_check$TMIN_gt_TMAX, "\n"))
}
cat("\n")

sink()

message(paste("✅ 处理日志保存至：", log_filename))

###############################################################################
########################## 8. 完成提示 #########################################
###############################################################################
message(paste0("\n", paste(rep("=", 50), collapse = "")))
message("🎉 处理完成！")
message(paste("📁 输出文件：", output_filename))
message(paste("📝 处理日志：", log_filename))
message(paste("📊 总列数：", ncol(processed_data)))
message(paste("⏰ 完成时间：", Sys.time()))
message(paste(rep("=", 50), collapse = ""))

# 显示输出文件的绝对路径
message("\n💡 输出文件的完整路径：")
message(normalizePath(output_filename))