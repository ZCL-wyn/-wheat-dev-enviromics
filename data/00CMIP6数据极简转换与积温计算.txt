###############################################################################
########################## CMIP6数据转换与日积温计算脚本 ########################
###############################################################################
# 作者：AI助手
# 版本：v1.2
# 日期：2024年
# 描述：读取CMIP6数据，进行单位转换，计算日积温相关变量
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
########################## 2. 站点经纬度数据 ##################################
###############################################################################
message("📍 设置站点经纬度数据...")

loc_data <- data.frame(
  site = c("NY", "PY", "YL", "ZMD"),
  latitude = c(32.9126, 35.7627, 34.2879, 33.014),
  longitude = c(112.4626, 115.0292, 108.0001, 114.0219),
  stringsAsFactors = FALSE
)

print(loc_data)

###############################################################################
########################## 3. 积温计算函数 ####################################
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

# 2) 普通积温（基温+封顶）
calc_GDD_simple <- function(Tmean, Tbase1 = 0, Tbase2 = 36) {
  Tmean <- as.numeric(Tmean)
  Tcap <- pmin(Tmean, Tbase2)  # 封顶温度
  pmax(0, Tcap - Tbase1)       # 减去基温，确保非负
}

# 3) 日长计算函数（FAO56天文公式）
calc_daylength_hours <- function(date, lat_deg) {
  # 年内日序 1..365/366
  J <- lubridate::yday(as.Date(date))
  
  # 纬度（弧度）
  phi <- lat_deg * pi / 180
  
  # 太阳赤纬（弧度）- FAO56公式
  delta <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  
  # 日落时角（弧度）
  ws_arg <- -tan(phi) * tan(delta)
  
  # 数值裁剪，避免 acos NaN
  ws_arg <- pmin(1, pmax(-1, ws_arg))
  
  # 计算日落时角
  ws <- acos(ws_arg)
  
  # 日长（小时）
  24 / pi * ws
}

###############################################################################
########################## 4. 主处理函数 ######################################
###############################################################################
process_cmip6_data <- function(input_file, loc_df, 
                               Tbase1 = 0, Topt1 = 15, 
                               Topt2 = 25, Tbase2 = 36) {
  
  message("\n🔄 开始处理CMIP6数据...")
  message("📂 输入文件：", input_file)
  
  # 检查输入文件是否存在
  if (!file.exists(input_file)) {
    stop(paste("❌ 输入文件不存在：", input_file))
  }
  
  # 读取数据
  message("📊 读取数据...")
  cmip6_data <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # 检查必要列
  required_cols <- c("site", "scenario", "time", "tas", "tasmax", "tasmin", 
                     "rsds", "pr", "hurs", "sfcwind", "huss", "ps")
  missing_cols <- setdiff(required_cols, names(cmip6_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("❌ 输入数据缺失必要列：", paste(missing_cols, collapse = ", ")))
  }
  
  # 转换时间格式
  cmip6_data$time <- as.Date(cmip6_data$time)
  
  # 检查站点 - 改为stop确保所有站点都有经纬度
  missing_sites <- setdiff(unique(cmip6_data$site), loc_df$site)
  if (length(missing_sites) > 0) {
    stop(paste("❌ 以下站点在经纬度表中缺失：", paste(missing_sites, collapse = ", ")))
  }
  
  message("✅ 数据读取完成")
  message(paste("  行数：", nrow(cmip6_data)))
  message(paste("  列数：", ncol(cmip6_data)))
  message(paste("  时间范围：", min(cmip6_data$time), "至", max(cmip6_data$time)))
  message(paste("  包含站点：", paste(unique(cmip6_data$site), collapse = ", ")))
  message(paste("  包含情景：", paste(unique(cmip6_data$scenario), collapse = ", ")))
  
  #############################################################################
  # 数据处理和转换
  #############################################################################
  message("\n🔄 开始单位转换和变量计算...")
  
  # 合并纬度信息
  result_data <- cmip6_data %>%
    left_join(loc_df %>% select(site, latitude), by = "site")
  
  # 温度转换：K → °C
  result_data <- result_data %>%
    mutate(
      TMEAN = round(tas - 273.15, 2),      # 平均温度
      TMAX = round(tasmax - 273.15, 2),    # 最高温度
      TMIN = round(tasmin - 273.15, 2)     # 最低温度
    )
  
  # 修复TMIN > TMAX的问题
  bad_records <- result_data$TMIN > result_data$TMAX
  n_bad <- sum(bad_records, na.rm = TRUE)
  
  if (n_bad > 0) {
    message(paste("⚠️ 发现", n_bad, "条TMIN>TMAX的记录，自动交换修复"))
    result_data <- result_data %>%
      mutate(
        temp = ifelse(TMIN > TMAX, TMIN, NA_real_),
        TMIN = ifelse(TMIN > TMAX, TMAX, TMIN),
        TMAX = ifelse(!is.na(temp), temp, TMAX)
      ) %>%
      select(-temp)
  }
  
  # 短波辐射转换：W/m² → MJ/m²/day
  result_data$ASSD <- round(result_data$rsds * 0.0864, 4)
  
  # 降水转换：kg m⁻² s⁻¹ → mm/day
  result_data$PR <- round(result_data$pr * 86400, 4)
  
  # 相对湿度：裁剪到[0,100]范围
  result_data$RH <- round(pmin(100, pmax(0, result_data$hurs)), 2)
  
  # 风速：直接使用（单位已是m/s）
  result_data$WS2M <- round(result_data$sfcwind, 3)
  
  # 比湿：直接使用（单位已是kg/kg）
  result_data$QV2M <- round(result_data$huss, 6)
  
  # 日温度较差
  result_data$DTR <- round(result_data$TMAX - result_data$TMIN, 2)
  
  # 光合有效辐射
  result_data$APAR <- round(result_data$ASSD * 0.48, 4)
  
  # 日长计算
  result_data$DL <- round(calc_daylength_hours(result_data$time, result_data$latitude), 3)
  
  # 露点温度计算
  result_data <- result_data %>%
    mutate(
      # 饱和水汽压 (kPa)
      es_kPa = 0.6108 * exp(17.27 * TMEAN / (TMEAN + 237.3)),
      # 实际水汽压 (kPa) - 使用裁剪后的RH
      ea_kPa = (RH / 100) * es_kPa,
      # 确保ea_kPa大于0（避免log(0)）
      ea_kPa_safe = pmax(ea_kPa, 1e-6),
      # 露点温度
      ln_term = log(ea_kPa_safe / 0.6108),
      T2MDEW = round((237.3 * ln_term) / (17.27 - ln_term), 2)
    ) %>%
    select(-es_kPa, -ea_kPa, -ea_kPa_safe, -ln_term)
  
  # 水汽压差 VPD (kPa)
  result_data <- result_data %>%
    mutate(
      es_kPa = 0.6108 * exp(17.27 * TMEAN / (TMEAN + 237.3)),
      ea_kPa = (RH / 100) * es_kPa,
      VPD = round(pmax(es_kPa - ea_kPa, 0), 4)
    ) %>%
    select(-es_kPa, -ea_kPa)
  
  #############################################################################
  # 日积温计算 - 统一使用(TMAX+TMIN)/2作为日均温
  #############################################################################
  message("\n🌡️ 计算日积温相关变量...")
  message(paste("  使用温度参数：Tbase1 =", Tbase1, "°C, Topt1 =", Topt1, 
                "°C, Topt2 =", Topt2, "°C, Tbase2 =", Tbase2, "°C"))
  
  # 计算日平均温度（用于积温计算）- 统一使用(TMAX+TMIN)/2
  result_data$Tmean_for_gdd <- (result_data$TMAX + result_data$TMIN) / 2
  
  # 生长度日 GDD (以0°C为基准) - 统一使用Tmean_for_gdd
  result_data$GDD_base0 <- round(pmax(0, result_data$Tmean_for_gdd - 0), 2)
  
  # 计算普通积温（有封顶）
  result_data$GDD_simple <- round(
    calc_GDD_simple(result_data$Tmean_for_gdd, Tbase1 = Tbase1, Tbase2 = Tbase2), 
    2
  )
  
  # 计算温度响应（0~1）
  result_data$FRUE <- round(
    calc_FRUE(result_data$Tmean_for_gdd, Tbase1 = Tbase1, Topt1 = Topt1, 
              Topt2 = Topt2, Tbase2 = Tbase2), 
    3
  )
  
  # 计算门限积温：有效温度 × 响应
  result_data$GDD_cardinal <- round(result_data$GDD_simple * result_data$FRUE, 2)
  
  # 移除中间变量
  result_data$Tmean_for_gdd <- NULL
  
  #############################################################################
  # 整理输出列
  #############################################################################
  message("\n📋 整理输出数据...")
  
  # 定义输出列顺序
  final_columns <- c(
    "site", "scenario", "time",
    "TMEAN", "TMAX", "TMIN", "T2MDEW", "ASSD", "APAR", "DL", "PR",
    "QV2M", "RH", "WS2M", "VPD", 
    "GDD_base0", "GDD_simple", "FRUE", "GDD_cardinal", "DTR"
  )
  
  # 选择并重排列
  final_data <- result_data %>%
    select(all_of(final_columns[final_columns %in% names(.)])) %>%
    arrange(site, scenario, time)
  
  message("✅ 数据处理完成")
  
  return(final_data)
}

###############################################################################
########################## 5. 执行数据处理 ####################################
###############################################################################

# 设置输入文件名（请根据实际情况修改）
input_filename <- "all_sites_ssp126_ssp585_2001_2050.csv"

# 检查文件是否存在
if (!file.exists(input_filename)) {
  # 列出当前目录下的文件，帮助用户找到正确的文件名
  message("当前目录下的文件：")
  print(list.files(pattern = "\\.csv$"))
  stop(paste("❌ 找不到输入文件：", input_filename))
}

# 设置温度参数
Tbase1 <- 0   # 基温
Topt1 <- 15   # 最适温度下限
Topt2 <- 25   # 最适温度上限
Tbase2 <- 36  # 上限温度

# 执行数据处理
processed_data <- process_cmip6_data(
  input_file = input_filename,
  loc_df = loc_data,
  Tbase1 = Tbase1,
  Topt1 = Topt1,
  Topt2 = Topt2,
  Tbase2 = Tbase2
)

###############################################################################
########################## 6. 保存结果 #########################################
###############################################################################
message("\n💾 保存处理结果...")

# 生成输出文件名（基于输入文件名）
output_filename <- gsub("\\.csv$", "_with_daily_thermal.csv", input_filename)
if (output_filename == input_filename) {
  output_filename <- paste0(tools::file_path_sans_ext(input_filename), "_with_daily_thermal.csv")
}

# 保存为CSV文件
write.csv(processed_data, output_filename, row.names = FALSE)

message(paste("✅ 结果保存至：", output_filename))
message(paste("📊 输出数据维度：", nrow(processed_data), "行 ×", ncol(processed_data), "列"))

###############################################################################
########################## 7. 数据质量检查 #####################################
###############################################################################
message("\n🔍 数据质量检查...")

# 1. 基本统计信息
message("📈 基本统计信息：")
cat(paste("  时间范围：", min(processed_data$time), "至", max(processed_data$time), "\n"))
cat(paste("  包含站点：", paste(unique(processed_data$site), collapse = ", "), "\n"))
cat(paste("  包含情景：", paste(unique(processed_data$scenario), collapse = ", "), "\n"))
cat(paste("  总记录数：", nrow(processed_data), "\n"))

# 2. 检查TMIN和TMAX
tmin_tmax_check <- processed_data %>%
  summarise(
    TMIN_gt_TMAX = sum(TMIN > TMAX, na.rm = TRUE),
    TMIN_min = round(min(TMIN, na.rm = TRUE), 1),
    TMAX_max = round(max(TMAX, na.rm = TRUE), 1),
    TMEAN_range = paste(round(min(TMEAN, na.rm = TRUE), 1), "-", round(max(TMEAN, na.rm = TRUE), 1))
  )

message("🌡️ 温度检查：")
cat(paste("  TMIN > TMAX的记录数：", tmin_tmax_check$TMIN_gt_TMAX, "\n"))
cat(paste("  TMIN最小值：", tmin_tmax_check$TMIN_min, "°C\n"))
cat(paste("  TMAX最大值：", tmin_tmax_check$TMAX_max, "°C\n"))
cat(paste("  TMEAN范围：", tmin_tmax_check$TMEAN_range, "°C\n"))

# 3. 积温相关统计
thermal_stats <- processed_data %>%
  summarise(
    GDD_simple_mean = round(mean(GDD_simple, na.rm = TRUE), 2),
    GDD_simple_min = round(min(GDD_simple, na.rm = TRUE), 2),
    GDD_simple_max = round(max(GDD_simple, na.rm = TRUE), 2),
    FRUE_mean = round(mean(FRUE, na.rm = TRUE), 3),
    FRUE_min = round(min(FRUE, na.rm = TRUE), 3),
    FRUE_max = round(max(FRUE, na.rm = TRUE), 3),
    FRUE_zero_days = sum(FRUE == 0, na.rm = TRUE),
    FRUE_one_days = sum(FRUE == 1, na.rm = TRUE),
    GDD_cardinal_mean = round(mean(GDD_cardinal, na.rm = TRUE), 2)
  )

message("🌡️ 日积温统计：")
cat(paste("  GDD_simple平均值：", thermal_stats$GDD_simple_mean, "°C-day\n"))
cat(paste("  GDD_simple范围：", thermal_stats$GDD_simple_min, "-", thermal_stats$GDD_simple_max, "°C-day\n"))
cat(paste("  FRUE平均值：", thermal_stats$FRUE_mean, "\n"))
cat(paste("  FRUE范围：", thermal_stats$FRUE_min, "-", thermal_stats$FRUE_max, "\n"))
cat(paste("  FRUE=0的天数：", thermal_stats$FRUE_zero_days, "\n"))
cat(paste("  FRUE=1的天数：", thermal_stats$FRUE_one_days, "\n"))
cat(paste("  GDD_cardinal平均值：", thermal_stats$GDD_cardinal_mean, "°C-day\n"))

# 4. 其他变量统计
other_stats <- processed_data %>%
  summarise(
    PR_total = round(sum(PR, na.rm = TRUE), 1),
    ASSD_mean = round(mean(ASSD, na.rm = TRUE), 2),
    RH_mean = round(mean(RH, na.rm = TRUE), 1),
    RH_min = round(min(RH, na.rm = TRUE), 1),
    RH_max = round(max(RH, na.rm = TRUE), 1),
    DL_min = round(min(DL, na.rm = TRUE), 2),
    DL_max = round(max(DL, na.rm = TRUE), 2),
    DL_mean = round(mean(DL, na.rm = TRUE), 2),
    VPD_mean = round(mean(VPD, na.rm = TRUE), 3)
  )

message("📊 其他变量统计：")
cat(paste("  总降水量：", other_stats$PR_total, "mm\n"))
cat(paste("  平均短波辐射：", other_stats$ASSD_mean, "MJ/m²/day\n"))
cat(paste("  平均相对湿度：", other_stats$RH_mean, "% (范围：", 
          other_stats$RH_min, "-", other_stats$RH_max, ")\n"))
cat(paste("  日长范围：", other_stats$DL_min, "-", other_stats$DL_max, 
          "小时 (平均：", other_stats$DL_mean, "小时)\n"))
cat(paste("  平均水汽压差：", other_stats$VPD_mean, "kPa\n"))

# 5. 分站点统计
message("\n🌍 分站点统计（平均值）：")
site_summary <- processed_data %>%
  group_by(site, scenario) %>%
  summarise(
    records = n(),
    avg_TMEAN = round(mean(TMEAN, na.rm = TRUE), 1),
    avg_GDD_simple = round(mean(GDD_simple, na.rm = TRUE), 2),
    avg_FRUE = round(mean(FRUE, na.rm = TRUE), 3),
    avg_GDD_cardinal = round(mean(GDD_cardinal, na.rm = TRUE), 2),
    total_PR = round(sum(PR, na.rm = TRUE), 1),
    avg_DL = round(mean(DL, na.rm = TRUE), 2),
    .groups = "drop"
  )

print(site_summary)

# 6. 数据预览
message("\n👀 数据预览（前10行）：")
print(head(processed_data, 10))

###############################################################################
########################## 8. 生成处理日志 #####################################
###############################################################################
message("\n📝 生成处理日志...")

log_filename <- gsub("\\.csv$", "_processing_log.txt", output_filename)
if (log_filename == output_filename) {
  log_filename <- paste0(tools::file_path_sans_ext(output_filename), "_log.txt")
}

sink(log_filename)

cat("CMIP6数据转换与日积温计算处理日志\n")
cat("==================================================\n\n")
cat("处理时间：", as.character(Sys.time()), "\n")
cat("工作目录：", work_dir, "\n")
cat("输入文件：", input_filename, "\n")
cat("输出文件：", output_filename, "\n\n")

cat("数据概览：\n")
cat("--------------------------------------------\n")
cat(paste("行数：", nrow(processed_data), "\n"))
cat(paste("列数：", ncol(processed_data), "\n"))
cat(paste("时间范围：", min(processed_data$time), "至", max(processed_data$time), "\n"))
cat(paste("包含站点：", paste(unique(processed_data$site), collapse = ", "), "\n"))
cat(paste("包含情景：", paste(unique(processed_data$scenario), collapse = ", "), "\n\n"))

cat("温度参数设置：\n")
cat("--------------------------------------------\n")
cat(paste("Tbase1（基温）：", Tbase1, "°C\n"))
cat(paste("Topt1（最适温度下限）：", Topt1, "°C\n"))
cat(paste("Topt2（最适温度上限）：", Topt2, "°C\n"))
cat(paste("Tbase2（上限温度）：", Tbase2, "°C\n\n"))

cat("输出变量说明（20个变量）：\n")
cat("--------------------------------------------\n")
cat("1. site          - 站点\n")
cat("2. scenario      - 情景（ssp126/ssp585等）\n")
cat("3. time          - 日期\n")
cat("4. TMEAN         - 日平均温度 (°C) [来自CMIP6 tas]\n")
cat("5. TMAX          - 日最高温度 (°C)\n")
cat("6. TMIN          - 日最低温度 (°C)\n")
cat("7. T2MDEW        - 露点温度 (°C)\n")
cat("8. ASSD          - 短波辐射 (MJ/m²/day)\n")
cat("9. APAR          - 光合有效辐射 (MJ/m²/day)\n")
cat("10. DL           - 日长 (小时)\n")
cat("11. PR           - 降水量 (mm/day)\n")
cat("12. QV2M         - 比湿 (kg/kg)\n")
cat("13. RH           - 相对湿度 (%)\n")
cat("14. WS2M         - 风速 (m/s)\n")
cat("15. VPD          - 水汽压差 (kPa)\n")
cat("16. GDD_base0    - 生长度日（0°C基准）(°C-day) [使用(TMAX+TMIN)/2]\n")
cat("17. GDD_simple   - 普通积温（有封顶）(°C-day) [使用(TMAX+TMIN)/2]\n")
cat("18. FRUE         - 温度响应函数 (0-1) [使用(TMAX+TMIN)/2]\n")
cat("19. GDD_cardinal - 门限积温 (°C-day)\n")
cat("20. DTR          - 日温度较差 (°C)\n\n")

cat("单位转换规则：\n")
cat("--------------------------------------------\n")
cat("• 温度：°C = K - 273.15\n")
cat("• 短波辐射：MJ/m²/day = W/m² × 0.0864\n")
cat("• 光合有效辐射：APAR = ASSD × 0.48\n")
cat("• 降水：mm/day = kg m⁻² s⁻¹ × 86400\n")
cat("• 日长：基于FAO56天文公式计算\n\n")

cat("积温计算规则：\n")
cat("--------------------------------------------\n")
cat("• 日均温计算：所有积温计算均使用 (TMAX+TMIN)/2\n")
cat("• GDD_base0：max(0, (TMAX+TMIN)/2 - 0°C)，无封顶\n")
cat("• GDD_simple：max(0, min((TMAX+TMIN)/2, Tbase2) - Tbase1)，有封顶\n")
cat("• FRUE：四参数温度响应函数（线性上升/平台/下降）\n")
cat("• GDD_cardinal：GDD_simple × FRUE，加权积温\n\n")

cat("数据质量检查结果：\n")
cat("--------------------------------------------\n")
cat(paste("TMIN > TMAX的记录数：", tmin_tmax_check$TMIN_gt_TMAX, "\n"))
cat(paste("温度范围：", tmin_tmax_check$TMEAN_range, "°C\n"))
cat(paste("FRUE=0的天数：", thermal_stats$FRUE_zero_days, "\n"))
cat(paste("FRUE=1的天数：", thermal_stats$FRUE_one_days, "\n\n"))

cat("分站点统计：\n")
cat("--------------------------------------------\n")
print(site_summary)

sink()

message(paste("✅ 处理日志保存至：", log_filename))

###############################################################################
########################## 9. 完成提示 #########################################
###############################################################################
message(paste0("\n", paste(rep("=", 50), collapse = "")))
message("🎉 处理完成！")
message(paste("📁 输出文件：", output_filename))
message(paste("📝 处理日志：", log_filename))
message("🌡️ 包含日积温变量：GDD_simple, FRUE, GDD_cardinal")
message(paste("📊 总变量数：", ncol(processed_data)))
message(paste("⏰ 完成时间：", Sys.time()))
message(paste(rep("=", 50), collapse = ""))

# 显示输出文件的绝对路径
message("\n💡 提示：输出文件的完整路径：")
message(normalizePath(output_filename))