########################## 1. 前置检查：自动安装并加载必要包 ##########################
if (!require("apsimx", quietly = TRUE)) {
  install.packages("apsimx", dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(apsimx)
}
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

# 清除环境残留与控制台输出
rm(list = ls(all = TRUE))
cat("\014")
message("✅ 前置环境清理完成")


########################## 2. 配置工作目录与APSIM路径 ##########################
work_dir <- "C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割"
if (!dir.exists(work_dir)) {
  stop(paste("❌ 工作目录不存在：", work_dir))
}
setwd(work_dir)
message(paste("✅ 工作目录：", work_dir))

# APSIM路径配置
apsimx_exe_path <- "C:/APSIM/bin/ApsimNG.exe"
apsimx_examples_path <- "C:/APSIM/Examples"

# 路径验证
if (!file.exists(apsimx_exe_path)) stop(paste("❌ APSIM可执行文件缺失：", apsimx_exe_path))
if (!file.exists(apsimx_examples_path)) stop(paste("❌ 示例目录缺失：", apsimx_examples_path))
if (!file.exists(file.path(apsimx_examples_path, "Wheat.apsimx"))) {
  stop(paste("❌ 小麦模板缺失：", file.path(apsimx_examples_path, "Wheat.apsimx")))
}

# 配置apsimx包
apsimx_options(
  exe.path = apsimx_exe_path,
  examples.path = apsimx_examples_path,
  warn.find.apsimx = FALSE
)
message("✅ APSIM路径配置正确")


########################## 3. 读取站点信息 ##########################
loc_data <- data.frame(
  site = c("NY", "PY", "YL", "ZMD"),
  latitude = c(32.9126, 35.7627, 34.2879, 33.014),
  longitude = c(112.4626, 115.0292, 108.0001, 114.0219),
  stringsAsFactors = FALSE
)
message("✅ 站点数据")
print(loc_data)


########################## 4. 创建未来气候模拟独立目录 ##########################
# 为未来气候模拟创建完全独立的目录结构
future_dirs <- list(
  weather = file.path(getwd(), "Future_APSIMmet"),
  script = file.path(getwd(), "Future_APSIMinput_apsimx"),
  output = file.path(getwd(), "Future_APSIMout"),
  config = file.path(getwd(), "Future_APSIMconfig")  # 新增：配置目录
)

# 创建所有目录
for (dir_name in names(future_dirs)) {
  if (!dir.exists(future_dirs[[dir_name]])) {
    dir.create(future_dirs[[dir_name]], recursive = TRUE)
    message(paste("✅ 新建未来模拟目录：", future_dirs[[dir_name]]))
  } else {
    # 清空目录（如果已存在）
    if (length(list.files(future_dirs[[dir_name]], all.files = TRUE, no.. = TRUE)) > 0) {
      unlink(list.files(future_dirs[[dir_name]], full.names = TRUE), recursive = TRUE)
      message(paste("⚠️ 清空目录：", future_dirs[[dir_name]]))
    }
    message(paste("✅ 未来模拟目录已存在：", future_dirs[[dir_name]]))
  }
}


########################## 5. 定义播种日期配置 ##########################
sowing_configs <- list(
  Z = list(
    name = "早播",
    sowing_date_offset = 0,  # 10月15日
    harvest_offset = 228     # 5月31日收获
  ),
  M = list(
    name = "中播",
    sowing_date_offset = 15,  # 10月30日
    harvest_offset = 227      # 6月14日收获
  ),
  W = list(
    name = "晚播",
    sowing_date_offset = 30,  # 11月14日
    harvest_offset = 227      # 6月29日收获
  )
)


########################## 6. 检查矫正后的气象数据是否存在 ##########################
corrected_file <- "all_sites_ssp126_ssp585_2001_2050_corrected.csv"
if (!file.exists(corrected_file)) {
  stop(paste("❌ 矫正后的气象数据缺失：", corrected_file))
}
message(paste("✅ 检测到矫正后的气象数据：", corrected_file))


########################## 7. 生成未来气候仿真配置（从2025-2026生长季开始）##########################
message("\n🔧 创建未来气候仿真配置...")

# 定义模拟年份：从2025年到2049年（模拟2025-2026到2049-2050生长季）
simulation_years <- 2025:2049
scenarios <- c("ssp126", "ssp585")
cultivars <- c("Keyu13")

# 创建仿真配置数据框
sim_configs <- expand.grid(
  site = loc_data$site,
  year = simulation_years,
  scenario = scenarios,
  sowing_code = names(sowing_configs),
  cultivar = cultivars,
  stringsAsFactors = FALSE
)

# 合并站点信息
sim_configs <- merge(sim_configs, loc_data, by = "site")

# 添加日期信息（播种前30天到收获后30天）
for (i in 1:nrow(sim_configs)) {
  config <- sowing_configs[[sim_configs$sowing_code[i]]]
  
  # 播种日期：year年10月15日 + offset
  sowing_date <- as.Date(paste0(sim_configs$year[i], "-10-15")) + days(config$sowing_date_offset)
  
  # 收获日期：下一年，根据播种日期计算
  harvest_date <- sowing_date + days(config$harvest_offset)
  
  # 气象数据时间范围：播种前30天到收获后30天
  start_date <- sowing_date - days(30)
  end_date <- harvest_date + days(30)
  
  # 生长季标识（例如：2025-2026）
  growing_season <- paste0(sim_configs$year[i], "-", sim_configs$year[i] + 1)
  
  sim_configs$growing_season[i] <- growing_season
  sim_configs$sowing_date[i] <- as.character(sowing_date)
  sim_configs$harvest_date[i] <- as.character(harvest_date)
  sim_configs$start_date[i] <- as.character(start_date)
  sim_configs$end_date[i] <- as.character(end_date)
  
  # 设置种植密度
  sim_configs$plant_density[i] <- 200
  
  # 计算实际模拟的天数
  sim_configs$simulation_days[i] <- as.integer(end_date - start_date)
}

# 保存配置到未来模拟配置目录
config_file <- file.path(future_dirs$config, "future_simulation_config_2025_2049.csv")
write.csv(sim_configs, config_file, row.names = FALSE)

# 生成配置摘要
config_summary <- sim_configs %>%
  group_by(site, scenario, sowing_code) %>%
  summarise(
    growing_seasons = n(),
    first_season = min(growing_season),
    last_season = max(growing_season),
    avg_simulation_days = round(mean(simulation_days), 1),
    .groups = "drop"
  )

summary_file <- file.path(future_dirs$config, "simulation_config_summary.txt")
sink(summary_file)
cat("未来气候小麦模拟配置摘要\n")
cat("=", paste(rep("=", 50), collapse = ""), "\n\n")
cat("模拟时间范围: 2025-2026 至 2049-2050 生长季\n")
cat("生长季数量:", length(simulation_years), "\n")
cat("总仿真配置数:", nrow(sim_configs), "\n\n")
cat("各站点配置详情:\n")
print(config_summary)
cat("\n各播期设置:\n")
for (code in names(sowing_configs)) {
  cat(paste0("  ", code, " (", sowing_configs[[code]]$name, "):\n"))
  cat(paste0("    播种日期: 10月15日 + ", sowing_configs[[code]]$sowing_date_offset, "天\n"))
  cat(paste0("    收获日期: 播种后 ", sowing_configs[[code]]$harvest_offset, "天\n"))
}
sink()

message(paste("✅ 仿真配置保存：", config_file))
message(paste("✅ 配置摘要保存：", summary_file))
message(paste("📊 总仿真配置数：", nrow(sim_configs), "（", length(simulation_years), "个生长季）"))


########################## 8. 批量仿真主循环 ##########################
message("\n🚀 开始批量仿真...")

# 读取矫正后的气象数据
corrected_weather <- read.csv(corrected_file, stringsAsFactors = FALSE)
corrected_weather$time <- as.Date(corrected_weather$time)

# 创建仿真进度记录文件
progress_file <- file.path(future_dirs$config, "simulation_progress.csv")
progress_data <- data.frame(
  config_id = integer(),
  env_id = character(),
  status = character(),
  weather_file_created = logical(),
  simulation_completed = logical(),
  output_saved = logical(),
  error_message = character(),
  timestamp = character(),
  stringsAsFactors = FALSE
)
write.csv(progress_data, progress_file, row.names = FALSE)

# 测试模式设置（TRUE=测试，FALSE=完整运行）
test_mode <- FALSE
if (test_mode) {
  max_simulations <- min(2, nrow(sim_configs))
  message(paste("🔧 测试模式：只运行前", max_simulations, "个仿真"))
} else {
  max_simulations <- nrow(sim_configs)
  message(paste("🔧 完整运行模式：运行全部", max_simulations, "个仿真"))
}

# 成功和失败计数
success_count <- 0
fail_count <- 0

# 主仿真循环
for (i in 1:max_simulations) {
  config <- sim_configs[i, ]
  
  # 1. 设置环境ID和文件名
  env_id <- paste0(config$growing_season, config$site, "_", config$scenario, "_", config$sowing_code)
  
  # 修正的显示方式
  separator_line <- paste(rep("-", 60), collapse = "")
  message("\n", separator_line)
  message(paste("🌾 处理环境", i, "/", max_simulations, "：", env_id))
  message(paste("📅 生长季：", config$growing_season))
  message(paste("📍 站点：", config$site, "(lat:", config$latitude, "lon:", config$longitude, ")"))
  message(paste("🌡️ 情景：", config$scenario, "| 播期：", config$sowing_code, "(", sowing_configs[[config$sowing_code]]$name, ")"))
  
  # 2. 气象文件设置
  weather_file <- paste0(env_id, ".met")
  weather_path <- file.path(future_dirs$weather, weather_file)
  
  # 更新进度记录
  progress_update <- data.frame(
    config_id = i,
    env_id = env_id,
    status = "开始处理",
    weather_file_created = FALSE,
    simulation_completed = FALSE,
    output_saved = FALSE,
    error_message = "",
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  
  # 3. 检查是否已生成气象文件
  weather_created <- FALSE
  if (file.exists(weather_path)) {
    message(paste("✅ 气象文件已存在：", weather_file))
    weather_created <- TRUE
  } else {
    tryCatch({
      # 从矫正后的数据中提取该站点的气象数据
      site_scenario_data <- corrected_weather %>%
        filter(site == config$site, scenario == config$scenario) %>%
        arrange(time)
      
      if (nrow(site_scenario_data) == 0) {
        stop(paste("没有气象数据：站点", config$site, "情景", config$scenario))
      }
      
      # 筛选时间范围内的数据
      met_data <- site_scenario_data %>%
        filter(time >= as.Date(config$start_date) & time <= as.Date(config$end_date)) %>%
        arrange(time)
      
      if (nrow(met_data) == 0) {
        stop(paste("没有气象数据在范围", config$start_date, "至", config$end_date))
      }
      
      # 检查数据完整性
      required_columns <- c("tasmax_corrected", "tasmin_corrected", "pr_corrected", 
                           "rsds_corrected", "hurs_corrected", "sfcwind_corrected", "vapr_corrected")
      missing_cols <- setdiff(required_columns, names(met_data))
      if (length(missing_cols) > 0) {
        stop(paste("缺失气象变量：", paste(missing_cols, collapse = ", ")))
      }
      
      # 检查缺失值
      missing_count <- sum(is.na(met_data[, required_columns]))
      if (missing_count > 0) {
        warning(paste("⚠️ 气象数据有", missing_count, "个缺失值"))
      }
      
      # 构建APSIM气象数据框
      apsim_met <- data.frame(
        date = met_data$time,
        year = year(met_data$time),
        day = yday(met_data$time),
        radn = met_data$rsds_corrected,       # MJ/m²/day
        maxt = met_data$tasmax_corrected,     # °C
        mint = met_data$tasmin_corrected,     # °C
        rain = met_data$pr_corrected,         # mm/day
        rh = met_data$hurs_corrected,         # %
        wind = met_data$sfcwind_corrected,    # m/s
        vapr = met_data$vapr_corrected        # kPa
      )
      
      # 计算tav和amp（使用该生长季数据计算）
      tavg <- (apsim_met$maxt + apsim_met$mint) / 2
      tav <- mean(tavg, na.rm = TRUE)
      
      # 计算月平均温度
      monthly_avg <- apsim_met %>%
        mutate(month = month(date)) %>%
        group_by(month) %>%
        summarise(avg_temp = mean((maxt + mint)/2, na.rm = TRUE))
      
      amp <- (max(monthly_avg$avg_temp, na.rm = TRUE) - 
              min(monthly_avg$avg_temp, na.rm = TRUE)) / 2
      
      # 创建met文件头
      met_header <- c(
        "[weather.met.weather]",
        paste("Station =", env_id),
        paste("latitude =", config$latitude, "(DECIMAL DEGREES)"),
        paste("longitude =", config$longitude, "(DECIMAL DEGREES)"),
        paste("tav =", round(tav, 3), "(oC) ! annual average ambient temperature"),
        paste("amp =", round(amp, 3), "(oC) ! annual amplitude in mean monthly temperature"),
        "year day radn maxt mint rain rh wind",
        "() () (MJ/m^2) (oC) (oC) (mm) (%) (m/s)"
      )
      
      # 写入文件
      writeLines(met_header, weather_path)
      write.table(apsim_met[, c("year", "day", "radn", "maxt", "mint", "rain", "rh", "wind")],
                  file = weather_path, append = TRUE,
                  row.names = FALSE, col.names = FALSE, sep = " ")
      
      weather_created <- TRUE
      message(paste("✅ 生成气象文件：", weather_file))
      message(paste("📅 气象时间范围：", config$start_date, "至", config$end_date, "（共", nrow(apsim_met), "天）"))
      
    }, error = function(e) {
      message(paste("❌ 生成气象文件失败：", e$message))
      progress_update$error_message <- paste("气象文件错误：", e$message)
      progress_update$status <- "失败-气象文件"
      return(NULL)
    })
  }
  
  # 如果气象文件创建失败，跳过该配置
  if (!weather_created) {
    fail_count <- fail_count + 1
    progress_data <- rbind(progress_data, progress_update)
    write.csv(progress_data, progress_file, row.names = FALSE)
    next
  }
  
  # 4. 准备仿真脚本
  sim_name <- paste0(env_id, "_", config$cultivar)
  sim_file <- paste0(sim_name, ".apsimx")
  sim_script_path <- file.path(future_dirs$script, sim_file)
  source_script <- file.path(apsimx_examples_path, "Wheat.apsimx")
  
  if (!file.exists(sim_script_path)) {
    file.copy(source_script, sim_script_path, overwrite = TRUE)
    message(paste("✅ 复制仿真脚本：", sim_file))
  }
  
  # 5. 修改仿真参数
  tryCatch({
    # 5.1 修改Clock节点（仿真时间）
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Clock",
      parm = c("Start", "End"),
      value = c(config$start_date, config$end_date),
      overwrite = TRUE
    )
    
    # 5.2 修改Weather节点（气象文件路径）
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Weather",
      parm = "FileName",
      value = normalizePath(weather_path),
      overwrite = TRUE
    )
    
    # 5.3 查找播种/施肥脚本路径
    sow_script <- inspect_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      node = "Manager",
      parm = list("Sow using a variable rule", NA),
      print.path = FALSE
    )
    fert_script <- inspect_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      node = "Manager",
      parm = list("Fertilise at sowing", NA),
      print.path = FALSE
    )
    
    # 5.4 修改播种参数
    original_locale <- Sys.getlocale("LC_TIME")
    Sys.setlocale("LC_TIME", "C")
    plant_date_formatted <- format(as.Date(config$sowing_date), "%d-%b")
    
    # 播种日期
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".StartDate"),
      value = plant_date_formatted,
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".EndDate"),
      value = plant_date_formatted,
      overwrite = TRUE
    )
    
    # 移除播种限制
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".MinRain"),
      value = 0,
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".MinESW"),
      value = 0,
      overwrite = TRUE
    )
    
    # 种植密度
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".Population"),
      value = config$plant_density,
      overwrite = TRUE
    )
    
    # 5.5 修改施肥参数
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Other",
      parm.path = paste0(fert_script, ".Amount"),
      value = 200,
      overwrite = TRUE
    )
    Sys.setlocale("LC_TIME", original_locale)  # 恢复区域设置
    
    # 5.6 修改Report节点（输出指标）
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Report",
      parm = "EventNames",
      value = c("[Clock].EndOfDay", "[Wheat].Harvesting"),
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      wrt.dir = future_dirs$script,
      node = "Report",
      parm = "VariableNames",
      value = c(
        "[Clock].Today", "[Wheat].LAI", "[Wheat].Phenology.Zadok.Stage",
        "[Wheat].Phenology.CurrentStageName", "[Wheat].AboveGround.Wt",
        "[Wheat].AboveGround.N", "[Wheat].Grain.Total.Wt * 10 as Yield",
        "[Wheat].Grain.Protein", "[Wheat].Grain.Size", "[Wheat].Grain.Number",
        "[Wheat].Grain.Total.Wt", "[Wheat].Grain.Total.N", "[Wheat].Total.Wt",
        "[Clock].Today", "[Wheat].Phenology.Stage", "[SurfaceOrganicMatter].Wt",
        "[Soil].Water.PAWmm", "[Soil].Water.PAW", "[Physical].PAWCmm", "[Physical].PAWC",
        "[ISoilWater].Eos", "[ISoilWater].Es", "[ISoilWater].Eo",
        "[ISoilWater].Runoff", "[ISoilWater].Infiltration", "[ISoilWater].Drainage",
        "[ISoilWater].LeachNO3", "[ISoilWater].LeachNH4", "[ISoilWater].LeachCl",
        "[NFlow].Natm", "[NFlow].N2Oatm", "[Hydrolysis].Value",
        "[Nutrient].MineralisedN", "[Nutrient].Catm"
      ),
      overwrite = TRUE
    )
    
    message("✅ 仿真脚本配置完成")
    
  }, error = function(e) {
    message(paste("❌ 仿真脚本配置失败：", e$message))
    progress_update$error_message <- paste("脚本配置错误：", e$message)
    progress_update$status <- "失败-脚本配置"
    fail_count <- fail_count + 1
    return(NULL)
  })
  
  # 6. 运行仿真
  message("▶️ 运行仿真...")
  result <- tryCatch({
    apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      value = "report",
      cleanup = FALSE,
      simplify = TRUE
    )
  }, error = function(e) {
    message(paste("❌ 仿真运行失败：", e$message))
    progress_update$error_message <- paste("仿真运行错误：", e$message)
    progress_update$status <- "失败-仿真运行"
    return(NULL)
  })
  
  # 7. 保存结果
  if (!is.null(result)) {
    output_file <- paste0(env_id, "_simulation_output.csv")
    output_path <- file.path(future_dirs$output, output_file)
    write.csv(result, output_path, row.names = FALSE)
    
    # 检查结果文件
    if (file.exists(output_path)) {
      file_info <- file.info(output_path)
      if (file_info$size > 1024) {  # 文件大小大于1KB
        message(paste("✅ 结果保存成功：", output_file, "（", round(file_info$size/1024, 1), "KB）"))
        progress_update$output_saved <- TRUE
        progress_update$simulation_completed <- TRUE
        progress_update$status <- "成功"
        success_count <- success_count + 1
      } else {
        message(paste("⚠️ 结果文件过小：", output_file))
        progress_update$error_message <- "结果文件过小"
        progress_update$status <- "警告-结果文件"
        fail_count <- fail_count + 1
      }
    } else {
      message(paste("❌ 结果文件未创建：", output_file))
      progress_update$error_message <- "结果文件未创建"
      progress_update$status <- "失败-文件保存"
      fail_count <- fail_count + 1
    }
  } else {
    fail_count <- fail_count + 1
  }
  
  # 更新进度记录
  progress_update$weather_file_created <- weather_created
  if (exists("progress_update$simulation_completed") && !progress_update$simulation_completed) {
    progress_update$simulation_completed <- FALSE
  }
  if (exists("progress_update$output_saved") && !progress_update$output_saved) {
    progress_update$output_saved <- FALSE
  }
  progress_update$timestamp <- as.character(Sys.time())
  
  # 保存进度
  progress_data <- rbind(progress_data, progress_update)
  write.csv(progress_data, progress_file, row.names = FALSE)
  
  # 显示进度
  progress_pct <- round(i / max_simulations * 100, 1)
  message(paste("📊 进度：", i, "/", max_simulations, "（", progress_pct, "%）| 成功：", success_count, "失败：", fail_count))
  
  # 避免过载，短暂延迟
  Sys.sleep(0.5)
}

# 修正这里的显示方式
separator_line_final <- paste(rep("=", 60), collapse = "")
message("\n", separator_line_final)


########################## 9. 生成仿真总结报告 ##########################
message("\n📊 生成仿真总结报告...")

# 读取进度数据
if (file.exists(progress_file)) {
  progress_summary <- read.csv(progress_file, stringsAsFactors = FALSE)
  
  # 统计状态
  status_summary <- progress_summary %>%
    group_by(status) %>%
    summarise(count = n(), .groups = "drop")
  
  # 按站点统计
  site_summary <- progress_summary %>%
    mutate(site = substr(env_id, 10, 11)) %>%  # 从env_id中提取站点代码
    group_by(site) %>%
    summarise(
      total = n(),
      success = sum(status == "成功"),
      success_rate = round(success / total * 100, 1),
      .groups = "drop"
    )
  
  # 按情景统计
  scenario_summary <- progress_summary %>%
    mutate(scenario = ifelse(grepl("ssp126", env_id), "ssp126", "ssp585")) %>%
    group_by(scenario) %>%
    summarise(
      total = n(),
      success = sum(status == "成功"),
      success_rate = round(success / total * 100, 1),
      .groups = "drop"
    )
} else {
  status_summary <- data.frame(status = "无进度数据", count = 0)
  site_summary <- data.frame(site = "无数据", total = 0, success = 0, success_rate = 0)
  scenario_summary <- data.frame(scenario = "无数据", total = 0, success = 0, success_rate = 0)
}

# 生成总结报告
summary_report <- file.path(future_dirs$config, "future_simulation_final_summary.txt")
sink(summary_report)

cat("未来气候小麦模拟完成总结报告\n")
cat("=", paste(rep("=", 50), collapse = ""), "\n\n")
cat("报告生成时间:", as.character(Sys.time()), "\n")
cat("工作目录:", work_dir, "\n\n")

cat("1. 模拟基本信息\n")
cat("   - 模拟类型: 未来气候情景\n")
cat("   - 时间范围: 2025-2026 至 2049-2050 生长季\n")
cat("   - 生长季数量:", length(simulation_years), "\n")
cat("   - 气候情景:", paste(scenarios, collapse = ", "), "\n")
cat("   - 站点数量:", nrow(loc_data), "\n")
cat("   - 播期设置:", paste(names(sowing_configs), collapse = ", "), "\n")
cat("   - 品种:", paste(cultivars, collapse = ", "), "\n")
cat("   - 总仿真配置数:", nrow(sim_configs), "\n")
cat("   - 实际运行数:", max_simulations, "\n\n")

cat("2. 运行结果统计\n")
cat("   - 成功仿真数:", success_count, "\n")
cat("   - 失败仿真数:", fail_count, "\n")
cat("   - 成功率:", ifelse(max_simulations > 0, round(success_count / max_simulations * 100, 1), 0), "%\n\n")

cat("3. 状态分布\n")
for (i in 1:nrow(status_summary)) {
  cat("   - ", status_summary$status[i], ": ", status_summary$count[i], "\n")
}
cat("\n")

cat("4. 按站点统计\n")
for (i in 1:nrow(site_summary)) {
  cat("   - ", site_summary$site[i], ": ", site_summary$success[i], "/", site_summary$total[i], 
      " (", site_summary$success_rate[i], "%)\n")
}
cat("\n")

cat("5. 按情景统计\n")
for (i in 1:nrow(scenario_summary)) {
  cat("   - ", scenario_summary$scenario[i], ": ", scenario_summary$success[i], "/", scenario_summary$total[i], 
      " (", scenario_summary$success_rate[i], "%)\n")
}
cat("\n")

cat("6. 目录结构\n")
cat("   - 气象文件目录:", future_dirs$weather, "\n")
cat("   - 仿真脚本目录:", future_dirs$script, "\n")
cat("   - 输出结果目录:", future_dirs$output, "\n")
cat("   - 配置文件目录:", future_dirs$config, "\n\n")

cat("7. 关键参数设置\n")
cat("   - 种植密度: 200 plants/m²\n")
cat("   - 施肥量: 200 kg/ha\n")
cat("   - 气象数据: CMIP6矫正数据 (2001-2050)\n")
cat("   - 气象时间范围: 播种前30天至收获后30天\n")
cat("   - 输出变量: 37个指标\n\n")

cat("8. 文件统计\n")
for (dir_name in names(future_dirs)) {
  dir_path <- future_dirs[[dir_name]]
  if (dir.exists(dir_path)) {
    file_count <- length(list.files(dir_path, pattern = "\\.(met|apsimx|csv|txt)$"))
    cat("   - ", dir_name, "目录文件数: ", file_count, "\n")
  }
}

sink()

message(paste("✅ 总结报告保存：", summary_report))


########################## 10. 最终输出 ##########################
# 创建最终分隔线
final_separator <- paste(rep("=", 80), collapse = "")
cat("\n", final_separator, "\n")
cat("🎉 未来气候小麦模拟任务完成！\n\n")
cat("📈 模拟统计:\n")
cat("   - 配置总数:", nrow(sim_configs), "\n")
cat("   - 实际运行:", max_simulations, "\n")
cat("   - 成功:", success_count, " (", round(success_count/max_simulations*100, 1), "%)\n")
cat("   - 失败:", fail_count, " (", round(fail_count/max_simulations*100, 1), "%)\n\n")

cat("📁 输出文件位置:\n")
for (dir_name in names(future_dirs)) {
  cat(paste0("   - ", dir_name, ": ", future_dirs[[dir_name]], "\n"))
}

cat("\n🔍 下一步建议:\n")
cat("   1. 检查 summary_report.txt 了解详细结果\n")
cat("   2. 查看 simulation_progress.csv 了解每个仿真的状态\n")
cat("   3. 分析 Future_APSIMout 目录中的结果文件\n")
cat("   4. 如有失败案例，检查错误信息进行调试\n")

cat("\n", final_separator, "\n")