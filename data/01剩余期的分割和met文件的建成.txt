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

# 清除环境残留与控制台输出
rm(list = ls(all = TRUE))
cat("\014")
message("✅ 前置环境清理完成")


########################## 2. 配置工作目录与APSIM路径 ##########################
work_dir <- "G:\\株高信息\\基础数据分析\\02APSIM\\TrainandTest\\planting"
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


########################## 3. 读取种植数据 ##########################
planting_file <- "plantingtrain.csv"
if (!file.exists(planting_file)) {
  stop(paste("❌ 种植数据缺失：", file.path(getwd(), planting_file)))
}
planting_info <- read.csv(planting_file, stringsAsFactors = FALSE)
message(paste("✅ 读取种植数据：", nrow(planting_info), "个环境"))
cultivars <- c("Keyu13")  # 模拟品种


########################## 4. 气象数据下载函数 ##########################
download_and_generate_weather <- function(lat, lon, start_date, end_date, weather_dir, weather_file) {
  max_attempts <- 5
  attempt <- 1
  start_time <- Sys.time()
  
  while (attempt <= max_attempts) {
    time_elapsed <- difftime(Sys.time(), start_time, units = "hours")
    if (as.numeric(time_elapsed) > 1) stop("❌ 气象下载超时（>1小时）")
    
    tryCatch({
      pwr <- get_power_apsim_met(lonlat = c(lon, lat), dates = c(start_date, end_date))
      write_apsim_met(pwr, wrt.dir = weather_dir, filename = weather_file)
      message(paste("✅ 气象文件保存：", file.path(weather_dir, weather_file)))
      break
    }, error = function(e) {
      if (attempt < max_attempts) {
        message(paste("❌ 下载失败（尝试", attempt, "/5）：", e$message, "→5秒后重试"))
        Sys.sleep(5)
      } else {
        stop(paste("❌ 下载失败（5次尝试）：", e$message))
      }
      attempt <- attempt + 1
    })
  }
}


########################## 5. 创建核心目录 ##########################
dirs <- list(
  weather = file.path(getwd(), "01APSIMmet"),
  script = file.path(getwd(), "01APSIMinput_apsimx"),
  output = file.path(getwd(), "01APSIMout")
)

for (dir_name in names(dirs)) {
  if (!dir.exists(dirs[[dir_name]])) {
    dir.create(dirs[[dir_name]], recursive = TRUE)
    message(paste("✅ 新建目录：", dirs[[dir_name]]))
  } else {
    message(paste("✅ 目录已存在：", dirs[[dir_name]]))
  }
}


########################## 6. 主循环：批量仿真 ##########################
for (job in 1:18) {
  for (cultivar in cultivars) {
    dat <- planting_info[job, ]
    env_id <- dat$year_loc
    message(paste("\n=== 处理环境", job, "/18：", env_id, "==="))
    
    # 1. 处理日期
    plant_date <- as.Date(dat$sowing_date, format = "%Y-%m-%d")
    harvest_date <- as.Date(dat$end_date, format = "%Y-%m-%d")
    if (is.na(harvest_date)) {
      harvest_date <- plant_date + 270
      message(paste("⚠️ 收获日期默认设为：", harvest_date))
    }
    
    # 2. 气象文件设置
    weather_file <- paste0(env_id, ".met")
    weather_path <- file.path(dirs$weather, weather_file)
    
    # 3. 气象时间范围
    start_date <- plant_date - 30
    end_date <- harvest_date + 30
    message(paste("📅 气象时间范围：", start_date, "至", end_date))
    
    # 4. 下载气象数据
    if (file.exists(weather_path)) {
      file.remove(weather_path)
      message(paste("ℹ️ 删除旧气象文件：", weather_file))
    }
    download_and_generate_weather(
      lat = dat$latitude,
      lon = dat$longitude,
      start_date = start_date,
      end_date = end_date,
      weather_dir = dirs$weather,
      weather_file = weather_file
    )
    
    # 5. 准备仿真脚本
    sim_name <- paste0(env_id, "_", cultivar)
    sim_file <- paste0(sim_name, ".apsimx")
    sim_script_path <- file.path(dirs$script, sim_file)
    source_script <- file.path(apsimx_examples_path, "Wheat.apsimx")
    
    if (!file.exists(sim_script_path)) {
      file.copy(source_script, sim_script_path, overwrite = TRUE)
      message(paste("✅ 复制仿真脚本：", sim_file))
    }
    
    # 6. 修改Clock节点（仿真时间）
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Clock",
      parm = c("Start", "End"),
      value = c(as.character(start_date), as.character(end_date)),
      overwrite = TRUE
    )
    
    # 7. 修改Weather节点（气象文件路径）
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Weather",
      parm = "FileName",
      value = normalizePath(weather_path),
      overwrite = TRUE
    )
    
    # 8. 查找播种/施肥脚本路径
    sow_script <- inspect_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      node = "Manager",
      parm = list("Sow using a variable rule", NA),
      print.path = FALSE
    )
    fert_script <- inspect_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      node = "Manager",
      parm = list("Fertilise at sowing", NA),
      print.path = FALSE
    )
    
    # 9. 修改播种参数
    original_locale <- Sys.getlocale("LC_TIME")
    Sys.setlocale("LC_TIME", "C")
    plant_date_formatted <- format(plant_date, "%d-%b")
    
    # 播种日期
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".StartDate"),
      value = plant_date_formatted,
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".EndDate"),
      value = plant_date_formatted,
      overwrite = TRUE
    )
    
    # 移除播种限制
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".MinRain"),
      value = 0,
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".MinESW"),
      value = 0,
      overwrite = TRUE
    )
    
    # 种植密度
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(sow_script, ".Population"),
      value = dat$plant_density,
      overwrite = TRUE
    )
    
    # 10. 修改施肥参数
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Other",
      parm.path = paste0(fert_script, ".Amount"),
      value = 200,
      overwrite = TRUE
    )
    Sys.setlocale("LC_TIME", original_locale)  # 恢复区域设置
    
    # 11. 修改Report节点（输出指标）
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Report",
      parm = "EventNames",
      value = c("[Clock].EndOfDay", "[Wheat].Harvesting"),
      overwrite = TRUE
    )
    edit_apsimx(
      file = sim_file,
      src.dir = dirs$script,
      wrt.dir = dirs$script,
      node = "Report",
      parm = "VariableNames",
      value = c(
        "[Clock].Today", "[Wheat].LAI", "[Wheat].Phenology.Zadok.Stage",
        "[Wheat].Phenology.CurrentStageName", "[Wheat].AboveGround.Wt",
        "[Wheat].AboveGround.N", "[Wheat].Grain.Total.Wt * 10 as Yield",  # 产量单位转换（t/ha）
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
    
    # 12. 运行仿真
    message(paste("▶️ 运行仿真：", sim_name))
    result <- tryCatch({
      apsimx(
        file = sim_file,
        src.dir = dirs$script,
        value = "report",
        cleanup = FALSE,
        simplify = TRUE
      )
    }, error = function(e) {
      stop(paste("❌ 仿真失败（环境", env_id, "）：", e$message))
    })
    
    # 13. 保存结果
    output_file <- paste0(env_id, "_simulation_output.csv")
    output_path <- file.path(dirs$output, output_file)
    write.csv(result, output_path, row.names = FALSE)
    message(paste("✅ 结果保存：", output_file))
    message(paste("=== 环境", job, "/18 处理完成 ==="))
  }
}


########################## 7. 运行完成提示 ##########################
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("✅ 所有18个环境仿真完成！\n")
cat("📂 文件位置：\n")
cat("- 气象文件：", dirs$weather, "\n")
cat("- 仿真脚本：", dirs$script, "\n")
cat("- 结果文件：", dirs$output, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
