# 文件名：run_apsim_future_linux_FULL_FINAL_FIX_v5_1_1_JSON_ALL_TT_FALLBACK_GDD.R
# 位置：/home/wangyanan/APSIM_Wheat_Project/
# 描述：未来气候小麦 APSIM 批量模拟（Linux物理机版，v5.1.1 完整版，不省略）
# v5.1.1 在 v5.1 基础上新增：
#   ✅ CLI 捕获日志（stdout+stderr），自动解析并剔除 Report 无效变量后重试（最多3次）
#   ✅ 仍保留 v5.1 所有修复：JSON 全强写 + TT 回退 + 无 TT 则 GDD(base0/base5) 气象兜底
#   ✅ 继续修复 Weather.ExcelWorkSheetName 被写坏成 {}：强制写回 null
#   ✅ wrapper 自动补 run：避免 Linux 上 CLI 调用兼容问题

###############################################################################
########################## 0. 初始化环境 ######################################
###############################################################################
cat("\014")
message("🎯 APSIM未来气候模拟 - Linux物理机版（FULL FINAL FIX v5.1.1，JSON全强写 + CLI智能剔除无效变量）")
message("====================================================================================")
message(paste("开始时间:", Sys.time()))
message(paste("R版本:", R.version$version.string))
message(paste("系统平台:", R.version$platform))

###############################################################################
########################## 1. 安装并加载必要包 #################################
###############################################################################
safe_require <- function(pkg, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/", quiet = TRUE) {
  if (!require(pkg, character.only = TRUE, quietly = quiet)) {
    install.packages(pkg, dependencies = TRUE, repos = repo)
    library(pkg, character.only = TRUE)
  }
  message(paste("✅", pkg, "加载成功"))
}

message("\n📦 加载R包...")
safe_require("apsimx")
safe_require("lubridate")
safe_require("dplyr")
safe_require("tidyr")
safe_require("jsonlite")
suppressWarnings(suppressMessages(library(parallel)))
message("✅ parallel 加载成功（基础包）")

###############################################################################
########################## 2. 配置工作目录与APSIM路径（Linux）###################
###############################################################################
message("\n⚙️ 配置工作目录与APSIM路径（Linux）...")

work_dir <- "/home/wangyanan/APSIM_Wheat_Project"
if (!dir.exists(work_dir)) stop(paste("❌ 工作目录不存在：", work_dir))
setwd(work_dir)
message(paste("✅ 工作目录：", work_dir))

apsimx_exe_path <- "/software/apsim_bin/ApsimX/bin/Release/net8.0/apsim_wrapper.sh"
if (!file.exists(apsimx_exe_path)) stop(paste("❌ APSIM wrapper 不存在：", apsimx_exe_path))
try(Sys.chmod(apsimx_exe_path, mode = "0755"), silent = TRUE)

apsimx_examples_path <- system.file("extdata", package = "apsimx")
if (apsimx_examples_path == "" || !dir.exists(apsimx_examples_path)) {
  apsimx_examples_path <- "/software/apsim_bin/ApsimX/Examples"
}
if (!dir.exists(apsimx_examples_path)) stop(paste("❌ 示例目录不存在：", apsimx_examples_path))
if (!file.exists(file.path(apsimx_examples_path, "Wheat.apsimx"))) {
  stop(paste("❌ 小麦模板缺失：", file.path(apsimx_examples_path, "Wheat.apsimx")))
}

# wrapper 自动补 run：避免 “Verb 'xxx.apsimx' is not recognized.”
read_wrapper_text <- function(path) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  paste(txt, collapse = "\n")
}
wrapper_text <- read_wrapper_text(apsimx_exe_path)
need_fix <- TRUE
if (nchar(wrapper_text) > 0) {
  if (grepl("apsim\\.dll\\s+run", wrapper_text) || grepl("run\\s+\"\\$@\"", wrapper_text)) need_fix <- FALSE
}
if (need_fix) {
  message("⚠️ wrapper 可能不兼容 apsimx 调用方式，写入兼容版（自动补 run）...")
  fixed_wrapper <- c(
    "#!/bin/bash",
    "set -e",
    "cd /software/apsim_bin/ApsimX",
    "if [[ \"$1\" == *.apsimx ]]; then",
    "  dotnet bin/Release/net8.0/apsim.dll run \"$@\"",
    "else",
    "  dotnet bin/Release/net8.0/apsim.dll \"$@\"",
    "fi"
  )
  writeLines(fixed_wrapper, apsimx_exe_path)
  Sys.chmod(apsimx_exe_path, mode = "0755")
  message("✅ wrapper 已修复并赋权")
} else {
  message("✅ wrapper 已包含 run 兼容逻辑")
}

apsimx::apsimx_options(
  exe.path = apsimx_exe_path,
  examples.path = apsimx_examples_path,
  warn.find.apsimx = FALSE
)
message("✅ apsimx_options 配置完成")
message(paste("✅ APSIM wrapper:", apsimx_exe_path))
message(paste("✅ APSIM examples:", apsimx_examples_path))

###############################################################################
########################## 3. 站点信息 ########################################
###############################################################################
message("\n📍 读取站点信息...")

loc_data <- data.frame(
  site = c("NY", "PY", "YL", "ZMD"),
  latitude = c(32.9126, 35.7627, 34.2879, 33.014),
  longitude = c(112.4626, 115.0292, 108.0001, 114.0219),
  stringsAsFactors = FALSE
)
message("✅ 站点数据：")
print(loc_data)

###############################################################################
########################## 4. 创建未来模拟目录（清空旧文件）#####################
###############################################################################
message("\n📁 创建未来模拟目录结构（清空旧文件）...")

future_dirs <- list(
  weather = file.path(getwd(), "Future_APSIMmet"),
  script  = file.path(getwd(), "Future_APSIMinput_apsimx"),
  output  = file.path(getwd(), "Future_APSIMout"),
  config  = file.path(getwd(), "Future_APSIMconfig"),
  logs    = file.path(getwd(), "logs")
)

for (dir_name in names(future_dirs)) {
  dir_path <- future_dirs[[dir_name]]
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    message(paste("✅ 新建目录：", dir_path))
  } else {
    files_in_dir <- list.files(dir_path, all.files = TRUE, no.. = TRUE)
    if (length(files_in_dir) > 0) {
      unlink(list.files(dir_path, full.names = TRUE), recursive = TRUE, force = TRUE)
      message(paste("⚠️ 清空目录：", dir_path))
    }
    message(paste("✅ 目录已存在：", dir_path))
  }
}

###############################################################################
########################## 5. 播种日期配置 #####################################
###############################################################################
message("\n🌱 定义播种日期配置...")

sowing_configs <- list(
  Z = list(name = "早播", sowing_date_offset = 0,  harvest_offset = 228),
  M = list(name = "中播", sowing_date_offset = 15, harvest_offset = 227),
  W = list(name = "晚播", sowing_date_offset = 30, harvest_offset = 227)
)
message("✅ 播期配置：")
print(sowing_configs)

###############################################################################
########################## 6. 气象数据检查 #####################################
###############################################################################
message("\n🌤️ 检查矫正气象数据...")

corrected_file <- file.path(work_dir, "all_sites_ssp126_ssp585_2001_2050_corrected.csv")
if (!file.exists(corrected_file)) stop(paste("❌ 矫正后的气象数据缺失：", corrected_file))
message(paste("✅ 检测到矫正后的气象数据：", corrected_file))

corrected_weather <- read.csv(corrected_file, stringsAsFactors = FALSE)
corrected_weather$time <- as.Date(corrected_weather$time)
message(paste("📊 气象数据维度：", nrow(corrected_weather), "行 ×", ncol(corrected_weather), "列"))

###############################################################################
########################## 7. 生成仿真配置（2025-2049）##########################
###############################################################################
message("\n🔧 创建未来气候仿真配置...")

simulation_years <- 2025:2049
scenarios <- c("ssp126", "ssp585")
cultivars <- c("Keyu13")

sim_configs <- expand.grid(
  site = loc_data$site,
  year = simulation_years,
  scenario = scenarios,
  sowing_code = names(sowing_configs),
  cultivar = cultivars,
  stringsAsFactors = FALSE
)

sim_configs <- merge(sim_configs, loc_data, by = "site")

for (i in 1:nrow(sim_configs)) {
  cfg <- sowing_configs[[sim_configs$sowing_code[i]]]
  sowing_date <- as.Date(paste0(sim_configs$year[i], "-10-15")) + lubridate::days(cfg$sowing_date_offset)
  harvest_date <- sowing_date + lubridate::days(cfg$harvest_offset)
  start_date <- sowing_date - lubridate::days(30)
  end_date <- harvest_date + lubridate::days(30)
  growing_season <- paste0(sim_configs$year[i], "-", sim_configs$year[i] + 1)

  sim_configs$growing_season[i] <- growing_season
  sim_configs$sowing_date[i] <- as.character(sowing_date)
  sim_configs$harvest_date[i] <- as.character(harvest_date)
  sim_configs$start_date[i] <- as.character(start_date)
  sim_configs$end_date[i] <- as.character(end_date)
  sim_configs$plant_density[i] <- 200
  sim_configs$simulation_days[i] <- as.integer(end_date - start_date)
}

config_file <- file.path(future_dirs$config, "future_simulation_config_2025_2049.csv")
write.csv(sim_configs, config_file, row.names = FALSE)
message(paste("✅ 仿真配置保存：", config_file))
message(paste("📊 总仿真配置数：", nrow(sim_configs)))

###############################################################################
########################## 8. 进度文件 #########################################
###############################################################################
progress_file <- file.path(future_dirs$config, "simulation_progress.csv")
progress_data <- data.frame(
  config_id = integer(),
  env_id = character(),
  status = character(),
  weather_file_created = logical(),
  script_configured = logical(),
  simulation_completed = logical(),
  output_saved = logical(),
  used_tt = logical(),
  error_message = character(),
  timestamp = character(),
  stringsAsFactors = FALSE
)
write.csv(progress_data, progress_file, row.names = FALSE)
message(paste("✅ 进度文件初始化：", progress_file))

###############################################################################
########################## 9. 工具函数（JSON强写 + TT回退 + GDD兜底）##############
###############################################################################

# APSIM d-mmm（英文月份小写）
format_apsim_d_mmm <- function(date_obj) {
  original_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C")
  out <- format(as.Date(date_obj), "%d-%b")
  Sys.setlocale("LC_TIME", original_locale)
  tolower(out)
}

# GDD 兜底
calc_gdd <- function(df_met, base_temp = 0) {
  tmean <- (df_met$maxt + df_met$mint) / 2
  gdd_daily <- pmax(0, tmean - base_temp)
  gdd_cum <- cumsum(gdd_daily)
  data.frame(date = as.Date(df_met$date), gdd_daily = gdd_daily, gdd_cum = gdd_cum, stringsAsFactors = FALSE)
}

# ---- JSON 读写基础 ----
read_apsimx_json <- function(sim_path) {
  txt <- readLines(sim_path, warn = FALSE, encoding = "UTF-8")
  jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = FALSE)
}

write_apsimx_json <- function(js, sim_path) {
  out <- jsonlite::toJSON(js, auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines(out, sim_path, useBytes = TRUE)
  TRUE
}

find_child_idx_by_name <- function(children, name) {
  if (is.null(children) || length(children) == 0) return(integer(0))
  which(vapply(children, function(x) !is.null(x$Name) && identical(x$Name, name), logical(1)))
}

get_simulation_node <- function(js) {
  sim_idx <- find_child_idx_by_name(js$Children, "Simulation")
  if (length(sim_idx) == 0) stop("找不到顶层 Simulation 节点（Name='Simulation'）")
  list(idx = sim_idx[1], node = js$Children[[sim_idx[1]]])
}

get_zone_node <- function(sim, zone_name = "Field") {
  zone_idx <- find_child_idx_by_name(sim$Children, zone_name)
  if (length(zone_idx) == 0) stop(paste0("找不到 Zone 节点：", zone_name))
  list(idx = zone_idx[1], node = sim$Children[[zone_idx[1]]])
}

# 修复：Weather.ExcelWorkSheetName 必须是 null 或 string；强制写回 null
fix_weather_excelworksheet_null <- function(sim) {
  w_idx <- find_child_idx_by_name(sim$Children, "Weather")
  if (length(w_idx) == 0) return(sim)
  w <- sim$Children[[w_idx[1]]]
  w$ExcelWorkSheetName <- NULL
  sim$Children[[w_idx[1]]] <- w
  sim
}

# 通用：在 Simulation.Children 里找一个节点（Clock/Weather/SummaryFile 等）并强写字段
edit_simnode_fields_json <- function(sim_path, node_name, fields_named_list) {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js)
  sim <- sim_info$node

  sim <- fix_weather_excelworksheet_null(sim)

  n_idx <- find_child_idx_by_name(sim$Children, node_name)
  if (length(n_idx) == 0) stop(paste0("找不到 Simulation 子节点：", node_name))
  n <- sim$Children[[n_idx[1]]]

  for (k in names(fields_named_list)) {
    n[[k]] <- fields_named_list[[k]]
  }

  sim$Children[[n_idx[1]]] <- n
  js$Children[[sim_info$idx]] <- sim
  write_apsimx_json(js, sim_path)
  TRUE
}

strict_edit_simnode_fields_json <- function(sim_path, node_name, fields_named_list) {
  ok <- tryCatch({
    edit_simnode_fields_json(sim_path, node_name, fields_named_list)
    TRUE
  }, error = function(e) {
    message(paste0("❌ JSON 强写节点失败：", node_name, " | ", e$message))
    FALSE
  })
  if (!ok) stop(paste0("严格写入失败：JSON 强写节点 ", node_name, " 失败"))
  TRUE
}

# ---- JSON 强写 Manager.Parameters(Key/Value) ----
edit_manager_param_json <- function(sim_path, zone_name = "Field", manager_name, key, value) {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js)
  sim <- sim_info$node

  sim <- fix_weather_excelworksheet_null(sim)

  zone_info <- get_zone_node(sim, zone_name = zone_name)
  zone <- zone_info$node

  mgr_idx <- find_child_idx_by_name(zone$Children, manager_name)
  if (length(mgr_idx) == 0) stop(paste0("找不到 Manager: ", manager_name))
  mgr <- zone$Children[[mgr_idx[1]]]

  if (is.null(mgr$Parameters)) mgr$Parameters <- list()

  hit <- which(vapply(mgr$Parameters, function(p) !is.null(p$Key) && identical(p$Key, key), logical(1)))
  if (length(hit) == 0) {
    mgr$Parameters[[length(mgr$Parameters) + 1]] <- list(Key = key, Value = as.character(value))
  } else {
    mgr$Parameters[[hit[1]]]$Value <- as.character(value)
  }

  zone$Children[[mgr_idx[1]]] <- mgr
  sim$Children[[zone_info$idx]] <- zone
  js$Children[[sim_info$idx]] <- sim
  write_apsimx_json(js, sim_path)
  TRUE
}

strict_edit_manager_json <- function(sim_path, manager_name, key, value, zone_name = "Field") {
  ok <- tryCatch({
    edit_manager_param_json(sim_path, zone_name = zone_name, manager_name = manager_name, key = key, value = value)
    TRUE
  }, error = function(e) {
    message(paste0("❌ JSON 强写 Manager 参数失败：", manager_name, ".", key, " | ", e$message))
    FALSE
  })
  if (!ok) stop(paste0("严格写入失败：JSON 强写 ", zone_name, "[", manager_name, "].", key, " 失败"))
  TRUE
}

# ---- JSON 强写 Report.EventNames / VariableNames ----
get_report_node <- function(zone, report_name = "Report") {
  rpt_idx <- find_child_idx_by_name(zone$Children, report_name)
  if (length(rpt_idx) == 0) {
    rpt_idx <- which(vapply(zone$Children, function(x) {
      !is.null(x$`$type`) && grepl("Models\\.Report", x$`$type`)
    }, logical(1)))
  }
  if (length(rpt_idx) == 0) stop("找不到 Report 节点（Name='Report' 或 $type=Models.Report）")
  list(idx = rpt_idx[1], node = zone$Children[[rpt_idx[1]]])
}

edit_report_json <- function(sim_path,
                             zone_name = "Field",
                             report_name = "Report",
                             event_names = NULL,
                             variable_names = NULL) {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js)
  sim <- sim_info$node

  sim <- fix_weather_excelworksheet_null(sim)

  zone_info <- get_zone_node(sim, zone_name = zone_name)
  zone <- zone_info$node

  rpt_info <- get_report_node(zone, report_name = report_name)
  rpt <- rpt_info$node

  if (!is.null(event_names)) {
    if (!is.character(event_names)) event_names <- as.character(event_names)
    rpt$EventNames <- unname(event_names)
  }
  if (!is.null(variable_names)) {
    if (!is.character(variable_names)) variable_names <- as.character(variable_names)
    rpt$VariableNames <- unname(variable_names)
  }

  zone$Children[[rpt_info$idx]] <- rpt
  sim$Children[[zone_info$idx]] <- zone
  js$Children[[sim_info$idx]] <- sim
  write_apsimx_json(js, sim_path)
  TRUE
}

strict_edit_report_json <- function(sim_path,
                                   zone_name = "Field",
                                   report_name = "Report",
                                   event_names = NULL,
                                   variable_names = NULL) {
  ok <- tryCatch({
    edit_report_json(sim_path,
                     zone_name = zone_name,
                     report_name = report_name,
                     event_names = event_names,
                     variable_names = variable_names)
    TRUE
  }, error = function(e) {
    message(paste0("❌ JSON 强写 Report 失败：", e$message))
    FALSE
  })
  if (!ok) stop("严格写入失败：JSON 强写 Report.EventNames / VariableNames 失败")
  TRUE
}

# TT 候选变量（v5.1.1 建议精简到最小集，最稳）
get_tt_candidates <- function() {
  c("[Wheat].Phenology.ThermalTime")
}

# 运行 apsimx 并抓错（失败返回 NULL + 错误文本）
run_apsim_once <- function(sim_file, src_dir) {
  err_txt <- ""
  res <- tryCatch({
    apsimx::apsimx(
      file = sim_file,
      src.dir = src_dir,
      value = "report",
      cleanup = FALSE,
      simplify = TRUE
    )
  }, error = function(e) {
    err_txt <<- e$message
    NULL
  })
  list(result = res, error = err_txt)
}

# 是否结果已经包含 TT/积温字段（宽松判断）
has_tt_in_result <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(FALSE)
  any(grepl("ThermalTime|AccumulatedThermalTime|CumulativeThermalTime|\\bTT\\b|DeltaTT", names(df), ignore.case = TRUE))
}

###############################################################################
# ========== >>>>> v5.1.1 新增：APSIM CLI 日志智能剔除无效变量 <<<<< ==========
###############################################################################

# ① 用 CLI 运行并完整捕获日志（stdout + stderr）
run_apsim_cli_capture <- function(apsim_exe, apsimx_path) {
  tryCatch({
    system2(
      apsim_exe,
      args = c("run", apsimx_path),
      stdout = TRUE, stderr = TRUE
    )
  }, error = function(e) {
    c(paste("SYSTEM2_ERROR:", e$message))
  })
}

# ② 从日志中解析“Invalid report variables found:”后面的变量名
extract_invalid_report_vars <- function(cli_lines) {
  k <- grep("Invalid report variables found:", cli_lines, fixed = TRUE)
  if (length(k) == 0) return(character(0))

  block <- cli_lines[(k[1] + 1):length(cli_lines)]

  # 遇到堆栈或空行即停止
  end_idx <- c(
    grep("^\\s*at\\s+", block),
    grep("^\\s*---", block),
    grep("^\\s*$", block)
  )
  if (length(end_idx) > 0) {
    block <- block[1:(min(end_idx) - 1)]
  }

  # 提取变量名（去掉冒号后内容）
  bad <- gsub(":.*$", "", block)
  bad <- bad[nchar(bad) > 0]

  # 转为 APSIM Report 路径（常见：Wheat.xxx -> [Wheat].xxx）
  bad <- gsub("^Wheat\\.", "[Wheat].", bad)

  unique(bad)
}

# ③ 读取当前 Report.VariableNames（JSON 层）
get_report_variablenames_json <- function(sim_path,
                                          zone_name = "Field",
                                          report_name = "Report") {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  zone_info <- get_zone_node(sim, zone_name); zone <- zone_info$node
  rpt_info <- get_report_node(zone, report_name); rpt <- rpt_info$node

  vn <- rpt$VariableNames
  if (is.null(vn)) character(0) else unlist(vn)
}

# ④ 删除无效变量并写回 apsimx
drop_report_vars_and_write <- function(sim_path,
                                       bad_vars,
                                       zone_name = "Field",
                                       report_name = "Report") {
  current <- get_report_variablenames_json(sim_path, zone_name, report_name)
  new_vars <- setdiff(current, bad_vars)

  strict_edit_report_json(
    sim_path = sim_path,
    zone_name = zone_name,
    report_name = report_name,
    event_names = NULL,
    variable_names = new_vars
  )

  new_vars
}

###############################################################################
########################## 10. 批量仿真主循环（v5.1.1）##########################
###############################################################################
message("\n🚀 开始批量仿真...")

test_mode <- FALSE
max_simulations <- if (test_mode) min(2, nrow(sim_configs)) else nrow(sim_configs)
message(if (test_mode) paste("🔧 测试模式：只运行前", max_simulations, "个仿真")
        else paste("🔧 完整运行模式：运行全部", max_simulations, "个仿真"))

success_count <- 0
fail_count <- 0

for (i in 1:max_simulations) {

  config <- sim_configs[i, ]
  env_id <- paste0(config$growing_season, config$site, "_", config$scenario, "_", config$sowing_code)

  message("\n", paste(rep("-", 75), collapse = ""))
  message(paste("🌾 处理环境", i, "/", max_simulations, "：", env_id))
  message(paste("📅 生长季：", config$growing_season))
  message(paste0("📍 站点： ", config$site, " (lat: ", config$latitude, " lon: ", config$longitude, " )"))
  message(paste0("🌡️ 情景： ", config$scenario, " | 播期： ", config$sowing_code, " ( ", sowing_configs[[config$sowing_code]]$name, " )"))

  status <- "开始处理"
  error_msg <- ""
  weather_created <- FALSE
  script_configured <- FALSE
  sim_completed <- FALSE
  output_saved <- FALSE
  used_tt <- FALSE

  ###########################################################################
  # 10.1 生成 met 气象文件 + GDD兜底
  ###########################################################################
  weather_file <- paste0(env_id, ".met")
  weather_path <- file.path(future_dirs$weather, weather_file)

  apsim_met <- NULL
  gdd0 <- NULL
  gdd5 <- NULL

  tryCatch({

    site_scenario_data <- corrected_weather %>%
      dplyr::filter(site == config$site, scenario == config$scenario) %>%
      dplyr::arrange(time)

    if (nrow(site_scenario_data) == 0) stop(paste("没有气象数据：站点", config$site, "情景", config$scenario))

    met_data <- site_scenario_data %>%
      dplyr::filter(time >= as.Date(config$start_date) & time <= as.Date(config$end_date)) %>%
      dplyr::arrange(time)

    if (nrow(met_data) == 0) stop(paste("没有气象数据在范围", config$start_date, "至", config$end_date))

    required_columns <- c("tasmax_corrected", "tasmin_corrected", "pr_corrected",
                          "rsds_corrected", "hurs_corrected", "sfcwind_corrected", "vapr_corrected")
    missing_cols <- setdiff(required_columns, names(met_data))
    if (length(missing_cols) > 0) stop(paste("缺失气象变量：", paste(missing_cols, collapse = ", ")))

    apsim_met <- data.frame(
      date = met_data$time,
      year = lubridate::year(met_data$time),
      day  = lubridate::yday(met_data$time),
      radn = met_data$rsds_corrected,
      maxt = met_data$tasmax_corrected,
      mint = met_data$tasmin_corrected,
      rain = met_data$pr_corrected,
      rh   = met_data$hurs_corrected,
      wind = met_data$sfcwind_corrected,
      vapr = met_data$vapr_corrected,
      stringsAsFactors = FALSE
    )

    tavg <- (apsim_met$maxt + apsim_met$mint) / 2
    tav <- mean(tavg, na.rm = TRUE)

    monthly_avg <- apsim_met %>%
      dplyr::mutate(month = lubridate::month(date)) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(avg_temp = mean((maxt + mint) / 2, na.rm = TRUE), .groups = "drop")
    amp <- (max(monthly_avg$avg_temp, na.rm = TRUE) - min(monthly_avg$avg_temp, na.rm = TRUE)) / 2

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

    writeLines(met_header, weather_path)
    write.table(
      apsim_met[, c("year", "day", "radn", "maxt", "mint", "rain", "rh", "wind")],
      file = weather_path, append = TRUE,
      row.names = FALSE, col.names = FALSE, sep = " "
    )

    gdd0 <- calc_gdd(apsim_met, base_temp = 0); names(gdd0) <- c("date", "GDD_base0_daily", "GDD_base0_cum")
    gdd5 <- calc_gdd(apsim_met, base_temp = 5); names(gdd5) <- c("date", "GDD_base5_daily", "GDD_base5_cum")

    weather_created <- TRUE
    message(paste("✅ 生成气象文件：", weather_file))
    message(paste0("📅 气象时间范围： ", config$start_date, " 至 ", config$end_date, "（共 ", nrow(apsim_met), " 天）"))

  }, error = function(e) {
    weather_created <<- FALSE
    status <<- "失败-气象文件"
    error_msg <<- paste("气象文件错误：", e$message)
    message(paste("❌ 生成气象文件失败：", e$message))
  })

  if (!weather_created) {
    fail_count <- fail_count + 1
    progress_data <- rbind(progress_data, data.frame(
      config_id = i, env_id = env_id, status = status,
      weather_file_created = FALSE, script_configured = FALSE, simulation_completed = FALSE, output_saved = FALSE,
      used_tt = FALSE, error_message = error_msg, timestamp = as.character(Sys.time()),
      stringsAsFactors = FALSE
    ))
    write.csv(progress_data, progress_file, row.names = FALSE)
    next
  }

  ###########################################################################
  # 10.2 复制模板脚本
  ###########################################################################
  sim_name <- paste0(env_id, "_", config$cultivar)
  sim_file <- paste0(sim_name, ".apsimx")
  sim_script_path <- file.path(future_dirs$script, sim_file)
  source_script <- file.path(apsimx_examples_path, "Wheat.apsimx")

  if (!file.exists(sim_script_path)) {
    file.copy(source_script, sim_script_path, overwrite = TRUE)
    message(paste("✅ 复制仿真脚本：", sim_file))
  }

  ###########################################################################
  # 10.3 JSON 修改脚本（Clock / Weather / Manager / Report（含TT候选））
  ###########################################################################
  script_ok <- TRUE

  # 常用变量清单（后续 TT 回退也要用）
  base_vars <- c(
    "[Clock].Today",
    "[Wheat].LAI",
    "[Wheat].Phenology.Zadok.Stage",
    "[Wheat].Phenology.CurrentStageName",
    "[Wheat].AboveGround.Wt",
    "[Wheat].AboveGround.N",
    "[Wheat].Grain.Total.Wt * 10 as Yield",
    "[Wheat].Grain.Protein",
    "[Wheat].Grain.Size",
    "[Wheat].Grain.Number",
    "[Wheat].Grain.Total.Wt",
    "[Wheat].Grain.Total.N",
    "[Wheat].Total.Wt"
  )

  tryCatch({

    # Clock：Start/End
    strict_edit_simnode_fields_json(sim_script_path, "Clock", list(
      Start = paste0(config$start_date, "T00:00:00"),
      End   = paste0(config$end_date,   "T00:00:00")
    ))
    message("✅ JSON 强写 Clock.Start/End 完成")

    # Weather：FileName + 修复 ExcelWorkSheetName 为空
    strict_edit_simnode_fields_json(sim_script_path, "Weather", list(
      FileName = normalizePath(weather_path, winslash = "/", mustWork = FALSE),
      ExcelWorkSheetName = NULL
    ))
    message("✅ JSON 强写 Weather.FileName 完成（并修复 ExcelWorkSheetName=null）")

    # 严格播种日（当天必播）
    plant_date_formatted <- format_apsim_d_mmm(as.Date(config$sowing_date))
    message(paste("📌 严格播种日（APSIM格式 d-mmm 小写）：", plant_date_formatted))

    strict_edit_manager_json(sim_script_path, "SowingRule1", "StartDate", plant_date_formatted)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "EndDate",   plant_date_formatted)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "MinESW",   -1)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "MinRain",  -1)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "RainDays",  1)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "Population", config$plant_density)
    strict_edit_manager_json(sim_script_path, "SowingRule1", "CultivarName", config$cultivar)

    strict_edit_manager_json(sim_script_path, "SowingFertiliser", "Amount", 200)
    strict_edit_manager_json(sim_script_path, "SowingFertiliser", "CropName", "wheat")

    # Report：EventNames/VariableNames（先写 TT 候选）
    event_names_to_set <- c("[Clock].EndOfDay", "[Wheat].Harvesting")
    tt_vars <- get_tt_candidates()

    strict_edit_report_json(
      sim_path = sim_script_path,
      zone_name = "Field",
      report_name = "Report",
      event_names = event_names_to_set,
      variable_names = c(base_vars, tt_vars)
    )
    used_tt <- TRUE
    message("✅ JSON 强写 Report.EventNames/VariableNames 完成（已包含 TT 候选字段）")

    script_configured <- TRUE
    status <- "脚本已配置"

  }, error = function(e) {
    script_ok <<- FALSE
    script_configured <<- FALSE
    status <<- "失败-脚本配置"
    error_msg <<- paste("脚本配置错误：", e$message)
    message(paste("❌ 仿真脚本配置失败：", e$message))
  })

  if (!script_ok) {
    fail_count <- fail_count + 1
    progress_data <- rbind(progress_data, data.frame(
      config_id = i, env_id = env_id, status = status,
      weather_file_created = TRUE, script_configured = FALSE, simulation_completed = FALSE, output_saved = FALSE,
      used_tt = FALSE, error_message = error_msg, timestamp = as.character(Sys.time()),
      stringsAsFactors = FALSE
    ))
    write.csv(progress_data, progress_file, row.names = FALSE)
    next
  }

  ###########################################################################
  # 10.4 运行仿真（v5.1.1：CLI 捕获 + 自动剔除无效变量 + 重试）
  ###########################################################################
  message("▶️ 运行仿真（v5.1.1：智能剔除无效变量）...")

  max_retry <- 3
  cli_log <- character(0)
  ran_ok <- FALSE

  for (rr in seq_len(max_retry)) {

    cli_log <- run_apsim_cli_capture(apsimx_exe_path, sim_script_path)

    # 若未出现 Invalid report variables found: 认为通过
    if (!any(grepl("Invalid report variables found:", cli_log, fixed = TRUE))) {
      ran_ok <- TRUE
      break
    }

    bad_vars <- extract_invalid_report_vars(cli_log)

    if (length(bad_vars) == 0) {
      message("⚠️ 出现 Invalid 但无法解析变量名，准备整体回退 TT…")
      break
    }

    message(paste0("⚠️ 第 ", rr, " 次剔除无效变量："))
    message(paste("   ", paste(bad_vars, collapse = ", ")))

    new_vars <- drop_report_vars_and_write(sim_script_path, bad_vars)

    message(paste0("✅ 剩余 Report 变量数：", length(new_vars)))

    if (length(new_vars) < 5) {
      message("⚠️ Report 变量过少，停止重试，后续走 TT 回退/兜底 GDD")
      break
    }
  }

  # ======= 读取仿真结果（不论 CLI 是否提示 Invalid，最终以 report 是否能读为准） =======
  message("▶️ 读取仿真结果…")

  result <- tryCatch({
    apsimx::apsimx(
      file = sim_file,
      src.dir = future_dirs$script,
      value = "report",
      cleanup = FALSE,
      simplify = TRUE
    )
  }, error = function(e) {
    NULL
  })

  # 若读不到 report：执行 v5.1 的 TT 一刀切回退（仅保留 base_vars）并用 apsimx::apsimx 再跑一次
  if (is.null(result)) {
    message("⚠️ APSIM 运行后无法读取 report：触发 TT 一刀切回退（仅保留基础变量）并重跑一次…")
    used_tt <- FALSE

    ok2 <- tryCatch({
      strict_edit_report_json(
        sim_path = sim_script_path,
        zone_name = "Field",
        report_name = "Report",
        event_names = c("[Clock].EndOfDay", "[Wheat].Harvesting"),
        variable_names = base_vars
      )
      TRUE
    }, error = function(e) FALSE)

    if (ok2) {
      run2 <- run_apsim_once(sim_file = sim_file, src_dir = future_dirs$script)
      result <- run2$result
      if (is.null(result)) {
        status <- "失败-仿真运行"
        error_msg <- paste0("回退后仍无法读取 report：", run2$error)
      }
    } else {
      status <- "失败-仿真运行"
      error_msg <- "TT 回退写入 Report 失败（JSON）"
      result <- NULL
    }
  }

  if (is.null(result)) {
    status <- "失败-仿真运行"
    if (error_msg == "") error_msg <- "APSIM 运行完成但无法读取 report"
    fail_count <- fail_count + 1

    progress_data <- rbind(progress_data, data.frame(
      config_id = i,
      env_id = env_id,
      status = status,
      weather_file_created = TRUE,
      script_configured = TRUE,
      simulation_completed = FALSE,
      output_saved = FALSE,
      used_tt = used_tt,
      error_message = error_msg,
      timestamp = as.character(Sys.time()),
      stringsAsFactors = FALSE
    ))
    write.csv(progress_data, progress_file, row.names = FALSE)

    message(paste("❌ 仿真读取失败：", error_msg))
    next
  }

  sim_completed <- TRUE
  status <- "仿真已完成"
  message("✅ 仿真运行成功")

  ###########################################################################
  # 10.5 合并积温（若没 TT，则合并 GDD 兜底）并保存结果
  ###########################################################################
  date_col <- NULL
  if ("Clock.Today" %in% names(result)) date_col <- "Clock.Today"
  if ("[Clock].Today" %in% names(result)) date_col <- "[Clock].Today"
  if ("Today" %in% names(result)) date_col <- "Today"

  if (!is.null(date_col)) {
    result$date <- as.Date(result[[date_col]])
  } else {
    result$date <- NA
    message("⚠️ 结果中未找到日期列（Clock.Today/Today），无法合并GDD，但仍会保存结果")
  }

  # 若结果里没有 TT，则合并气象 GDD 兜底
  if (!all(is.na(result$date)) && !is.null(gdd0) && !is.null(gdd5)) {
    if (!has_tt_in_result(result)) {
      result <- result %>%
        dplyr::left_join(gdd0, by = "date") %>%
        dplyr::left_join(gdd5, by = "date")
      message("✅ 结果未检测到 TT 字段：已合并气象 GDD(base0/base5) 兜底输出")
    } else {
      message("✅ 检测到 APSIM 内部 TT/积温字段：不再合并 GDD（避免重复）")
    }
  }

  output_file <- paste0(env_id, "_simulation_output.csv")
  output_path <- file.path(future_dirs$output, output_file)
  write.csv(result, output_path, row.names = FALSE)

  if (file.exists(output_path)) {
    file_info <- file.info(output_path)
    if (!is.na(file_info$size) && file_info$size > 1024) {
      message(paste("✅ 结果保存成功：", output_file, "（", round(file_info$size / 1024, 1), "KB）"))
      output_saved <- TRUE
      status <- "成功"
      success_count <- success_count + 1
    } else {
      status <- "警告-结果文件"
      error_msg <- "结果文件过小"
      message(paste("⚠️ 结果文件过小：", output_file))
      fail_count <- fail_count + 1
    }
  } else {
    status <- "失败-文件保存"
    error_msg <- "结果文件未创建"
    message(paste("❌ 结果文件未创建：", output_file))
    fail_count <- fail_count + 1
  }

  ###########################################################################
  # 10.6 写入进度
  ###########################################################################
  progress_data <- rbind(progress_data, data.frame(
    config_id = i,
    env_id = env_id,
    status = status,
    weather_file_created = TRUE,
    script_configured = TRUE,
    simulation_completed = sim_completed,
    output_saved = output_saved,
    used_tt = used_tt,
    error_message = error_msg,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  ))
  write.csv(progress_data, progress_file, row.names = FALSE)

  message(paste("📊 进度：", i, "/", max_simulations, "| 成功：", success_count, "失败：", fail_count))
  Sys.sleep(0.05)
}

message("\n", paste(rep("=", 75), collapse = ""))
message("✅ 批量仿真循环结束")
message(paste("🎯 最终：成功", success_count, "失败", fail_count))
message(paste("进度文件：", progress_file))

###############################################################################
########################## 11. 总结报告 ########################################
###############################################################################
message("\n📊 生成总结报告...")

progress_summary <- if (file.exists(progress_file)) read.csv(progress_file, stringsAsFactors = FALSE) else data.frame()
status_summary <- if (nrow(progress_summary) > 0) {
  progress_summary %>% dplyr::group_by(status) %>% dplyr::summarise(count = dplyr::n(), .groups = "drop")
} else {
  data.frame(status = "无数据", count = 0)
}

summary_report <- file.path(future_dirs$config, "future_simulation_final_summary_v5_1_1.txt")
sink(summary_report)

cat("未来气候小麦模拟完成总结报告（Linux FULL FINAL FIX v5.1.1：JSON全强写 + CLI智能剔除无效变量 + TT回退 + GDD兜底）\n")
cat(paste0(rep("=", 110), collapse = ""), "\n\n")
cat("报告生成时间:", as.character(Sys.time()), "\n")
cat("工作目录:", work_dir, "\n\n")
cat("总配置数:", nrow(sim_configs), "\n")
cat("实际运行数:", max_simulations, "\n")
cat("成功数:", success_count, "\n")
cat("失败数:", fail_count, "\n")
cat("成功率:", ifelse(max_simulations > 0, round(success_count / max_simulations * 100, 2), 0), "%\n\n")
cat("状态分布：\n")
print(status_summary)

cat("\n\n输出目录：\n")
for (nm in names(future_dirs)) cat(" - ", nm, ": ", future_dirs[[nm]], "\n", sep = "")

cat("\n严格播种说明（必播逻辑）：\n")
cat(" - StartDate = EndDate = 播种日（播种窗口=1天）\n")
cat(" - MinESW = -1（ESW > -1 恒成立）\n")
cat(" - MinRain = -1（accumulatedRain.Sum > -1 恒成立，即使 0 雨也能播）\n")
cat(" - RainDays = 1\n")

cat("\n积温/TT 输出说明：\n")
cat(" - 先尝试在 Report.VariableNames 写入 TT 候选字段（默认仅 ThermalTime，最稳）\n")
cat(" - v5.1.1：若 CLI 提示 Invalid report variables found，则自动解析并精准删除无效变量后重试（最多3次）\n")
cat(" - 若最终 report 仍无法读取，则回退：删除 TT 字段，仅保留基础变量并重跑一次\n")
cat(" - 若最终结果未检测到 TT 字段，则合并气象计算的 GDD(base0/base5) 兜底输出\n")

cat("\nJSON 强写说明：\n")
cat(" - 本 v5.1.1 完全不依赖 edit_apsimx()，Clock/Weather/Manager/Report 全部走 JSON 强写\n")
cat(" - 同时强制修复 Weather.ExcelWorkSheetName = null，避免被写坏为 {} 导致 APSIM 解析崩溃\n")

sink()
message(paste("✅ 总结报告保存：", summary_report))

final_sep <- paste(rep("=", 80), collapse = "")
cat("\n", final_sep, "\n", sep = "")
cat("🎉 未来气候小麦模拟任务完成（v5.1.1）！\n\n")
cat("成功:", success_count, " 失败:", fail_count, "\n")
cat("进度文件:", progress_file, "\n")
cat("总结报告:", summary_report, "\n")
cat("\n", final_sep, "\n", sep = "")
message(paste("结束时间:", Sys.time()))
