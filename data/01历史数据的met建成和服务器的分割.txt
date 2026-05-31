###############################################################################
# 文件名：run_apsim_history_power_linux_PARALLEL_FINAL_v1_2_FULL.R
# 位置：/home/wangyanan/APSIM_Wheat_Project_history/
# 描述：历史种植信息 + NASA POWER 自动下载气象（函数完全不动） + APSIM 批量并行（Linux）
# v1.2 FULL 关键能力：
#   ✅ 并行：parallel::mclapply（Linux），支持 8核/16核
#   ✅ POWER 下载限流：文件锁控制最大同时下载数（默认 1，防止 SSL reset）
#   ✅ JSON 全强写：Clock / Weather / Manager / Report 全部 JSON 写入（避免并行 edit_apsimx 写坏）
#   ✅ wrapper 自动补 run：避免 Linux 上 CLI 调用 “Verb xxx.apsimx not recognized”
#   ✅ CLI 捕获日志：自动解析并剔除 Report 无效变量，最多重试 3 次
#   ✅ report 读取失败兜底：自动降级 Report 变量集 + CLI 再跑一次，再读 report
#   ✅ 失败不中断：每个环境独立返回 ok/err，最终汇总失败清单并保存
#   ✅ 播种/施肥 Manager 自动探测（兼容不同 Wheat.apsimx 模板）
#   ✅ 并行前模板预检（提前报错，避免全跑全挂）
#   ✅ 复制模板提前（先确保模板/Manager OK，再去请求 POWER，避免无意义请求）
###############################################################################

cat("\014")
message("🎯 APSIM 历史气象（POWER）并行版 - Linux Server（FINAL v1.2 FULL）")
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
safe_require("jsonlite")
suppressWarnings(suppressMessages(library(parallel)))
message("✅ parallel 加载成功（基础包）")

# 网络/下载类的总体超时（不改下载函数体，只做全局设置）
options(timeout = 300)

###############################################################################
########################## 2. 配置工作目录与 APSIM（Linux）######################
###############################################################################
message("\n⚙️ 配置工作目录与 APSIM（Linux）...")

# 服务器工作目录（按你的服务器实际路径改）
work_dir <- "/home/wangyanan/APSIM_Wheat_Project_history"
if (!dir.exists(work_dir)) stop(paste("❌ 工作目录不存在：", work_dir))
setwd(work_dir)
message(paste("✅ 工作目录：", work_dir))

# APSIM Linux wrapper（按你的服务器实际安装改）
apsimx_exe_path <- "/software/apsim_bin/ApsimX/bin/Release/net8.0/apsim_wrapper.sh"
if (!file.exists(apsimx_exe_path)) stop(paste("❌ APSIM wrapper 不存在：", apsimx_exe_path))
try(Sys.chmod(apsimx_exe_path, mode = "0755"), silent = TRUE)

apsimx_examples_path <- "/software/apsim_bin/ApsimX/Examples"
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
    "  echo \"自动补上 'run' 命令: $1\"",
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
########################## 3. 读取种植信息 #####################################
###############################################################################
message("\n📍 读取种植信息...")

planting_file <- "planting2000_2024.csv"
if (!file.exists(planting_file)) stop(paste("❌ 种植数据缺失：", file.path(getwd(), planting_file)))
planting_info <- read.csv(planting_file, stringsAsFactors = FALSE)
message(paste("✅ 读取种植数据：", nrow(planting_info), "个环境"))

cultivars <- c("Keyu13")

###############################################################################
########################## 4. 目录结构 ########################################
###############################################################################
message("\n📁 创建目录结构...")

dirs <- list(
  weather = file.path(getwd(), "01APSIMmet"),
  script  = file.path(getwd(), "01APSIMinput_apsimx"),
  output  = file.path(getwd(), "01APSIMout"),
  logs    = file.path(getwd(), "logs"),
  config  = file.path(getwd(), "config")
)
for (nm in names(dirs)) {
  if (!dir.exists(dirs[[nm]])) dir.create(dirs[[nm]], recursive = TRUE, showWarnings = FALSE)
}
message("✅ 目录准备完成")

###############################################################################
########################## 5. 你指定的气象函数（完全不动）#######################
###############################################################################
# 这是从你的Windows代码中提取的有效函数（绝对不改动）
download_and_generate_weather <- function(lat, lon, start_date, end_date, weather_dir, weather_file) {
  max_attempts <- 5
  attempt <- 1
  start_time <- Sys.time()
  
  while (attempt <= max_attempts) {
    time_elapsed <- difftime(Sys.time(), start_time, units = "hours")
    if (as.numeric(time_elapsed) > 1) stop("❌ 气象下载超时（>1小时）")
    
    tryCatch({
      pwr <- apsimx::get_power_apsim_met(lonlat = c(lon, lat), dates = c(start_date, end_date))
      apsimx::write_apsim_met(pwr, wrt.dir = weather_dir, filename = weather_file)
      message(paste("✅ 气象文件保存：", file.path(weather_dir, weather_file)))
      return(TRUE)
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

###############################################################################
########################## 6. POWER 下载限流（文件锁）###########################
###############################################################################
acquire_download_slot <- function(lock_dir, max_slots = 1, wait_s = 2) {
  dir.create(lock_dir, showWarnings = FALSE, recursive = TRUE)
  repeat {
    locks <- list.files(lock_dir, full.names = TRUE)

    # 清理过期锁（>30min）
    if (length(locks) > 0) {
      info <- file.info(locks)
      old <- which(difftime(Sys.time(), info$mtime, units = "mins") > 30)
      if (length(old) > 0) unlink(locks[old], force = TRUE)
    }

    locks <- list.files(lock_dir, full.names = TRUE)
    if (length(locks) < max_slots) {
      f <- file.path(lock_dir, paste0("slot_", Sys.getpid(), "_", as.integer(runif(1, 1, 1e9))))
      file.create(f)
      return(f)
    }
    Sys.sleep(wait_s + runif(1, 0, 1))
  }
}

release_download_slot <- function(lock_file) {
  if (!is.null(lock_file) && file.exists(lock_file)) unlink(lock_file, force = TRUE)
}

# 外层包锁：不改 download_and_generate_weather 本体
download_weather_with_lock <- function(lat, lon, start_date, end_date, weather_dir, weather_file,
                                       max_parallel_download = 1) {
  lock_dir <- file.path(dirs$config, "power_download_locks")
  lock_file <- acquire_download_slot(lock_dir, max_slots = max_parallel_download)
  on.exit(release_download_slot(lock_file), add = TRUE)

  # 可选：轻微抖动，减少 TLS 同步冲击（不改函数本体）
  Sys.sleep(runif(1, 0.2, 1.5))

  message(sprintf("🌤️ [pid=%s] got lock, downloading %s", Sys.getpid(), weather_file))

  download_and_generate_weather(
    lat = lat, lon = lon,
    start_date = start_date, end_date = end_date,
    weather_dir = weather_dir, weather_file = weather_file
  )
}

###############################################################################
########################## 7. 工具函数：日期/JSON/Manager/Report/CLI ###########
###############################################################################
format_apsim_d_mmm <- function(date_obj) {
  original_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C")
  out <- format(as.Date(date_obj), "%d-%b")
  Sys.setlocale("LC_TIME", original_locale)
  tolower(out)
}

# ---- JSON 读写 ----
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

# Weather.ExcelWorkSheetName 强制 null
fix_weather_excelworksheet_null <- function(sim) {
  w_idx <- find_child_idx_by_name(sim$Children, "Weather")
  if (length(w_idx) == 0) return(sim)
  w <- sim$Children[[w_idx[1]]]
  w$ExcelWorkSheetName <- NULL
  sim$Children[[w_idx[1]]] <- w
  sim
}

# 写 Simulation 子节点字段（Clock/Weather）
strict_edit_simnode_fields_json <- function(sim_path, node_name, fields_named_list) {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  sim <- fix_weather_excelworksheet_null(sim)

  n_idx <- find_child_idx_by_name(sim$Children, node_name)
  if (length(n_idx) == 0) stop(paste0("找不到 Simulation 子节点：", node_name))
  n <- sim$Children[[n_idx[1]]]
  for (k in names(fields_named_list)) n[[k]] <- fields_named_list[[k]]

  sim$Children[[n_idx[1]]] <- n
  js$Children[[sim_info$idx]] <- sim
  write_apsimx_json(js, sim_path)
  TRUE
}

# 写 Manager 参数（Key/Value）
strict_edit_manager_json <- function(sim_path, manager_name, key, value, zone_name = "Field") {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  sim <- fix_weather_excelworksheet_null(sim)

  zone_info <- get_zone_node(sim, zone_name); zone <- zone_info$node
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

# 写 Report
strict_edit_report_json <- function(sim_path, event_names = NULL, variable_names = NULL,
                                   zone_name = "Field", report_name = "Report") {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  sim <- fix_weather_excelworksheet_null(sim)

  zone_info <- get_zone_node(sim, zone_name); zone <- zone_info$node
  rpt_info <- get_report_node(zone, report_name); rpt <- rpt_info$node

  if (!is.null(event_names)) rpt$EventNames <- unname(as.character(event_names))
  if (!is.null(variable_names)) rpt$VariableNames <- unname(as.character(variable_names))

  zone$Children[[rpt_info$idx]] <- rpt
  sim$Children[[zone_info$idx]] <- zone
  js$Children[[sim_info$idx]] <- sim
  write_apsimx_json(js, sim_path)
  TRUE
}

###############################################################################
########################## 8. Manager 自动探测（兼容 Wheat 模板）###############
###############################################################################
is_manager_node <- function(x) {
  !is.null(x$`$type`) && grepl("Models\\.Manager", x$`$type`)
}
find_manager_idx_fuzzy <- function(zone, patterns, prefer_names = character(0)) {
  if (is.null(zone$Children) || length(zone$Children) == 0) return(integer(0))
  idx_mgr <- which(vapply(zone$Children, is_manager_node, logical(1)))
  if (length(idx_mgr) == 0) return(integer(0))

  if (length(prefer_names) > 0) {
    hit <- idx_mgr[vapply(idx_mgr, function(i) {
      nm <- zone$Children[[i]]$Name
      !is.null(nm) && nm %in% prefer_names
    }, logical(1))]
    if (length(hit) > 0) return(hit[1])
  }

  hit <- idx_mgr[vapply(idx_mgr, function(i) {
    node <- zone$Children[[i]]
    nm <- if (!is.null(node$Name)) tolower(node$Name) else ""
    sc <- ""
    if (!is.null(node$Script)) sc <- tolower(paste(node$Script, collapse = "\n"))
    any(vapply(patterns, function(p) grepl(p, nm) || grepl(p, sc), logical(1)))
  }, logical(1))]

  if (length(hit) == 0) integer(0) else hit[1]
}

auto_get_manager_name <- function(sim_path, zone_name = "Field", purpose = c("sowing", "fertiliser")) {
  purpose <- match.arg(purpose)

  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  zone_info <- get_zone_node(sim, zone_name); zone <- zone_info$node

  if (purpose == "sowing") {
    prefer <- c("SowingRule1", "SowingRule", "Sowing", "Sow", "Plant", "PlantingRule", "SowingManager")
    patterns <- c("sowing", "\\bsow\\b", "plant", "cultivar", "population", "startdate", "enddate")
    idx <- find_manager_idx_fuzzy(zone, patterns = patterns, prefer_names = prefer)
  } else {
    prefer <- c("SowingFertiliser", "Fertiliser", "Fertilizer", "Fertilise", "Fertilize", "ApplyFertiliser")
    patterns <- c("fert", "urea", "nitrate", "apply", "amount", "cropname")
    idx <- find_manager_idx_fuzzy(zone, patterns = patterns, prefer_names = prefer)
  }

  if (length(idx) == 0) return(NA_character_)
  as.character(zone$Children[[idx]]$Name)
}

###############################################################################
########################## 9. CLI 捕获 + invalid 变量解析/剔除 ##################
###############################################################################
run_apsim_cli_capture <- function(apsim_exe, apsimx_path) {
  tryCatch({
    system2(apsim_exe, args = c("run", apsimx_path), stdout = TRUE, stderr = TRUE)
  }, error = function(e) c(paste("SYSTEM2_ERROR:", e$message)))
}

# 关键：把 CLI 中的 "ISoilWater.LeachCl" 映射为 "[ISoilWater].LeachCl"
extract_invalid_report_vars <- function(cli_lines) {
  k <- grep("Invalid report variables found:", cli_lines, fixed = TRUE)
  if (length(k) == 0) return(character(0))

  block <- cli_lines[(k[1] + 1):length(cli_lines)]

  # 截断到堆栈/空行/分隔线前
  end_idx <- c(grep("^\\s*at\\s+", block), grep("^\\s*---", block), grep("^\\s*$", block))
  if (length(end_idx) > 0) block <- block[1:(min(end_idx) - 1)]

  x <- gsub(":.*$", "", block)
  x <- trimws(x)
  x <- x[nchar(x) > 0]

  x <- vapply(x, function(s) {
    if (grepl("^\\[.+\\]\\.\\w+", s)) return(s)
    if (grepl("^[A-Za-z0-9_]+\\.[A-Za-z0-9_\\.]+$", s)) {
      parts <- strsplit(s, "\\.", fixed = FALSE)[[1]]
      return(paste0("[", parts[1], "].", paste(parts[-1], collapse = ".")))
    }
    s
  }, character(1))

  unique(x)
}

get_report_variablenames_json <- function(sim_path, zone_name = "Field", report_name = "Report") {
  js <- read_apsimx_json(sim_path)
  sim_info <- get_simulation_node(js); sim <- sim_info$node
  zone_info <- get_zone_node(sim, zone_name); zone <- zone_info$node
  rpt_info <- get_report_node(zone, report_name); rpt <- rpt_info$node
  vn <- rpt$VariableNames
  if (is.null(vn)) character(0) else unlist(vn)
}

drop_report_vars_and_write <- function(sim_path, bad_vars, zone_name = "Field", report_name = "Report") {
  current <- get_report_variablenames_json(sim_path, zone_name, report_name)
  new_vars <- setdiff(current, bad_vars)
  strict_edit_report_json(sim_path, event_names = NULL, variable_names = new_vars,
                         zone_name = zone_name, report_name = report_name)
  new_vars
}

read_report_safe <- function(sim_file, src_dir) {
  tryCatch({
    apsimx::apsimx(
      file = sim_file,
      src.dir = src_dir,
      value = "report",
      cleanup = FALSE,
      simplify = TRUE
    )
  }, error = function(e) NULL)
}

###############################################################################
########################## 10. 并行参数 ########################################
###############################################################################
message("\n🧵 配置并行参数...")

n_cores <- 16
avail <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
if (!is.na(avail)) n_cores <- max(1, min(n_cores, max(1, avail - 1)))
message(paste0("✅ 并行核数 n_cores = ", n_cores, "（detectCores=", avail, "）"))

# 关键：下载必须串行或极低并发（你遇到 SSL reset 时必须为 1）
max_parallel_download <- 1
message(paste0("✅ POWER 最大同时下载数 = ", max_parallel_download))

###############################################################################
########################## 11. 模板预检（避免全跑全挂）##########################
###############################################################################
message("\n🔎 模板预检：探测播种/施肥 Manager ...")

tmp_probe <- file.path(dirs$script, "__manager_probe__.apsimx")
file.copy(file.path(apsimx_examples_path, "Wheat.apsimx"), tmp_probe, overwrite = TRUE)

probe_sow  <- auto_get_manager_name(tmp_probe, purpose = "sowing")
probe_fert <- auto_get_manager_name(tmp_probe, purpose = "fertiliser")

message(paste0("🔎 模板播种 Manager 探测结果: ", probe_sow))
message(paste0("🔎 模板施肥 Manager 探测结果: ", probe_fert))

if (is.na(probe_sow) || is.na(probe_fert)) {
  stop("❌ 模板预检失败：无法识别播种/施肥 Manager。请先打开 Wheat.apsimx 确认 Field 下 Manager 名称/脚本。")
}

###############################################################################
########################## 12. 单任务函数（每个环境）############################
###############################################################################
run_one_env <- function(job, cultivar = "Keyu13") {

  tryCatch({

    dat <- planting_info[job, ]
    env_id <- as.character(dat$year_loc)

    # 输出文件/日志
    weather_file <- paste0(env_id, ".met")
    weather_path <- file.path(dirs$weather, weather_file)

    cli_log_path <- file.path(dirs$logs, paste0(env_id, "_cli.log"))

    sim_name <- paste0(env_id, "_", cultivar)
    sim_file <- paste0(sim_name, ".apsimx")
    sim_script_path <- file.path(dirs$script, sim_file)

    # 1) 日期
    plant_date <- as.Date(dat$sowing_date)
    harvest_date <- as.Date(dat$end_date)
    if (is.na(harvest_date)) harvest_date <- plant_date + 270

    start_date <- plant_date - 30
    end_date <- harvest_date + 30

    # 2) 复制模板（提前：先确保模板/Manager OK，再请求 POWER）
    source_script <- file.path(apsimx_examples_path, "Wheat.apsimx")
    file.copy(source_script, sim_script_path, overwrite = TRUE)

    # 3) 下载气象（每个 env 独立 met；必须用你的函数；外层加锁）
    if (file.exists(weather_path)) file.remove(weather_path)

    download_weather_with_lock(
      lat = dat$latitude,
      lon = dat$longitude,
      start_date = start_date,
      end_date = end_date,
      weather_dir = dirs$weather,
      weather_file = weather_file,
      max_parallel_download = max_parallel_download
    )

    if (!file.exists(weather_path)) stop("met 写出失败（文件不存在）")

    # 4) JSON 强写 Clock / Weather
    strict_edit_simnode_fields_json(sim_script_path, "Clock", list(
      Start = paste0(as.character(start_date), "T00:00:00"),
      End   = paste0(as.character(end_date),   "T00:00:00")
    ))
    strict_edit_simnode_fields_json(sim_script_path, "Weather", list(
      FileName = normalizePath(weather_path, winslash = "/", mustWork = FALSE),
      ExcelWorkSheetName = NULL
    ))

    # 5) 播种/施肥（自动探测 Manager）
    plant_date_formatted <- format_apsim_d_mmm(plant_date)

    sowing_mgr <- auto_get_manager_name(sim_script_path, zone_name = "Field", purpose = "sowing")
    fert_mgr   <- auto_get_manager_name(sim_script_path, zone_name = "Field", purpose = "fertiliser")

    if (is.na(sowing_mgr) || sowing_mgr == "") {
      stop("找不到播种 Manager（模板中无 sow/sowing/plant 相关 Manager）。请检查 Wheat.apsimx 的 Field 下 Manager 名称。")
    }
    if (is.na(fert_mgr) || fert_mgr == "") {
      stop("找不到施肥 Manager（模板中无 fert/fertiliser 相关 Manager）。请检查 Wheat.apsimx 的 Field 下 Manager 名称。")
    }

    strict_edit_manager_json(sim_script_path, sowing_mgr, "StartDate", plant_date_formatted)
    strict_edit_manager_json(sim_script_path, sowing_mgr, "EndDate",   plant_date_formatted)
    strict_edit_manager_json(sim_script_path, sowing_mgr, "MinRain",  -1)
    strict_edit_manager_json(sim_script_path, sowing_mgr, "MinESW",   -1)
    try(strict_edit_manager_json(sim_script_path, sowing_mgr, "RainDays", 1), silent = TRUE)
    strict_edit_manager_json(sim_script_path, sowing_mgr, "Population", dat$plant_density)
    strict_edit_manager_json(sim_script_path, sowing_mgr, "CultivarName", cultivar)

    strict_edit_manager_json(sim_script_path, fert_mgr, "Amount", 200)
    try(strict_edit_manager_json(sim_script_path, fert_mgr, "CropName", "wheat"), silent = TRUE)

    # 6) Report：先写完整变量集
    report_vars <- c(
      "[Clock].Today", "[Wheat].LAI", "[Wheat].Phenology.Zadok.Stage",
      "[Wheat].Phenology.CurrentStageName", "[Wheat].AboveGround.Wt",
      "[Wheat].AboveGround.N", "[Wheat].Grain.Total.Wt * 10 as Yield",
      "[Wheat].Grain.Protein", "[Wheat].Grain.Size", "[Wheat].Grain.Number",
      "[Wheat].Grain.Total.Wt", "[Wheat].Grain.Total.N", "[Wheat].Total.Wt",
      "[Wheat].Phenology.Stage", "[SurfaceOrganicMatter].Wt",
      "[Soil].Water.PAWmm", "[Soil].Water.PAW", "[Physical].PAWCmm", "[Physical].PAWC",
      "[ISoilWater].Eos", "[ISoilWater].Es", "[ISoilWater].Eo",
      "[ISoilWater].Runoff", "[ISoilWater].Infiltration", "[ISoilWater].Drainage",
      "[ISoilWater].LeachNO3", "[ISoilWater].LeachNH4", "[ISoilWater].LeachCl",
      "[NFlow].Natm", "[NFlow].N2Oatm", "[Hydrolysis].Value",
      "[Nutrient].MineralisedN", "[Nutrient].Catm"
    )

    strict_edit_report_json(
      sim_script_path,
      event_names = c("[Clock].EndOfDay", "[Wheat].Harvesting"),
      variable_names = report_vars
    )

    # 7) CLI：invalid report vars 自动剔除，最多3次
    for (rr in 1:3) {

      cli_log <- run_apsim_cli_capture(apsimx_exe_path, sim_script_path)

      try(write(paste0("\n\n===== TRY ", rr, " @ ", Sys.time(), " =====\n"), file = cli_log_path, append = TRUE), silent = TRUE)
      try(write(cli_log, file = cli_log_path, append = TRUE), silent = TRUE)

      if (!any(grepl("Invalid report variables", cli_log, fixed = TRUE))) break

      bad_vars <- extract_invalid_report_vars(cli_log)
      if (length(bad_vars) == 0) break

      new_vars <- drop_report_vars_and_write(sim_script_path, bad_vars)
      if (length(new_vars) < 5) break
    }

    # 8) 读取 report（R 端）
    result <- read_report_safe(sim_file, dirs$script)

    # 9) report 读取失败兜底：降级 Report 变量集 + CLI 再跑一次 + 再读
    if (is.null(result)) {

      fallback_vars <- c(
        "[Clock].Today",
        "[Wheat].LAI",
        "[Wheat].AboveGround.Wt",
        "[Wheat].Grain.Total.Wt * 10 as Yield"
      )

      ok_fb <- tryCatch({
        strict_edit_report_json(
          sim_script_path,
          event_names = c("[Clock].EndOfDay", "[Wheat].Harvesting"),
          variable_names = fallback_vars
        )
        TRUE
      }, error = function(e) FALSE)

      if (ok_fb) {
        cli_log2 <- run_apsim_cli_capture(apsimx_exe_path, sim_script_path)
        try(write(paste0("\n\n===== FALLBACK RUN @ ", Sys.time(), " =====\n"), file = cli_log_path, append = TRUE), silent = TRUE)
        try(write(cli_log2, file = cli_log_path, append = TRUE), silent = TRUE)

        result <- read_report_safe(sim_file, dirs$script)
      }
    }

    if (is.null(result)) stop("APSIM 运行完成但无法读取 report（已触发兜底仍失败）")

    # 10) 保存结果
    output_file <- paste0(env_id, "_simulation_output.csv")
    output_path <- file.path(dirs$output, output_file)
    write.csv(result, output_path, row.names = FALSE)

    list(ok = TRUE, env_id = env_id, output = output_path, err = "")

  }, error = function(e) {
    env_id <- tryCatch(as.character(planting_info$year_loc[job]), error = function(x) paste0("job_", job))
    list(ok = FALSE, env_id = env_id, output = NA, err = e$message)
  })
}

###############################################################################
########################## 13. 并行运行 #######################################
###############################################################################
# 默认：先跑一小段测试；确认稳定后改为 1:nrow(planting_info)
jobs <- 1:20

message("\n🚀 开始并行运行：")
message(paste0(" - 环境数：", length(jobs)))
message(paste0(" - 核数：", n_cores))
message(paste0(" - POWER 同时下载：", max_parallel_download))

res <- parallel::mclapply(
  jobs,
  FUN = function(j) run_one_env(j, cultivar = cultivars[1]),
  mc.cores = n_cores,
  mc.preschedule = TRUE
)

ok <- vapply(res, function(x) isTRUE(x$ok), logical(1))

message("====================================================================================")
message(paste0("✅ 完成：成功 ", sum(ok), " 失败 ", length(ok) - sum(ok)))
message(paste0("📂 输出目录：", dirs$output))
message(paste0("📂 气象目录：", dirs$weather))
message(paste0("📂 脚本目录：", dirs$script))
message(paste0("📂 日志目录：", dirs$logs))

# 失败清单
bad <- Filter(function(x) !isTRUE(x$ok), res)
if (length(bad) > 0) {
  message("\n❌ 失败环境清单：")
  for (x in bad) {
    message(paste0(" - ", x$env_id, " | ", x$err))
  }

  # 把失败清单写文件，便于二次补跑
  fail_csv <- file.path(dirs$logs, paste0("failed_envs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  fail_df <- data.frame(
    env_id = vapply(bad, `[[`, character(1), "env_id"),
    err    = vapply(bad, `[[`, character(1), "err"),
    stringsAsFactors = FALSE
  )
  try(write.csv(fail_df, fail_csv, row.names = FALSE), silent = TRUE)
  message(paste0("🧾 失败清单已保存：", fail_csv))
}

message(paste("结束时间:", Sys.time()))
