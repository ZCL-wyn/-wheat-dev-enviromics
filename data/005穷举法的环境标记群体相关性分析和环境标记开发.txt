# ============================================================
# 全流程脚本（从读取数据开始）- env_id 升级版 + 加速版
# 路径：/mnt/7t_storage/zhangcl/weather/
#
# 主要输出（你要的矩阵文件名已修正）：
# tables/【环境×开花偏移窗口矩阵env_id版】env_id_flowering_offset_window_matrix_all.csv
#
# 说明：
# - 第 1~15 步：和你原流程一致（env_id、开花/成熟、共同环境、相关性穷举）
# - 第 16 步：重写为 chunk + 并行 foreach + 分块写盘（显著加速）
# - 所有输入数据默认都在 /mnt/7t_storage/zhangcl/weather/
# ============================================================

rm(list = ls()); gc()

# ===================== 1. 工作目录 =====================
currentWorkingDir <- "/mnt/7t_storage/zhangcl/weather"
setwd(currentWorkingDir)

# ===================== 2. 加载必要包 =====================
cat("加载必要包...\n")
packages_needed <- c(
  "data.table", "dplyr", "tidyr", "lubridate",
  "ggplot2", "readr", "stringr", "zoo",
  "future.apply"
)
for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ===================== 3. 创建输出目录 =====================
cat("创建输出目录...\n")
output_dir <- file.path(currentWorkingDir, "FloweringOffset_EnvID_Analysis_Complete_V2")
dirs_to_create <- c(
  file.path(output_dir, "tables"),
  file.path(output_dir, "plots"),
  file.path(output_dir, "debug_info"),
  file.path(output_dir, "logs"),
  file.path(output_dir, "flowering_offset_analysis"),
  file.path(output_dir, "window_values"),
  file.path(output_dir, "window_values", "mapping"),
  file.path(output_dir, "all_environments"),
  file.path(output_dir, "data_quality")
)
for (dir_path in c(output_dir, dirs_to_create)) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

start_time <- Sys.time()
cat("\n分析开始时间:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# ============================================================
# 4. 读取数据（全部来自 /mnt/7t_storage/zhangcl/weather/）
# ============================================================
cat("\n=== Step 1: 读取数据 ===\n")

# 4.1 planting
planting_file <- file.path(currentWorkingDir, "planting.txt")
if (!file.exists(planting_file)) stop("planting.txt 不存在: ", planting_file)
planting_info <- tryCatch(
  {
    read.table(planting_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               quote = "", fill = TRUE, na.strings = c("", "NA", "null", "NULL", "NaN", "#N/A", "-999"))
  },
  error = function(e) {
    read.table(planting_file, header = TRUE, sep = ",", stringsAsFactors = FALSE,
               quote = "", fill = TRUE, na.strings = c("", "NA", "null", "NULL", "NaN", "#N/A", "-999"))
  }
)
cat("planting_info:", dim(planting_info), "\n")

# 4.2 phenotype
phenotype_file <- file.path(currentWorkingDir, "TKW_mean_table.txt")
if (!file.exists(phenotype_file)) stop("TKW_mean_table.txt 不存在: ", phenotype_file)
phenotype_data <- read.table(
  phenotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  na.strings = c("", "NA", "null", "NULL", "NaN", "#N/A", "-999")
)
cat("phenotype_data:", dim(phenotype_data), "\n")

# 4.3 merged weather
merged_weather_file <- file.path(currentWorkingDir, "merged_historical_future_growthstage_data.csv")
if (!file.exists(merged_weather_file)) stop("merged_historical_future_growthstage_data.csv 不存在: ", merged_weather_file)
weather_data <- data.table::fread(
  merged_weather_file, stringsAsFactors = FALSE,
  na.strings = c("", "NA", "null", "NULL", "NaN", "#N/A", "-999", ".")
)
cat("weather_data:", dim(weather_data), "\n")

# ============================================================
# 5. env_id 升级（env_id = ENV + '__' + Data_Source）
# ============================================================
cat("\n=== Step 2: env_id 升级 ===\n")

if (!"Data_Source" %in% names(weather_data)) {
  weather_data[, Data_Source := "UNKNOWN"]
}
if (!"env" %in% names(weather_data)) stop("weather_data 缺少 env 列")

weather_data[, env_id := ifelse(!is.na(Data_Source) & Data_Source != "",
                               paste(toupper(trimws(env)), toupper(trimws(Data_Source)), sep = "__"),
                               toupper(trimws(env)))]

# env 映射表（env_upper + Data_Source -> env_id）
env_mapping <- weather_data %>%
  dplyr::select(env_raw = env, Data_Source, env_id) %>%
  dplyr::mutate(env_upper = toupper(trimws(env_raw))) %>%
  dplyr::distinct(env_upper, Data_Source, env_id) %>%
  dplyr::arrange(env_upper, Data_Source)

multi_source_envs <- env_mapping %>%
  dplyr::group_by(env_upper) %>%
  dplyr::summarise(n_sources = n_distinct(Data_Source), .groups = "drop") %>%
  dplyr::filter(n_sources > 1) %>%
  dplyr::arrange(desc(n_sources))

if (nrow(multi_source_envs) > 0) {
  write.csv(multi_source_envs, file.path(output_dir, "tables", "multi_source_environments.csv"), row.names = FALSE)

  default_sources <- weather_data %>%
    dplyr::filter(!is.na(Data_Source)) %>%
    dplyr::mutate(env_upper = toupper(trimws(env))) %>%
    dplyr::group_by(env_upper, Data_Source) %>%
    dplyr::summarise(n_records = n(), .groups = "drop") %>%
    dplyr::group_by(env_upper) %>%
    dplyr::mutate(
      is_historical = grepl("HISTORICAL|HIST|历史", Data_Source, ignore.case = TRUE),
      priority = ifelse(is_historical, 1, 2)
    ) %>%
    dplyr::arrange(env_upper, priority, desc(n_records)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(env_upper, default_Data_Source = Data_Source) %>%
    dplyr::mutate(default_env_id = paste(env_upper, default_Data_Source, sep = "__"))

  default_mapping <- default_sources
}

# planting env_id
if (!"env" %in% names(planting_info)) stop("planting_info 缺少 env 列")
planting_info$env_upper <- toupper(trimws(planting_info$env))
planting_info$env_id <- NA_character_
for (i in seq_len(nrow(planting_info))) {
  eu <- planting_info$env_upper[i]
  mids <- env_mapping$env_id[env_mapping$env_upper == eu]
  if (length(mids) == 0) {
    planting_info$env_id[i] <- eu
  } else if (length(mids) == 1) {
    planting_info$env_id[i] <- mids[1]
  } else {
    if (exists("default_mapping")) {
      did <- default_mapping$default_env_id[default_mapping$env_upper == eu]
      planting_info$env_id[i] <- ifelse(length(did) > 0, did[1], mids[1])
    } else {
      planting_info$env_id[i] <- mids[1]
    }
  }
}

# phenotype env_id
if ("env_code" %in% names(phenotype_data)) {
  phenotype_data$env_upper <- toupper(trimws(phenotype_data$env_code))
} else if ("env" %in% names(phenotype_data)) {
  phenotype_data$env_upper <- toupper(trimws(phenotype_data$env))
} else {
  env_col_candidates <- grep("env|Env|ENV", names(phenotype_data), value = TRUE, ignore.case = TRUE)
  if (length(env_col_candidates) == 0) stop("phenotype_data 找不到环境列")
  phenotype_data$env_upper <- toupper(trimws(phenotype_data[[env_col_candidates[1]]]))
}
phenotype_data$env_id <- NA_character_
for (i in seq_len(nrow(phenotype_data))) {
  eu <- phenotype_data$env_upper[i]
  mids <- env_mapping$env_id[env_mapping$env_upper == eu]
  if (length(mids) == 0) {
    phenotype_data$env_id[i] <- eu
  } else if (length(mids) == 1) {
    phenotype_data$env_id[i] <- mids[1]
  } else {
    if (exists("default_mapping")) {
      did <- default_mapping$default_env_id[default_mapping$env_upper == eu]
      phenotype_data$env_id[i] <- ifelse(length(did) > 0, did[1], mids[1])
    } else {
      phenotype_data$env_id[i] <- mids[1]
    }
  }
}

cat("唯一 env_id 数：weather=", length(unique(weather_data$env_id)),
    " planting=", length(unique(planting_info$env_id)),
    " phenotype=", length(unique(phenotype_data$env_id)), "\n")

write.csv(
  data.frame(
    data_source = c("weather", "planting", "phenotype"),
    unique_env_ids = c(length(unique(weather_data$env_id)),
                       length(unique(planting_info$env_id)),
                       length(unique(phenotype_data$env_id))),
    total_records = c(nrow(weather_data), nrow(planting_info), nrow(phenotype_data))
  ),
  file.path(output_dir, "tables", "env_id_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 6. 关键列名与 Zadok 提取
# ============================================================
cat("\n=== Step 3: 关键列 + Zadok ===\n")

phenology_stage_col <- if ("Wheat.Phenology.CurrentStageName" %in% names(weather_data)) {
  "Wheat.Phenology.CurrentStageName"
} else if ("CurrentStageName" %in% names(weather_data)) {
  "CurrentStageName"
} else if ("StageName" %in% names(weather_data)) {
  "StageName"
} else {
  NA_character_
}
if (is.na(phenology_stage_col)) stop("weather_data 找不到物候阶段列")

zadok_col <- if ("Wheat.Phenology.Zadok.Stage" %in% names(weather_data)) {
  "Wheat.Phenology.Zadok.Stage"
} else if ("zadok" %in% names(weather_data)) {
  "zadok"
} else if ("Zadok" %in% names(weather_data)) {
  "Zadok"
} else {
  NA_character_
}

# 使用 env_id 作为环境列
weather_data$env_standard <- trimws(as.character(weather_data$env_id))
planting_info$env_standard <- trimws(as.character(planting_info$env_id))
phenotype_data$env_standard <- trimws(as.character(phenotype_data$env_id))

# 日期列（你的数据里是 YYYYMMDD，但你原代码用 %Y/%m/%d；这里做双保险）
date_col_weather <- "YYYYMMDD"
if (!date_col_weather %in% names(weather_data)) stop("weather_data 缺少 YYYYMMDD 列")

parse_date_safe <- function(x) {
  x <- as.character(x)
  # 尝试多种格式
  d <- suppressWarnings(as.Date(x, format = "%Y/%m/%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%Y%m%d"))
  d
}
weather_data$date_clean <- parse_date_safe(weather_data[[date_col_weather]])
if (all(is.na(weather_data$date_clean))) stop("天气日期解析失败：检查 YYYYMMDD 格式")

# Zadok 数值
if (!is.na(zadok_col) && zadok_col %in% names(weather_data)) {
  weather_data$zadok_numeric <- suppressWarnings(as.numeric(stringr::str_extract(
    as.character(weather_data[[zadok_col]]), "[-+]?[0-9]*\\.?[0-9]+"
  )))
  weather_data$zadok_int <- round(weather_data$zadok_numeric)
} else {
  # 退化：从 stage name 推断
  weather_data <- weather_data %>%
    dplyr::mutate(
      zadok_int = dplyr::case_when(
        grepl("Emergence", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 10,
        grepl("Tillering", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 20,
        grepl("Stem", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 30,
        grepl("Booting", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 40,
        grepl("Heading", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 50,
        grepl("Anthesis|Flowering", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 65,
        grepl("GrainFill", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 70,
        grepl("Maturity|EndGrainFill", .data[[phenology_stage_col]], ignore.case = TRUE) ~ 90,
        TRUE ~ NA_real_
      ),
      zadok_numeric = as.numeric(zadok_int)
    )
}

# ============================================================
# 7. 提取开花/成熟（按 env_id）
# ============================================================
cat("\n=== Step 4: 提取开花/成熟 ===\n")

flowering_zadok_by_env <- weather_data %>%
  dplyr::filter(grepl("Anthesis|Flowering", .data[[phenology_stage_col]], ignore.case = TRUE)) %>%
  dplyr::filter(!is.na(zadok_numeric)) %>%
  dplyr::group_by(env_standard) %>%
  dplyr::arrange(date_clean) %>%
  dplyr::summarise(
    flowering_zadok = round(first(zadok_numeric)),
    flowering_date = first(date_clean),
    flowering_n_days = dplyr::n(),
    .groups = "drop"
  )

maturity_zadok_by_env <- weather_data %>%
  dplyr::filter(grepl("EndGrainFill|Maturity", .data[[phenology_stage_col]], ignore.case = TRUE)) %>%
  dplyr::filter(!is.na(zadok_numeric)) %>%
  dplyr::group_by(env_standard) %>%
  dplyr::arrange(date_clean) %>%
  dplyr::summarise(
    maturity_zadok = round(last(zadok_numeric)),
    maturity_date = last(date_clean),
    maturity_n_days = dplyr::n(),
    .groups = "drop"
  )

zadok_phenology_data <- dplyr::full_join(flowering_zadok_by_env, maturity_zadok_by_env, by = "env_standard") %>%
  dplyr::mutate(
    zadok_range = maturity_zadok - flowering_zadok,
    days_range = as.numeric(difftime(maturity_date, flowering_date, units = "days")),
    has_flowering = !is.na(flowering_zadok),
    has_maturity  = !is.na(maturity_zadok)
  )

write.csv(
  zadok_phenology_data,
  file.path(output_dir, "flowering_offset_analysis", "zadok_based_phenology_data_envid.csv"),
  row.names = FALSE
)

# ============================================================
# 8. 环境集合（共同环境 & 矩阵环境）
# ============================================================
cat("\n=== Step 5: 环境集合 ===\n")

weather_envs <- unique(weather_data$env_standard[!is.na(weather_data$env_standard)])
planting_envs <- unique(planting_info$env_standard[!is.na(planting_info$env_standard)])
phenotype_envs <- unique(phenotype_data$env_standard[!is.na(phenotype_data$env_standard)])

complete_zadok_envs <- zadok_phenology_data %>%
  dplyr::filter(has_flowering & has_maturity) %>%
  dplyr::pull(env_standard) %>% unique()

all_weather_envs_with_flowering <- unique(zadok_phenology_data$env_standard[zadok_phenology_data$has_flowering])
all_envs_for_matrix <- all_weather_envs_with_flowering

common_envs <- Reduce(intersect, list(weather_envs, planting_envs, phenotype_envs, complete_zadok_envs))
if (length(common_envs) == 0) common_envs <- Reduce(intersect, list(weather_envs, phenotype_envs, complete_zadok_envs))
if (length(common_envs) == 0) common_envs <- Reduce(intersect, list(weather_envs, phenotype_envs))
if (length(common_envs) == 0) stop("没有共同环境，无法继续")

write.csv(
  data.frame(environment_codes = common_envs),
  file.path(output_dir, "tables", "common_environment_codes_envid.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(environment_codes = all_envs_for_matrix),
  file.path(output_dir, "tables", "all_environments_for_matrix_envid.csv"),
  row.names = FALSE
)

# ============================================================
# 9. 清洗数据 + 表型均值
# ============================================================
cat("\n=== Step 6: 数据清洗 + 表型均值 ===\n")

weather_data_clean <- weather_data[weather_data$env_standard %in% common_envs, ]
zadok_phenology_data_clean <- zadok_phenology_data[zadok_phenology_data$env_standard %in% common_envs, ]

weather_data_all <- weather_data[weather_data$env_standard %in% all_envs_for_matrix, ]
zadok_phenology_data_all <- zadok_phenology_data[zadok_phenology_data$env_standard %in% all_envs_for_matrix, ]

# 表型均值（共同环境）
if (!"TKW" %in% names(phenotype_data)) stop("phenotype_data 缺少 TKW 列")
phenotype_data_clean <- phenotype_data[phenotype_data$env_standard %in% common_envs, ]

env_phenotype_mean <- phenotype_data_clean %>%
  dplyr::group_by(env_standard) %>%
  dplyr::summarise(
    TKW_mean = mean(TKW, na.rm = TRUE),
    TKW_sd = sd(TKW, na.rm = TRUE),
    TKW_n = dplyr::n(),
    TKW_min = min(TKW, na.rm = TRUE),
    TKW_max = max(TKW, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::rename(env_code = env_standard)

# 所有环境表型
all_env_phenotype <- phenotype_data %>%
  dplyr::group_by(env_standard) %>%
  dplyr::summarise(
    TKW_mean = mean(TKW, na.rm = TRUE),
    TKW_sd = sd(TKW, na.rm = TRUE),
    TKW_n = dplyr::n(),
    TKW_min = min(TKW, na.rm = TRUE),
    TKW_max = max(TKW, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::rename(env_code = env_standard)

write.csv(env_phenotype_mean, file.path(output_dir, "tables", "environment_phenotype_means_common_envid.csv"), row.names = FALSE)
write.csv(all_env_phenotype, file.path(output_dir, "tables", "environment_phenotype_means_all_envid.csv"), row.names = FALSE)

# ============================================================
# 10. 可用气象因子 + factor_types（包含所有 GDD 列）
# ============================================================
cat("\n=== Step 7: 因子准备 ===\n")

# 预定义可能因子（不含 GDD，后面自动识别）
new_possible_factors <- c(
  "TMEAN", "TMAX", "TMIN", "T2MDEW", "ASSD", "APAR", "DL",
  "PR", "QV2M", "RH", "WS2M", "VPD", "DTR"
)

# 从全量数据中识别所有 GDD 列（不区分大小写）
gdd_columns <- grep("^GDD", names(weather_data_all), value = TRUE, ignore.case = TRUE)
if (length(gdd_columns) > 0) {
  cat("检测到 GDD 列（基于全量数据）:", paste(gdd_columns, collapse = ", "), "\n")
  new_possible_factors <- unique(c(new_possible_factors, gdd_columns))
}

# 基于 weather_data_clean 确定可用因子（用于相关性分析，但矩阵构建会基于全量数据）
available_factors <- new_possible_factors[new_possible_factors %in% names(weather_data_clean)]
# 缺失因子尝试按相似列补齐（仅在 weather_data_clean 中）
missing_factors <- setdiff(new_possible_factors, available_factors)
if (length(missing_factors) > 0) {
  for (factor in missing_factors) {
    similar_cols <- grep(factor, names(weather_data_clean), ignore.case = TRUE, value = TRUE)
    if (length(similar_cols) > 0) {
      weather_data_clean[[factor]] <- weather_data_clean[[similar_cols[1]]]
      weather_data_all[[factor]]  <- weather_data_all[[similar_cols[1]]]
      available_factors <- c(available_factors, factor)
    }
  }
}
available_factors <- unique(available_factors)
cat("可用因子（用于相关性分析）:", paste(available_factors, collapse = ", "), "\n")

# 构建因子类型表（所有 GDD 因子同时计算 mean 和 sum）
factor_types <- data.frame(factor = character(), type = character(), stringsAsFactors = FALSE)
for (factor in available_factors) {
  if (grepl("^GDD", factor, ignore.case = TRUE)) {
    # 所有 GDD 因子都同时计算 mean 和 sum
    factor_types <- rbind(factor_types, data.frame(factor = factor, type = "mean", stringsAsFactors = FALSE))
    factor_types <- rbind(factor_types, data.frame(factor = factor, type = "sum", stringsAsFactors = FALSE))
  } else if (factor %in% c("TMEAN", "TMAX", "TMIN", "T2MDEW", "RH", "WS2M", "VPD",
                           "DTR", "APAR", "DL", "QV2M")) {
    factor_types <- rbind(factor_types, data.frame(factor = factor, type = "mean", stringsAsFactors = FALSE))
  } else if (factor %in% c("PR", "ASSD")) {
    factor_types <- rbind(factor_types, data.frame(factor = factor, type = "sum", stringsAsFactors = FALSE))
  } else {
    # 默认 mean
    factor_types <- rbind(factor_types, data.frame(factor = factor, type = "mean", stringsAsFactors = FALSE))
  }
}
# 去重（避免重复添加）
factor_types <- unique(factor_types)
cat("factor_types 行数:", nrow(factor_types), "\n")

# ============================================================
# 11. 参数 + 生成 offset_window_combinations
# ============================================================
cat("\n=== Step 8: 参数 + 窗口组合 ===\n")

min_window_size <- 5
max_window_size <- 25
pre_flowering   <- 10
post_flowering  <- 25
min_envs_for_correlation <- 5
min_days_per_window <- 3

offset_window_combinations <- expand.grid(
  start_offset = seq(-pre_flowering, post_flowering - min_window_size + 1),
  end_offset   = seq(-pre_flowering + min_window_size - 1, post_flowering)
)
offset_window_combinations$window_size <- offset_window_combinations$end_offset - offset_window_combinations$start_offset + 1
offset_window_combinations <- offset_window_combinations[
  offset_window_combinations$window_size >= min_window_size &
    offset_window_combinations$window_size <= max_window_size, ]

offset_window_combinations <- offset_window_combinations %>%
  dplyr::mutate(
    offset_label = ifelse(start_offset < 0 & end_offset < 0,
                          paste0("PreFlowering_", abs(start_offset), "_", abs(end_offset)),
                          ifelse(start_offset >= 0 & end_offset >= 0,
                                 paste0("PostFlowering", start_offset, "_", end_offset),
                                 paste0("CrossFlowering", abs(start_offset), "_", end_offset))),
    window_type = ifelse(start_offset < 0 & end_offset < 0, "pre_flowering",
                         ifelse(start_offset >= 0 & end_offset >= 0, "post_flowering", "cross_flowering")),
    window_description = ifelse(start_offset < 0 & end_offset < 0,
                                paste0("开花前", abs(start_offset), "至", abs(end_offset), "阶段"),
                                ifelse(start_offset >= 0 & end_offset >= 0,
                                       paste0("开花后", start_offset, "至", end_offset, "阶段"),
                                       paste0("开花前", abs(start_offset), "阶段至开花后", end_offset, "阶段")))
  )
offset_window_combinations$window_id <- seq_len(nrow(offset_window_combinations))

write.csv(offset_window_combinations,
          file.path(output_dir, "tables", "flowering_offset_window_combinations_envid.csv"),
          row.names = FALSE)

cat("窗口总数:", nrow(offset_window_combinations), "\n")

# ============================================================
# 12. 核心函数（与你版本一致：窗口因子 & 相关性）
# ============================================================
cat("\n=== Step 9: 核心函数 ===\n")

calculate_env_factors_flowering_offset_enhanced <- function(
  env_code, start_offset, end_offset, zadok_phenology_df, weather_df, factors
) {
  env_zadok_phenology <- zadok_phenology_df[zadok_phenology_df$env_standard == env_code, ]
  if (nrow(env_zadok_phenology) == 0 || is.na(env_zadok_phenology$flowering_zadok[1])) {
    return(list(summary=NULL, daily_data=NULL,
                quality_metrics=list(env_code=env_code, has_flowering=FALSE, n_days=NA,
                                     actual_window_size=NA, data_completeness=NA,
                                     window_too_small=FALSE, no_data_in_window=FALSE,
                                     insufficient_days=FALSE, has_valid_data=FALSE)))
  }

  flowering_zadok <- env_zadok_phenology$flowering_zadok[1]
  target_start_zadok <- flowering_zadok + start_offset
  target_end_zadok <- flowering_zadok + end_offset

  global_min_z <- min(weather_df$zadok_int, na.rm = TRUE)
  global_max_z <- max(weather_df$zadok_int, na.rm = TRUE)

  actual_start_zadok <- max(floor(target_start_zadok), global_min_z)
  actual_end_zadok <- min(ceiling(target_end_zadok), global_max_z)

  actual_window_size <- actual_end_zadok - actual_start_zadok + 1
  if (actual_window_size < min_window_size) {
    return(list(summary=NULL, daily_data=NULL,
                quality_metrics=list(env_code=env_code, has_flowering=TRUE, n_days=NA,
                                     actual_window_size=actual_window_size, data_completeness=NA,
                                     window_too_small=TRUE, no_data_in_window=FALSE,
                                     insufficient_days=FALSE, has_valid_data=FALSE)))
  }

  env_weather <- weather_df[
    weather_df$env_standard == env_code &
      !is.na(weather_df$zadok_int) &
      weather_df$zadok_int >= actual_start_zadok &
      weather_df$zadok_int <= actual_end_zadok, ]

  if (nrow(env_weather) == 0) {
    return(list(summary=NULL, daily_data=NULL,
                quality_metrics=list(env_code=env_code, has_flowering=TRUE, n_days=0,
                                     actual_window_size=actual_window_size, data_completeness=0,
                                     window_too_small=FALSE, no_data_in_window=TRUE,
                                     insufficient_days=FALSE, has_valid_data=FALSE)))
  }

  if (nrow(env_weather) < min_days_per_window) {
    return(list(summary=NULL, daily_data=NULL,
                quality_metrics=list(env_code=env_code, has_flowering=TRUE, n_days=nrow(env_weather),
                                     actual_window_size=actual_window_size,
                                     data_completeness=nrow(env_weather)/actual_window_size,
                                     window_too_small=FALSE, no_data_in_window=FALSE,
                                     insufficient_days=TRUE, has_valid_data=FALSE)))
  }

  window_size_planned <- end_offset - start_offset + 1
  window_size_actual <- actual_window_size

  if (start_offset < 0 && end_offset < 0) {
    offset_label <- paste0("PreFlowering_", abs(start_offset), "_", abs(end_offset))
    window_type <- "pre_flowering"
  } else if (start_offset >= 0 && end_offset >= 0) {
    offset_label <- paste0("PostFlowering", start_offset, "_", end_offset)
    window_type <- "post_flowering"
  } else {
    offset_label <- paste0("CrossFlowering", abs(start_offset), "_", end_offset)
    window_type <- "cross_flowering"
  }

  actual_zadok_label <- paste0("Z", actual_start_zadok, "_", actual_end_zadok)

  result_summary <- list(
    env_code = env_code,
    env_standard = env_code,
    start_offset = start_offset,
    end_offset = end_offset,
    offset_label = offset_label,
    window_type = window_type,
    window_size_planned = window_size_planned,
    flowering_zadok = flowering_zadok,
    target_start_zadok = target_start_zadok,
    target_end_zadok = target_end_zadok,
    actual_start_zadok = actual_start_zadok,
    actual_end_zadok = actual_end_zadok,
    actual_zadok_label = actual_zadok_label,
    window_size_actual = window_size_actual,
    n_days = nrow(env_weather),
    min_zadok_in_window = min(env_weather$zadok_int, na.rm = TRUE),
    max_zadok_in_window = max(env_weather$zadok_int, na.rm = TRUE),
    start_date = min(env_weather$date_clean, na.rm = TRUE),
    end_date = max(env_weather$date_clean, na.rm = TRUE),
    data_completeness = nrow(env_weather) / window_size_actual
  )

  for (factor in factors) {
    if (!factor %in% names(env_weather)) {
      result_summary[[paste0(factor, "_mean")]] <- NA
      result_summary[[paste0(factor, "_sum")]] <- NA
      result_summary[[paste0(factor, "_sd")]] <- NA
      result_summary[[paste0(factor, "_n")]] <- 0
      result_summary[[paste0(factor, "_completeness")]] <- 0
      next
    }

    values <- env_weather[[factor]]
    valid_values <- values[!is.na(values)]
    n_valid <- length(valid_values)

    result_summary[[paste0(factor, "_n")]] <- n_valid
    result_summary[[paste0(factor, "_completeness")]] <- n_valid / nrow(env_weather)

    if (n_valid < min_days_per_window) {
      result_summary[[paste0(factor, "_mean")]] <- NA
      result_summary[[paste0(factor, "_sum")]] <- NA
      result_summary[[paste0(factor, "_sd")]] <- NA
    } else {
      result_summary[[paste0(factor, "_mean")]] <- mean(valid_values)
      result_summary[[paste0(factor, "_sum")]] <- sum(valid_values)
      result_summary[[paste0(factor, "_sd")]] <- sd(valid_values)
    }
  }

  quality_metrics <- list(
    env_code = env_code, has_flowering = TRUE, n_days = nrow(env_weather),
    actual_window_size = actual_window_size,
    data_completeness = nrow(env_weather) / window_size_actual,
    window_too_small = FALSE, no_data_in_window = FALSE,
    insufficient_days = FALSE, has_valid_data = TRUE
  )

  return(list(summary = as.data.frame(result_summary, stringsAsFactors = FALSE),
              daily_data = NULL,
              quality_metrics = quality_metrics))
}

calculate_correlations_flowering_offset_enhanced <- function(env_factors_df, phenotype_mean_df, factor_types_df,
                                                            min_envs = min_envs_for_correlation) {
  combined <- merge(env_factors_df, phenotype_mean_df, by = "env_code", all.x = TRUE)
  combined <- combined[!is.na(combined$TKW_mean), ]
  if (nrow(combined) < min_envs) return(NULL)

  start_offset_u <- unique(env_factors_df$start_offset)
  end_offset_u <- unique(env_factors_df$end_offset)
  offset_label_u <- unique(env_factors_df$offset_label)
  window_type_u <- unique(env_factors_df$window_type)
  window_size_planned_u <- unique(env_factors_df$window_size_planned)
  if (length(start_offset_u) != 1 || length(end_offset_u) != 1 || length(offset_label_u) != 1) return(NULL)

  actual_start_min <- min(env_factors_df$actual_start_zadok, na.rm = TRUE)
  actual_start_max <- max(env_factors_df$actual_start_zadok, na.rm = TRUE)
  actual_end_min <- min(env_factors_df$actual_end_zadok, na.rm = TRUE)
  actual_end_max <- max(env_factors_df$actual_end_zadok, na.rm = TRUE)
  n_actual_labels <- length(unique(env_factors_df$actual_zadok_label))
  label_table <- sort(table(env_factors_df$actual_zadok_label), decreasing = TRUE)
  actual_label_mode <- ifelse(length(label_table) > 0, names(label_table)[1], NA)
  actual_label_mode_n <- ifelse(length(label_table) > 0, as.integer(label_table[1]), NA)

  mean_flowering_zadok <- mean(env_factors_df$flowering_zadok, na.rm = TRUE)
  sd_flowering_zadok <- sd(env_factors_df$flowering_zadok, na.rm = TRUE)

  results <- list()

  for (i in seq_len(nrow(factor_types_df))) {
    factor <- factor_types_df$factor[i]
    ftype <- factor_types_df$type[i]
    col_name <- if (ftype == "mean") paste0(factor, "_mean") else paste0(factor, "_sum")
    if (!col_name %in% names(combined)) next

    completeness_col <- paste0(factor, "_completeness")
    if (completeness_col %in% names(combined)) {
      valid_data <- combined[!is.na(combined[[col_name]]) & !is.na(combined$TKW_mean) &
                               combined[[completeness_col]] > 0.5, ]
    } else {
      valid_data <- combined[!is.na(combined[[col_name]]) & !is.na(combined$TKW_mean), ]
    }
    if (nrow(valid_data) < min_envs) next

    # 去 3-sigma
    factor_values <- valid_data[[col_name]]
    mu <- mean(factor_values, na.rm = TRUE)
    sdv <- sd(factor_values, na.rm = TRUE)
    if (!is.na(sdv) && sdv > 0) {
      z <- abs((factor_values - mu) / sdv)
      valid_data <- valid_data[z <= 3, ]
    }
    if (nrow(valid_data) < min_envs) next

    tryCatch({
      cor_test <- cor.test(valid_data[[col_name]], valid_data$TKW_mean,
                           method = "pearson", use = "complete.obs")
      lm_fit <- lm(TKW_mean ~ value, data = data.frame(
        TKW_mean = valid_data$TKW_mean,
        value = valid_data[[col_name]]
      ))
      ci <- confint(lm_fit, "value", level = 0.95)
      n_envs_original <- nrow(combined[!is.na(combined[[col_name]]) & !is.na(combined$TKW_mean), ])
      n_envs_after_qc <- nrow(valid_data)

      results[[length(results) + 1]] <- data.frame(
        start_offset = start_offset_u,
        end_offset = end_offset_u,
        offset_label = offset_label_u,
        window_type = window_type_u,
        window_size_planned = window_size_planned_u,
        actual_start_min = actual_start_min,
        actual_start_max = actual_start_max,
        actual_end_min = actual_end_min,
        actual_end_max = actual_end_max,
        n_unique_actual_labels = n_actual_labels,
        actual_label_mode = actual_label_mode,
        actual_label_mode_n = actual_label_mode_n,
        mean_flowering_zadok = mean_flowering_zadok,
        sd_flowering_zadok = sd_flowering_zadok,
        factor = factor,
        value_type = ftype,
        correlation = as.numeric(cor_test$estimate),
        p_value = cor_test$p.value,
        n_envs_original = n_envs_original,
        n_envs_after_qc = n_envs_after_qc,
        r_squared = summary(lm_fit)$r.squared,
        slope = coef(lm_fit)[2],
        slope_ci_lower = ci[1],
        slope_ci_upper = ci[2],
        intercept = coef(lm_fit)[1],
        mean_data_completeness = if (completeness_col %in% names(valid_data)) {
          mean(valid_data[[completeness_col]], na.rm = TRUE)
        } else NA_real_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {})
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# ============================================================
# 13. 穷举窗口相关性（共同环境）
# ============================================================
cat("\n=== Step 10: 穷举窗口相关性（共同环境）===\n")

all_results_offset <- list()
all_window_factor_values_offset <- list()
all_quality_metrics_offset <- list()

total_windows <- nrow(offset_window_combinations)
progress_interval <- max(1, floor(total_windows / 20))
cat("开始分析窗口数:", total_windows, "\n")

for (i in seq_len(total_windows)) {
  win <- offset_window_combinations[i, ]
  start_offset <- win$start_offset
  end_offset   <- win$end_offset
  window_id    <- win$window_id
  offset_label <- win$offset_label
  window_type  <- win$window_type
  window_desc  <- win$window_description

  window_summaries_list <- list()
  quality_metrics_list <- list()

  for (env_code in common_envs) {
    factor_data <- calculate_env_factors_flowering_offset_enhanced(
      env_code = env_code,
      start_offset = start_offset,
      end_offset = end_offset,
      zadok_phenology_df = zadok_phenology_data_clean,
      weather_df = weather_data_clean,
      factors = available_factors
    )

    if (!is.null(factor_data$summary)) {
      factor_data$summary$window_id <- window_id
      factor_data$summary$offset_label <- offset_label
      factor_data$summary$window_type <- window_type
      factor_data$summary$window_description <- window_desc
      window_summaries_list[[env_code]] <- factor_data$summary
    }

    if (!is.null(factor_data$quality_metrics)) {
      qm <- factor_data$quality_metrics
      qm$window_id <- window_id
      qm$offset_label <- offset_label
      qm$window_type <- window_type
      qm$start_offset <- start_offset
      qm$end_offset <- end_offset
      quality_metrics_list[[env_code]] <- qm
    }
  }

  if (length(window_summaries_list) > 0) {
    window_summaries_df <- dplyr::bind_rows(window_summaries_list)

    correlations <- calculate_correlations_flowering_offset_enhanced(
      env_factors_df = window_summaries_df,
      phenotype_mean_df = env_phenotype_mean,
      factor_types_df = factor_types
    )
    if (!is.null(correlations)) {
      correlations$window_id <- window_id
      correlations$window_description <- window_desc
      all_results_offset[[i]] <- correlations
    }

    window_factor_values <- window_summaries_df %>%
      dplyr::left_join(env_phenotype_mean, by = "env_code") %>%
      dplyr::mutate(
        window_id = window_id,
        offset_label = offset_label,
        window_type = window_type,
        window_description = window_desc
      )
    all_window_factor_values_offset[[i]] <- window_factor_values
  }

  if (length(quality_metrics_list) > 0) {
    quality_metrics_df <- data.table::rbindlist(
      lapply(quality_metrics_list, function(x) as.data.frame(as.list(x), stringsAsFactors = FALSE)),
      fill = TRUE
    )
    all_quality_metrics_offset[[i]] <- as.data.frame(quality_metrics_df)
  }

  if (i %% progress_interval == 0 || i == total_windows) {
    progress <- round(i / total_windows * 100, 1)
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("进度: %d/%d (%.1f%%) | 已运行: %.1f分钟 | 当前窗口: %s\n",
                i, total_windows, progress, elapsed, offset_label))
  }
}

# 保存相关性结果
if (length(all_results_offset) > 0) {
  results_offset_df <- dplyr::bind_rows(all_results_offset) %>%
    dplyr::mutate(
      abs_correlation = abs(correlation),
      log10_pvalue = -log10(p_value),
      fdr = p.adjust(p_value, method = "fdr"),
      bonferroni = p.adjust(p_value, method = "bonferroni"),
      significance = dplyr::case_when(
        bonferroni < 0.001 ~ "***",
        bonferroni < 0.01 ~ "**",
        bonferroni < 0.05 ~ "*",
        TRUE ~ ""
      ),
      window_label_composite = paste0(
        offset_label, " | actual(Zstart_min-max=",
        actual_start_min, "-", actual_start_max,
        ", Zend_min-max=", actual_end_min, "-", actual_end_max,
        ", mode=", actual_label_mode, "(", actual_label_mode_n, " envs), ",
        "n_unique_actual=", n_unique_actual_labels, ")"
      )
    ) %>% dplyr::arrange(desc(abs_correlation))

  write.csv(results_offset_df,
            file.path(output_dir, "tables", "【完整相关性结果】flowering_offset_correlation_results_common_envid.csv"),
            row.names = FALSE)

  write.csv(head(results_offset_df[order(-results_offset_df$abs_correlation), ], 100),
            file.path(output_dir, "tables", "【完整相关性结果】top_100_flowering_offset_correlations_common_envid.csv"),
            row.names = FALSE)

  cat("相关性条数:", nrow(results_offset_df), "\n")
}

# 保存窗口汇总数据
if (length(all_window_factor_values_offset) > 0) {
  all_factor_values_offset_complete <- dplyr::bind_rows(all_window_factor_values_offset)
  write.csv(all_factor_values_offset_complete,
            file.path(output_dir, "window_values", "【所有窗口汇总数据】all_flowering_offset_window_factor_values_common_envid.csv"),
            row.names = FALSE)
}

# 保存质量指标
if (length(all_quality_metrics_offset) > 0) {
  all_quality_metrics_complete <- data.table::rbindlist(all_quality_metrics_offset, fill = TRUE)
  write.csv(all_quality_metrics_complete,
            file.path(output_dir, "data_quality", "window_quality_metrics_common_envid.csv"),
            row.names = FALSE)
}

# ============================================================
# 14. 关键加速：构建全环境矩阵（chunk + 并行 foreach + 分块输出写盘）
# ============================================================
cat("\n=== Step 11: 构建【环境×开花偏移窗口矩阵】（加速版）===\n")

# 输出文件名（你指定的“正确写法”）
out_matrix_file <- file.path(
  output_dir, "tables",
  "【环境×开花偏移窗口矩阵env_id版】env_id_flowering_offset_window_matrix_all.csv"
)
if (file.exists(out_matrix_file)) file.remove(out_matrix_file)

# 并行参数
n_workers <- max(1, parallel::detectCores() - 1)
chunk_size <- 20
cat("并行 workers:", n_workers, " | chunk_size:", chunk_size, "\n")

# data.table 预处理 + 按 env 分组缓存
wdt <- as.data.table(weather_data_all)
zdt <- as.data.table(zadok_phenology_data_all)

# 仅保留必要列（减少内存与拷贝）
need_cols <- unique(c("env_standard", "zadok_int", "date_clean", available_factors))
wdt <- wdt[, ..need_cols]

setkeyv(wdt, c("env_standard", "zadok_int"))
setkey(zdt, env_standard)

weather_by_env <- split(wdt, by = "env_standard", keep.by = FALSE, drop = TRUE)
phen_by_env    <- split(zdt, by = "env_standard", keep.by = FALSE, drop = TRUE)

# 预生成全部矩阵列名（稳定顺序）
all_matrix_cols <- character(0)
for (win_idx in seq_len(nrow(offset_window_combinations))) {
  olab <- offset_window_combinations$offset_label[win_idx]
  for (ii in seq_len(nrow(factor_types))) {
    f  <- factor_types$factor[ii]
    vt <- factor_types$type[ii]
    all_matrix_cols <- c(all_matrix_cols, paste0(f, "&", vt, "&", olab))
  }
}
cat("矩阵因子列数:", length(all_matrix_cols), "\n")

# 表型映射（允许 NA）
ph <- as.data.table(all_env_phenotype)
if ("env_code" %in% names(ph)) setnames(ph, "env_code", "env")
if (!("env" %in% names(ph) && "TKW_mean" %in% names(ph))) stop("all_env_phenotype 需要 env(或 env_code) 与 TKW_mean")
ph <- unique(ph[, .(env, TKW_mean)])

# 写 header
header_dt <- data.table(env = character(0))
for (cn in all_matrix_cols) header_dt[[cn]] <- numeric(0)
header_dt[, TKW_mean := numeric(0)]
fwrite(header_dt, out_matrix_file)

# 单环境计算一行（fast）
calc_one_env_row_fast <- function(env_code) {
  env_weather <- weather_by_env[[env_code]]
  env_pheno  <- phen_by_env[[env_code]]

  # 默认全 NA
  out_vec <- rep(NA_real_, length(all_matrix_cols))
  names(out_vec) <- all_matrix_cols

  if (is.null(env_weather) || is.null(env_pheno) || nrow(env_pheno) == 0 || is.na(env_pheno$flowering_zadok[1])) {
    row <- as.list(out_vec); row$env <- env_code
    return(as.data.table(row)[, c("env", all_matrix_cols), with = FALSE])
  }

  flowering_z <- env_pheno$flowering_zadok[1]
  zmin <- suppressWarnings(min(env_weather$zadok_int, na.rm = TRUE))
  zmax <- suppressWarnings(max(env_weather$zadok_int, na.rm = TRUE))
  if (!is.finite(zmin) || !is.finite(zmax)) {
    row <- as.list(out_vec); row$env <- env_code
    return(as.data.table(row)[, c("env", all_matrix_cols), with = FALSE])
  }

  for (win_idx in seq_len(nrow(offset_window_combinations))) {
    st_off <- offset_window_combinations$start_offset[win_idx]
    ed_off <- offset_window_combinations$end_offset[win_idx]
    olab   <- offset_window_combinations$offset_label[win_idx]

    target_start <- flowering_z + st_off
    target_end   <- flowering_z + ed_off
    actual_start <- max(floor(target_start), zmin)
    actual_end   <- min(ceiling(target_end), zmax)
    actual_size  <- actual_end - actual_start + 1
    if (actual_size < min_window_size) next

    wsub <- env_weather[zadok_int >= actual_start & zadok_int <= actual_end]
    if (nrow(wsub) < min_days_per_window) next

    for (ii in seq_len(nrow(factor_types))) {
      f  <- factor_types$factor[ii]
      vt <- factor_types$type[ii]
      if (!f %in% names(wsub)) next

      vals <- wsub[[f]]
      vals <- vals[!is.na(vals)]
      if (length(vals) < min_days_per_window) next

      dst <- paste0(f, "&", vt, "&", olab)
      if (vt == "mean") out_vec[dst] <- mean(vals) else out_vec[dst] <- sum(vals)
    }
  }

  row <- as.list(out_vec); row$env <- env_code
  as.data.table(row)[, c("env", all_matrix_cols), with = FALSE]
}

# chunk + 并行 + 分块写盘
envs <- as.character(all_envs_for_matrix)
envs <- envs[!is.na(envs) & envs != ""]
env_chunks <- split(envs, ceiling(seq_along(envs) / chunk_size))
cat("总环境:", length(envs), " | chunks:", length(env_chunks), "\n")

future::plan(future::multisession, workers = n_workers)

matrix_start <- Sys.time()
for (ck in seq_along(env_chunks)) {
  chunk_envs <- env_chunks[[ck]]

  chunk_rows <- future.apply::future_lapply(
    chunk_envs,
    FUN = calc_one_env_row_fast,
    future.seed = TRUE
  )

  chunk_dt <- data.table::rbindlist(chunk_rows, fill = TRUE)

  # 补缺列、顺序固定
  miss_cols <- setdiff(c("env", all_matrix_cols), names(chunk_dt))
  if (length(miss_cols) > 0) for (m in miss_cols) chunk_dt[[m]] <- NA_real_
  setcolorder(chunk_dt, c("env", all_matrix_cols))

  # 加表型
  chunk_dt <- merge(chunk_dt, ph, by = "env", all.x = TRUE, sort = FALSE)

  # 分块追加写盘
  fwrite(chunk_dt, out_matrix_file, append = TRUE, col.names = FALSE)

  elapsed <- round(as.numeric(difftime(Sys.time(), matrix_start, units = "mins")), 2)
  cat(sprintf("chunk %d/%d 完成 | env=%d | 已用 %.2f 分钟 | 写入 %s\n",
              ck, length(env_chunks), length(chunk_envs), elapsed, basename(out_matrix_file)))
}
future::plan(future::sequential)

cat("矩阵构建完成：", out_matrix_file, "\n")

# ============================================================
# 15. 列名映射 + env_id 分解（可选但建议保留）
# ============================================================
cat("\n=== Step 12: 列名映射 + env_id 分解 ===\n")

# 列名映射
env_matrix_cols <- c("env", all_matrix_cols, "TKW_mean")
column_mapping_all <- data.frame(
  column_name = env_matrix_cols,
  description = "",
  factor = "",
  value_type = "",
  offset_label = "",
  stringsAsFactors = FALSE
)

for (k in seq_len(nrow(column_mapping_all))) {
  cn <- column_mapping_all$column_name[k]
  if (cn == "env") {
    column_mapping_all$description[k] <- "环境ID（env__Data_Source）"
  } else if (cn == "TKW_mean") {
    column_mapping_all$description[k] <- "千粒重表型均值（可能为NA）"
  } else {
    parts <- strsplit(cn, "&", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      column_mapping_all$factor[k] <- parts[1]
      column_mapping_all$value_type[k] <- parts[2]
      column_mapping_all$offset_label[k] <- parts[3]
      column_mapping_all$description[k] <- if (parts[2] == "mean") {
        paste0(parts[1], " 在窗口 ", parts[3], " 的平均值")
      } else {
        paste0(parts[1], " 在窗口 ", parts[3], " 的总和")
      }
    }
  }
}
write.csv(column_mapping_all,
          file.path(output_dir, "tables", "column_name_mapping_for_flowering_offset_matrix_all_envid.csv"),
          row.names = FALSE)

# env_id 组成分解
env_id_components <- data.frame(
  env_id = all_envs_for_matrix,
  original_env = "",
  data_source = "",
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(env_id_components))) {
  eid <- env_id_components$env_id[i]
  parts <- strsplit(eid, "__", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    env_id_components$original_env[i] <- parts[1]
    env_id_components$data_source[i] <- paste(parts[-1], collapse = "__")
  } else {
    env_id_components$original_env[i] <- eid
    env_id_components$data_source[i] <- "UNKNOWN"
  }
}
write.csv(env_id_components, file.path(output_dir, "tables", "env_id_components_breakdown.csv"), row.names = FALSE)

# ============================================================
# 16. 报告 + 保存工作空间
# ============================================================
cat("\n=== Step 13: 报告 + 保存 ===\n")

n_total_correlations <- if (exists("results_offset_df")) nrow(results_offset_df) else 0
max_correlation <- if (exists("results_offset_df")) max(results_offset_df$abs_correlation, na.rm = TRUE) else NA
best_factor <- if (exists("results_offset_df")) results_offset_df$factor[which.max(results_offset_df$abs_correlation)] else NA
best_window <- if (exists("results_offset_df")) results_offset_df$offset_label[which.max(results_offset_df$abs_correlation)] else NA

analysis_report_offset <- data.frame(
  Category = c(
    "分析开始时间", "分析结束时间", "总运行时间（分钟）",
    "所有天气环境数", "有开花的环境数", "种植信息环境数",
    "表型数据环境数", "完整物候环境数", "共同环境数（相关性分析）",
    "所有用于矩阵的环境数",
    "总开花偏移窗口组合数",
    "计算的相关性数",
    "最强相关系数(abs)", "最强相关因子", "最强相关窗口(offset_label)",
    "矩阵输出文件名"
  ),
  Value = c(
    format(start_time, "%Y-%m-%d %H:%M:%S"),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2),
    length(weather_envs),
    length(all_weather_envs_with_flowering),
    length(planting_envs),
    length(phenotype_envs),
    length(complete_zadok_envs),
    length(common_envs),
    length(all_envs_for_matrix),
    nrow(offset_window_combinations),
    n_total_correlations,
    round(max_correlation, 4),
    best_factor,
    best_window,
    basename(out_matrix_file)
  ),
  stringsAsFactors = FALSE
)

write.csv(analysis_report_offset,
          file.path(output_dir, "tables", "flowering_offset_analysis_report_complete_envid.csv"),
          row.names = FALSE)

save.image(file.path(output_dir, "flowering_offset_analysis_workspace_complete_envid.RData"))

cat("\n=== 全流程完成 ===\n")
cat("输出目录:", output_dir, "\n")
cat("矩阵文件（正确文件名）:\n  ", out_matrix_file, "\n")
cat("总用时(分钟):", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "\n")
