# 设置工作目录
setwd("C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\001历史气象数据和未来气象数据的整理和met文件的整理和apsim对气象的分割\\01APSIMout")

# 加载必要的库
library(dplyr)

# 第一步：对每个CSV文件进行去重处理
cat("=== 第一步：处理重复数据 ===\n")

# 获取目录下所有的.csv 文件
csv_files <- list.files(pattern = "\\.csv$")

if (length(csv_files) == 0) {
  cat("没有找到CSV文件。\n")
} else {
  cat(paste0("找到 ", length(csv_files), " 个CSV文件:\n"))
  for (file in csv_files) {
    cat("  - ", file, "\n")
  }
  cat("\n")
}

# 遍历每个.csv 文件进行去重处理
for (file in csv_files) {
  cat("处理文件:", file, "...\n")
  
  tryCatch({
    # 读取文件
    data <- read.csv(file)
    original_rows <- nrow(data)
    
    # 删除所有重复行
    data <- unique(data)
    unique_rows <- nrow(data)
    removed_complete_duplicates <- original_rows - unique_rows
    
    # 检查是否有 Date 列
    if ("Date" %in% names(data)) {
      # 找出 Date 列重复的日期
      duplicate_dates <- unique(data$Date[duplicated(data$Date)])
      
      if (length(duplicate_dates) > 0) {
        cat("  发现", length(duplicate_dates), "个重复日期\n")
        
        rows_removed <- 0
        # 处理每个重复日期
        for (date in duplicate_dates) {
          # 提取该日期对应的所有行
          date_subset <- data[data$Date == date, ]
          
          # 计算每行的非空数据个数
          non_na_counts <- rowSums(!is.na(date_subset))
          
          # 找出非空数据个数最多的行的索引
          max_count_index <- which.max(non_na_counts)
          
          # 如果有多个行的非空数据个数相同且都是最多的，取第一个
          if (length(max_count_index) > 1) {
            max_count_index <- max_count_index[1]
          }
          
          # 提取要保留的行
          row_to_keep <- date_subset[max_count_index, ]
          
          # 从原数据中删除该日期对应的所有行
          data <- data[data$Date != date, ]
          
          # 将保留的行添加回原数据
          data <- rbind(data, row_to_keep)
          
          rows_removed <- rows_removed + (nrow(date_subset) - 1)
        }
        cat("  移除了", rows_removed, "个重复日期行\n")
      } else {
        cat("  没有发现重复日期\n")
      }
    } else {
      cat("  警告：文件中没有Date列\n")
    }
    
    # 保存处理后的数据到原文件
    write.csv(data, file = file, row.names = FALSE)
    
    final_rows <- nrow(data)
    cat("  处理完成: 原始", original_rows, "行 -> 最终", final_rows, "行\n")
    cat("  共移除重复行:", original_rows - final_rows, "\n\n")
    
  }, error = function(e) {
    cat("  错误处理文件 ", file, ": ", e$message, "\n", sep = "")
  })
}

# 第二步：合并所有处理后的CSV文件
cat("=== 第二步：合并CSV文件 ===\n")

# 定义函数用于合并 CSV 文件
combine_csv_files <- function(pattern, output_dir = "02APSIMcombined") {
  # 查找所有匹配模式的文件
  file_paths <- list.files(pattern = pattern)
  
  if (length(file_paths) == 0) {
    cat("没有找到匹配的文件。\n")
    return(NULL)
  }
  
  cat(paste0("找到 ", length(file_paths), " 个文件需要合并:\n"))
  for (file_path in file_paths) {
    cat("  - ", file_path, "\n")
  }
  
  # 初始化一个空列表来存储数据框
  dfs <- list()
  
  # 循环读取每个文件并添加到列表中
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    
    tryCatch({
      df <- read.csv(file_path)
      
      # 从文件名中提取信息
      file_name <- basename(file_path)
      base_name <- sub("\\.csv$", "", file_name)
      
      # 添加新列 'file_source'
      df$file_source <- base_name
      
      # 如果文件名符合特定模式，提取更多信息
      if (grepl("_simulation_output", base_name)) {
        # 提取位置和年份信息
        location_year <- sub("_simulation_output", "", base_name)
        df$location_year <- location_year
        
        # 尝试提取年份（假设年份在末尾或开头）
        year_match <- regmatches(location_year, regexpr("\\d{4}", location_year))
        if (length(year_match) > 0) {
          df$year <- as.numeric(year_match[1])
        }
      }
      
      # 将数据框添加到列表中
      dfs[[i]] <- df
      
      cat("  已读取文件: ", file_path, " (", nrow(df), " 行, ", ncol(df), " 列)\n", sep = "")
    }, error = function(e) {
      cat("  错误: 无法读取文件 ", file_path, ": ", e$message, "\n", sep = "")
    })
  }
  
  # 检查是否有成功读取的文件
  valid_dfs <- dfs[!sapply(dfs, is.null)]
  if (length(valid_dfs) == 0) {
    cat("没有成功读取任何文件。\n")
    return(NULL)
  }
  
  # 合并所有数据框
  combined_df <- bind_rows(valid_dfs)
  
  # 删除完全重复的行
  combined_df <- distinct(combined_df)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("已创建输出目录: ", output_dir, "\n", sep = "")
  }
  
  # 保存合并后的数据框到CSV文件
  output_file <- file.path(output_dir, "combined_simulation_output.csv")
  write.csv(combined_df, output_file, row.names = FALSE)
  
  # 同时保存为TXT文件（制表符分隔）
  txt_file <- file.path(output_dir, "combined_simulation_output.txt")
  write.table(combined_df, txt_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n合并完成:\n")
  cat("  输出CSV文件: ", output_file, "\n", sep = "")
  cat("  输出TXT文件: ", txt_file, "\n", sep = "")
  cat("  总行数: ", nrow(combined_df), "\n", sep = "")
  cat("  总列数: ", ncol(combined_df), "\n", sep = "")
  
  # 统计信息
  if ("file_source" %in% names(combined_df)) {
    cat("  包含的文件来源: ", length(unique(combined_df$file_source)), " 个\n", sep = "")
    cat("  各文件行数统计:\n")
    source_counts <- table(combined_df$file_source)
    for (source in names(source_counts)) {
      cat("    - ", source, ": ", source_counts[source], " 行\n", sep = "")
    }
  }
  
  # 返回合并后的数据框
  return(combined_df)
}

# 调用合并函数（使用特定模式）
cat("\n开始合并文件...\n")
result <- combine_csv_files("*_simulation_output.csv", output_dir = "02APSIMcombined")

# 生成详细统计报告
if (!is.null(result)) {
  cat("\n=== 第三步：生成统计报告 ===\n")
  
  # 创建详细统计文件
  stats_file <- file.path("02APSIMcombined", "detailed_statistics.txt")
  sink(stats_file)
  
  cat("APSIM模拟数据合并统计报告\n")
  cat("==========================\n")
  cat("报告生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("工作目录: ", getwd(), "\n")
  cat("输出目录: 02APSIMcombined\n")
  
  cat("\n数据概览:\n")
  cat("  总行数: ", nrow(result), "\n")
  cat("  总列数: ", ncol(result), "\n")
  
  cat("\n列信息:\n")
  for (i in seq_along(names(result))) {
    col_name <- names(result)[i]
    col_type <- class(result[[col_name]])[1]
    na_count <- sum(is.na(result[[col_name]]))
    na_percent <- round(na_count / nrow(result) * 100, 2)
    
    cat(sprintf("  %2d. %-20s %-10s NA: %4d (%5.1f%%)\n", 
                i, col_name, col_type, na_count, na_percent))
  }
  
  if ("file_source" %in% names(result)) {
    cat("\n文件来源统计:\n")
    sources <- unique(result$file_source)
    for (source in sort(sources)) {
      source_rows <- sum(result$file_source == source)
      source_percent <- round(source_rows / nrow(result) * 100, 2)
      cat("  - ", source, ": ", source_rows, " 行 (", source_percent, "%)\n", sep = "")
    }
  }
  
  if ("Date" %in% names(result)) {
    cat("\n日期范围:\n")
    dates <- as.Date(result$Date)
    cat("  最早日期: ", min(dates, na.rm = TRUE), "\n")
    cat("  最晚日期: ", max(dates, na.rm = TRUE), "\n")
    cat("  日期数量: ", length(unique(dates)), "\n")
  }
  
  # 数值列的统计摘要
  numeric_cols <- names(result)[sapply(result, is.numeric)]
  if (length(numeric_cols) > 0) {
    cat("\n数值列统计摘要:\n")
    for (col in numeric_cols) {
      if (col != "year") {  # 跳过年份列
        cat("  ", col, ":\n", sep = "")
        cat("    最小值: ", round(min(result[[col]], na.rm = TRUE), 2), "\n")
        cat("    最大值: ", round(max(result[[col]], na.rm = TRUE), 2), "\n")
        cat("    平均值: ", round(mean(result[[col]], na.rm = TRUE), 2), "\n")
        cat("    中位数: ", round(median(result[[col]], na.rm = TRUE), 2), "\n")
        cat("    标准差: ", round(sd(result[[col]], na.rm = TRUE), 2), "\n")
      }
    }
  }
  
  sink()
  cat("详细统计报告已保存到: ", stats_file, "\n", sep = "")
  
  # 显示前几行数据预览
  cat("\n合并数据预览（前3行）:\n")
  print(head(result, 3))
}

cat("\n=== 处理完成 ===\n")
cat("所有CSV文件已去重并合并到 '02APSIMcombined' 目录。\n")