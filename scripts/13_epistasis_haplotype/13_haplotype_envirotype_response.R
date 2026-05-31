#!/usr/bin/env Rscript

# ==================== 0. 初始设置 ====================
cat("========== Thousand Kernel Weight Haplotype and Environment Response Analysis ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# 设置基础目录
base_dir <- "/mnt/7t_storage/zhangcl/TKW"

# 设置标签参数（用于图形标注）
label_params <- list(
  label_size = 6,
  label_fontface = "bold",
  label_hjust = -0.5
)

# ==================== 1. 加载必需的R包 ====================
cat("\n1. Loading required R packages...\n")

required_packages <- c(
  "vcfR",      # 用于读取VCF文件
  "ggplot2",   # 用于绘图
  "RColorBrewer", # 用于颜色调色板
  "cowplot",   # 用于图形组合
  "dplyr",     # 用于数据处理
  "tidyr",     # 用于数据整理
  "viridis",   # 用于更多颜色选项
  "scales"     # 用于更好的图形比例
)

# 安装缺失的包
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

# 加载所有包
cat("Loading packages:\n")
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
  cat(sprintf("  ✓ %s (version: %s)\n", pkg, packageVersion(pkg)))
}

# ==================== 2. 定义辅助函数 ====================
cat("\n2. Defining helper functions...\n")

# 函数: 将基因型转换为数值 (0/2)
convert_gt_to_numeric <- function(gt_mat) {
  cat("  Converting genotypes to numeric (0/1/2)...\n")
  out <- matrix(NA_real_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  
  for (i in 1:nrow(gt_mat)) {
    for (j in 1:ncol(gt_mat)) {
      gt_str <- gt_mat[i, j]
      
      # 处理缺失值
      if (is.na(gt_str) || gt_str %in% c("./.", ".|.", ".", "NA", "")) {
        out[i, j] <- NA_real_
      } else {
        # 分割等位基因
        alleles <- strsplit(gt_str, "[|/]")[[1]]
        
        if (length(alleles) == 2) {
          # 尝试转换为数值
          a1 <- suppressWarnings(as.numeric(alleles[1]))
          a2 <- suppressWarnings(as.numeric(alleles[2]))
          
          if (is.finite(a1) && is.finite(a2)) {
            out[i, j] <- a1 + a2
          } else {
            # 如果是字符型等位基因（如A/T/C/G），转换为数值
            if (alleles[1] == alleles[2]) {
              # 纯合子
              out[i, j] <- 0  # 默认设为0，代表参考等位基因
            } else {
              # 杂合子
              out[i, j] <- 1
            }
          }
        } else {
          out[i, j] <- NA_real_
        }
      }
    }
  }
  
  rownames(out) <- rownames(gt_mat)
  colnames(out) <- colnames(gt_mat)
  
  cat(sprintf("  Conversion completed: %d rows × %d columns\n", nrow(out), ncol(out)))
  return(out)
}

# 函数: 创建颜色映射
create_color_palette <- function(n_colors) {
  cat(sprintf("  Creating color palette (%d colors)...\n", n_colors))
  
  if (n_colors <= 8) {
    palette <- brewer.pal(max(3, n_colors), "Set2")
    if (n_colors < 3) {
      palette <- palette[1:n_colors]
    }
  } else if (n_colors <= 12) {
    palette <- brewer.pal(n_colors, "Paired")
  } else {
    # 使用viridis颜色方案
    palette <- viridis(n_colors)
  }
  
  return(palette)
}

# 函数: 创建形状映射
create_shape_palette <- function(n_shapes) {
  cat(sprintf("  Creating shape palette (%d shapes)...\n", n_shapes))
  
  # 基础形状集（0-25是ggplot2的标准形状）
  basic_shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10, 11, 12, 13, 14)
  
  if (n_shapes <= length(basic_shapes)) {
    return(basic_shapes[1:n_shapes])
  } else {
    # 如果需要的形状多于基础形状，重复使用
    warning(sprintf("Only %d unique shapes available, will reuse shapes for %d categories", 
                   length(basic_shapes), n_shapes))
    return(rep(basic_shapes, length.out = n_shapes))
  }
}

# ==================== 3. 定义主分析函数 ====================
cat("\n3. Defining main analysis function...\n")

plot_e_corrected <- function() {
  cat("\n========== Starting Figure e (Haplotype Fitted Lines) Analysis ==========\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # ==================== 1. 路径设置 ====================
  cat("\n[Step 1] Setting up paths...\n")
  output_dir <- file.path(base_dir, "Haplotype_Environment_Analysis_Corrected")
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    cat(sprintf("  Creating output directory: %s\n", output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
  } else {
    cat(sprintf("  Output directory already exists: %s\n", output_dir))
  }
  
  # 定义输入文件路径
  vcf_file <- file.path(base_dir, "extracted_G1_to_G339.vcf")
  qtl_file <- file.path(base_dir, "QTL4.csv")
  phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")
  environment_file <- file.path(base_dir, "ECs_results.csv")
  
  target_factor <- "PAR_TEMP&GS67"
  
  cat(sprintf("  VCF file: %s\n", vcf_file))
  cat(sprintf("  QTL file: %s\n", qtl_file))
  cat(sprintf("  Phenotype file: %s\n", phenotype_file))
  cat(sprintf("  Environment factor file: %s\n", environment_file))
  cat(sprintf("  Target environment factor: %s\n", target_factor))
  
  # ==================== 2. 读取QTL标记 ====================
  cat("\n[Step 2] Reading QTL markers...\n")
  
  # 检查文件是否存在
  if (!file.exists(qtl_file)) {
    stop(sprintf("Error: QTL file does not exist: %s", qtl_file))
  }
  
  # 读取QTL数据
  qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  Read QTL data: %d rows × %d columns\n", nrow(qtl_df), ncol(qtl_df)))
  
  # 显示列名
  cat("  Column names:", paste(colnames(qtl_df), collapse = ", "), "\n")
  
  # 自动识别SNP/Marker列
  snp_col <- NULL
  if ("SNP" %in% colnames(qtl_df)) {
    snp_col <- "SNP"
  } else if ("Marker" %in% colnames(qtl_df)) {
    snp_col <- "Marker"
  } else if ("snp" %in% colnames(qtl_df)) {
    snp_col <- "snp"
  } else {
    candidate_cols <- grep("SNP|snp|Marker|marker|ID|id", colnames(qtl_df),
                           value = TRUE, ignore.case = TRUE)
    if (length(candidate_cols) > 0) {
      snp_col <- candidate_cols[1]
      cat(sprintf("  Using column '%s' as SNP identifier column\n", snp_col))
    }
  }
  
  if (is.null(snp_col)) {
    stop("Error: Cannot find SNP/Marker column in QTL file")
  }
  
  cat(sprintf("  Using column '%s' as SNP identifier\n", snp_col))
  
  # 提取有效的QTL SNPs
  qtl_snps <- unique(qtl_df[[snp_col]])
  qtl_snps <- qtl_snps[!is.na(qtl_snps) & qtl_snps != ""]
  
  if (length(qtl_snps) == 0) {
    stop("Error: No valid QTL SNPs found in QTL file")
  }
  
  cat(sprintf("  Found %d unique QTL markers\n", length(qtl_snps)))
  cat("  First 10 QTL markers:", paste(head(qtl_snps, 10), collapse = ", "), "\n")
  
  # ==================== 3. 读取VCF并提取目标SNPs ====================
  cat("\n[Step 3] Reading VCF file...\n")
  
  if (!file.exists(vcf_file)) {
    stop(sprintf("Error: VCF file does not exist: %s", vcf_file))
  }
  
  # 读取VCF文件
  vcf_data <- read.vcfR(vcf_file, verbose = FALSE)
  cat(sprintf("  Read VCF: %d variants × %d samples\n", 
              nrow(vcf_data@fix), ncol(vcf_data@gt)))
  
  # 匹配SNP IDs
  vcf_snp_ids <- vcf_data@fix[, "ID"]
  matched_indices <- which(vcf_snp_ids %in% qtl_snps)
  
  cat(sprintf("  Matched %d markers by ID\n", length(matched_indices)))
  
  # 如果通过ID匹配失败，尝试通过染色体位置匹配
  if (length(matched_indices) == 0) {
    cat("  Warning: No matches by ID, trying chromosome position matching...\n")
    
    # 识别染色体和位置列
    chr_col <- grep("CHR|Chr|chr|chromosome|Chromosome", colnames(qtl_df),
                    value = TRUE, ignore.case = TRUE)[1]
    pos_col <- grep("POS|Pos|pos|position|Position", colnames(qtl_df),
                    value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.null(chr_col) && !is.null(pos_col)) {
      cat(sprintf("  Using columns '%s' and '%s' for position matching\n", chr_col, pos_col))
      
      # 创建位置ID
      qtl_pos_ids <- paste(qtl_df[[chr_col]], qtl_df[[pos_col]], sep = "_")
      vcf_pos_ids <- paste(vcf_data@fix[, "CHROM"], vcf_data@fix[, "POS"], sep = "_")
      
      matched_indices <- which(vcf_pos_ids %in% qtl_pos_ids)
      cat(sprintf("  Matched %d markers by position\n", length(matched_indices)))
    }
  }
  
  if (length(matched_indices) == 0) {
    stop("Error: No matching SNP markers found in VCF")
  }
  
  # 提取目标SNPs
  extracted_vcf <- vcf_data[matched_indices, ]
  cat(sprintf("  Extracted %d matching SNP markers\n", length(matched_indices)))
  
  # ==================== 4. 提取基因型并转换为0/1/2 ====================
  cat("\n[Step 4] Extracting and converting genotypes...\n")
  
  # 提取基因型
  gt_matrix <- extract.gt(extracted_vcf, element = "GT")
  cat(sprintf("  Genotype matrix: %d rows × %d columns\n", nrow(gt_matrix), ncol(gt_matrix)))
  
  # 转换为数值型
  gt_numeric <- convert_gt_to_numeric(gt_matrix)
  
  # 显示基因型统计
  genotype_counts <- table(gt_numeric, useNA = "always")
  cat("  Genotype distribution:\n")
  for (val in names(genotype_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", 
                val, genotype_counts[val], 
                genotype_counts[val]/sum(!is.na(gt_numeric))*100))
  }
  
  # ==================== 5. 构建单倍型 ====================
  cat("\n[Step 5] Constructing haplotypes...\n")
  
  # 转换为数据框
  haplotype_df <- as.data.frame(t(gt_numeric), stringsAsFactors = FALSE)
  haplotype_df$Sample <- rownames(haplotype_df)
  
  snp_cols <- setdiff(colnames(haplotype_df), "Sample")
  cat(sprintf("  Number of SNP columns: %d\n", length(snp_cols)))
  cat("  SNP columns:", paste(head(snp_cols, 5), collapse = ", "), 
      if(length(snp_cols) > 5) "...\n" else "\n")
  
  # 过滤杂合子样本（基因型为1的）
  cat("  Filtering heterozygous samples...\n")
  haplotype_df$has_heterozygous <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
    any(x == 1, na.rm = TRUE)
  })
  
  # 统计杂合子数量
  n_heterozygous <- sum(haplotype_df$has_heterozygous)
  cat(sprintf("  Found %d heterozygous samples (%.1f%%)\n", 
              n_heterozygous, n_heterozygous/nrow(haplotype_df)*100))
  
  # 移除杂合子
  haplotype_df <- haplotype_df[!haplotype_df$has_heterozygous, ]
  cat(sprintf("  Number of samples after filtering heterozygotes: %d\n", nrow(haplotype_df)))
  
  # 构建单倍型字符串
  cat("  Constructing haplotype strings...\n")
  haplotype_df$Haplotype <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
    paste(ifelse(is.na(x), "N", x), collapse = "-")
  })
  
  # 统计单倍型频率
  haplotype_counts <- as.data.frame(table(haplotype_df$Haplotype), stringsAsFactors = FALSE)
  colnames(haplotype_counts) <- c("Haplotype", "Frequency")
  haplotype_counts <- haplotype_counts[order(haplotype_counts$Frequency, decreasing = TRUE), ]
  
  # 获取所有单倍型
  all_haps <- haplotype_counts$Haplotype
  n_haps <- length(all_haps)
  
  # 生成单倍型名称
  roman_nums <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
  if (n_haps <= length(roman_nums)) {
    haps_named <- paste0("Hap", roman_nums[1:n_haps])
  } else {
    haps_named <- paste0("Hap", 1:n_haps)
  }
  names(haps_named) <- all_haps
  
  cat(sprintf("  Detected %d haplotypes\n", n_haps))
  cat(sprintf("  Haplotype names: %s\n", paste(haps_named, collapse = ", ")))
  
  # ==================== 6. 创建详细的单倍型信息表 ====================
  cat("\n[Step 6] Creating detailed haplotype information tables...\n")
  
  # 6.1 获取SNP位点信息
  cat("  Extracting SNP marker information...\n")
  snp_info <- extracted_vcf@fix
  snp_info_df <- as.data.frame(snp_info, stringsAsFactors = FALSE)
  snp_info_df$SNP_ID <- snp_info_df$ID
  
  # 如果ID为空，使用CHROM_POS作为ID
  snp_info_df$SNP_ID <- ifelse(snp_info_df$SNP_ID == "" | is.na(snp_info_df$SNP_ID),
                               paste(snp_info_df$CHROM, snp_info_df$POS, sep = "_"),
                               snp_info_df$SNP_ID)
  
  # 选择需要的列
  snp_info_df <- snp_info_df[, c("CHROM", "POS", "SNP_ID", "REF", "ALT")]
  colnames(snp_info_df) <- c("Chromosome", "Position", "SNP_ID", "Reference_Allele", "Alternative_Allele")
  
  # 添加更多信息
  if ("QUAL" %in% colnames(snp_info)) {
    snp_info_df$Quality <- snp_info[, "QUAL"]
  }
  if ("FILTER" %in% colnames(snp_info)) {
    snp_info_df$Filter <- snp_info[, "FILTER"]
  }
  
  # 保存SNP位点信息
  snp_info_file <- file.path(output_dir, "SNP_Marker_Information.csv")
  write.csv(snp_info_df, snp_info_file, row.names = FALSE)
  cat(sprintf("  ✓ SNP marker information saved: %s\n", snp_info_file))
  cat(sprintf("    Contains %d SNP markers\n", nrow(snp_info_df)))
  
  # 6.2 添加单倍型命名到主数据框
  haplotype_df$Haplotype_Named <- haps_named[match(haplotype_df$Haplotype, names(haps_named))]
  
  # 6.3 创建详细的单倍型-株系对应表
  cat("  Creating detailed haplotype-line correspondence table...\n")
  haplotype_detail_df <- haplotype_df[, c("Sample", "Haplotype_Named", "Haplotype", snp_cols)]
  
  # 重命名SNP列为实际的SNP ID
  if (length(snp_cols) == nrow(snp_info_df)) {
    colnames(haplotype_detail_df)[4:ncol(haplotype_detail_df)] <- snp_info_df$SNP_ID
  } else {
    warning("Number of SNP columns doesn't match, keeping original column names")
  }
  
  # 添加株系类型信息
  haplotype_detail_df$Line_Type <- "Unknown"
  
  # 保存详细的单倍型-株系对应表
  haplotype_detail_file <- file.path(output_dir, "Haplotype_Line_Correspondence_Detailed.csv")
  write.csv(haplotype_detail_df, haplotype_detail_file, row.names = FALSE)
  cat(sprintf("  ✓ Detailed haplotype-line correspondence table saved: %s\n", haplotype_detail_file))
  cat(sprintf("    Contains %d lines × %d columns\n", nrow(haplotype_detail_df), ncol(haplotype_detail_df)))
  
  # 6.4 计算单倍型频率和占比
  cat("  Calculating haplotype frequencies and proportions...\n")
  total_samples <- nrow(haplotype_df)
  haplotype_frequency_df <- haplotype_counts
  haplotype_frequency_df$Haplotype_Named <- haps_named[match(haplotype_frequency_df$Haplotype, names(haps_named))]
  haplotype_frequency_df$Proportion <- haplotype_frequency_df$Frequency / total_samples * 100
  haplotype_frequency_df$Percentage <- sprintf("%.2f%%", haplotype_frequency_df$Proportion)
  
  # 重新排序列
  haplotype_frequency_df <- haplotype_frequency_df[, c("Haplotype_Named", "Haplotype", "Frequency", "Proportion", "Percentage")]
  
  # 保存单倍型频率表
  frequency_file <- file.path(output_dir, "Haplotype_Frequency_Table.csv")
  write.csv(haplotype_frequency_df, frequency_file, row.names = FALSE)
  cat(sprintf("  ✓ Haplotype frequency table saved: %s\n", frequency_file))
  
  # 6.5 创建每个单倍型的基因型模式表
  cat("  Creating haplotype genotype patterns table...\n")
  haplotype_patterns <- data.frame()
  
  for (hap in all_haps) {
    # 获取该单倍型的一个样本
    hap_sample <- haplotype_df[haplotype_df$Haplotype == hap, "Sample"][1]
    hap_genotypes <- haplotype_df[haplotype_df$Sample == hap_sample, snp_cols]
    
    # 创建基因型模式行
    pattern_row <- data.frame(
      Haplotype_Named = haps_named[hap],
      Haplotype = hap,
      stringsAsFactors = FALSE
    )
    
    # 添加每个SNP的基因型
    for (i in 1:length(snp_cols)) {
      snp_name <- snp_cols[i]
      if (i <= nrow(snp_info_df)) {
        snp_id <- snp_info_df$SNP_ID[i]
      } else {
        snp_id <- paste0("SNP_", i)
      }
      pattern_row[[snp_id]] <- as.character(hap_genotypes[[snp_name]])
    }
    
    haplotype_patterns <- rbind(haplotype_patterns, pattern_row)
  }
  
  # 保存单倍型基因型模式表
  patterns_file <- file.path(output_dir, "Haplotype_Genotype_Patterns.csv")
  write.csv(haplotype_patterns, patterns_file, row.names = FALSE)
  cat(sprintf("  ✓ Haplotype genotype patterns table saved: %s\n", patterns_file))
  
  # 6.6 输出单倍型统计摘要
  cat("\n  === Haplotype Statistics Summary ===\n")
  cat(sprintf("  Total number of samples: %d\n", total_samples))
  cat(sprintf("  Total number of haplotypes: %d\n", n_haps))
  cat("  Haplotype frequency distribution:\n")
  
  for (i in 1:min(10, nrow(haplotype_frequency_df))) {
    hap_info <- haplotype_frequency_df[i, ]
    cat(sprintf("    %s: %d samples (%.2f%%)\n", 
                hap_info$Haplotype_Named, 
                hap_info$Frequency, 
                hap_info$Proportion))
  }
  
  if (n_haps > 10) {
    cat(sprintf("    ... and %d more haplotypes\n", n_haps - 10))
  }
  
  # ==================== 7. 读取表型和环境因子数据 ====================
  cat("\n[Step 7] Reading phenotype and environment factor data...\n")
  
  # 检查文件是否存在
  if (!file.exists(phenotype_file)) {
    stop(sprintf("Error: Phenotype file does not exist: %s", phenotype_file))
  }
  
  if (!file.exists(environment_file)) {
    stop(sprintf("Error: Environment factor file does not exist: %s", environment_file))
  }
  
  # 读取表型数据
  cat("  Reading phenotype data...\n")
  phenotype_data <- read.table(
    phenotype_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, na.strings = c("", "NA")
  )
  
  cat(sprintf("  Phenotype data: %d rows × %d columns\n", nrow(phenotype_data), ncol(phenotype_data)))
  cat("  Phenotype data column names:", paste(colnames(phenotype_data), collapse = ", "), "\n")
  
  # 数据清理和重命名
  phenotype_data$genotype <- as.character(phenotype_data$genotype)
  phenotype_data$env_code <- as.character(phenotype_data$env_code)
  phenotype_data$TKW <- as.numeric(phenotype_data$TKW)
  
  # 重命名列
  colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
  colnames(phenotype_data)[colnames(phenotype_data) == "TKW"] <- "PH"
  
  # 移除缺失值
  n_before <- nrow(phenotype_data)
  phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]
  n_after <- nrow(phenotype_data)
  
  cat(sprintf("  Removing missing values: %d → %d rows (removed %d rows)\n", 
              n_before, n_after, n_before - n_after))
  
  # 读取环境因子数据
  cat("  Reading environment factor data...\n")
  ecs_data <- read.csv(environment_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  
  # 重命名第一列
  colnames(ecs_data)[1] <- "env_code"
  ecs_data$env_code <- as.character(ecs_data$env_code)
  
  cat(sprintf("  Environment factor data: %d rows × %d columns\n", nrow(ecs_data), ncol(ecs_data)))
  cat("  Environment factor column names:", paste(head(colnames(ecs_data), 10), collapse = ", "), 
      if(ncol(ecs_data) > 10) "...\n" else "\n")
  
  # 确认目标环境因子存在
  if (!target_factor %in% colnames(ecs_data)) {
    cat("  Warning: Target factor not found, trying fuzzy matching...\n")
    possible <- grep("PAR.*TEMP", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
    if (length(possible) > 0) {
      target_factor <- possible[1]
      cat(sprintf("  Using: %s\n", target_factor))
    } else {
      stop("Error: Cannot find matching environment factor")
    }
  }
  
  cat(sprintf("  Target environment factor: %s\n", target_factor))
  
  # 提取环境因子数据
  env_factor_data <- ecs_data[, c("env_code", target_factor)]
  colnames(env_factor_data)[2] <- "env_factor_value"
  env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)
  
  # 移除缺失值
  n_before_env <- nrow(env_factor_data)
  env_factor_data <- env_factor_data[!is.na(env_factor_data$env_factor_value), ]
  n_after_env <- nrow(env_factor_data)
  
  cat(sprintf("  Environment factor data: %d → %d rows (removed %d missing values)\n", 
              n_before_env, n_after_env, n_before_env - n_after_env))
  
  # ==================== 8. 合并数据 ====================
  cat("\n[Step 8] Merging data...\n")
  
  # 合并表型和单倍型数据
  cat("  Merging phenotype and haplotype data...\n")
  pheno_haplo <- merge(
    phenotype_data,
    haplotype_df[, c("Sample", "Haplotype", "Haplotype_Named")],
    by.x = "line_code", by.y = "Sample", 
    all.x = FALSE, all.y = FALSE
  )
  
  cat(sprintf("  Phenotype-haplotype merge: %d rows\n", nrow(pheno_haplo)))
  
  # 合并环境因子
  cat("  Merging environment factor data...\n")
  pheno_haplo_env <- merge(
    pheno_haplo,
    env_factor_data,
    by = "env_code", all.x = FALSE, all.y = FALSE
  )
  
  # 移除缺失值
  pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$env_factor_value), ]
  
  cat(sprintf("  Final merged data: %d rows × %d columns\n", 
              nrow(pheno_haplo_env), ncol(pheno_haplo_env)))
  cat(sprintf("  Number of lines: %d\n", length(unique(pheno_haplo_env$line_code))))
  cat(sprintf("  Number of haplotypes: %d\n", length(unique(pheno_haplo_env$Haplotype_Named))))
  cat(sprintf("  Number of environments: %d\n", length(unique(pheno_haplo_env$env_code))))
  
  # 保存合并的数据
  merged_data_file <- file.path(output_dir, "Merged_Pheno_Haplo_Env_Data.csv")
  write.csv(pheno_haplo_env, merged_data_file, row.names = FALSE)
  cat(sprintf("  ✓ Merged data saved: %s\n", merged_data_file))
  
  # ==================== 9. 株系拟合（按单倍型分组） ====================
  cat("\n[Step 9] Line fitting by haplotype group...\n")
  
  haplotype_line_models <- data.frame()
  model_count <- 0
  
  for (hap_name in unique(pheno_haplo_env$Haplotype_Named)) {
    # 获取该单倍型的所有数据
    hap_data <- pheno_haplo_env[pheno_haplo_env$Haplotype_Named == hap_name, ]
    
    # 获取该单倍型的所有株系
    lines_in_hap <- unique(hap_data$line_code)
    
    for (line in lines_in_hap) {
      line_data <- hap_data[hap_data$line_code == line, ]
      
      # 确保株系有足够的数据点
      if (nrow(line_data) >= 2) {
        m <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
        if (!inherits(m, "try-error")) {
          fit_summary <- summary(m)
          
          haplotype_line_models <- rbind(haplotype_line_models, data.frame(
            Haplotype_Named = hap_name,
            Haplotype = unique(line_data$Haplotype)[1],
            line_code = line,
            Slope = fit_summary$coefficients["env_factor_value", "Estimate"],
            Intercept = fit_summary$coefficients["(Intercept)", "Estimate"],
            Slope_SE = fit_summary$coefficients["env_factor_value", "Std. Error"],
            Intercept_SE = fit_summary$coefficients["(Intercept)", "Std. Error"],
            Slope_p = fit_summary$coefficients["env_factor_value", "Pr(>|t|)"],
            Intercept_p = fit_summary$coefficients["(Intercept)", "Pr(>|t|)"],
            R_squared = fit_summary$r.squared,
            Adj_R_squared = fit_summary$adj.r.squared,
            n_obs = nrow(line_data),
            stringsAsFactors = FALSE
          ))
          model_count <- model_count + 1
        }
      }
    }
  }
  
  cat(sprintf("  Number of line fitting models: %d\n", model_count))
  
  if (model_count == 0) {
    warning("Warning: No line fitting models generated")
  } else {
    # 保存株系拟合模型结果
    line_models_file <- file.path(output_dir, "Haplotype_Line_Models.csv")
    write.csv(haplotype_line_models, line_models_file, row.names = FALSE)
    cat(sprintf("  ✓ Line fitting model results saved: %s\n", line_models_file))
  }
  
  # ==================== 10. 单倍型平均拟合 ====================
  cat("\n[Step 10] Performing haplotype mean fitting...\n")
  
  haplotype_models <- data.frame()
  hap_model_count <- 0
  
  for (hap_name in unique(pheno_haplo_env$Haplotype_Named)) {
    # 获取该单倍型的所有数据
    hap_data <- pheno_haplo_env[pheno_haplo_env$Haplotype_Named == hap_name, ]
    
    # 确保有足够的数据点
    if (nrow(hap_data) >= 4) {
      m <- try(lm(PH ~ env_factor_value, data = hap_data), silent = TRUE)
      if (!inherits(m, "try-error")) {
        fit_summary <- summary(m)
        
        # 计算置信区间
        conf_int <- try(confint(m, level = 0.95), silent = TRUE)
        
        if (inherits(conf_int, "try-error")) {
          conf_int <- matrix(NA, nrow = 2, ncol = 2)
          rownames(conf_int) <- c("(Intercept)", "env_factor_value")
          colnames(conf_int) <- c("2.5 %", "97.5 %")
        }
        
        haplotype_models <- rbind(haplotype_models, data.frame(
          Haplotype_Named = hap_name,
          Haplotype = unique(hap_data$Haplotype)[1],
          Slope = fit_summary$coefficients["env_factor_value", "Estimate"],
          Intercept = fit_summary$coefficients["(Intercept)", "Estimate"],
          Slope_SE = fit_summary$coefficients["env_factor_value", "Std. Error"],
          Intercept_SE = fit_summary$coefficients["(Intercept)", "Std. Error"],
          Slope_p = fit_summary$coefficients["env_factor_value", "Pr(>|t|)"],
          Intercept_p = fit_summary$coefficients["(Intercept)", "Pr(>|t|)"],
          Slope_CI_lower = conf_int["env_factor_value", 1],
          Slope_CI_upper = conf_int["env_factor_value", 2],
          Intercept_CI_lower = conf_int["(Intercept)", 1],
          Intercept_CI_upper = conf_int["(Intercept)", 2],
          R_squared = fit_summary$r.squared,
          Adj_R_squared = fit_summary$adj.r.squared,
          F_statistic = ifelse(!is.null(fit_summary$fstatistic[1]), fit_summary$fstatistic[1], NA),
          F_p_value = ifelse(!is.null(fit_summary$fstatistic), 
                            pf(fit_summary$fstatistic[1], fit_summary$fstatistic[2], 
                               fit_summary$fstatistic[3], lower.tail = FALSE), NA),
          n_lines = length(unique(hap_data$line_code)),
          n_obs = nrow(hap_data),
          stringsAsFactors = FALSE
        ))
        hap_model_count <- hap_model_count + 1
      }
    }
  }
  
  cat(sprintf("  Number of haplotype mean fitting models: %d\n", hap_model_count))
  
  if (hap_model_count == 0) {
    warning("Warning: No haplotype mean fitting models generated")
  } else {
    # 保存单倍型平均拟合模型结果
    mean_models_file <- file.path(output_dir, "Haplotype_Mean_Models.csv")
    write.csv(haplotype_models, mean_models_file, row.names = FALSE)
    cat(sprintf("  ✓ Haplotype mean fitting model results saved: %s\n", mean_models_file))
  }
  
  # ==================== 11. 统计检验 ====================
  cat("\n[Step 11] Performing statistical tests...\n")
  
  # 检验不同单倍型之间株系斜率的差异
  if (nrow(haplotype_line_models) >= 5 && 
      length(unique(haplotype_line_models$Haplotype_Named)) >= 2) {
    
    cat("  Performing slope ANOVA test...\n")
    # 斜率ANOVA检验
    slope_anova <- aov(Slope ~ Haplotype_Named, data = haplotype_line_models)
    slope_anova_summary <- summary(slope_anova)
    
    # 斜率Tukey HSD检验
    cat("  Performing slope Tukey HSD test...\n")
    slope_tukey <- TukeyHSD(slope_anova)
    slope_tukey_df <- as.data.frame(slope_tukey$Haplotype_Named)
    slope_tukey_df$Comparison <- rownames(slope_tukey_df)
    
    # 保存斜率Tukey检验结果
    slope_tukey_file <- file.path(output_dir, "Slope_Tukey_HSD_Results.csv")
    write.csv(slope_tukey_df, slope_tukey_file, row.names = FALSE)
    cat(sprintf("  ✓ Slope Tukey HSD test results saved: %s\n", slope_tukey_file))
    
    cat("  Performing intercept ANOVA test...\n")
    # 截距ANOVA检验
    intercept_anova <- aov(Intercept ~ Haplotype_Named, data = haplotype_line_models)
    intercept_anova_summary <- summary(intercept_anova)
    
    # 截距Tukey HSD检验
    cat("  Performing intercept Tukey HSD test...\n")
    intercept_tukey <- TukeyHSD(intercept_anova)
    intercept_tukey_df <- as.data.frame(intercept_tukey$Haplotype_Named)
    intercept_tukey_df$Comparison <- rownames(intercept_tukey_df)
    
    # 保存截距Tukey检验结果
    intercept_tukey_file <- file.path(output_dir, "Intercept_Tukey_HSD_Results.csv")
    write.csv(intercept_tukey_df, intercept_tukey_file, row.names = FALSE)
    cat(sprintf("  ✓ Intercept Tukey HSD test results saved: %s\n", intercept_tukey_file))
    
    # 输出统计结果
    cat("\n  === Slope ANOVA Results ===\n")
    print(slope_anova_summary)
    
    cat("\n  === Intercept ANOVA Results ===\n")
    print(intercept_anova_summary)
    
  } else {
    cat("  Warning: Insufficient data for ANOVA/Tukey tests\n")
    cat(sprintf("    Need at least 5 models and 2 haplotypes, current: %d models, %d haplotypes\n", 
                nrow(haplotype_line_models), 
                length(unique(haplotype_line_models$Haplotype_Named))))
  }
  
  # ==================== 12. 颜色和形状映射 ====================
  cat("\n[Step 12] Creating color and shape mappings...\n")
  
  # 颜色映射
  n_colors <- length(unique(pheno_haplo_env$Haplotype_Named))
  hap_cols <- create_color_palette(n_colors)
  names(hap_cols) <- unique(pheno_haplo_env$Haplotype_Named)
  
  cat(sprintf("  Created color mapping for %d haplotypes\n", n_colors))
  
  # 形状映射
  env_codes <- unique(pheno_haplo_env$env_code)
  n_shapes <- length(env_codes)
  env_shapes <- create_shape_palette(n_shapes)
  names(env_shapes) <- env_codes
  
  cat(sprintf("  Created shape mapping for %d environments\n", n_shapes))
  cat("  Environment codes:", paste(env_codes, collapse = ", "), "\n")
  
  # ==================== 13. 创建预测数据 ====================
  cat("\n[Step 13] Creating prediction data...\n")
  
  env_range <- range(pheno_haplo_env$env_factor_value, na.rm = TRUE)
  env_seq <- seq(env_range[1], env_range[2], length.out = 100)
  
  cat(sprintf("  Environment factor range: %.2f to %.2f\n", env_range[1], env_range[2]))
  cat(sprintf("  Creating prediction sequence with %d points\n", length(env_seq)))
  
  # 株系预测数据
  cat("  Creating line prediction data...\n")
  line_pred_list <- list()
  if (nrow(haplotype_line_models) > 0) {
    for (i in 1:nrow(haplotype_line_models)) {
      model <- haplotype_line_models[i, ]
      line_pred_list[[i]] <- data.frame(
        Haplotype_Named = model$Haplotype_Named,
        line_code = model$line_code,
        env_factor_value = env_seq,
        Predicted_PH = model$Intercept + model$Slope * env_seq,
        stringsAsFactors = FALSE
      )
    }
    line_pred <- do.call(rbind, line_pred_list)
    cat(sprintf("  Created %d line prediction curves\n", length(line_pred_list)))
  } else {
    line_pred <- data.frame()
    cat("  Warning: No line prediction data\n")
  }
  
  # 单倍型平均预测数据
  cat("  Creating haplotype mean prediction data...\n")
  hap_pred_list <- list()
  if (nrow(haplotype_models) > 0) {
    for (i in 1:nrow(haplotype_models)) {
      model <- haplotype_models[i, ]
      hap_pred_list[[i]] <- data.frame(
        Haplotype_Named = model$Haplotype_Named,
        env_factor_value = env_seq,
        Predicted_PH = model$Intercept + model$Slope * env_seq,
        stringsAsFactors = FALSE
      )
    }
    hap_pred <- do.call(rbind, hap_pred_list)
    cat(sprintf("  Created %d haplotype mean prediction curves\n", length(hap_pred_list)))
  } else {
    hap_pred <- data.frame()
    cat("  Warning: No haplotype mean prediction data\n")
  }
  
  # ==================== 14. 创建包含统计信息的标签 ====================
  cat("\n[Step 14] Creating statistical labels...\n")
  
  haplotype_labels <- character(length(unique(pheno_haplo_env$Haplotype_Named)))
  names(haplotype_labels) <- unique(pheno_haplo_env$Haplotype_Named)
  
  for (hap_name in unique(pheno_haplo_env$Haplotype_Named)) {
    if (hap_name %in% haplotype_models$Haplotype_Named) {
      model <- haplotype_models[haplotype_models$Haplotype_Named == hap_name, ]
      
      # 格式化斜率和截距
      slope_fmt <- sprintf("%.2f", model$Slope)
      intercept_fmt <- sprintf("%.1f", model$Intercept)
      
      # 创建标签格式: HapI(截距，斜率，样本数)
      haplotype_labels[hap_name] <- sprintf("%s(%.1f, %s, %d)", 
                                            hap_name, 
                                            model$Intercept,
                                            slope_fmt,
                                            model$n_lines)
    } else {
      haplotype_labels[hap_name] <- hap_name
    }
  }
  
  cat("  Created labels:\n")
  for (label in haplotype_labels) {
    cat(sprintf("    %s\n", label))
  }
  
  # ==================== 15. 创建图例 ====================
  cat("\n[Step 15] Creating legends...\n")
  
  # 创建单倍型颜色图例
  legend_data_hap <- data.frame(
    Haplotype_Named = names(haplotype_labels),
    Label = haplotype_labels,
    stringsAsFactors = FALSE
  )
  
  p_legend_hap <- ggplot(legend_data_hap, aes(x = 0, y = 0, color = Haplotype_Named)) +
    geom_point(size = 3) +
    scale_color_manual(values = hap_cols, 
                       labels = haplotype_labels,
                       name = "Haplotype\n(intercept, slope, n)") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9, lineheight = 1.1),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(1.2, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  legend_grob_hap <- get_legend(p_legend_hap)
  cat("  Haplotype legend created\n")
  
  # 创建环境形状图例
  legend_data_env <- data.frame(
    env_code = env_codes,
    Label = env_codes,
    stringsAsFactors = FALSE
  )
  
  p_legend_env <- ggplot(legend_data_env, aes(x = 0, y = 0, shape = env_code)) +
    geom_point(size = 3) +
    scale_shape_manual(values = env_shapes,
                       name = "Environment") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(1.2, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(shape = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  legend_grob_env <- get_legend(p_legend_env)
  cat("  Environment legend created\n")
  
  # ==================== 16. 绘制主图 ====================
  cat("\n[Step 16] Drawing main plot...\n")
  
  # 计算y轴范围
  y_min <- min(pheno_haplo_env$PH, na.rm = TRUE)
  y_max <- max(pheno_haplo_env$PH, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 扩展y轴范围，为图例留出空间
  y_expand_factor <- 0.25
  
  cat(sprintf("  Phenotype range: %.2f to %.2f (range: %.2f)\n", y_min, y_max, y_range))
  cat(sprintf("  Y-axis expansion factor: %.2f\n", y_expand_factor))
  
  # 创建主图
  p_main <- ggplot() +
    # 株系实际观测点（使用颜色区分单倍型，形状区分环境）
    geom_point(
      data = pheno_haplo_env,
      aes(x = env_factor_value, y = PH, color = Haplotype_Named, shape = env_code),
      alpha = 0.7, size = 2.0
    ) +
    # 株系拟合线（虚线）
    geom_line(
      data = line_pred,
      aes(x = env_factor_value, y = Predicted_PH, group = line_code, color = Haplotype_Named),
      linetype = "dashed", alpha = 0.3, linewidth = 0.5
    ) +
    # 单倍型平均拟合线（实线）
    geom_line(
      data = hap_pred,
      aes(x = env_factor_value, y = Predicted_PH, color = Haplotype_Named),
      linetype = "solid", alpha = 0.9, linewidth = 1.5
    ) +
    # 颜色映射
    scale_color_manual(values = hap_cols, guide = "none") +
    # 形状映射
    scale_shape_manual(values = env_shapes, guide = "none") +
    # 坐标轴标签（使用英文）
    labs(
      x = "PAR_TEMP&GS67 (PAR to temperature ratio during early grain filling)",
      y = "Thousand Kernel Weight (g)"
    ) +
    # 子图标签
    annotate("text",
             x = -Inf,
             y = y_max + y_range * y_expand_factor,
             label = "e",
             size = label_params$label_size,
             fontface = label_params$label_fontface,
             hjust = label_params$label_hjust,
             vjust = 0.5) +
    # 主题设置
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.0),
      plot.margin = margin(10, 10, 10, 20),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    # 坐标轴范围
    coord_cartesian(
      xlim = env_range,
      ylim = c(y_min - y_range * 0.05, y_max + y_range * y_expand_factor),
      clip = "off"
    )
  
  cat("  Main plot drawn\n")
  
  # ==================== 17. 组合图形和图例 ====================
  cat("\n[Step 17] Combining plot and legends...\n")
  
  # 创建图例组合
  p_combined_legends <- plot_grid(
    legend_grob_hap,
    legend_grob_env,
    ncol = 1,
    align = "v",
    rel_heights = c(1, 0.6)
  )
  
  # 组合主图和图例
  p_final <- plot_grid(
    p_main,
    p_combined_legends,
    ncol = 2,
    rel_widths = c(3, 1),
    align = "h",
    axis = "tb"
  )
  
  cat("  Plot and legends combined\n")
  
  # ==================== 18. 保存图形 ====================
  cat("\n[Step 18] Saving figures...\n")
  
  # 保存PNG格式
  png_file <- file.path(output_dir, "Figure_e_Haplotype_Response_Corrected.png")
  ggsave(png_file, p_final, width = 16, height = 10, dpi = 300, bg = "white")
  cat(sprintf("  ✓ PNG figure saved: %s\n", png_file))
  
  # 保存PDF格式
  pdf_file <- file.path(output_dir, "Figure_e_Haplotype_Response_Corrected.pdf")
  ggsave(pdf_file, p_final, width = 16, height = 10, bg = "white")
  cat(sprintf("  ✓ PDF figure saved: %s\n", pdf_file))
  
  # ==================== 19. 创建汇总报告 ====================
  cat("\n[Step 19] Creating summary report...\n")
  
  # 创建汇总报告文件
  summary_file <- file.path(output_dir, "Analysis_Summary.txt")
  
  summary_text <- c(
    "==================================================",
    "       Thousand Kernel Weight Haplotype Analysis",
    "==================================================",
    sprintf("Analysis time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Output directory: %s", output_dir),
    "",
    "1. Input Files:",
    sprintf("  VCF file: %s", vcf_file),
    sprintf("  QTL file: %s", qtl_file),
    sprintf("  Phenotype file: %s", phenotype_file),
    sprintf("  Environment factor file: %s", environment_file),
    sprintf("  Target environment factor: %s", target_factor),
    "",
    "2. Data Statistics:",
    sprintf("  Total samples: %d", total_samples),
    sprintf("  SNP markers: %d", length(snp_cols)),
    sprintf("  Haplotypes: %d", n_haps),
    sprintf("  Haplotype sample distribution: %s", 
            paste(sapply(1:min(5, nrow(haplotype_frequency_df)), function(i) {
              hap <- haplotype_frequency_df[i, ]
              sprintf("%s(%d, %.1f%%)", hap$Haplotype_Named, hap$Frequency, hap$Proportion)
            }), collapse = ", ")),
    "",
    "3. Analysis Results:",
    sprintf("  Merged data rows: %d", nrow(pheno_haplo_env)),
    sprintf("  Number of lines: %d", length(unique(pheno_haplo_env$line_code))),
    sprintf("  Number of environments: %d", length(unique(pheno_haplo_env$env_code))),
    sprintf("  Line fitting models: %d", nrow(haplotype_line_models)),
    sprintf("  Haplotype mean fitting models: %d", nrow(haplotype_models)),
    "",
    "4. Generated Files:",
    "  1. SNP_Marker_Information.csv - SNP marker details",
    "  2. Haplotype_Line_Correspondence_Detailed.csv - Detailed haplotype-line correspondence",
    "  3. Haplotype_Frequency_Table.csv - Haplotype frequency statistics",
    "  4. Haplotype_Genotype_Patterns.csv - Haplotype genotype patterns",
    "  5. Merged_Pheno_Haplo_Env_Data.csv - Merged data table",
    "  6. Haplotype_Line_Models.csv - Line fitting model results",
    "  7. Haplotype_Mean_Models.csv - Haplotype mean fitting model results",
    "  8. Figure_e_Haplotype_Response_Corrected.png/pdf - Figure results"
  )
  
  # 如果有统计检验结果，添加到报告中
  if (nrow(haplotype_line_models) >= 5 && 
      length(unique(haplotype_line_models$Haplotype_Named)) >= 2) {
    summary_text <- c(summary_text,
                      "  9. Slope_Tukey_HSD_Results.csv - Slope Tukey test results",
                      "  10. Intercept_Tukey_HSD_Results.csv - Intercept Tukey test results")
  }
  
  summary_text <- c(summary_text,
                    "",
                    "5. Haplotype Response Characteristics:",
                    sprintf("  Environment factor range: %.2f to %.2f", env_range[1], env_range[2]),
                    sprintf("  Phenotype range: %.2f to %.2f g", y_min, y_max))
  
  # 添加单倍型具体信息
  if (nrow(haplotype_models) > 0) {
    summary_text <- c(summary_text, "")
    for (i in 1:nrow(haplotype_models)) {
      model <- haplotype_models[i, ]
      summary_text <- c(summary_text,
                       sprintf("  %s: intercept=%.1f, slope=%.3f, R²=%.3f, n=%d", 
                               model$Haplotype_Named, 
                               model$Intercept,
                               model$Slope, 
                               model$R_squared, 
                               model$n_lines))
    }
  }
  
  summary_text <- c(summary_text,
                    "",
                    "==================================================",
                    "Analysis completed!",
                    "==================================================")
  
  # 写入汇总报告
  writeLines(summary_text, summary_file, useBytes = TRUE)
  cat(sprintf("  ✓ Summary report saved: %s\n", summary_file))
  
  # ==================== 20. 完成 ====================
  cat("\n========== Analysis Completed ==========\n")
  cat(sprintf("✅ All results saved to: %s\n", output_dir))
  cat("Generated files:\n")
  
  # 列出生成的文件
  generated_files <- list.files(output_dir, full.names = TRUE)
  for (i in 1:length(generated_files)) {
    file_size <- file.info(generated_files[i])$size / 1024  # KB
    cat(sprintf("  %d. %s (%.1f KB)\n", i, basename(generated_files[i]), file_size))
  }
  
  cat(sprintf("\nTotal files: %d\n", length(generated_files)))
  cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("Total runtime:", format(Sys.time() - start_time, digits = 2), "\n")
  cat("========================================\n")
  
  return(p_final)
}

# ==================== 4. 主程序 ====================
cat("\n4. Starting main program...\n")

# 记录开始时间
start_time <- Sys.time()

# 执行主函数
tryCatch({
  result_plot <- plot_e_corrected()
  cat("\n✅ Program executed successfully!\n")
}, error = function(e) {
  cat(sprintf("\n❌ Program execution failed: %s\n", e$message))
  cat("\nError details:\n")
  print(e)
  cat("\nPlease check if the input file paths and formats are correct.\n")
  quit(status = 1)
})

# ==================== 5. 程序结束 ====================
cat("\n5. Program finished\n")
cat("Total runtime:", format(Sys.time() - start_time, digits = 2), "\n")

# 结束程序
cat("\n========== Thousand Kernel Weight Haplotype Analysis Completed ==========\n")