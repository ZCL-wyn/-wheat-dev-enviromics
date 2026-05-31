#!/usr/bin/env Rscript

# ============================================================================
# TKW基因组选择 - 独立外部验证（最终定稿版：MLM + GAPIT/BGLR 双框架）
# 训练集：TKW_mean_table.txt
# 验证集：testTKW.txt
# 基因型：myGD.csv
# 图谱：myGM.csv
# 环境参数：EC8.csv (使用 DTR&mean&PostFlowering1_14)
# 染色体映射：严格映射表，只接受预定义的22种染色体名
#
# 配置总数：
#   2个框架 × 3种方法 × 2种预测方式 = 12种配置
#
# 框架：
#   - GAPIT
#   - BGLR
#
# 方法：
#   - MAS
#   - GBLUP
#   - MAS+GBLUP
#
# 预测方式：
#   - 截距
#   - 截距+斜率
#
# 核心规则（最终版）：
#   1. GWAS 使用 GAPIT MLM
#   2. 显著标记筛选：
#        P.value < 0.0001  或  H&B.P.Value < 0.05
#   3. 标记去冗余：
#        同染色体内按效应绝对值优先，LD > 0.7 的候选标记只保留一个
#   4. MAS+GBLUP 背景标记过滤：
#        删除与显著标记同染色体、3Mb内、且 LD >= 0.7 的背景标记
# ============================================================================

# ------------------------------ 1. 配置区 ------------------------------
config <- list(
  workdir      = "/mnt/7t_storage/zhangcl/TKW",
  out_dir      = "TKW_External_Validation_MLM_Final_All12Configs",

  # 文件
  train_pheno  = "TKW_mean_table.txt",
  test_pheno   = "testTKW.txt",
  geno_file    = "myGD.csv",
  map_file     = "myGM.csv",
  env_file     = "EC8.csv",
  env_para     = "DTR&mean&PostFlowering1_14",

  # 质量控制
  maf_min      = 0.03,
  pNA_max      = 0.10,

  # 显著标记筛选参数
  pval_thresh    = 0.0001,
  hb_pval_thresh = 0.05,
  ld_threshold   = 0.7,
  physical_dist  = 3e6,   # 3 Mb

  # BGLR参数
  nIter        = 8000,
  burnIn       = 2000,
  thin         = 5,
  seed         = 195021,

  # 并行核心数（检测核心的一半）
  ncores       = max(1, floor(parallel::detectCores() * 0.5))
)

# ------------------------------ 2. 加载包 ------------------------------
required_packages <- c(
  "data.table", "BGLR", "glmnet", "ggplot2",
  "rrBLUP", "parallel", "doParallel", "foreach",
  "Matrix", "dplyr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少R包：", pkg, "。请先安装后再运行。")
  }
}
invisible(lapply(required_packages, library, character.only = TRUE))

setwd(config$workdir)
cat("工作目录已设置为:", getwd(), "\n")

# GAPIT文件路径（使用绝对路径）
gapit_lib <- file.path(config$workdir, "GAPIT.library.R")
gapit_fun <- file.path(config$workdir, "gapit_functions.txt")

if (!file.exists(gapit_lib) || !file.exists(gapit_fun)) {
  stop("请先将 GAPIT.library.R 和 gapit_functions.txt 下载到工作目录: ", config$workdir)
}

source(gapit_lib)
source(gapit_fun)
cat("GAPIT本地加载成功\n")

# 创建输出目录
dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ 3. 工具函数 ------------------------------

# 3.1 染色体严格映射函数（只接受预定义名称）
map_chr_to_numeric <- function(chr_vec) {
  chr_map <- c(
    "Chr1A" = 1,  "Chr1B" = 2,  "Chr1D" = 3,
    "Chr2A" = 4,  "Chr2B" = 5,  "Chr2D" = 6,
    "Chr3A" = 7,  "Chr3B" = 8,  "Chr3D" = 9,
    "Chr4A" = 10, "Chr4B" = 11, "Chr4D" = 12,
    "Chr5A" = 13, "Chr5B" = 14, "Chr5D" = 15,
    "Chr6A" = 16, "Chr6B" = 17, "Chr6D" = 18,
    "Chr7A" = 19, "Chr7B" = 20, "Chr7D" = 21,
    "ChrUnknown" = 22
  )

  chr_vec <- trimws(as.character(chr_vec))
  out <- unname(chr_map[chr_vec])

  bad <- unique(chr_vec[is.na(out)])
  if (length(bad) > 0) {
    stop(
      "存在未定义染色体名: ", paste(bad, collapse = ", "),
      "\n请检查 myGM.csv 中染色体列，只允许以下名称:\n",
      paste(names(chr_map), collapse = " ")
    )
  }

  as.integer(out)
}

# 3.2 SNP名称标准化
standardize_snp_names <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub(";", "_", x)
  x <- gsub(":", "_", x)
  x <- gsub("-", "_", x)
  x <- gsub("/", "_", x)
  x <- gsub("\\|", "_", x)
  x <- gsub("\\(", "_", x)
  x <- gsub("\\)", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

# 3.3 计算评估指标
calc_metrics <- function(actual, pred) {
  ok <- complete.cases(actual, pred)
  if (sum(ok) < 3) {
    return(data.frame(
      n = sum(ok),
      Correlation = NA_real_,
      RMSE = NA_real_,
      MAE = NA_real_,
      Bias = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  actual <- actual[ok]
  pred <- pred[ok]

  data.frame(
    n = length(actual),
    Correlation = suppressWarnings(cor(actual, pred)),
    RMSE = sqrt(mean((actual - pred)^2)),
    MAE = mean(abs(actual - pred)),
    Bias = mean(pred - actual),
    stringsAsFactors = FALSE
  )
}

# 3.4 均值预测后备函数
safe_mean_pred <- function(y_train, test_lines) {
  pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
  names(pred) <- test_lines
  pred
}

# 3.5 计算亲缘关系矩阵（VanRaden，带稳健分母检查）
calc_kinship_matrix <- function(geno_matrix) {
  cat("      计算亲缘关系矩阵（VanRaden方法）...\n")

  if (!is.matrix(geno_matrix) && !is.data.frame(geno_matrix)) {
    stop("geno_matrix 必须是 matrix 或 data.frame")
  }

  geno_numeric <- data.matrix(geno_matrix)

  if (nrow(geno_numeric) < 2 || ncol(geno_numeric) < 2) {
    stop("用于计算亲缘矩阵的基因型矩阵维度不足")
  }

  p <- colMeans(geno_numeric, na.rm = TRUE) / 2
  Z <- sweep(geno_numeric, 2, 2 * p, "-")
  denom <- sum(2 * p * (1 - p), na.rm = TRUE)

  if (!is.finite(denom) || denom <= 0) {
    stop("VanRaden分母非正或非有限，无法计算亲缘矩阵")
  }

  K <- tcrossprod(Z) / denom
  K <- (K + t(K)) / 2
  diag(K) <- diag(K) + 1e-6

  cat("      亲缘关系矩阵维度:", dim(K), "\n")
  K
}

# 3.6 GAPIT预测结果提取函数
extract_gapit_prediction <- function(pred_df, test_lines) {
  if (is.null(pred_df) || nrow(pred_df) == 0) return(NULL)

  pred_col <- NULL
  possible_pred_cols <- c("Prediction", "pred", "Pred", "BLUP", "u", "predict", "Predicted_value")

  for (col in possible_pred_cols) {
    if (col %in% colnames(pred_df)) {
      pred_col <- col
      break
    }
  }

  if (is.null(pred_col)) {
    numeric_cols <- colnames(pred_df)[sapply(pred_df, is.numeric)]
    if (length(numeric_cols) > 0) pred_col <- numeric_cols[1]
  }

  if (is.null(pred_col)) return(NULL)

  taxa_col <- if ("Taxa" %in% colnames(pred_df)) {
    "Taxa"
  } else if ("taxa" %in% colnames(pred_df)) {
    "taxa"
  } else {
    NULL
  }

  if (is.null(taxa_col)) return(NULL)

  pred_test <- pred_df[pred_df[[taxa_col]] %in% test_lines, c(taxa_col, pred_col), drop = FALSE]
  if (nrow(pred_test) == 0) return(NULL)

  colnames(pred_test) <- c("Taxa", "Prediction")

  pred_values <- pred_test$Prediction
  names(pred_values) <- pred_test$Taxa
  pred_values
}

# 3.7 从GWAS结果中提取显著标记（最终规则版）
select_markers_by_pvalue_ld <- function(gwas_dir, geno_all, GM_all,
                                        pval_thresh = 0.0001,
                                        hb_pval_thresh = 0.05,
                                        ld_threshold = 0.7) {
  files <- list.files(gwas_dir, pattern = "\\.csv$", full.names = TRUE)

  gwas_files <- files[
    grepl("GWAS[._ ]?Results", basename(files), ignore.case = TRUE) |
      (grepl("GWAS", basename(files), ignore.case = TRUE) &
         grepl("Results", basename(files), ignore.case = TRUE))
  ]

  if (length(gwas_files) == 0) {
    gwas_files <- files[grepl("GWAS", basename(files), ignore.case = TRUE)]
  }

  if (length(gwas_files) == 0) {
    cat("    未找到GWAS结果文件\n")
    return(character(0))
  }

  gwas_file <- gwas_files[order(file.info(gwas_files)$mtime, decreasing = TRUE)][1]
  gwas_data <- read.csv(gwas_file, stringsAsFactors = FALSE, check.names = FALSE)

  needed_cols <- c("SNP", "Chr", "Pos", "P.value", "H&B.P.Value", "Effect")
  if (!all(needed_cols %in% colnames(gwas_data))) {
    missing <- setdiff(needed_cols, colnames(gwas_data))
    cat("    GWAS结果缺少列:", paste(missing, collapse = ", "), "\n")
    return(character(0))
  }

  gwas_data$SNP <- standardize_snp_names(gwas_data$SNP)

  pval <- suppressWarnings(as.numeric(gwas_data[["P.value"]]))
  hb   <- suppressWarnings(as.numeric(gwas_data[["H&B.P.Value"]]))
  eff  <- suppressWarnings(as.numeric(gwas_data[["Effect"]]))
  chrv <- suppressWarnings(as.numeric(gwas_data[["Chr"]]))
  posv <- suppressWarnings(as.numeric(gwas_data[["Pos"]]))

  # 最终规则：P.value < 0.0001 或 H&B.P.Value < 0.05
  idx <- which((!is.na(pval) & pval < pval_thresh) |
                 (!is.na(hb) & hb < hb_pval_thresh))

  if (length(idx) == 0) {
    cat("    没有标记满足阈值：P.value < ", pval_thresh,
        " 或 H&B.P.Value < ", hb_pval_thresh, "\n", sep = "")
    return(character(0))
  }

  candidate <- data.frame(
    SNP = gwas_data$SNP[idx],
    Chr = chrv[idx],
    Pos = posv[idx],
    Effect = eff[idx],
    stringsAsFactors = FALSE
  )

  candidate <- candidate[complete.cases(candidate[, c("SNP", "Chr", "Pos", "Effect")]), , drop = FALSE]
  candidate <- candidate[candidate$SNP %in% colnames(geno_all), , drop = FALSE]
  candidate <- candidate[candidate$SNP %in% GM_all$SNP, , drop = FALSE]

  if (nrow(candidate) == 0) {
    cat("    候选标记经清理后为空\n")
    return(character(0))
  }

  candidate$Effect_abs <- abs(candidate$Effect)
  candidate <- candidate[order(candidate$Chr, candidate$Pos), , drop = FALSE]

  cat("    初选候选标记数:", nrow(candidate), "\n")

  keep_markers <- character(0)
  chromosomes <- unique(candidate$Chr)

  for (chr in chromosomes) {
    chr_markers <- candidate[candidate$Chr == chr, , drop = FALSE]

    if (nrow(chr_markers) <= 1) {
      keep_markers <- c(keep_markers, chr_markers$SNP)
      next
    }

    # 同染色体按效应绝对值降序，LD > 0.7 的只保留一个
    chr_markers <- chr_markers[order(chr_markers$Effect_abs, decreasing = TRUE), , drop = FALSE]
    keep <- rep(TRUE, nrow(chr_markers))

    for (i in seq_len(nrow(chr_markers) - 1)) {
      if (!keep[i]) next

      snp_i <- chr_markers$SNP[i]
      if (!(snp_i %in% colnames(geno_all))) {
        keep[i] <- FALSE
        next
      }

      for (j in (i + 1):nrow(chr_markers)) {
        if (!keep[j]) next

        snp_j <- chr_markers$SNP[j]
        if (!(snp_j %in% colnames(geno_all))) {
          keep[j] <- FALSE
          next
        }

        ld <- suppressWarnings(cor(geno_all[, snp_i], geno_all[, snp_j], use = "complete.obs")^2)
        if (!is.na(ld) && ld > ld_threshold) {
          keep[j] <- FALSE
        }
      }
    }

    keep_markers <- c(keep_markers, chr_markers$SNP[keep])
  }

  keep_markers <- unique(keep_markers)

  cat("    LD过滤后保留标记数:", length(keep_markers), "\n")
  if (length(keep_markers) > 0) {
    cat("    前几个标记:", paste(head(keep_markers, 5), collapse = ", "), "\n")
  }

  keep_markers
}

# 3.8 背景标记过滤（最终规则版）
# 删除与显著标记同染色体、3Mb内、且 LD >= 0.7 的背景标记
filter_background_markers <- function(significant_markers, geno_all, GM_all,
                                      physical_distance = 3e6,
                                      ld_threshold = 0.7) {
  cat("        过滤背景标记（删除显著标记3Mb内且LD>=0.7的背景标记）...\n")

  if (length(significant_markers) == 0) {
    return(colnames(geno_all))
  }

  significant_markers <- intersect(significant_markers, colnames(geno_all))
  if (length(significant_markers) == 0) {
    return(colnames(geno_all))
  }

  sig_info <- GM_all[GM_all$SNP %in% significant_markers, , drop = FALSE]
  if (nrow(sig_info) == 0) {
    return(colnames(geno_all))
  }

  all_markers <- colnames(geno_all)
  all_info <- GM_all[GM_all$SNP %in% all_markers, , drop = FALSE]
  markers_to_exclude <- character(0)

  for (i in seq_len(nrow(sig_info))) {
    sig_snp <- sig_info$SNP[i]
    sig_chr <- sig_info$Chromosome[i]
    sig_pos <- sig_info$Position[i]

    chr_markers <- all_info[all_info$Chromosome == sig_chr, , drop = FALSE]
    if (nrow(chr_markers) == 0) next

    nearby_markers <- chr_markers[
      abs(chr_markers$Position - sig_pos) <= physical_distance &
        chr_markers$SNP != sig_snp,
      , drop = FALSE
    ]

    if (nrow(nearby_markers) == 0) next

    for (bg_snp in nearby_markers$SNP) {
      if (bg_snp %in% colnames(geno_all) && sig_snp %in% colnames(geno_all)) {
        ld <- suppressWarnings(cor(geno_all[, sig_snp], geno_all[, bg_snp], use = "complete.obs")^2)
        if (!is.na(ld) && ld >= ld_threshold) {
          markers_to_exclude <- c(markers_to_exclude, bg_snp)
        }
      }
    }
  }

  markers_to_exclude <- unique(markers_to_exclude)
  background_markers <- setdiff(all_markers, c(significant_markers, markers_to_exclude))

  cat("        需排除背景标记数:", length(markers_to_exclude), "\n")
  cat("        保留背景标记数:", length(background_markers), "\n")

  background_markers
}

# ------------------------------ 4. FW参数计算 ------------------------------
calculate_FW_parameters <- function(train_pheno, env_factors, env_means) {
  cat("    计算FW模型参数...\n")

  genotypes <- unique(train_pheno$genotype)
  train_envs <- unique(train_pheno$env_code)

  train_env_factors <- env_factors[env_factors$env_code %in% train_envs, , drop = FALSE]
  mean_kPara_train <- mean(train_env_factors$kPara, na.rm = TRUE)

  cat("    训练集环境参数平均值:", mean_kPara_train, "\n")

  results <- data.frame()

  for (geno in genotypes) {
    geno_data <- subset(train_pheno, genotype == geno)
    if (nrow(geno_data) < 3) next

    geno_data_merged <- merge(geno_data, env_factors, by = "env_code")
    geno_data_merged <- merge(geno_data_merged, env_means, by = "env_code")

    geno_data_clean <- geno_data_merged[
      complete.cases(geno_data_merged$TKW, geno_data_merged$kPara),
      , drop = FALSE
    ]

    if (nrow(geno_data_clean) < 2) next

    fit_res <- tryCatch({
      geno_data_clean$kPara_centered <- geno_data_clean$kPara - mean_kPara_train
      lm_para <- lm(TKW ~ kPara_centered, data = geno_data_clean)
      coefs <- coef(lm_para)

      if (length(coefs) < 2) return(NULL)

      data.frame(
        genotype = geno,
        Intcp_para_adj = as.numeric(coefs[1]),
        Slope_para = as.numeric(coefs[2]),
        R2_para = summary(lm_para)$r.squared,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("      品种", geno, "的FW参数计算失败:", e$message, "\n")
      NULL
    })

    if (!is.null(fit_res)) {
      results <- rbind(results, fit_res)
    }
  }

  cat("    成功计算FW参数的品种数量:", nrow(results), "\n")

  list(
    FW_params = results,
    mean_kPara_train = mean_kPara_train
  )
}

# ------------------------------ 5. GAPIT GWAS函数（MLM版） ------------------------------
perform_gapit_gwas_mlm <- function(train_FW_params, GD_all_gapit, GM_all, trait, output_dir) {
  cat("    在训练集上进行", trait, "的 GAPIT MLM GWAS 分析...\n")

  if (!is.data.frame(train_FW_params) || nrow(train_FW_params) == 0) {
    return(NULL)
  }

  pheno_df <- data.frame(
    Taxa = train_FW_params$genotype,
    Trait = train_FW_params[[trait]],
    stringsAsFactors = FALSE
  )
  pheno_df <- pheno_df[complete.cases(pheno_df$Trait), , drop = FALSE]

  if (nrow(pheno_df) < 10) {
    cat("      有效样本不足，跳过", trait, "的 GWAS\n")
    return(NULL)
  }

  common_genotypes <- intersect(pheno_df$Taxa, GD_all_gapit$taxa)
  if (length(common_genotypes) < 10) {
    cat("      共同基因型不足，跳过", trait, "的 GWAS\n")
    return(NULL)
  }

  pheno_gapit <- pheno_df[match(common_genotypes, pheno_df$Taxa), , drop = FALSE]
  GD_gapit <- GD_all_gapit[match(common_genotypes, GD_all_gapit$taxa), , drop = FALSE]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  original_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_wd), add = TRUE)

  tryCatch({
    GAPIT(
      Y = pheno_gapit,
      G = NULL,
      GD = GD_gapit,
      GM = GM_all,
      model = "MLM",
      PCA.total = 2,
      Multiple_analysis = FALSE,
      file.output = TRUE,
      Major.allele.zero = FALSE,
      SNP.MAF = 0.01,
      cutOff = 0.05,
      kinship.algorithm = "VanRaden",
      memo = paste("GWAS", trait, sep = "_")
    )

    list(output_dir = output_dir)
  }, error = function(e) {
    cat("    GAPIT分析失败:", e$message, "\n")
    NULL
  })
}

# ------------------------------ 6. GAPIT 预测函数 ------------------------------

# 6.1 GAPIT MAS
GAPIT_MAS_prediction_simple <- function(y_train, train_lines, test_lines,
                                        GD_all, GM_all,
                                        significant_markers) {
  cat("      执行 GAPIT MAS 预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = NA_real_,
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  available_markers <- intersect(significant_markers, colnames(GD_all))
  if (length(available_markers) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  markers_geno <- GD_all[all_lines, c("taxa", available_markers), drop = FALSE]
  colnames(markers_geno)[1] <- "Taxa"

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_MAS_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  on.exit({
    setwd(original_wd)
    unlink(temp_dir, recursive = TRUE)
  }, add = TRUE)

  setwd(temp_dir)

  tryCatch({
    gapit_result <- GAPIT(
      Y = Y_all,
      CV = markers_geno,
      GD = NULL,
      GM = NULL,
      model = "GLM",
      SNP.test = FALSE,
      file.output = TRUE,
      memo = "MAS_Prediction_NoPCA"
    )

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        return(pred_values)
      }
    }
    safe_mean_pred(y_train, test_lines)
  }, error = function(e) {
    cat("        GAPIT MAS 预测失败:", e$message, "\n")
    safe_mean_pred(y_train, test_lines)
  })
}

# 6.2 GAPIT GBLUP
GAPIT_GBLUP_prediction_simple <- function(y_train, train_lines, test_lines,
                                          GD_all, GM_all) {
  cat("      执行 GAPIT GBLUP 预测 (无PCA)...\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = NA_real_,
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_GBLUP_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  on.exit({
    setwd(original_wd)
    unlink(temp_dir, recursive = TRUE)
  }, add = TRUE)

  setwd(temp_dir)

  tryCatch({
    gapit_result <- GAPIT(
      Y = Y_all,
      G = NULL,
      GD = GD_all[all_lines, , drop = FALSE],
      GM = GM_all,
      model = "gBLUP",
      PCA.total = 0,
      file.output = FALSE,
      Major.allele.zero = FALSE,
      memo = "GBLUP_NoPCA"
    )

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        return(pred_values)
      }
    }
    safe_mean_pred(y_train, test_lines)
  }, error = function(e) {
    cat("        GAPIT GBLUP 预测失败:", e$message, "\n")
    safe_mean_pred(y_train, test_lines)
  })
}

# 6.3 GAPIT MAS+GBLUP
GAPIT_MAS_gBLUP_prediction_simple <- function(y_train, train_lines, test_lines,
                                              GD_all, GM_all,
                                              significant_markers,
                                              physical_distance = 3e6,
                                              ld_threshold = 0.7) {
  cat("      执行 GAPIT MAS+GBLUP 预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    return(GAPIT_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all))
  }

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  background_markers <- filter_background_markers(
    significant_markers = significant_markers,
    geno_all = GD_all[, -1, drop = FALSE],
    GM_all = GM_all,
    physical_distance = physical_distance,
    ld_threshold = ld_threshold
  )

  if (length(background_markers) == 0) {
    return(GAPIT_MAS_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all,
                                       significant_markers))
  }

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = NA_real_,
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  available_sig <- intersect(significant_markers, colnames(GD_all))
  available_bg  <- intersect(background_markers, colnames(GD_all))

  if (length(available_sig) == 0) {
    return(GAPIT_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all))
  }

  cv_matrix <- GD_all[all_lines, c("taxa", available_sig), drop = FALSE]
  colnames(cv_matrix)[1] <- "Taxa"

  if (length(available_bg) > 0) {
    GD_bg <- GD_all[all_lines, c("taxa", available_bg), drop = FALSE]
    GM_bg <- GM_all[GM_all$SNP %in% available_bg, , drop = FALSE]
    GM_bg <- GM_bg[match(colnames(GD_bg)[-1], GM_bg$SNP), , drop = FALSE]
  } else {
    return(GAPIT_MAS_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all,
                                       significant_markers))
  }

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_MASgBLUP_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  on.exit({
    setwd(original_wd)
    unlink(temp_dir, recursive = TRUE)
  }, add = TRUE)

  setwd(temp_dir)

  tryCatch({
    gapit_result <- GAPIT(
      Y = Y_all,
      G = NULL,
      GD = GD_bg,
      GM = GM_bg,
      CV = cv_matrix,
      model = "gBLUP",
      PCA.total = 0,
      file.output = TRUE,
      Major.allele.zero = FALSE,
      memo = "MAS+gBLUP_NoPCA"
    )

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        return(pred_values)
      }
    }
    safe_mean_pred(y_train, test_lines)
  }, error = function(e) {
    cat("        GAPIT MAS+GBLUP 预测失败:", e$message, "\n")
    safe_mean_pred(y_train, test_lines)
  })
}

# ------------------------------ 7. BGLR 预测函数 ------------------------------

# 7.1 BGLR MAS
MAS_prediction <- function(y_train, train_lines, test_lines, geno_all_markers,
                           significant_markers,
                           nIter = config$nIter,
                           burnIn = config$burnIn,
                           thin = config$thin,
                           seed = config$seed) {
  cat("      执行 BGLR MAS 预测...\n")

  if (length(significant_markers) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, rownames(geno_all_markers))

  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  available_markers <- intersect(significant_markers, colnames(geno_all_markers))
  if (length(available_markers) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  test_geno <- geno_all_markers[test_lines, available_markers, drop = FALSE]
  test_var <- apply(test_geno, 2, var, na.rm = TRUE)
  test_var[is.na(test_var)] <- 0
  if (length(test_var) > 0 && all(test_var == 0)) {
    cat("        警告：所有显著标记在测试集中均为单态，预测可能效果不佳\n")
  }

  X <- data.matrix(geno_all_markers[all_lines, available_markers, drop = FALSE])
  rownames(X) <- all_lines

  y_all <- rep(NA_real_, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx  <- which(all_lines %in% test_lines)

  tryCatch({
    set.seed(seed)

    ETA <- list(
      markers = list(X = X, model = "BRR", saveEffects = FALSE)
    )

    fm <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      verbose = FALSE
    )

    pred <- fm$yHat[test_idx]
    names(pred) <- all_lines[test_idx]
    pred
  }, error = function(e) {
    cat("        BGLR BRR失败:", e$message, "\n尝试岭回归...\n")

    tryCatch({
      X_train <- X[train_idx, , drop = FALSE]
      X_test  <- X[test_idx, , drop = FALSE]
      y_tr <- as.numeric(y_train[train_lines])

      cv_fit <- cv.glmnet(
        x = X_train,
        y = y_tr,
        alpha = 0,
        nfolds = min(5, length(y_tr))
      )

      pred <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.min"))
      names(pred) <- all_lines[test_idx]
      pred
    }, error = function(e2) {
      cat("        岭回归也失败:", e2$message, "\n")
      safe_mean_pred(y_train, test_lines)
    })
  })
}

# 7.2 BGLR GBLUP
GBLUP_prediction <- function(y_train, train_lines, test_lines, geno_all_markers,
                             nIter = config$nIter,
                             burnIn = config$burnIn,
                             thin = config$thin,
                             seed = config$seed) {
  cat("      执行 BGLR GBLUP 预测...\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, rownames(geno_all_markers))

  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  geno_sub <- geno_all_markers[all_lines, , drop = FALSE]

  K <- tryCatch(
    calc_kinship_matrix(geno_sub),
    error = function(e) {
      cat("        亲缘矩阵计算失败:", e$message, "\n")
      NULL
    }
  )

  if (is.null(K)) {
    return(safe_mean_pred(y_train, test_lines))
  }

  y_all <- rep(NA_real_, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx  <- which(all_lines %in% test_lines)

  tryCatch({
    set.seed(seed)

    ETA <- list(
      K = list(K = K, model = "RKHS", saveEffects = FALSE)
    )

    fm <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      verbose = FALSE
    )

    pred <- fm$yHat[test_idx]
    names(pred) <- all_lines[test_idx]
    pred
  }, error = function(e) {
    cat("        BGLR RKHS失败:", e$message, "\n尝试rrBLUP备用...\n")

    tryCatch({
      X_train <- data.matrix(geno_sub[train_lines, , drop = FALSE])
      X_test  <- data.matrix(geno_sub[test_lines, , drop = FALSE])

      p_train <- colMeans(X_train, na.rm = TRUE) / 2
      X_train_center <- sweep(X_train, 2, 2 * p_train, "-")
      X_test_center  <- sweep(X_test, 2, 2 * p_train, "-")

      y_tr <- as.numeric(y_train[train_lines])

      mix <- rrBLUP::mixed.solve(y = y_tr, Z = X_train_center, method = "REML")
      pred <- as.numeric(mix$beta + X_test_center %*% mix$u)
      names(pred) <- test_lines
      pred
    }, error = function(e2) {
      cat("        rrBLUP也失败:", e2$message, "\n")
      safe_mean_pred(y_train, test_lines)
    })
  })
}

# 7.3 BGLR MAS+GBLUP
MAS_GBLUP_prediction <- function(y_train, train_lines, test_lines,
                                 geno_all_markers, GM_all,
                                 significant_markers,
                                 physical_distance = 3e6,
                                 ld_threshold = 0.7,
                                 nIter = config$nIter,
                                 burnIn = config$burnIn,
                                 thin = config$thin,
                                 seed = config$seed) {
  cat("      执行 BGLR MAS+GBLUP 预测...\n")

  if (length(significant_markers) == 0) {
    return(GBLUP_prediction(
      y_train = y_train,
      train_lines = train_lines,
      test_lines = test_lines,
      geno_all_markers = geno_all_markers,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      seed = seed
    ))
  }

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, rownames(geno_all_markers))

  train_lines <- intersect(train_lines, all_lines)
  test_lines  <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  if (length(train_lines) < 5 || length(test_lines) == 0) {
    return(safe_mean_pred(y_train, test_lines))
  }

  available_sig <- intersect(significant_markers, colnames(geno_all_markers))
  if (length(available_sig) == 0) {
    return(GBLUP_prediction(
      y_train = y_train,
      train_lines = train_lines,
      test_lines = test_lines,
      geno_all_markers = geno_all_markers,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      seed = seed
    ))
  }

  random_markers <- filter_background_markers(
    significant_markers = available_sig,
    geno_all = geno_all_markers,
    GM_all = GM_all,
    physical_distance = physical_distance,
    ld_threshold = ld_threshold
  )

  random_markers <- setdiff(random_markers, available_sig)

  if (length(random_markers) == 0) {
    return(MAS_prediction(
      y_train = y_train,
      train_lines = train_lines,
      test_lines = test_lines,
      geno_all_markers = geno_all_markers,
      significant_markers = available_sig,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      seed = seed
    ))
  }

  X_fixed <- data.matrix(geno_all_markers[all_lines, available_sig, drop = FALSE])
  rownames(X_fixed) <- all_lines

  geno_rand <- geno_all_markers[all_lines, random_markers, drop = FALSE]
  K_rand <- tryCatch(
    calc_kinship_matrix(geno_rand),
    error = function(e) {
      cat("        随机背景亲缘矩阵计算失败:", e$message, "\n")
      NULL
    }
  )

  if (is.null(K_rand)) {
    return(MAS_prediction(
      y_train = y_train,
      train_lines = train_lines,
      test_lines = test_lines,
      geno_all_markers = geno_all_markers,
      significant_markers = available_sig,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      seed = seed
    ))
  }

  y_all <- rep(NA_real_, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx  <- which(all_lines %in% test_lines)

  tryCatch({
    set.seed(seed)

    ETA <- list(
      fixed = list(X = X_fixed, model = "BRR", saveEffects = FALSE),
      random = list(K = K_rand, model = "RKHS", saveEffects = FALSE)
    )

    fm <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      verbose = FALSE
    )

    pred <- fm$yHat[test_idx]
    names(pred) <- all_lines[test_idx]
    pred
  }, error = function(e) {
    cat("        BGLR混合模型失败:", e$message, "\n尝试BayesB备用...\n")

    tryCatch({
      all_markers <- unique(c(available_sig, random_markers))
      X_all <- data.matrix(geno_all_markers[all_lines, all_markers, drop = FALSE])
      rownames(X_all) <- all_lines

      ETA2 <- list(
        all = list(X = X_all, model = "BayesB", saveEffects = FALSE)
      )

      fm2 <- BGLR(
        y = y_all,
        ETA = ETA2,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = FALSE
      )

      pred <- fm2$yHat[test_idx]
      names(pred) <- all_lines[test_idx]
      pred
    }, error = function(e2) {
      cat("        BayesB也失败:", e2$message, "\n退化为GBLUP\n")
      GBLUP_prediction(
        y_train = y_train,
        train_lines = train_lines,
        test_lines = test_lines,
        geno_all_markers = geno_all_markers,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        seed = seed
      )
    })
  })
}

# ------------------------------ 8. 数据读取与预处理 ------------------------------
cat("\n=== 数据准备 ===\n")

# 8.1 训练表型
train_raw <- data.table::fread(config$train_pheno, header = TRUE, data.table = FALSE)
if (!all(c("genotype", "env_code", "TKW") %in% colnames(train_raw))) {
  if (ncol(train_raw) == 3) {
    colnames(train_raw) <- c("genotype", "env_code", "TKW")
  } else {
    stop("训练表型列名不符，且列数不为3")
  }
}
train_raw$genotype <- as.character(train_raw$genotype)
train_raw$env_code <- as.character(train_raw$env_code)
train_raw$TKW <- as.numeric(train_raw$TKW)
train_raw <- train_raw[!is.na(train_raw$TKW), , drop = FALSE]

# 8.2 验证表型
test_raw <- data.table::fread(config$test_pheno, header = TRUE, data.table = FALSE)
if (!all(c("genotype", "env_code", "TKW") %in% colnames(test_raw))) {
  if (ncol(test_raw) == 3) {
    colnames(test_raw) <- c("genotype", "env_code", "TKW")
  } else {
    stop("验证表型列名不符，且列数不为3")
  }
}
test_raw$genotype <- as.character(test_raw$genotype)
test_raw$env_code <- as.character(test_raw$env_code)
test_raw$TKW <- as.numeric(test_raw$TKW)
test_raw <- test_raw[!is.na(test_raw$TKW), , drop = FALSE]

# 8.3 环境参数
env_data <- read.csv(
  config$env_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(env_data)[1] <- "env_code"

kPara_Name <- config$env_para
if (!kPara_Name %in% colnames(env_data)) {
  kPara_Name <- gsub("&", ".", kPara_Name)
}
if (!kPara_Name %in% colnames(env_data)) {
  stop("无法匹配环境参数列：", config$env_para)
}

env_factors <- env_data[, c("env_code", kPara_Name), drop = FALSE]
colnames(env_factors)[2] <- "kPara"
env_factors$kPara <- as.numeric(env_factors$kPara)

env_mean <- aggregate(TKW ~ env_code, data = train_raw, mean, na.rm = TRUE)
colnames(env_mean)[2] <- "meanY"

# 8.4 基因型
geno_dt <- data.table::fread(config$geno_file, header = TRUE, data.table = FALSE)
if (ncol(geno_dt) < 2) stop("基因型文件格式错误：至少应包含样本列和1列标记")

sample_names <- as.character(geno_dt[[1]])
geno_matrix <- as.matrix(geno_dt[, -1, drop = FALSE])
rownames(geno_matrix) <- sample_names
storage.mode(geno_matrix) <- "double"

# 8.5 图谱
map_data <- data.table::fread(config$map_file, header = TRUE, data.table = FALSE)
if (ncol(map_data) < 3) stop("myGM.csv 至少需要 3 列")
colnames(map_data)[1:3] <- c("SNP", "Chromosome", "Position")

# 8.6 标准化SNP名
colnames(geno_matrix) <- standardize_snp_names(colnames(geno_matrix))
map_data$SNP <- standardize_snp_names(map_data$SNP)

# 8.7 染色体映射
map_data$Chromosome <- map_chr_to_numeric(map_data$Chromosome)

# 8.8 对齐基因型与图谱
common_snps <- intersect(colnames(geno_matrix), map_data$SNP)
cat("共同SNP数量:", length(common_snps), "\n")
if (length(common_snps) < 10) stop("共同SNP过少，无法继续分析")

geno_matrix <- geno_matrix[, common_snps, drop = FALSE]
map_data <- map_data[map_data$SNP %in% common_snps, , drop = FALSE]
map_data <- map_data[match(colnames(geno_matrix), map_data$SNP), , drop = FALSE]

# 8.9 去除重复标记
dup <- duplicated(colnames(geno_matrix))
if (any(dup)) {
  cat("检测到重复SNP标记数:", sum(dup), "，将保留首次出现者\n")
  geno_matrix <- geno_matrix[, !dup, drop = FALSE]
  map_data <- map_data[!dup, , drop = FALSE]
}

# 8.10 表型与基因型对齐
train_geno_ok <- intersect(unique(train_raw$genotype), rownames(geno_matrix))
test_geno_ok  <- intersect(unique(test_raw$genotype), rownames(geno_matrix))

if (length(train_geno_ok) < 5 || length(test_geno_ok) == 0) {
  stop("有效基因型不足：训练或验证集中与基因型文件对齐的样本太少")
}

train_raw <- train_raw[train_raw$genotype %in% train_geno_ok, , drop = FALSE]
test_raw  <- test_raw[test_raw$genotype %in% test_geno_ok, , drop = FALSE]

# 8.11 基因型质量控制
cat("\n=== 基因型质量控制 ===\n")

maf_vals <- colMeans(geno_matrix, na.rm = TRUE) / 2
maf <- pmin(maf_vals, 1 - maf_vals)
keep_maf <- which(maf >= config$maf_min)

geno_matrix <- geno_matrix[, keep_maf, drop = FALSE]
map_data <- map_data[keep_maf, , drop = FALSE]
cat("MAF过滤后剩余标记:", ncol(geno_matrix), "\n")

miss_rate <- apply(geno_matrix, 2, function(x) mean(is.na(x)))
keep_miss <- which(miss_rate <= config$pNA_max)

geno_matrix <- geno_matrix[, keep_miss, drop = FALSE]
map_data <- map_data[keep_miss, , drop = FALSE]
cat("缺失率过滤后剩余标记:", ncol(geno_matrix), "\n")

if (ncol(geno_matrix) < 10) {
  stop("QC后剩余标记过少，无法继续分析")
}

# 8.12 缺失值插补
for (j in seq_len(ncol(geno_matrix))) {
  if (any(is.na(geno_matrix[, j]))) {
    marker_mean <- mean(geno_matrix[, j], na.rm = TRUE)
    geno_matrix[, j] <- ifelse(is.na(geno_matrix[, j]), marker_mean, geno_matrix[, j])
  }
}

# 8.13 最终GM
GM_all <- map_data[, 1:3, drop = FALSE]
colnames(GM_all) <- c("SNP", "Chromosome", "Position")
GM_all$Chromosome <- as.integer(GM_all$Chromosome)
GM_all$Position <- as.integer(as.numeric(GM_all$Position))

# 8.14 最终纯marker矩阵
geno_all_markers <- geno_matrix[, , drop = FALSE]

# 8.15 GAPIT输入矩阵
GD_all_gapit <- data.frame(
  taxa = rownames(geno_all_markers),
  geno_all_markers,
  check.names = FALSE
)
rownames(GD_all_gapit) <- GD_all_gapit$taxa

stopifnot(identical(colnames(geno_all_markers), GM_all$SNP))
cat("最终纯marker矩阵维度:", dim(geno_all_markers), "\n")
cat("最终GAPIT输入GD维度:", dim(GD_all_gapit), "\n")

write.csv(geno_all_markers, file.path(config$out_dir, "geno_all_markers.csv"), row.names = TRUE)
write.csv(GM_all, file.path(config$out_dir, "GM_all.csv"), row.names = FALSE)

# ------------------------------ 9. 训练集FW参数计算 ------------------------------
cat("\n=== 训练集FW参数计算 ===\n")
fw_res <- calculate_FW_parameters(train_raw, env_factors, env_mean)
train_FW <- fw_res$FW_params
mean_kPara_train <- fw_res$mean_kPara_train

if (nrow(train_FW) == 0) {
  stop("未能从训练集中成功计算任何FW参数")
}

write.csv(train_FW, file.path(config$out_dir, "训练集FW参数.csv"), row.names = FALSE)

# ------------------------------ 10. GWAS分析（MLM） ------------------------------
cat("\n=== GWAS分析（MLM） ===\n")

traits <- c("Intcp_para_adj", "Slope_para")
gwas_dir <- file.path(config$out_dir, "GWAS_results")
dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)

gwas_results_list <- list()

for (trait in traits) {
  trait_dir <- file.path(gwas_dir, trait)
  res <- perform_gapit_gwas_mlm(
    train_FW_params = train_FW,
    GD_all_gapit = GD_all_gapit,
    GM_all = GM_all,
    trait = trait,
    output_dir = trait_dir
  )
  gwas_results_list[[trait]] <- res
}

# ------------------------------ 11. 从GWAS结果中提取显著标记 ------------------------------
cat("\n=== 提取显著标记 ===\n")
significant_markers_list <- list()

for (trait in traits) {
  if (!is.null(gwas_results_list[[trait]])) {
    trait_dir <- gwas_results_list[[trait]]$output_dir
    markers <- select_markers_by_pvalue_ld(
      gwas_dir = trait_dir,
      geno_all = geno_all_markers,
      GM_all = GM_all,
      pval_thresh = config$pval_thresh,
      hb_pval_thresh = config$hb_pval_thresh,
      ld_threshold = config$ld_threshold
    )
    significant_markers_list[[trait]] <- markers
  } else {
    significant_markers_list[[trait]] <- character(0)
  }
}

cat("Intcp_para_adj 显著标记数:", length(significant_markers_list[["Intcp_para_adj"]]), "\n")
cat("Slope_para 显著标记数:", length(significant_markers_list[["Slope_para"]]), "\n")

writeLines(
  text = c(
    paste0("Intcp_para_adj: ", paste(significant_markers_list[["Intcp_para_adj"]], collapse = ",")),
    paste0("Slope_para: ", paste(significant_markers_list[["Slope_para"]], collapse = ","))
  ),
  con = file.path(config$out_dir, "显著标记汇总.txt")
)

# ------------------------------ 12. 准备验证集环境参数 ------------------------------
cat("\n=== 验证集环境参数准备 ===\n")

test_env_para <- env_factors[env_factors$env_code %in% unique(test_raw$env_code), , drop = FALSE]
test_env_para$kPara_centered <- test_env_para$kPara - mean_kPara_train

test_raw <- merge(
  test_raw,
  test_env_para[, c("env_code", "kPara", "kPara_centered"), drop = FALSE],
  by = "env_code",
  all.x = TRUE
)

if (any(is.na(test_raw$kPara))) {
  warning("部分验证环境无环境参数，将剔除这些记录")
  test_raw <- test_raw[!is.na(test_raw$kPara), , drop = FALSE]
}

test_genotypes_all <- unique(test_raw$genotype)
cat("验证集参与预测基因型数:", length(test_genotypes_all), "\n")

if (length(test_genotypes_all) == 0) {
  stop("验证集经环境参数匹配后无可用基因型")
}

# ------------------------------ 13. 并行设置 ------------------------------
cat("\n=== 启动并行 ===\n")

cl <- parallel::makeCluster(config$ncores)
doParallel::registerDoParallel(cl)

parallel::clusterExport(
  cl,
  varlist = c(
    "config",
    "gapit_lib",
    "gapit_fun",
    "train_FW",
    "test_raw",
    "test_genotypes_all",
    "geno_all_markers",
    "GD_all_gapit",
    "GM_all",
    "significant_markers_list",
    "safe_mean_pred",
    "calc_kinship_matrix",
    "filter_background_markers",
    "extract_gapit_prediction",
    "GAPIT_MAS_prediction_simple",
    "GAPIT_GBLUP_prediction_simple",
    "GAPIT_MAS_gBLUP_prediction_simple",
    "MAS_prediction",
    "GBLUP_prediction",
    "MAS_GBLUP_prediction"
  ),
  envir = environment()
)

parallel::clusterEvalQ(cl, {
  setwd(config$workdir)
  library(data.table)
  library(BGLR)
  library(glmnet)
  library(rrBLUP)
  library(Matrix)
  library(dplyr)
  source(gapit_lib)
  source(gapit_fun)
  NULL
})

on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
}, add = TRUE)

# ------------------------------ 14. 定义方法组合 ------------------------------
methods_list <- c("MAS", "GBLUP", "MAS+GBLUP")
pred_types <- c("截距", "截距+斜率")
frameworks <- c("BGLR", "GAPIT")

param_grid <- expand.grid(
  framework = frameworks,
  method = methods_list,
  pred_type = pred_types,
  stringsAsFactors = FALSE
)

cat("\n总计", nrow(param_grid), "种配置 (",
    length(frameworks), "框架 × ",
    length(methods_list), "方法 × ",
    length(pred_types), "预测类型)\n", sep = "")

# ------------------------------ 15. 并行预测 ------------------------------
results_list <- foreach::foreach(
  k = seq_len(nrow(param_grid)),
  .packages = c("data.table", "BGLR", "glmnet", "rrBLUP", "Matrix")
) %dopar% {

  setwd(config$workdir)
  source(gapit_lib)
  source(gapit_fun)

  framework <- param_grid$framework[k]
  method    <- param_grid$method[k]
  pred_type <- param_grid$pred_type[k]

  cat("  开始配置:", framework, "-", method, "-", pred_type, "\n")

  if (pred_type == "截距") {
    trait <- "Intcp_para_adj"
    y_tr <- train_FW[[trait]]
    names(y_tr) <- train_FW$genotype
    y_tr <- y_tr[!is.na(y_tr)]

    if (length(y_tr) < 5) return(NULL)
    train_geno_valid <- names(y_tr)

    if (framework == "BGLR") {
      if (method == "MAS") {
        pred_intcp <- MAS_prediction(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          significant_markers = significant_markers_list[[trait]]
        )
      } else if (method == "GBLUP") {
        pred_intcp <- GBLUP_prediction(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers
        )
      } else if (method == "MAS+GBLUP") {
        pred_intcp <- MAS_GBLUP_prediction(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          GM_all = GM_all,
          significant_markers = significant_markers_list[[trait]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }
    } else if (framework == "GAPIT") {
      if (method == "MAS") {
        pred_intcp <- GAPIT_MAS_prediction_simple(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[[trait]]
        )
      } else if (method == "GBLUP") {
        pred_intcp <- GAPIT_GBLUP_prediction_simple(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all
        )
      } else if (method == "MAS+GBLUP") {
        pred_intcp <- GAPIT_MAS_gBLUP_prediction_simple(
          y_train = y_tr,
          train_lines = train_geno_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[[trait]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }

    pred_slope <- setNames(rep(NA_real_, length(test_genotypes_all)), test_genotypes_all)

  } else {  # 截距+斜率
    y_intcp <- train_FW$Intcp_para_adj
    names(y_intcp) <- train_FW$genotype
    y_intcp <- y_intcp[!is.na(y_intcp)]
    if (length(y_intcp) < 5) return(NULL)
    train_intcp_valid <- names(y_intcp)

    y_slope <- train_FW$Slope_para
    names(y_slope) <- train_FW$genotype
    y_slope <- y_slope[!is.na(y_slope)]
    if (length(y_slope) < 5) return(NULL)
    train_slope_valid <- names(y_slope)

    if (framework == "BGLR") {
      # 截距
      if (method == "MAS") {
        pred_intcp <- MAS_prediction(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          significant_markers = significant_markers_list[["Intcp_para_adj"]]
        )
      } else if (method == "GBLUP") {
        pred_intcp <- GBLUP_prediction(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers
        )
      } else if (method == "MAS+GBLUP") {
        pred_intcp <- MAS_GBLUP_prediction(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Intcp_para_adj"]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }

      # 斜率
      if (method == "MAS") {
        pred_slope <- MAS_prediction(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          significant_markers = significant_markers_list[["Slope_para"]]
        )
      } else if (method == "GBLUP") {
        pred_slope <- GBLUP_prediction(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers
        )
      } else if (method == "MAS+GBLUP") {
        pred_slope <- MAS_GBLUP_prediction(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          geno_all_markers = geno_all_markers,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Slope_para"]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }

    } else if (framework == "GAPIT") {
      # 截距
      if (method == "MAS") {
        pred_intcp <- GAPIT_MAS_prediction_simple(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Intcp_para_adj"]]
        )
      } else if (method == "GBLUP") {
        pred_intcp <- GAPIT_GBLUP_prediction_simple(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all
        )
      } else if (method == "MAS+GBLUP") {
        pred_intcp <- GAPIT_MAS_gBLUP_prediction_simple(
          y_train = y_intcp,
          train_lines = train_intcp_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Intcp_para_adj"]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }

      # 斜率
      if (method == "MAS") {
        pred_slope <- GAPIT_MAS_prediction_simple(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Slope_para"]]
        )
      } else if (method == "GBLUP") {
        pred_slope <- GAPIT_GBLUP_prediction_simple(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all
        )
      } else if (method == "MAS+GBLUP") {
        pred_slope <- GAPIT_MAS_gBLUP_prediction_simple(
          y_train = y_slope,
          train_lines = train_slope_valid,
          test_lines = test_genotypes_all,
          GD_all = GD_all_gapit,
          GM_all = GM_all,
          significant_markers = significant_markers_list[["Slope_para"]],
          physical_distance = config$physical_dist,
          ld_threshold = config$ld_threshold
        )
      } else {
        return(NULL)
      }

    } else {
      return(NULL)
    }
  }

  pred_df <- data.frame(
    genotype = test_genotypes_all,
    stringsAsFactors = FALSE
  )
  pred_df$FW_矫正截距 <- pred_intcp[match(test_genotypes_all, names(pred_intcp))]
  pred_df$FW_斜率 <- pred_slope[match(test_genotypes_all, names(pred_slope))]

  merged <- merge(test_raw, pred_df, by = "genotype", all.x = TRUE)

  if (pred_type == "截距") {
    merged$预测值 <- merged$FW_矫正截距
  } else {
    merged$预测值 <- merged$FW_矫正截距 +
      ifelse(is.na(merged$FW_斜率), 0, merged$FW_斜率) * merged$kPara_centered
  }

  merged <- merged[!is.na(merged$预测值), , drop = FALSE]
  if (nrow(merged) == 0) return(NULL)

  data.frame(
    框架 = framework,
    基因型ID = merged$genotype,
    环境ID = merged$env_code,
    实际值 = merged$TKW,
    预测值 = merged$预测值,
    FW_矫正截距 = merged$FW_矫正截距,
    FW_斜率 = merged$FW_斜率,
    方法 = method,
    预测类型 = pred_type,
    环境参数 = merged$kPara,
    中心化环境参数 = merged$kPara_centered,
    计算时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

all_results <- data.table::rbindlist(results_list, fill = TRUE)

# ------------------------------ 16. 汇总统计与可视化 ------------------------------
cat("\n=== 汇总统计与可视化 ===\n")

if (nrow(all_results) > 0) {
  data.table::fwrite(all_results, file.path(config$out_dir, "所有预测结果.csv"))

  summary_stats <- data.frame()

  for (fw in unique(all_results$框架)) {
    for (meth in unique(all_results$方法)) {
      for (pt in unique(all_results$预测类型)) {
        sub <- all_results[
          all_results$框架 == fw &
            all_results$方法 == meth &
            all_results$预测类型 == pt,
          , drop = FALSE
        ]

        if (nrow(sub) < 5) next

        valid <- complete.cases(sub$实际值, sub$预测值)
        if (sum(valid) < 5) next

        act <- sub$实际值[valid]
        pred <- sub$预测值[valid]

        mt <- calc_metrics(act, pred)

        summary_stats <- rbind(summary_stats, data.frame(
          框架 = fw,
          方法 = meth,
          预测类型 = pt,
          记录数 = nrow(sub),
          有效记录数 = mt$n,
          相关性 = round(mt$Correlation, 4),
          RMSE = round(mt$RMSE, 4),
          MAE = round(mt$MAE, 4),
          偏差 = round(mt$Bias, 4),
          唯一预测值个数 = length(unique(round(pred, 8))),
          备注 = "独立验证_MLM版_12配置",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  data.table::fwrite(summary_stats, file.path(config$out_dir, "验证集性能汇总.csv"))

  if (nrow(summary_stats) > 0) {
    cat("\n性能汇总：\n")
    print(summary_stats)

    best_idx <- which.max(summary_stats$相关性)
    if (length(best_idx) > 0 && !is.na(best_idx)) {
      best_row <- summary_stats[best_idx, , drop = FALSE]
      cat("\n最佳方法：\n")
      print(best_row)
      write.csv(best_row, file.path(config$out_dir, "最佳方法.csv"), row.names = FALSE)
    }
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(summary_stats) > 0) {
    p1 <- ggplot2::ggplot(summary_stats, aes(x = 方法, y = 相关性, fill = 框架)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      facet_wrap(~预测类型) +
      labs(title = "独立验证 - 各方法相关性比较", x = "方法", y = "相关性") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(config$out_dir, "相关性比较.png"), p1, width = 14, height = 7, dpi = 300)

    p2 <- ggplot2::ggplot(summary_stats, aes(x = 方法, y = RMSE, fill = 框架)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      facet_wrap(~预测类型) +
      labs(title = "独立验证 - 各方法RMSE比较", x = "方法", y = "RMSE") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(config$out_dir, "RMSE比较.png"), p2, width = 14, height = 7, dpi = 300)

    p3 <- ggplot2::ggplot(summary_stats, aes(x = 方法, y = MAE, fill = 框架)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      facet_wrap(~预测类型) +
      labs(title = "独立验证 - 各方法MAE比较", x = "方法", y = "MAE") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(config$out_dir, "MAE比较.png"), p3, width = 14, height = 7, dpi = 300)
  }

  cat("结果已保存至:", config$out_dir, "\n")
} else {
  cat("警告：没有生成任何有效预测结果！\n")
}

cat("\n=== 分析完成 ===\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")