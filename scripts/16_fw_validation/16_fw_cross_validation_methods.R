#!/usr/bin/env Rscript

# ===========================================
# 基因组选择分析完整代码 - 最终定稿版
# 预测阶段无PCA，GWAS阶段含PCA - 针对TKW性状
#
# 修正内容：
# 1. 修复字符串分隔线输出错误（使用 strrep）
# 2. 修复 BGLR_MAS_gBLUP_prediction_simple 对 GM_all 的隐式全局依赖
# 3. GAPIT 优先本地加载，网络失败时再尝试远程
# 4. 统一 GAPIT 预测结果列提取逻辑
# 5. 质量控制阶段按 SNP 名称同步过滤，降低错位风险
# 6. 修正 MAS+GBLUP 背景标记过滤逻辑：
#    删除与显著标记同染色体、3Mb以内、且 LD >= 0.7 的背景标记
# 7. 修正汇总统计中最佳方法筛选逻辑：
#    使用精确匹配，不再把 MAS 和 MAS+GBLUP 混在一起
# 8. 保留原始分析流程、输出结构和执行顺序
# ===========================================

# 设置工作目录
setwd("/mnt/7t_storage/zhangcl/TKW")

cat("=== TKW基因组选择分析开始（最终定稿版）===\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("执行顺序：1.环境留一基因型5折 → 2.基因型5折 → 3.环境留一\n")

# ===========================================
# 1. 初始设置和包加载
# ===========================================

cat("\n=== 初始设置和包加载 ===\n")

num_cores <- parallel::detectCores()
cat("检测到CPU核心数:", num_cores, "\n")
used_cores <- max(1, floor(num_cores * 0.5))
cat("将使用", used_cores, "个CPU核心进行并行计算\n")

cat("加载必要的R包...\n")
required_packages <- c(
  "data.table", "BGLR", "glmnet", "caret", "ggplot2",
  "vcfR", "rrBLUP", "parallel", "doParallel", "foreach",
  "Matrix", "dplyr", "tidyr", "reshape2"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装包:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
    cat("成功加载包:", pkg, "\n")
  } else {
    cat("包已加载:", pkg, "\n")
  }
}

if (!require("BGLR", character.only = TRUE, quietly = TRUE)) {
  cat("安装BGLR包...\n")
  install.packages("BGLR", dependencies = TRUE, repos = "https://cloud.r-project.org")
  library("BGLR")
  cat("成功加载BGLR包\n")
} else {
  cat("BGLR包已加载\n")
}

cat("注册并行后端...\n")
cl <- makeCluster(used_cores)
doParallel::registerDoParallel(cl)

# 在并行worker上预先加载GAPIT和所需包
cat("在并行worker上加载GAPIT函数和包...\n")
clusterEvalQ(cl, {
  if (file.exists("GAPIT.library.R") && file.exists("gapit_functions.txt")) {
    source("GAPIT.library.R")
    source("gapit_functions.txt")
  } else {
    tryCatch({
      source("http://zzlab.net/GAPIT/GAPIT.library.R")
      source("http://zzlab.net/GAPIT/gapit_functions.txt")
    }, error = function(e) {
      stop("并行worker无法加载GAPIT，请将 GAPIT.library.R 和 gapit_functions.txt 放到当前目录。")
    })
  }
  library(data.table)
  library(BGLR)
  library(glmnet)
  library(rrBLUP)
  library(dplyr)
  NULL
})
cat("并行后端注册完成\n")

# ===========================================
# 2. 加载GAPIT函数（主进程也需要）
# ===========================================

cat("\n=== 加载GAPIT函数 ===\n")

if (file.exists("GAPIT.library.R") && file.exists("gapit_functions.txt")) {
  source("GAPIT.library.R")
  source("gapit_functions.txt")
  cat("从本地加载GAPIT成功\n")
} else {
  cat("本地未发现GAPIT文件，尝试从网络加载GAPIT...\n")
  tryCatch({
    source("http://zzlab.net/GAPIT/GAPIT.library.R")
    source("http://zzlab.net/GAPIT/gapit_functions.txt")
    cat("网络GAPIT函数加载成功\n")
  }, error = function(e) {
    cat("网络GAPIT加载失败:", e$message, "\n")
    stop("无法加载GAPIT函数。请将 GAPIT.library.R 和 gapit_functions.txt 放置在当前目录。")
  })
}

# ===========================================
# 3. 数据准备 - 针对TKW性状
# ===========================================

cat("\n=== 数据准备 ===\n")

set.seed(195021)

# 创建输出目录
gs_output_dir <- "TKW_GS_Complete_NoPCA"
if (!dir.exists(gs_output_dir)) {
  dir.create(gs_output_dir, recursive = TRUE)
  cat("创建输出目录:", gs_output_dir, "\n")
}

# 读取表型数据
cat("读取TKW表型数据...\n")
phenotype_data <- read.table("TKW_mean_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("原始TKW表型数据维度:", dim(phenotype_data), "\n")
cat("列名:", colnames(phenotype_data), "\n")

# 数据清洗 - 检查列名
expected_cols <- c("genotype", "env_code", "TKW")
if (!all(expected_cols %in% colnames(phenotype_data))) {
  cat("警告：表型数据列名不匹配。实际列名:", colnames(phenotype_data), "\n")
  cat("期望的列名:", expected_cols, "\n")

  if (length(colnames(phenotype_data)) == 3) {
    colnames(phenotype_data) <- expected_cols
    cat("已重命名列名\n")
  } else {
    stop("表型数据列数与期望不符")
  }
}

# 去除NA值
phenotype_data_clean <- phenotype_data[!is.na(phenotype_data$TKW), ]
cat("去除NA值后的TKW表型数据维度:", dim(phenotype_data_clean), "\n")

# 读取环境参数
cat("读取环境参数数据...\n")
ecs_data <- read.csv("EC8.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cat("环境参数数据维度:", dim(ecs_data), "\n")
cat("环境参数列名:", colnames(ecs_data), "\n")

# 重命名第一列为env_code
colnames(ecs_data)[1] <- "env_code"

# 环境过滤
env_to_remove <- c()
if (length(env_to_remove) > 0) {
  phenotype_data_filtered <- phenotype_data_clean[!phenotype_data_clean$env_code %in% env_to_remove, ]
  cat("环境过滤后的TKW表型数据维度:", dim(phenotype_data_filtered), "\n")
} else {
  phenotype_data_filtered <- phenotype_data_clean
  cat("没有需要过滤的环境\n")
}

# 读取基因型数据
cat("读取基因型数据...\n")
geno_dt <- fread("myGD2.csv", header = TRUE, data.table = FALSE)
cat("基因型数据维度:", dim(geno_dt), "\n")

sample_names <- geno_dt[[1]]
geno_matrix <- as.matrix(geno_dt[, -1])
rownames(geno_matrix) <- sample_names
cat("基因型矩阵维度:", dim(geno_matrix), "\n")
cat("前5个品种名:", head(sample_names, 5), "\n")

# 读取SNP图谱
cat("读取SNP图谱...\n")
map_data <- fread("myGM2.csv", header = TRUE, data.table = FALSE)
cat("SNP图谱数据维度:", dim(map_data), "\n")
cat("SNP图谱列名:", colnames(map_data), "\n")

if (ncol(map_data) >= 3) {
  colnames(map_data)[1:3] <- c("SNP", "Chromosome", "Position")
  cat("标准化后的SNP图谱列名:", colnames(map_data)[1:3], "\n")
} else {
  stop("SNP图谱数据至少需要3列")
}

# 环境参数
kPara_Name <- "DTR&mean&PostFlowering1_14"
cat("选择的环境参数:", kPara_Name, "\n")

if (!kPara_Name %in% colnames(ecs_data)) {
  cat("错误：环境参数", kPara_Name, "不在数据中\n")
  cat("可用的环境参数:", colnames(ecs_data), "\n")
  stop("请选择正确的环境参数")
}

envMeanPara <- ecs_data[, c("env_code", kPara_Name)]
colnames(envMeanPara)[2] <- "kPara"
cat("环境参数数据维度:", dim(envMeanPara), "\n")

# 计算环境均值
env_mean_trait <- aggregate(TKW ~ env_code, data = phenotype_data_filtered, mean, na.rm = TRUE)
colnames(env_mean_trait)[2] <- "meanY"
cat("环境均值数据维度:", dim(env_mean_trait), "\n")

# 数据对齐
pheno_lines <- unique(phenotype_data_filtered$genotype)
geno_lines <- rownames(geno_matrix)
cat("表型数据中的品种数量:", length(pheno_lines), "\n")
cat("基因型数据中的品种数量:", length(geno_lines), "\n")

common_lines <- intersect(pheno_lines, geno_lines)
cat("共同品种数量:", length(common_lines), "\n")

if (length(common_lines) < 10) {
  cat("警告：共同品种数量不足\n")
  cat("表型数据中的前10个品种:", head(pheno_lines, 10), "\n")
  cat("基因型数据中的前10个品种:", head(geno_lines, 10), "\n")
}

phenotype_data_final <- phenotype_data_filtered[phenotype_data_filtered$genotype %in% common_lines, ]
geno_matrix_final <- geno_matrix[common_lines, , drop = FALSE]
cat("最终表型数据维度:", dim(phenotype_data_final), "\n")
cat("最终基因型数据维度:", dim(geno_matrix_final), "\n")

# 检查SNP一致性
geno_snps <- colnames(geno_matrix_final)
map_snps <- map_data$SNP
common_snps <- intersect(geno_snps, map_snps)
cat("共同SNP数量:", length(common_snps), "\n")

if (length(common_snps) < 100) {
  cat("警告：共同SNP数量较少\n")
  cat("基因型数据中的前10个SNP:", head(geno_snps, 10), "\n")
  cat("图谱数据中的前10个SNP:", head(map_snps, 10), "\n")
}

geno_matrix_final <- geno_matrix_final[, common_snps, drop = FALSE]
map_data_final <- map_data[map_data$SNP %in% common_snps, , drop = FALSE]
map_data_final <- map_data_final[match(colnames(geno_matrix_final), map_data_final$SNP), , drop = FALSE]
cat("过滤后基因型矩阵维度:", dim(geno_matrix_final), "\n")
cat("过滤后图谱数据维度:", dim(map_data_final), "\n")

# 质量控制
cat("\n=== 基因型数据质量控制 ===\n")
MAF.min <- 0.03
pNA.max <- 0.1

geno_matrix_qc <- geno_matrix_final

# 尝试转成数值矩阵
geno_matrix_qc <- apply(geno_matrix_qc, 2, as.numeric)
geno_matrix_qc <- as.matrix(geno_matrix_qc)
rownames(geno_matrix_qc) <- rownames(geno_matrix_final)
colnames(geno_matrix_qc) <- colnames(geno_matrix_final)

# MAF过滤（默认基因型编码为0/1/2）
MAF_values <- apply(geno_matrix_qc, 2, function(x) mean(x, na.rm = TRUE) / 2)
keep_maf_snps <- names(MAF_values)[!(MAF_values < MAF.min | MAF_values > 1 - MAF.min)]
cat("MAF过滤前位点数:", ncol(geno_matrix_qc), "\n")
cat("MAF过滤后保留位点数:", length(keep_maf_snps), "\n")

geno_matrix_qc <- geno_matrix_qc[, keep_maf_snps, drop = FALSE]
map_data_final <- map_data_final[map_data_final$SNP %in% keep_maf_snps, , drop = FALSE]
map_data_final <- map_data_final[match(colnames(geno_matrix_qc), map_data_final$SNP), , drop = FALSE]

# 缺失率过滤
missing_rates <- apply(geno_matrix_qc, 2, function(x) mean(is.na(x)))
keep_na_snps <- names(missing_rates)[missing_rates <= pNA.max]
cat("缺失率过滤前位点数:", ncol(geno_matrix_qc), "\n")
cat("缺失率过滤后保留位点数:", length(keep_na_snps), "\n")

geno_matrix_qc <- geno_matrix_qc[, keep_na_snps, drop = FALSE]
map_data_final <- map_data_final[map_data_final$SNP %in% keep_na_snps, , drop = FALSE]
map_data_final <- map_data_final[match(colnames(geno_matrix_qc), map_data_final$SNP), , drop = FALSE]

cat("质量控制后基因型数据维度:", dim(geno_matrix_qc), "\n")

# 缺失值插补
cat("进行缺失值插补...\n")
for (j in seq_len(ncol(geno_matrix_qc))) {
  tmp <- geno_matrix_qc[, j]
  if (any(is.na(tmp))) {
    geno_matrix_qc[, j] <- ifelse(is.na(tmp), mean(tmp, na.rm = TRUE), tmp)
  }
}
cat("缺失值插补完成\n")

# ===========================================
# 4. 准备GAPIT格式数据
# ===========================================

cat("\n=== 准备GAPIT格式数据 ===\n")

standardize_snp_names <- function(names_vec) {
  names_vec <- gsub(";", "_", names_vec)
  names_vec <- gsub("-", "_", names_vec)
  names_vec <- gsub(":", "_", names_vec)
  names_vec <- gsub("/", "_", names_vec)
  names_vec <- gsub("\\|", "_", names_vec)
  names_vec <- gsub("\\(", "_", names_vec)
  names_vec <- gsub("\\)", "_", names_vec)
  names_vec
}

new_geno_snp_names <- standardize_snp_names(colnames(geno_matrix_qc))
colnames(geno_matrix_qc) <- new_geno_snp_names
new_map_snp_names <- standardize_snp_names(map_data_final$SNP)
map_data_final$SNP <- new_map_snp_names

# 去除重复标记
keep_indices <- !duplicated(colnames(geno_matrix_qc))
geno_matrix_qc <- geno_matrix_qc[, keep_indices, drop = FALSE]
map_data_final <- map_data_final[keep_indices, , drop = FALSE]

# GAPIT 所需格式
GD_all <- data.frame(taxa = rownames(geno_matrix_qc), geno_matrix_qc, check.names = FALSE)
rownames(GD_all) <- GD_all$taxa

GM_all <- map_data_final
colnames(GM_all)[1:3] <- c("SNP", "Chromosome", "Position")
GM_all$Chromosome <- as.numeric(as.character(GM_all$Chromosome))
GM_all$Position <- as.numeric(as.character(GM_all$Position))

cat("GAPIT GD数据维度:", dim(GD_all), "\n")
cat("GAPIT GM数据维度:", dim(GM_all), "\n")
cat("染色体分布:\n")
print(table(GM_all$Chromosome, useNA = "ifany"))

# ===========================================
# 5. 核心函数定义
# ===========================================

cat("\n=== 定义核心函数 ===\n")

# 5.1 FW参数计算函数
calculate_FW_parameters <- function(train_pheno, env_factors, env_means) {
  cat("    计算FW模型参数...\n")
  genotypes <- unique(train_pheno$genotype)
  cat("    训练集中品种数量:", length(genotypes), "\n")

  train_envs <- unique(train_pheno$env_code)
  train_env_factors <- env_factors[env_factors$env_code %in% train_envs, , drop = FALSE]
  mean_kPara_train <- mean(train_env_factors$kPara, na.rm = TRUE)
  cat("    训练集环境参数平均值:", mean_kPara_train, "\n")

  results <- data.frame()
  for (geno in genotypes) {
    geno_data <- subset(train_pheno, genotype == geno)

    if (nrow(geno_data) >= 3) {
      geno_data_merged <- merge(geno_data, env_factors, by = "env_code")
      geno_data_merged <- merge(geno_data_merged, env_means, by = "env_code")
      geno_data_clean <- geno_data_merged[complete.cases(geno_data_merged$TKW, geno_data_merged$kPara), , drop = FALSE]

      if (nrow(geno_data_clean) >= 2) {
        tryCatch({
          geno_data_clean$kPara_centered <- geno_data_clean$kPara - mean_kPara_train
          lm_para <- lm(TKW ~ kPara_centered, data = geno_data_clean)
          coefficients <- coef(lm_para)
          r_squared <- summary(lm_para)$r.squared

          intcp_para_adj <- as.numeric(coefficients[1])
          slope_para <- as.numeric(coefficients[2])

          results <- rbind(
            results,
            data.frame(
              genotype = geno,
              Intcp_para_adj = intcp_para_adj,
              Slope_para = slope_para,
              R2_para = r_squared,
              stringsAsFactors = FALSE
            )
          )
        }, error = function(e) {
          cat("      品种", geno, "的FW参数计算失败:", e$message, "\n")
        })
      }
    }
  }

  cat("    成功计算FW参数的品种数量:", nrow(results), "\n")

  return(list(
    FW_params = results,
    mean_kPara_train = mean_kPara_train
  ))
}

# 统一提取GAPIT预测列
extract_gapit_prediction <- function(pred_df, test_lines) {
  if (is.null(pred_df) || nrow(pred_df) == 0) {
    return(NULL)
  }

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
    if (length(numeric_cols) > 0) {
      pred_col <- numeric_cols[1]
    }
  }

  if (is.null(pred_col)) {
    return(NULL)
  }

  taxa_col <- if ("Taxa" %in% colnames(pred_df)) "Taxa" else if ("taxa" %in% colnames(pred_df)) "taxa" else NULL
  if (is.null(taxa_col)) {
    return(NULL)
  }

  pred_test <- pred_df[pred_df[[taxa_col]] %in% test_lines, c(taxa_col, pred_col), drop = FALSE]
  colnames(pred_test) <- c("Taxa", "Prediction")

  pred_values <- pred_test$Prediction
  names(pred_values) <- pred_test$Taxa

  pred_values
}

# 5.2 GAPIT GWAS函数
perform_gapit_gwas <- function(train_FW_params, GD_all, GM_all, trait, output_dir) {
  cat("    在训练集上进行", trait, "的GAPIT GWAS分析...\n")

  if (!is.data.frame(train_FW_params) || nrow(train_FW_params) == 0) {
    cat("    错误：train_FW_params 无效\n")
    return(NULL)
  }

  pheno_df <- data.frame(
    Taxa = train_FW_params$genotype,
    Trait = train_FW_params[[trait]],
    stringsAsFactors = FALSE
  )
  pheno_df <- pheno_df[complete.cases(pheno_df$Trait), , drop = FALSE]

  cat("    ", trait, "有效样本数量:", nrow(pheno_df), "\n")

  if (nrow(pheno_df) < 10) {
    cat("      有效样本不足，跳过", trait, "的GWAS\n")
    return(NULL)
  }

  common_genotypes <- intersect(pheno_df$Taxa, GD_all$taxa)
  cat("    共同基因型数量:", length(common_genotypes), "\n")

  if (length(common_genotypes) < 10) {
    cat("      共同基因型不足，跳过", trait, "的GWAS\n")
    return(NULL)
  }

  pheno_gapit <- pheno_df[pheno_df$Taxa %in% common_genotypes, , drop = FALSE]
  GD_gapit <- GD_all[GD_all$taxa %in% common_genotypes, , drop = FALSE]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  original_wd <- getwd()

  tryCatch({
    setwd(output_dir)

    cat("    运行GAPIT MLM分析...\n")
    gapit_result <- GAPIT(
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
      cutOff = 1,
      memo = paste("GWAS", trait, sep = "_")
    )

    setwd(original_wd)
    return(list(output_dir = output_dir, result = gapit_result))
  }, error = function(e) {
    setwd(original_wd)
    cat("    GAPIT分析失败:", e$message, "\n")
    return(NULL)
  })
}

# 5.3 从GWAS结果中提取显著标记（双阈值+LD过滤）
select_markers_by_pvalue_ld <- function(gwas_dir, geno_all, GM_all,
                                        pval_thresh = 0.0001,
                                        hb_pval_thresh = 0.05,
                                        ld_threshold = 0.7) {
  files <- list.files(gwas_dir, pattern = "\\.csv$", full.names = TRUE)
  gwas_files <- files[
    grepl("GWAS[._ ]?Results", basename(files), ignore.case = TRUE) |
      (grepl("GWAS", basename(files), ignore.case = TRUE) & grepl("Results", basename(files), ignore.case = TRUE))
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
    missing_cols <- setdiff(needed_cols, colnames(gwas_data))
    cat("    GWAS结果缺少列:", paste(missing_cols, collapse = ", "), "\n")
    return(character(0))
  }

  gwas_data$SNP <- standardize_snp_names(gwas_data$SNP)

  pval <- suppressWarnings(as.numeric(gwas_data[["P.value"]]))
  hb <- suppressWarnings(as.numeric(gwas_data[["H&B.P.Value"]]))
  eff <- suppressWarnings(as.numeric(gwas_data[["Effect"]]))
  chrv <- suppressWarnings(as.numeric(gwas_data[["Chr"]]))
  posv <- suppressWarnings(as.numeric(gwas_data[["Pos"]]))

  idx <- which((!is.na(pval) & pval < pval_thresh) |
                 (!is.na(hb) & hb < hb_pval_thresh))

  if (length(idx) == 0) {
    cat("    没有标记满足P值阈值\n")
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

  if (nrow(candidate) == 0) {
    cat("    候选标记经清理后为空\n")
    return(character(0))
  }

  candidate <- candidate[candidate$SNP %in% GM_all$SNP, , drop = FALSE]
  if (nrow(candidate) == 0) {
    cat("    候选标记不在GM_all中\n")
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

# 5.4 背景标记过滤（用于MAS+GBLUP）
# 修正版：删除同染色体、3Mb内、且 LD >= 0.7 的背景标记
filter_background_markers <- function(significant_markers, GD_all, GM_all,
                                      physical_distance = 3e6, ld_threshold = 0.7) {
  cat("        过滤背景标记（去除与显著标记相关的标记）...\n")

  if (length(significant_markers) == 0) {
    cat("        没有显著标记，返回所有标记\n")
    return(colnames(GD_all)[-1])
  }

  sig_info <- GM_all[GM_all$SNP %in% significant_markers, , drop = FALSE]

  if (nrow(sig_info) == 0) {
    cat("        显著标记不在图谱中，返回所有标记\n")
    return(colnames(GD_all)[-1])
  }

  all_markers <- colnames(GD_all)[-1]
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
        chr_markers$SNP != sig_snp, , drop = FALSE
    ]

    if (nrow(nearby_markers) == 0) next

    for (j in seq_len(nrow(nearby_markers))) {
      bg_snp <- nearby_markers$SNP[j]

      if (sig_snp %in% colnames(GD_all) && bg_snp %in% colnames(GD_all)) {
        geno_sig <- GD_all[, sig_snp]
        geno_bg <- GD_all[, bg_snp]

        ld_value <- suppressWarnings(cor(geno_sig, geno_bg, use = "complete.obs")^2)

        # 关键修正：由 > 改为 >=
        if (!is.na(ld_value) && ld_value >= ld_threshold) {
          markers_to_exclude <- c(markers_to_exclude, bg_snp)
        }
      }
    }
  }

  markers_to_exclude <- unique(markers_to_exclude)
  cat("        需要排除的标记数量:", length(markers_to_exclude), "\n")

  background_markers <- setdiff(all_markers, c(significant_markers, markers_to_exclude))
  cat("        背景标记数量:", length(background_markers), "\n")

  return(background_markers)
}

# ===========================================
# 6. GAPIT预测函数定义 - 预测无PCA
# ===========================================

cat("\n=== GAPIT预测函数定义 - 简化版本（预测无PCA） ===\n")

# 6.1 GAPIT MAS预测函数 - 无PCA
GAPIT_MAS_prediction_simple <- function(y_train, train_lines, test_lines, GD_all, GM_all,
                                        significant_markers, n_top_markers = 15) {
  cat("      执行GAPIT MAS预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    cat("        没有显著标记，MAS使用均值预测\n")
    mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
    names(mean_pred) <- test_lines
    return(mean_pred)
  }

  cat("        使用", length(significant_markers), "个显著标记进行MAS预测\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = rep(NA, length(all_lines)),
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  available_markers <- intersect(significant_markers, colnames(GD_all))
  cat("        在基因型数据中可用的标记数量:", length(available_markers), "\n")

  if (length(available_markers) == 0) {
    cat("        没有可用的标记，MAS使用均值预测\n")
    mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
    names(mean_pred) <- test_lines
    return(mean_pred)
  }

  markers_geno <- GD_all[all_lines, c("taxa", available_markers), drop = FALSE]
  colnames(markers_geno)[1] <- "Taxa"

  cat("        标记基因型矩阵维度:", dim(markers_geno), "\n")

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_MAS_")
  dir.create(temp_dir, recursive = TRUE)

  tryCatch({
    setwd(temp_dir)

    cat("        运行GAPIT GLM模型进行MAS预测...\n")
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

    setwd(original_wd)

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        cat("        MAS预测完成，测试集样本数:", length(pred_values), "\n")
        return(pred_values)
      }
    }

    cat("        GAPIT没有返回可识别的预测结果\n")
  }, error = function(e) {
    setwd(original_wd)
    cat("        GAPIT MAS预测失败:", e$message, "\n")
  }, finally = {
    unlink(temp_dir, recursive = TRUE)
  })

  cat("        MAS预测失败，使用均值预测\n")
  mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
  names(mean_pred) <- test_lines
  return(mean_pred)
}

# 6.2 GAPIT GBLUP预测函数 - 无PCA
GAPIT_GBLUP_prediction_simple <- function(y_train, train_lines, test_lines, GD_all, GM_all) {
  cat("      执行GAPIT GBLUP预测 (无PCA)...\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = rep(NA, length(all_lines)),
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_GBLUP_")
  dir.create(temp_dir, recursive = TRUE)

  tryCatch({
    setwd(temp_dir)

    cat("        运行GAPIT gBLUP模型...\n")
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

    setwd(original_wd)

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        cat("        GAPIT GBLUP预测完成，测试集样本数:", length(pred_values), "\n")
        return(pred_values)
      } else {
        cat("        GAPIT GBLUP没有返回可识别预测结果\n")
      }
    } else {
      cat("        GAPIT GBLUP没有返回预测结果\n")
    }
  }, error = function(e) {
    setwd(original_wd)
    cat("        GAPIT GBLUP预测失败:", e$message, "\n")
  }, finally = {
    unlink(temp_dir, recursive = TRUE)
  })

  mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
  names(mean_pred) <- test_lines
  return(mean_pred)
}

# 6.3 GAPIT MAS+GBLUP预测函数 - 无PCA
GAPIT_MAS_gBLUP_prediction_simple <- function(y_train, train_lines, test_lines, GD_all, GM_all,
                                              significant_markers, n_top_markers = 15,
                                              filter_distance = 3e6, filter_ld = 0.7) {
  cat("      执行GAPIT MAS+GBLUP预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    cat("        没有显著标记，MAS+GBLUP使用GBLUP预测\n")
    return(GAPIT_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all))
  }

  cat("        使用", length(significant_markers), "个显著标记进行MAS+GBLUP预测\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  background_markers <- filter_background_markers(
    significant_markers, GD_all, GM_all,
    physical_distance = filter_distance,
    ld_threshold = filter_ld
  )

  if (length(background_markers) == 0) {
    cat("        没有背景标记，MAS+GBLUP使用MAS预测\n")
    return(GAPIT_MAS_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all,
                                       significant_markers, n_top_markers))
  }

  cat("        背景标记数量:", length(background_markers), "\n")

  Y_all <- data.frame(
    Taxa = all_lines,
    Trait = rep(NA, length(all_lines)),
    stringsAsFactors = FALSE
  )
  Y_all$Trait[Y_all$Taxa %in% train_lines] <- y_train[train_lines]

  available_sig_markers <- intersect(significant_markers, colnames(GD_all))
  available_bg_markers <- intersect(background_markers, colnames(GD_all))

  cat("        可用的显著标记数量:", length(available_sig_markers), "\n")
  cat("        可用的背景标记数量:", length(available_bg_markers), "\n")

  if (length(available_sig_markers) == 0 && length(available_bg_markers) == 0) {
    cat("        没有可用的标记，使用GBLUP预测\n")
    return(GAPIT_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all))
  }

  if (length(available_sig_markers) > 0) {
    sig_geno <- GD_all[all_lines, c("taxa", available_sig_markers), drop = FALSE]
  } else {
    sig_geno <- data.frame(taxa = all_lines, stringsAsFactors = FALSE)
  }

  if (length(available_bg_markers) > 0) {
    GD_background <- GD_all[, c("taxa", available_bg_markers), drop = FALSE]
  } else {
    GD_background <- GD_all[, "taxa", drop = FALSE]
  }

  cv_matrix <- sig_geno
  colnames(cv_matrix)[1] <- "Taxa"

  cat("        CV矩阵维度:", dim(cv_matrix), "\n")

  original_wd <- getwd()
  temp_dir <- tempfile("GAPIT_MAS_gBLUP_")
  dir.create(temp_dir, recursive = TRUE)

  tryCatch({
    setwd(temp_dir)

    cat("        运行GAPIT gBLUP模型进行MAS+GBLUP预测...\n")
    gapit_result <- GAPIT(
      Y = Y_all,
      G = NULL,
      GD = GD_background[all_lines, , drop = FALSE],
      GM = GM_all[GM_all$SNP %in% available_bg_markers, , drop = FALSE],
      CV = cv_matrix,
      model = "gBLUP",
      SNP.test = FALSE,
      PCA.total = 0,
      file.output = TRUE,
      Major.allele.zero = FALSE,
      memo = "MAS+gBLUP_NoPCA"
    )

    setwd(original_wd)

    if (!is.null(gapit_result$Pred)) {
      pred_values <- extract_gapit_prediction(gapit_result$Pred, test_lines)
      if (!is.null(pred_values) && length(pred_values) > 0) {
        cat("        MAS+GBLUP预测完成，测试集样本数:", length(pred_values), "\n")
        return(pred_values)
      }
    }

    cat("        GAPIT MAS+gBLUP没有返回可识别预测结果\n")
  }, error = function(e) {
    setwd(original_wd)
    cat("        GAPIT MAS+gBLUP预测失败:", e$message, "\n")
  }, finally = {
    unlink(temp_dir, recursive = TRUE)
  })

  cat("        MAS+GBLUP预测失败，使用GBLUP预测\n")
  return(GAPIT_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all, GM_all))
}

# ===========================================
# 7. BGLR预测函数定义 - 简化版本（无PCA）
# ===========================================

cat("\n=== BGLR预测函数定义 - 简化版本（无PCA） ===\n")

# 7.1 计算亲缘关系矩阵函数
calculate_kinship_matrix <- function(geno_matrix, method = "VanRaden") {
  cat("      计算亲缘关系矩阵（方法：", method, ")...\n", sep = "")

  geno_numeric <- apply(geno_matrix, 2, as.numeric)
  geno_numeric <- as.matrix(geno_numeric)
  rownames(geno_numeric) <- rownames(geno_matrix)
  colnames(geno_numeric) <- colnames(geno_matrix)

  p <- colMeans(geno_numeric, na.rm = TRUE) / 2

  if (method == "VanRaden") {
    W <- geno_numeric - 2 * matrix(p, nrow = nrow(geno_numeric), ncol = ncol(geno_numeric), byrow = TRUE)
    denom <- sum(2 * p * (1 - p))
    if (denom <= 0 || is.na(denom)) {
      denom <- 1
    }
    K <- tcrossprod(W) / denom
  } else if (method == "linear") {
    K <- tcrossprod(geno_numeric) / ncol(geno_numeric)
  } else if (method == "gaussian") {
    dist_mat <- as.matrix(dist(geno_numeric, method = "euclidean"))
    sigma <- median(dist_mat)
    if (is.na(sigma) || sigma == 0) sigma <- 1
    K <- exp(-dist_mat^2 / (2 * sigma^2))
  } else {
    W <- geno_numeric - 2 * matrix(p, nrow = nrow(geno_numeric), ncol = ncol(geno_numeric), byrow = TRUE)
    denom <- sum(2 * p * (1 - p))
    if (denom <= 0 || is.na(denom)) {
      denom <- 1
    }
    K <- tcrossprod(W) / denom
  }

  K <- (K + t(K)) / 2
  diag(K) <- diag(K) + 0.01

  cat("      亲缘关系矩阵维度:", dim(K), "\n")
  return(K)
}

# 7.2 BGLR MAS预测函数
BGLR_MAS_prediction_simple <- function(y_train, train_lines, test_lines, GD_all,
                                       significant_markers, nIter = 25000, burnIn = 5000, thin = 10) {
  cat("      执行BGLR MAS预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    cat("        没有显著标记，BGLR MAS使用均值预测\n")
    mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
    names(mean_pred) <- test_lines
    return(mean_pred)
  }

  cat("        使用", length(significant_markers), "个显著标记进行BGLR MAS预测\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  available_markers <- intersect(significant_markers, colnames(GD_all))
  cat("        在基因型数据中可用的标记数量:", length(available_markers), "\n")

  if (length(available_markers) == 0) {
    cat("        没有可用的标记，BGLR MAS使用均值预测\n")
    mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
    names(mean_pred) <- test_lines
    return(mean_pred)
  }

  geno_markers <- GD_all[all_lines, available_markers, drop = FALSE]
  X_markers <- apply(geno_markers, 2, as.numeric)
  X_markers <- as.matrix(X_markers)
  rownames(X_markers) <- all_lines
  colnames(X_markers) <- available_markers

  y_all <- rep(NA, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx <- which(all_lines %in% test_lines)

  cat("        训练样本索引:", length(train_idx), "\n")
  cat("        测试样本索引:", length(test_idx), "\n")
  cat("        运行BGLR BRR模型...\n")
  cat("        迭代参数: nIter =", nIter, ", burnIn =", burnIn, ", thin =", thin, "\n")

  tryCatch({
    ETA <- list()
    ETA$markers <- list(
      X = X_markers,
      model = "BRR",
      saveEffects = FALSE
    )

    set.seed(195021)
    bglr_fit <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      saveAt = paste0(tempfile(), "_"),
      verbose = FALSE
    )

    pred_values <- bglr_fit$yHat[test_idx]
    names(pred_values) <- all_lines[test_idx]

    cat("        BGLR MAS预测完成，测试集样本数:", length(pred_values), "\n")
    cat("        预测值范围:", range(pred_values, na.rm = TRUE), "\n")
    return(pred_values)
  }, error = function(e) {
    cat("        BGLR MAS预测失败:", e$message, "\n")

    cat("        尝试使用岭回归作为备选方案...\n")
    tryCatch({
      X_train <- X_markers[train_idx, , drop = FALSE]
      X_test <- X_markers[test_idx, , drop = FALSE]

      y_train_numeric <- as.numeric(y_train[train_lines])

      cv_fit <- cv.glmnet(X_train, y_train_numeric, alpha = 0, nfolds = 5)
      pred_values <- predict(cv_fit, newx = X_test, s = "lambda.min")
      pred_values <- as.numeric(pred_values)
      names(pred_values) <- all_lines[test_idx]

      cat("        岭回归预测完成\n")
      return(pred_values)
    }, error = function(e2) {
      cat("        岭回归也失败:", e2$message, "\n")

      mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
      names(mean_pred) <- test_lines
      return(mean_pred)
    })
  })
}

# 7.3 BGLR GBLUP预测函数
BGLR_GBLUP_prediction_simple <- function(y_train, train_lines, test_lines, GD_all,
                                         nIter = 25000, burnIn = 5000, thin = 10) {
  cat("      执行BGLR GBLUP预测 (无PCA)...\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  geno_all <- GD_all[all_lines, -1, drop = FALSE]

  cat("        计算亲缘关系矩阵...\n")
  K <- calculate_kinship_matrix(geno_all, method = "VanRaden")

  y_all <- rep(NA, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx <- which(all_lines %in% test_lines)

  cat("        训练样本索引:", length(train_idx), "\n")
  cat("        测试样本索引:", length(test_idx), "\n")
  cat("        运行BGLR RKHS模型...\n")
  cat("        迭代参数: nIter =", nIter, ", burnIn =", burnIn, ", thin =", thin, "\n")

  tryCatch({
    ETA <- list()
    ETA$K <- list(
      K = K,
      model = "RKHS",
      saveEffects = FALSE
    )

    set.seed(195021)
    bglr_fit <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      saveAt = paste0(tempfile(), "_"),
      verbose = FALSE
    )

    pred_values <- bglr_fit$yHat[test_idx]
    names(pred_values) <- all_lines[test_idx]

    cat("        BGLR GBLUP预测完成，测试集样本数:", length(pred_values), "\n")
    cat("        预测值范围:", range(pred_values, na.rm = TRUE), "\n")
    return(pred_values)
  }, error = function(e) {
    cat("        BGLR GBLUP预测失败:", e$message, "\n")

    cat("        尝试使用rrBLUP作为备选方案...\n")
    tryCatch({
      library(rrBLUP)

      pheno_train <- data.frame(
        line = train_lines,
        y = y_train[train_lines],
        stringsAsFactors = FALSE
      )

      geno_train <- geno_all[train_idx, , drop = FALSE]
      geno_test <- geno_all[test_idx, , drop = FALSE]

      mix_model <- mixed.solve(
        y = pheno_train$y,
        Z = geno_train,
        K = NULL,
        X = NULL,
        method = "REML"
      )

      pred_values <- as.numeric(geno_test %*% mix_model$u)
      names(pred_values) <- all_lines[test_idx]

      cat("        rrBLUP预测完成\n")
      return(pred_values)
    }, error = function(e2) {
      cat("        rrBLUP也失败:", e2$message, "\n")

      mean_pred <- rep(mean(y_train, na.rm = TRUE), length(test_lines))
      names(mean_pred) <- test_lines
      return(mean_pred)
    })
  })
}

# 7.4 BGLR MAS+GBLUP预测函数
BGLR_MAS_gBLUP_prediction_simple <- function(y_train, train_lines, test_lines, GD_all, GM_all,
                                             significant_markers, nIter = 25000, burnIn = 5000,
                                             thin = 10, filter_distance = 3e6, filter_ld = 0.7) {
  cat("      执行BGLR MAS+GBLUP预测 (无PCA)...\n")

  if (length(significant_markers) == 0) {
    cat("        没有显著标记，BGLR MAS+gBLUP使用GBLUP预测\n")
    return(BGLR_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all,
                                        nIter, burnIn, thin))
  }

  cat("        使用", length(significant_markers), "个显著标记进行BGLR MAS+gBLUP预测\n")

  all_lines <- unique(c(train_lines, test_lines))
  all_lines <- intersect(all_lines, GD_all$taxa)
  train_lines <- intersect(train_lines, all_lines)
  test_lines <- intersect(test_lines, all_lines)

  y_train <- y_train[train_lines]

  cat("        所有品种数量:", length(all_lines), "\n")
  cat("        训练品种数量:", length(train_lines), "\n")
  cat("        测试品种数量:", length(test_lines), "\n")

  background_markers <- filter_background_markers(
    significant_markers, GD_all, GM_all,
    physical_distance = filter_distance,
    ld_threshold = filter_ld
  )

  if (length(background_markers) == 0) {
    cat("        没有背景标记，BGLR MAS+gBLUP使用MAS预测\n")
    return(BGLR_MAS_prediction_simple(y_train, train_lines, test_lines, GD_all,
                                      significant_markers, nIter, burnIn, thin))
  }

  cat("        背景标记数量:", length(background_markers), "\n")

  available_sig_markers <- intersect(significant_markers, colnames(GD_all))
  cat("        在基因型数据中可用的显著标记数量:", length(available_sig_markers), "\n")

  if (length(available_sig_markers) == 0) {
    cat("        没有可用的显著标记，使用GBLUP预测\n")
    return(BGLR_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all,
                                        nIter, burnIn, thin))
  }

  geno_sig <- GD_all[all_lines, available_sig_markers, drop = FALSE]
  X_sig <- apply(geno_sig, 2, as.numeric)
  X_sig <- as.matrix(X_sig)
  rownames(X_sig) <- all_lines
  colnames(X_sig) <- available_sig_markers

  available_bg_markers <- intersect(background_markers, colnames(GD_all))
  cat("        在基因型数据中可用的背景标记数量:", length(available_bg_markers), "\n")

  if (length(available_bg_markers) == 0) {
    cat("        没有可用的背景标记，使用MAS预测\n")
    return(BGLR_MAS_prediction_simple(y_train, train_lines, test_lines, GD_all,
                                      significant_markers, nIter, burnIn, thin))
  }

  geno_bg <- GD_all[all_lines, available_bg_markers, drop = FALSE]

  cat("        计算背景标记的亲缘关系矩阵...\n")
  K_bg <- calculate_kinship_matrix(geno_bg, method = "VanRaden")

  y_all <- rep(NA, length(all_lines))
  names(y_all) <- all_lines
  y_all[train_lines] <- y_train[train_lines]

  train_idx <- which(all_lines %in% train_lines)
  test_idx <- which(all_lines %in% test_lines)

  cat("        训练样本索引:", length(train_idx), "\n")
  cat("        测试样本索引:", length(test_idx), "\n")
  cat("        运行BGLR组合模型（BRR + RKHS）...\n")
  cat("        迭代参数: nIter =", nIter, ", burnIn =", burnIn, ", thin =", thin, "\n")

  tryCatch({
    ETA <- list()
    ETA$markers_sig <- list(
      X = X_sig,
      model = "BRR",
      saveEffects = FALSE
    )
    ETA$kinship_bg <- list(
      K = K_bg,
      model = "RKHS",
      saveEffects = FALSE
    )

    set.seed(195021)
    bglr_fit <- BGLR(
      y = y_all,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      saveAt = paste0(tempfile(), "_"),
      verbose = FALSE
    )

    pred_values <- bglr_fit$yHat[test_idx]
    names(pred_values) <- all_lines[test_idx]

    cat("        BGLR MAS+gBLUP预测完成，测试集样本数:", length(pred_values), "\n")
    cat("        预测值范围:", range(pred_values, na.rm = TRUE), "\n")
    return(pred_values)
  }, error = function(e) {
    cat("        BGLR MAS+gBLUP预测失败:", e$message, "\n")

    cat("        尝试使用BayesB模型作为备选方案...\n")
    tryCatch({
      all_available_markers <- unique(c(available_sig_markers, available_bg_markers))
      geno_all_markers <- GD_all[all_lines, all_available_markers, drop = FALSE]
      X_all_markers <- apply(geno_all_markers, 2, as.numeric)
      X_all_markers <- as.matrix(X_all_markers)
      rownames(X_all_markers) <- all_lines
      colnames(X_all_markers) <- all_available_markers

      ETA_simple <- list()
      ETA_simple$all_markers <- list(
        X = X_all_markers,
        model = "BayesB",
        saveEffects = FALSE
      )

      set.seed(195021)
      bglr_fit2 <- BGLR(
        y = y_all,
        ETA = ETA_simple,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        saveAt = paste0(tempfile(), "_"),
        verbose = FALSE
      )

      pred_values <- bglr_fit2$yHat[test_idx]
      names(pred_values) <- all_lines[test_idx]

      cat("        BGLR BayesB预测完成\n")
      return(pred_values)
    }, error = function(e2) {
      cat("        BGLR BayesB也失败:", e2$message, "\n")

      cat("        使用BGLR GBLUP作为最终备用\n")
      return(BGLR_GBLUP_prediction_simple(y_train, train_lines, test_lines, GD_all,
                                          nIter, burnIn, thin))
    })
  })
}

# ===========================================
# 8. 通用交叉验证函数定义（修正版）
# ===========================================

cat("\n=== 通用交叉验证函数定义（修正版）===\n")

# 8.1 环境留一基因型5折组合交叉验证通用函数
run_geno5_envloo_cv_simple <- function(pheno_data, GD_all, GM_all, env_factors, env_means, nfolds = 5,
                                       methods = c("GAPIT_MAS", "GAPIT_GBLUP", "GAPIT_MAS+GBLUP",
                                                   "BGLR_MAS", "BGLR_GBLUP", "BGLR_MAS+GBLUP"),
                                       prediction_types = c("截距", "截距+斜率")) {
  cat("执行环境留一基因型5折组合交叉验证...\n")
  cat("重要：测试集 = 测试基因型 × 测试环境（cold × cold）\n")
  cat("训练集 = 训练基因型 × 训练环境\n")
  cat("预测方式:", paste(prediction_types, collapse = ", "), "\n")

  genotypes <- unique(pheno_data$genotype)
  genotypes <- intersect(genotypes, GD_all$taxa)
  cat("总品种数量:", length(genotypes), "\n")

  env_codes <- unique(pheno_data$env_code)
  cat("总环境数量:", length(env_codes), "\n")

  set.seed(195021)
  if (length(genotypes) < nfolds) {
    nfolds <- length(genotypes)
    cat("调整折数为品种数量:", nfolds, "\n")
  }

  geno_folds <- caret::createFolds(1:length(genotypes), k = nfolds, list = TRUE, returnTrain = FALSE)

  final_results_file <- file.path(gs_output_dir, "Geno5_EnvLOO", "环境留一基因型5折组合交叉验证_最终结果.csv")

  if (!dir.exists(dirname(final_results_file))) {
    dir.create(dirname(final_results_file), recursive = TRUE)
  }

  if (file.exists(final_results_file)) {
    cat("发现已存在的最终结果文件，尝试加载...\n")
    tryCatch({
      previous_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      cat("成功加载", nrow(previous_results), "条已有结果\n")
      results_all <- previous_results
      completed_combinations <- unique(paste(previous_results$基因型折数, previous_results$环境ID,
                                             previous_results$方法, previous_results$预测类型, sep = "_"))
      cat("已有", length(completed_combinations), "个已完成的计算组合\n")
    }, error = function(e) {
      cat("无法加载已有结果文件，将创建新文件:", e$message, "\n")
      results_all <- data.frame()
    })
  } else {
    cat("没有找到已有结果文件，将创建新文件\n")
    results_all <- data.frame()
  }

  total_iterations <- nfolds * length(env_codes)
  current_iteration <- 0

  all_combinations <- expand.grid(fold = 1:nfolds, env_index = 1:length(env_codes))
  total_combinations <- nrow(all_combinations)

  cat("总交叉验证组合数量:", total_combinations, "\n")
  cat("每个组合包含", length(methods), "种方法 ×", length(prediction_types), "种预测方式 =",
      length(methods) * length(prediction_types), "个配置\n")

  for (combo_index in seq_len(total_combinations)) {
    fold <- all_combinations$fold[combo_index]
    env_index <- all_combinations$env_index[combo_index]

    current_iteration <- current_iteration + 1
    cat("\n>>> 组合", combo_index, "/", total_combinations,
        " (总体进度: ", round(current_iteration / total_iterations * 100, 1), "%)\n", sep = "")
    cat("    基因型第", fold, "/", nfolds, "折，环境留一第", env_index, "/", length(env_codes), "\n", sep = "")

    test_env <- env_codes[env_index]
    combo_id <- paste(fold, test_env, sep = "_")

    fold_dir <- file.path(gs_output_dir, "Geno5_EnvLOO", paste0("Fold_", fold))
    if (!dir.exists(fold_dir)) {
      dir.create(fold_dir, recursive = TRUE)
    }

    combo_marker_file <- file.path(fold_dir, paste0("组合_", combo_id, "_完成标记.txt"))

    if (file.exists(combo_marker_file)) {
      cat("    发现完成标记文件，跳过这个组合:", combo_id, "\n")
      next
    }

    test_genotypes_fold <- genotypes[geno_folds[[fold]]]
    train_genotypes_fold <- setdiff(genotypes, test_genotypes_fold)

    cat("    训练基因型数量:", length(train_genotypes_fold), "\n")
    cat("    测试基因型数量:", length(test_genotypes_fold), "\n")

    test_env <- env_codes[env_index]
    train_envs <- setdiff(env_codes, test_env)

    cat("    测试环境:", test_env, "\n")
    cat("    训练环境:", paste(train_envs, collapse = ", "), "\n")

    train_pheno <- pheno_data[pheno_data$genotype %in% train_genotypes_fold &
                                pheno_data$env_code %in% train_envs, , drop = FALSE]

    test_pheno_env <- pheno_data[pheno_data$genotype %in% test_genotypes_fold &
                                   pheno_data$env_code == test_env, , drop = FALSE]

    cat("    训练集记录数（训练基因型×训练环境）:", nrow(train_pheno), "\n")
    cat("    测试集记录数（测试基因型×测试环境）:", nrow(test_pheno_env), "\n")

    if (nrow(train_pheno) < 10 || nrow(test_pheno_env) < 1) {
      cat("    数据不足，跳过这个组合\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 数据不足跳过"), combo_marker_file)
      next
    }

    cat("    步骤1: 在训练集上计算FW参数...\n")
    fw_result <- calculate_FW_parameters(train_pheno, env_factors, env_means)
    train_FW <- fw_result$FW_params
    mean_kPara_train <- fw_result$mean_kPara_train

    if (nrow(train_FW) < 5) {
      cat("    FW参数不足，跳过这个组合\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: FW参数不足跳过"), combo_marker_file)
      next
    }

    test_env_factor <- env_factors$kPara[env_factors$env_code == test_env]
    if (length(test_env_factor) == 0) {
      cat("    找不到测试环境参数，跳过这个组合\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 找不到环境参数跳过"), combo_marker_file)
      next
    }

    test_env_factor_centered <- test_env_factor[1] - mean_kPara_train
    cat("    测试环境参数:", test_env_factor[1], "\n")
    cat("    中心化的测试环境参数:", test_env_factor_centered, "\n")

    cat("    步骤2: 在训练集上进行GWAS分析...\n")
    traits <- c("Intcp_para_adj", "Slope_para")
    significant_markers_list <- list()

    gwas_base_dir <- file.path(fold_dir, paste0("GWAS_Fold", fold, "_Env", test_env))
    if (!dir.exists(gwas_base_dir)) {
      dir.create(gwas_base_dir, recursive = TRUE)
    }

    for (trait in traits) {
      cat("      在训练集上进行", trait, "的GWAS分析...\n")
      gwas_dir <- file.path(gwas_base_dir, paste0("GWAS_", trait))
      gwas_res <- perform_gapit_gwas(train_FW, GD_all, GM_all, trait, gwas_dir)

      if (!is.null(gwas_res)) {
        cat("      从", trait, "的GWAS结果中选择显著标记（双阈值+LD过滤）...\n")
        significant_markers <- select_markers_by_pvalue_ld(
          gwas_dir = gwas_dir,
          geno_all = GD_all[, -1, drop = FALSE],
          GM_all = GM_all,
          pval_thresh = 0.0001,
          hb_pval_thresh = 0.05,
          ld_threshold = 0.7
        )
        significant_markers_list[[trait]] <- significant_markers
        cat("      ", trait, "选择的显著标记数量:", length(significant_markers), "\n")

        write.table(
          data.frame(SNP = significant_markers, stringsAsFactors = FALSE),
          file.path(gwas_base_dir, paste0("Significant_Markers_", trait, ".txt")),
          row.names = FALSE, col.names = FALSE, quote = FALSE
        )
      } else {
        cat("      ", trait, "GWAS失败\n")
        significant_markers_list[[trait]] <- character(0)
      }
    }

    cat("    步骤4: 进行GS预测（", length(methods), "种方法 × ", length(prediction_types), "种预测方式）...\n", sep = "")
    combo_results <- data.frame()

    for (method_name in methods) {
      cat("      处理方法:", method_name, "...\n")

      method_parts <- strsplit(method_name, "_")[[1]]
      framework <- method_parts[1]
      method_type <- paste(method_parts[-1], collapse = "_")

      for (pred_type in prediction_types) {
        cat("        预测方式:", pred_type, "...\n")
        trait_predictions <- list()

        if (pred_type == "截距") {
          trait <- "Intcp_para_adj"
          y_train <- train_FW[[trait]]
          names(y_train) <- train_FW$genotype
          valid_train <- !is.na(y_train)

          if (sum(valid_train) < 5) {
            cat("        ", trait, "有效训练数据不足，跳过\n")
            next
          }

          y_train_valid <- y_train[valid_train]
          train_genotypes_valid <- names(y_train_valid)
          test_genotypes_valid <- unique(test_pheno_env$genotype)

          if (framework == "GAPIT") {
            if (method_type == "MAS") {
              pred_test <- GAPIT_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- GAPIT_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          } else if (framework == "BGLR") {
            if (method_type == "MAS") {
              pred_test <- BGLR_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- BGLR_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          }

          trait_predictions[[trait]] <- pred_test
          cat("        ", trait, method_name, pred_type, "预测完成\n")

          trait_predictions[["Slope_para"]] <- rep(NA, length(pred_test))
          if (!is.null(names(pred_test))) {
            names(trait_predictions[["Slope_para"]]) <- names(pred_test)
          }

        } else if (pred_type == "截距+斜率") {
          for (trait in traits) {
            y_train <- train_FW[[trait]]
            names(y_train) <- train_FW$genotype
            valid_train <- !is.na(y_train)

            if (sum(valid_train) < 5) {
              cat("        ", trait, "有效训练数据不足，跳过\n")
              next
            }

            y_train_valid <- y_train[valid_train]
            train_genotypes_valid <- names(y_train_valid)
            test_genotypes_valid <- unique(test_pheno_env$genotype)

            if (framework == "GAPIT") {
              if (method_type == "MAS") {
                pred_test <- GAPIT_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- GAPIT_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            } else if (framework == "BGLR") {
              if (method_type == "MAS") {
                pred_test <- BGLR_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- BGLR_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            }

            trait_predictions[[trait]] <- pred_test
            cat("        ", trait, method_name, pred_type, "预测完成\n")
          }
        }

        if (length(trait_predictions) >= 1 &&
            !is.null(trait_predictions[["Intcp_para_adj"]])) {

          intcp_pred <- trait_predictions[["Intcp_para_adj"]]
          slope_pred <- trait_predictions[["Slope_para"]]

          for (geno in unique(test_pheno_env$genotype)) {
            if (geno %in% names(intcp_pred)) {
              intcp_value <- intcp_pred[geno]

              if (pred_type == "截距") {
                pheno_pred <- intcp_value
                slope_value <- NA
              } else if (pred_type == "截距+斜率") {
                if (!is.null(slope_pred) && geno %in% names(slope_pred) && !is.na(slope_pred[geno])) {
                  slope_value <- slope_pred[geno]
                  pheno_pred <- intcp_value + slope_value * test_env_factor_centered
                } else {
                  pheno_pred <- intcp_value
                  slope_value <- NA
                }
              }

              actual_rows <- test_pheno_env[test_pheno_env$genotype == geno, , drop = FALSE]
              if (nrow(actual_rows) > 0) {
                actual_value <- actual_rows$TKW[1]
                if (!is.na(actual_value)) {
                  result_row <- data.frame(
                    基因型ID = geno,
                    环境ID = test_env,
                    基因型折数 = fold,
                    环境留一折数 = env_index,
                    实际值 = actual_value,
                    预测值 = pheno_pred,
                    FW_矫正截距 = intcp_value,
                    FW_斜率 = ifelse(is.na(slope_value), NA, slope_value),
                    方法 = method_name,
                    框架 = framework,
                    方法类型 = method_type,
                    预测类型 = pred_type,
                    使用PCA = "否",
                    交叉验证类型 = "环境留一基因型5折组合交叉验证",
                    环境参数 = test_env_factor[1],
                    中心化环境参数 = test_env_factor_centered,
                    备注 = "测试集 = 测试基因型 × 测试环境（cold × cold）",
                    计算时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    stringsAsFactors = FALSE
                  )

                  combo_results <- rbind(combo_results, result_row)
                }
              }
            }
          }

          cat("        ", method_name, pred_type, "预测记录生成完成\n")
        }
      }
    }

    cat("    步骤5: 立即保存当前组合的结果...\n")

    if (nrow(combo_results) > 0) {
      combo_results_file <- file.path(fold_dir, paste0("组合_", combo_id, "_结果.csv"))
      write.csv(combo_results, combo_results_file, row.names = FALSE)
      cat("      保存到组合文件:", combo_results_file, " (", nrow(combo_results), "条记录)\n", sep = "")

      if (file.exists(final_results_file)) {
        write.table(combo_results, final_results_file, sep = ",",
                    col.names = FALSE, row.names = FALSE, append = TRUE)
      } else {
        write.csv(combo_results, final_results_file, row.names = FALSE)
      }
      cat("      追加到总结果文件:", final_results_file, "\n")

      results_all <- rbind(results_all, combo_results)

      writeLines(
        paste(
          "处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "\n状态: 成功完成",
          "\n记录数:", nrow(combo_results),
          "\n方法数:", length(unique(combo_results$方法)),
          "\n预测类型:", paste(unique(combo_results$预测类型), collapse = ", "),
          "\n测试环境:", test_env,
          "\n测试基因型数:", length(unique(combo_results$基因型ID)),
          "\nCV类型: cold × cold（测试基因型×测试环境）"
        ),
        combo_marker_file
      )

      cat("      创建完成标记文件:", combo_marker_file, "\n")
    } else {
      cat("      当前组合没有生成有效结果\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 无有效结果"), combo_marker_file)
    }

    rm(combo_results, trait_predictions, significant_markers_list, train_FW)
    gc()

    cat("    当前组合处理完成!\n")
  }

  cat("\n环境留一基因型5折组合交叉验证完成!\n")

  if (nrow(results_all) > 0) {
    if (file.exists(final_results_file)) {
      cat("重新读取并整理最终结果文件...\n")
      final_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      final_results <- final_results[!duplicated(final_results), , drop = FALSE]
      write.csv(final_results, final_results_file, row.names = FALSE)

      cat("最终结果文件已整理:", final_results_file, "\n")
      cat("总记录数:", nrow(final_results), "\n")

      if (nrow(final_results) > 0) {
        cat("\n基本统计信息:\n")
        cat("总组合数:", length(unique(paste(final_results$基因型折数, final_results$环境ID))), "\n")
        cat("总方法数:", length(unique(final_results$方法)), "\n")
        cat("总预测类型数:", length(unique(final_results$预测类型)), "\n")
        cat("总品种数:", length(unique(final_results$基因型ID)), "\n")
        cat("总环境数:", length(unique(final_results$环境ID)), "\n")
      }
    }
  } else {
    cat("没有生成有效结果\n")
  }

  return(results_all)
}

# 8.2 基因型5折交叉验证通用函数
run_geno5_fold_cv_simple <- function(pheno_data, GD_all, GM_all, env_factors, env_means, nfolds = 5,
                                     methods = c("GAPIT_MAS", "GAPIT_GBLUP", "GAPIT_MAS+GBLUP",
                                                 "BGLR_MAS", "BGLR_GBLUP", "BGLR_MAS+GBLUP"),
                                     prediction_types = c("截距", "截距+斜率")) {
  cat("执行基因型5折交叉验证...\n")
  cat("预测方式:", paste(prediction_types, collapse = ", "), "\n")
  cat("注意：基因型5折的测试集 = 测试基因型 × 所有环境\n")

  genotypes <- unique(pheno_data$genotype)
  genotypes <- intersect(genotypes, GD_all$taxa)
  cat("总品种数量:", length(genotypes), "\n")

  set.seed(195021)
  if (length(genotypes) < nfolds) {
    nfolds <- length(genotypes)
    cat("调整折数为品种数量:", nfolds, "\n")
  }

  geno_folds <- caret::createFolds(1:length(genotypes), k = nfolds, list = TRUE, returnTrain = FALSE)

  final_results_file <- file.path(gs_output_dir, "Geno5_fold", "基因型5折交叉验证_最终结果.csv")

  if (!dir.exists(dirname(final_results_file))) {
    dir.create(dirname(final_results_file), recursive = TRUE)
  }

  if (file.exists(final_results_file)) {
    cat("发现已存在的最终结果文件，尝试加载...\n")
    tryCatch({
      previous_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      cat("成功加载", nrow(previous_results), "条已有结果\n")
      results_all <- previous_results
      completed_folds <- unique(previous_results$折数)
      cat("已有折:", paste(completed_folds, collapse = ", "), "\n")
    }, error = function(e) {
      cat("无法加载已有结果文件，将创建新文件:", e$message, "\n")
      results_all <- data.frame()
    })
  } else {
    cat("没有找到已有结果文件，将创建新文件\n")
    results_all <- data.frame()
  }

  for (fold in seq_len(nfolds)) {
    cat("\n>>> 第", fold, "/", nfolds, "折：基因型5折\n")

    fold_marker_file <- file.path(gs_output_dir, "Geno5_fold", paste0("Fold_", fold), "完成标记.txt")

    if (file.exists(fold_marker_file)) {
      cat("    发现完成标记文件，跳过第", fold, "折\n", sep = "")
      next
    }

    fold_dir <- file.path(gs_output_dir, "Geno5_fold", paste0("Fold_", fold))
    if (!dir.exists(fold_dir)) {
      dir.create(fold_dir, recursive = TRUE)
    }

    test_genotypes <- genotypes[geno_folds[[fold]]]
    train_genotypes <- setdiff(genotypes, test_genotypes)

    cat("    训练基因型数量:", length(train_genotypes), "\n")
    cat("    测试基因型数量:", length(test_genotypes), "\n")

    train_pheno <- pheno_data[pheno_data$genotype %in% train_genotypes, , drop = FALSE]
    test_pheno <- pheno_data[pheno_data$genotype %in% test_genotypes, , drop = FALSE]

    cat("    训练数据记录数（训练基因型×所有环境）:", nrow(train_pheno), "\n")
    cat("    测试数据记录数（测试基因型×所有环境）:", nrow(test_pheno), "\n")

    if (nrow(train_pheno) < 10 || nrow(test_pheno) < 1) {
      cat("    数据不足，跳过\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 数据不足跳过"), fold_marker_file)
      next
    }

    cat("    步骤1: 在训练集上计算FW参数...\n")
    fw_result <- calculate_FW_parameters(train_pheno, env_factors, env_means)
    train_FW <- fw_result$FW_params
    mean_kPara_train <- fw_result$mean_kPara_train

    if (nrow(train_FW) < 5) {
      cat("    FW参数不足，跳过\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: FW参数不足跳过"), fold_marker_file)
      next
    }

    write.csv(train_FW, file.path(fold_dir, "训练集FW参数.csv"), row.names = FALSE)
    cat("    训练集FW参数计算完成，品种数量:", nrow(train_FW), "\n")

    cat("    步骤2: 在训练集上进行GWAS分析...\n")
    traits <- c("Intcp_para_adj", "Slope_para")
    significant_markers_list <- list()

    for (trait in traits) {
      cat("      在训练集上进行", trait, "的GWAS分析...\n")
      gwas_dir <- file.path(fold_dir, paste0("GWAS_", trait))
      gwas_res <- perform_gapit_gwas(train_FW, GD_all, GM_all, trait, gwas_dir)

      if (!is.null(gwas_res)) {
        cat("      从", trait, "的GWAS结果中选择显著标记（双阈值+LD过滤）...\n")
        significant_markers <- select_markers_by_pvalue_ld(
          gwas_dir = gwas_dir,
          geno_all = GD_all[, -1, drop = FALSE],
          GM_all = GM_all,
          pval_thresh = 0.0001,
          hb_pval_thresh = 0.05,
          ld_threshold = 0.7
        )
        significant_markers_list[[trait]] <- significant_markers
        cat("      ", trait, "选择的显著标记数量:", length(significant_markers), "\n")

        write.table(
          data.frame(SNP = significant_markers, stringsAsFactors = FALSE),
          file.path(fold_dir, paste0("Significant_Markers_", trait, ".txt")),
          row.names = FALSE, col.names = FALSE, quote = FALSE
        )
      } else {
        cat("      ", trait, "GWAS失败\n")
        significant_markers_list[[trait]] <- character(0)
      }
    }

    cat("    步骤4: 进行GS预测（", length(methods), "种方法 × ", length(prediction_types), "种预测方式）...\n", sep = "")

    fold_results_df <- data.frame()

    for (method_name in methods) {
      cat("      处理方法:", method_name, "...\n")

      method_parts <- strsplit(method_name, "_")[[1]]
      framework <- method_parts[1]
      method_type <- paste(method_parts[-1], collapse = "_")

      for (pred_type in prediction_types) {
        cat("        预测方式:", pred_type, "...\n")
        trait_predictions <- list()

        if (pred_type == "截距") {
          trait <- "Intcp_para_adj"
          y_train <- train_FW[[trait]]
          names(y_train) <- train_FW$genotype
          valid_train <- !is.na(y_train)

          if (sum(valid_train) < 5) {
            cat("        ", trait, "有效训练数据不足，跳过\n")
            next
          }

          y_train_valid <- y_train[valid_train]
          train_genotypes_valid <- names(y_train_valid)
          test_genotypes_valid <- unique(test_pheno$genotype)

          if (framework == "GAPIT") {
            if (method_type == "MAS") {
              pred_test <- GAPIT_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- GAPIT_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          } else if (framework == "BGLR") {
            if (method_type == "MAS") {
              pred_test <- BGLR_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- BGLR_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          }

          trait_predictions[[trait]] <- pred_test
          cat("        ", trait, method_name, pred_type, "预测完成\n")

          trait_predictions[["Slope_para"]] <- rep(NA, length(pred_test))
          if (!is.null(names(pred_test))) {
            names(trait_predictions[["Slope_para"]]) <- names(pred_test)
          }

        } else if (pred_type == "截距+斜率") {
          for (trait in traits) {
            y_train <- train_FW[[trait]]
            names(y_train) <- train_FW$genotype
            valid_train <- !is.na(y_train)

            if (sum(valid_train) < 5) {
              cat("        ", trait, "有效训练数据不足，跳过\n")
              next
            }

            y_train_valid <- y_train[valid_train]
            train_genotypes_valid <- names(y_train_valid)
            test_genotypes_valid <- unique(test_pheno$genotype)

            if (framework == "GAPIT") {
              if (method_type == "MAS") {
                pred_test <- GAPIT_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- GAPIT_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            } else if (framework == "BGLR") {
              if (method_type == "MAS") {
                pred_test <- BGLR_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- BGLR_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            }

            trait_predictions[[trait]] <- pred_test
            cat("        ", trait, method_name, pred_type, "预测完成\n")
          }
        }

        if (length(trait_predictions) >= 1 &&
            !is.null(trait_predictions[["Intcp_para_adj"]])) {

          intcp_pred <- trait_predictions[["Intcp_para_adj"]]
          slope_pred <- trait_predictions[["Slope_para"]]

          for (env in unique(test_pheno$env_code)) {
            env_factor <- env_factors$kPara[env_factors$env_code == env]
            if (length(env_factor) > 0) {
              env_factor_centered <- env_factor[1] - mean_kPara_train

              for (geno in unique(test_pheno$genotype)) {
                if (geno %in% names(intcp_pred)) {
                  intcp_value <- intcp_pred[geno]

                  if (pred_type == "截距") {
                    pheno_pred <- intcp_value
                    slope_value <- NA
                  } else if (pred_type == "截距+斜率") {
                    if (!is.null(slope_pred) && geno %in% names(slope_pred) && !is.na(slope_pred[geno])) {
                      slope_value <- slope_pred[geno]
                      pheno_pred <- intcp_value + slope_value * env_factor_centered
                    } else {
                      pheno_pred <- intcp_value
                      slope_value <- NA
                    }
                  }

                  actual_rows <- test_pheno[test_pheno$genotype == geno & test_pheno$env_code == env, , drop = FALSE]
                  if (nrow(actual_rows) > 0) {
                    actual_value <- actual_rows$TKW[1]
                    if (!is.na(actual_value)) {
                      result_row <- data.frame(
                        基因型ID = geno,
                        环境ID = env,
                        折数 = fold,
                        实际值 = actual_value,
                        预测值 = pheno_pred,
                        FW_矫正截距 = intcp_value,
                        FW_斜率 = ifelse(is.na(slope_value), NA, slope_value),
                        方法 = method_name,
                        框架 = framework,
                        方法类型 = method_type,
                        预测类型 = pred_type,
                        使用PCA = "否",
                        交叉验证类型 = "基因型5折交叉验证",
                        基因型折数 = fold,
                        环境参数 = env_factor[1],
                        中心化环境参数 = env_factor_centered,
                        备注 = "测试集 = 测试基因型 × 所有环境",
                        计算时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        stringsAsFactors = FALSE
                      )

                      fold_results_df <- rbind(fold_results_df, result_row)
                    }
                  }
                }
              }
            }
          }

          cat("        ", method_name, pred_type, "预测记录生成完成\n")
        }
      }
    }

    cat("    步骤5: 立即保存当前折的结果...\n")

    if (nrow(fold_results_df) > 0) {
      fold_results_file <- file.path(fold_dir, paste0("Fold", fold, "_结果.csv"))
      write.csv(fold_results_df, fold_results_file, row.names = FALSE)
      cat("      保存到折文件:", fold_results_file, " (", nrow(fold_results_df), "条记录)\n", sep = "")

      if (file.exists(final_results_file)) {
        write.table(fold_results_df, final_results_file, sep = ",",
                    col.names = FALSE, row.names = FALSE, append = TRUE)
      } else {
        write.csv(fold_results_df, final_results_file, row.names = FALSE)
      }
      cat("      追加到总结果文件:", final_results_file, "\n")

      results_all <- rbind(results_all, fold_results_df)

      writeLines(
        paste(
          "处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "\n状态: 成功完成",
          "\n记录数:", nrow(fold_results_df),
          "\n方法数:", length(unique(fold_results_df$方法)),
          "\n预测类型:", paste(unique(fold_results_df$预测类型), collapse = ", "),
          "\n测试环境数:", length(unique(fold_results_df$环境ID)),
          "\n测试基因型数:", length(unique(fold_results_df$基因型ID))
        ),
        fold_marker_file
      )

      cat("      创建完成标记文件:", fold_marker_file, "\n")
    } else {
      cat("      当前折没有生成有效结果\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 无有效结果"), fold_marker_file)
    }

    rm(fold_results_df, trait_predictions, significant_markers_list, train_FW)
    gc()

    progress <- round(fold / nfolds * 100, 1)
    cat("    总体进度:", progress, "%\n")
  }

  cat("\n基因型5折交叉验证完成!\n")

  if (nrow(results_all) > 0) {
    if (file.exists(final_results_file)) {
      cat("重新读取并整理最终结果文件...\n")
      final_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      final_results <- final_results[!duplicated(final_results), , drop = FALSE]
      write.csv(final_results, final_results_file, row.names = FALSE)
      cat("最终结果文件已整理:", final_results_file, "\n")
      cat("总记录数:", nrow(final_results), "\n")
    }
  } else {
    cat("没有生成有效结果\n")
  }

  return(results_all)
}

# 8.3 环境留一交叉验证通用函数
run_env_loo_cv_simple <- function(pheno_data, GD_all, GM_all, env_factors, env_means,
                                  methods = c("GAPIT_MAS", "GAPIT_GBLUP", "GAPIT_MAS+GBLUP",
                                              "BGLR_MAS", "BGLR_GBLUP", "BGLR_MAS+GBLUP"),
                                  prediction_types = c("截距", "截距+斜率")) {
  cat("执行环境留一交叉验证...\n")
  cat("预测方式:", paste(prediction_types, collapse = ", "), "\n")
  cat("注意：环境留一的测试集 = 所有基因型 × 测试环境\n")

  env_codes <- unique(pheno_data$env_code)
  cat("总环境数量:", length(env_codes), "\n")

  final_results_file <- file.path(gs_output_dir, "EnvLOO", "环境留一交叉验证_最终结果.csv")

  if (!dir.exists(dirname(final_results_file))) {
    dir.create(dirname(final_results_file), recursive = TRUE)
  }

  if (file.exists(final_results_file)) {
    cat("发现已存在的最终结果文件，尝试加载...\n")
    tryCatch({
      previous_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      cat("成功加载", nrow(previous_results), "条已有结果\n")
      results_all <- previous_results
      completed_envs <- unique(previous_results$环境ID)
      cat("已有环境:", paste(completed_envs, collapse = ", "), "\n")
    }, error = function(e) {
      cat("无法加载已有结果文件，将创建新文件:", e$message, "\n")
      results_all <- data.frame()
    })
  } else {
    cat("没有找到已有结果文件，将创建新文件\n")
    results_all <- data.frame()
  }

  for (env_index in seq_along(env_codes)) {
    cat("\n>>> 第", env_index, "/", length(env_codes), "折：环境留一\n")

    test_env <- env_codes[env_index]
    env_marker_file <- file.path(gs_output_dir, "EnvLOO", paste0("Env_", test_env), "完成标记.txt")

    if (file.exists(env_marker_file)) {
      cat("    发现完成标记文件，跳过环境:", test_env, "\n")
      next
    }

    train_envs <- setdiff(env_codes, test_env)

    cat("    测试环境:", test_env, "\n")
    cat("    训练环境:", paste(train_envs, collapse = ", "), "\n")

    env_dir <- file.path(gs_output_dir, "EnvLOO", paste0("Env_", test_env))
    if (!dir.exists(env_dir)) {
      dir.create(env_dir, recursive = TRUE)
    }

    train_pheno <- pheno_data[pheno_data$env_code %in% train_envs, , drop = FALSE]
    test_pheno <- pheno_data[pheno_data$env_code == test_env, , drop = FALSE]

    cat("    训练数据记录数（所有基因型×训练环境）:", nrow(train_pheno), "\n")
    cat("    测试数据记录数（所有基因型×测试环境）:", nrow(test_pheno), "\n")

    if (nrow(train_pheno) < 10 || nrow(test_pheno) < 1) {
      cat("    数据不足，跳过\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 数据不足跳过"), env_marker_file)
      next
    }

    cat("    步骤1: 在训练集上计算FW参数...\n")
    fw_result <- calculate_FW_parameters(train_pheno, env_factors, env_means)
    train_FW <- fw_result$FW_params
    mean_kPara_train <- fw_result$mean_kPara_train

    if (nrow(train_FW) < 5) {
      cat("    FW参数不足，跳过\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: FW参数不足跳过"), env_marker_file)
      next
    }

    write.csv(train_FW, file.path(env_dir, "训练集FW参数.csv"), row.names = FALSE)
    cat("    训练集FW参数计算完成，品种数量:", nrow(train_FW), "\n")

    test_env_factor <- env_factors$kPara[env_factors$env_code == test_env]
    if (length(test_env_factor) == 0) {
      cat("    找不到测试环境参数，跳过\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 找不到环境参数跳过"), env_marker_file)
      next
    }

    test_env_factor_centered <- test_env_factor[1] - mean_kPara_train
    cat("    测试环境参数:", test_env_factor[1], "\n")
    cat("    中心化的测试环境参数:", test_env_factor_centered, "\n")

    cat("    步骤2: 在训练集上进行GWAS分析...\n")
    traits <- c("Intcp_para_adj", "Slope_para")
    significant_markers_list <- list()

    for (trait in traits) {
      cat("      在训练集上进行", trait, "的GWAS分析...\n")
      gwas_dir <- file.path(env_dir, paste0("GWAS_", trait))
      gwas_res <- perform_gapit_gwas(train_FW, GD_all, GM_all, trait, gwas_dir)

      if (!is.null(gwas_res)) {
        cat("      从", trait, "的GWAS结果中选择显著标记（双阈值+LD过滤）...\n")
        significant_markers <- select_markers_by_pvalue_ld(
          gwas_dir = gwas_dir,
          geno_all = GD_all[, -1, drop = FALSE],
          GM_all = GM_all,
          pval_thresh = 0.0001,
          hb_pval_thresh = 0.05,
          ld_threshold = 0.7
        )
        significant_markers_list[[trait]] <- significant_markers
        cat("      ", trait, "选择的显著标记数量:", length(significant_markers), "\n")

        write.table(
          data.frame(SNP = significant_markers, stringsAsFactors = FALSE),
          file.path(env_dir, paste0("Significant_Markers_", trait, ".txt")),
          row.names = FALSE, col.names = FALSE, quote = FALSE
        )
      } else {
        cat("      ", trait, "GWAS失败\n")
        significant_markers_list[[trait]] <- character(0)
      }
    }

    cat("    步骤4: 进行GS预测（", length(methods), "种方法 × ", length(prediction_types), "种预测方式）...\n", sep = "")

    env_results_df <- data.frame()

    for (method_name in methods) {
      cat("      处理方法:", method_name, "...\n")

      method_parts <- strsplit(method_name, "_")[[1]]
      framework <- method_parts[1]
      method_type <- paste(method_parts[-1], collapse = "_")

      for (pred_type in prediction_types) {
        cat("        预测方式:", pred_type, "...\n")
        trait_predictions <- list()

        if (pred_type == "截距") {
          trait <- "Intcp_para_adj"
          y_train <- train_FW[[trait]]
          names(y_train) <- train_FW$genotype
          valid_train <- !is.na(y_train)

          if (sum(valid_train) < 5) {
            cat("        ", trait, "有效训练数据不足，跳过\n")
            next
          }

          y_train_valid <- y_train[valid_train]
          train_genotypes_valid <- names(y_train_valid)
          test_genotypes_valid <- unique(test_pheno$genotype)

          if (framework == "GAPIT") {
            if (method_type == "MAS") {
              pred_test <- GAPIT_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- GAPIT_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          } else if (framework == "BGLR") {
            if (method_type == "MAS") {
              pred_test <- BGLR_MAS_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all,
                significant_markers_list[[trait]]
              )
            } else if (method_type == "GBLUP") {
              pred_test <- BGLR_GBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all
              )
            } else if (method_type == "MAS+GBLUP") {
              pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                y_train_valid, train_genotypes_valid, test_genotypes_valid,
                GD_all, GM_all,
                significant_markers_list[[trait]]
              )
            }
          }

          trait_predictions[[trait]] <- pred_test
          cat("        ", trait, method_name, pred_type, "预测完成\n")

          trait_predictions[["Slope_para"]] <- rep(NA, length(pred_test))
          if (!is.null(names(pred_test))) {
            names(trait_predictions[["Slope_para"]]) <- names(pred_test)
          }

        } else if (pred_type == "截距+斜率") {
          for (trait in traits) {
            y_train <- train_FW[[trait]]
            names(y_train) <- train_FW$genotype
            valid_train <- !is.na(y_train)

            if (sum(valid_train) < 5) {
              cat("        ", trait, "有效训练数据不足，跳过\n")
              next
            }

            y_train_valid <- y_train[valid_train]
            train_genotypes_valid <- names(y_train_valid)
            test_genotypes_valid <- unique(test_pheno$genotype)

            if (framework == "GAPIT") {
              if (method_type == "MAS") {
                pred_test <- GAPIT_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- GAPIT_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- GAPIT_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            } else if (framework == "BGLR") {
              if (method_type == "MAS") {
                pred_test <- BGLR_MAS_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all,
                  significant_markers_list[[trait]]
                )
              } else if (method_type == "GBLUP") {
                pred_test <- BGLR_GBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all
                )
              } else if (method_type == "MAS+GBLUP") {
                pred_test <- BGLR_MAS_gBLUP_prediction_simple(
                  y_train_valid, train_genotypes_valid, test_genotypes_valid,
                  GD_all, GM_all,
                  significant_markers_list[[trait]]
                )
              }
            }

            trait_predictions[[trait]] <- pred_test
            cat("        ", trait, method_name, pred_type, "预测完成\n")
          }
        }

        if (length(trait_predictions) >= 1 &&
            !is.null(trait_predictions[["Intcp_para_adj"]])) {

          intcp_pred <- trait_predictions[["Intcp_para_adj"]]
          slope_pred <- trait_predictions[["Slope_para"]]

          for (geno in unique(test_pheno$genotype)) {
            if (geno %in% names(intcp_pred)) {
              intcp_value <- intcp_pred[geno]

              if (pred_type == "截距") {
                pheno_pred <- intcp_value
                slope_value <- NA
              } else if (pred_type == "截距+斜率") {
                if (!is.null(slope_pred) && geno %in% names(slope_pred) && !is.na(slope_pred[geno])) {
                  slope_value <- slope_pred[geno]
                  pheno_pred <- intcp_value + slope_value * test_env_factor_centered
                } else {
                  pheno_pred <- intcp_value
                  slope_value <- NA
                }
              }

              actual_rows <- test_pheno[test_pheno$genotype == geno, , drop = FALSE]
              if (nrow(actual_rows) > 0) {
                actual_value <- actual_rows$TKW[1]
                if (!is.na(actual_value)) {
                  result_row <- data.frame(
                    基因型ID = geno,
                    环境ID = test_env,
                    折数 = env_index,
                    实际值 = actual_value,
                    预测值 = pheno_pred,
                    FW_矫正截距 = intcp_value,
                    FW_斜率 = ifelse(is.na(slope_value), NA, slope_value),
                    方法 = method_name,
                    框架 = framework,
                    方法类型 = method_type,
                    预测类型 = pred_type,
                    使用PCA = "否",
                    交叉验证类型 = "环境留一交叉验证",
                    环境参数 = test_env_factor[1],
                    中心化环境参数 = test_env_factor_centered,
                    备注 = "测试集 = 所有基因型 × 测试环境",
                    计算时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    stringsAsFactors = FALSE
                  )

                  env_results_df <- rbind(env_results_df, result_row)
                }
              }
            }
          }

          cat("        ", method_name, pred_type, "预测记录生成完成\n")
        }
      }
    }

    cat("    步骤5: 立即保存当前环境的结果...\n")

    if (nrow(env_results_df) > 0) {
      env_results_file <- file.path(env_dir, paste0("Env_", test_env, "_结果.csv"))
      write.csv(env_results_df, env_results_file, row.names = FALSE)
      cat("      保存到环境文件:", env_results_file, " (", nrow(env_results_df), "条记录)\n", sep = "")

      if (file.exists(final_results_file)) {
        write.table(env_results_df, final_results_file, sep = ",",
                    col.names = FALSE, row.names = FALSE, append = TRUE)
      } else {
        write.csv(env_results_df, final_results_file, row.names = FALSE)
      }
      cat("      追加到总结果文件:", final_results_file, "\n")

      results_all <- rbind(results_all, env_results_df)

      writeLines(
        paste(
          "处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "\n状态: 成功完成",
          "\n记录数:", nrow(env_results_df),
          "\n方法数:", length(unique(env_results_df$方法)),
          "\n预测类型:", paste(unique(env_results_df$预测类型), collapse = ", "),
          "\n测试基因型数:", length(unique(env_results_df$基因型ID))
        ),
        env_marker_file
      )

      cat("      创建完成标记文件:", env_marker_file, "\n")
    } else {
      cat("      当前环境没有生成有效结果\n")
      writeLines(paste("处理时间:", Sys.time(), "状态: 无有效结果"), env_marker_file)
    }

    rm(env_results_df, trait_predictions, significant_markers_list, train_FW)
    gc()

    progress <- round(env_index / length(env_codes) * 100, 1)
    cat("    总体进度:", progress, "%\n")
  }

  cat("\n环境留一交叉验证完成!\n")

  if (nrow(results_all) > 0) {
    if (file.exists(final_results_file)) {
      cat("重新读取并整理最终结果文件...\n")
      final_results <- read.csv(final_results_file, stringsAsFactors = FALSE)
      final_results <- final_results[!duplicated(final_results), , drop = FALSE]
      write.csv(final_results, final_results_file, row.names = FALSE)
      cat("最终结果文件已整理:", final_results_file, "\n")
      cat("总记录数:", nrow(final_results), "\n")
    }
  } else {
    cat("没有生成有效结果\n")
  }

  return(results_all)
}

# ===========================================
# 9. 执行所有交叉验证
# ===========================================

cat("\n=== 执行所有交叉验证 ===\n")

methods_no_pca <- c(
  "GAPIT_MAS", "GAPIT_GBLUP", "GAPIT_MAS+GBLUP",
  "BGLR_MAS", "BGLR_GBLUP", "BGLR_MAS+GBLUP"
)

prediction_types <- c("截距", "截距+斜率")

cat("将执行以下配置:\n")
cat("方法 (6种):", paste(methods_no_pca, collapse = ", "), "\n")
cat("预测方式 (2种):", paste(prediction_types, collapse = ", "), "\n")
cat("总计: 6种方法 × 2种预测方式 = 12种配置\n")

# 9.1 第一步：环境留一基因型5折组合交叉验证
cat("\n=== 第一步：环境留一基因型5折组合交叉验证 ===\n")
geno5_envloo_results <- run_geno5_envloo_cv_simple(
  phenotype_data_final, GD_all, GM_all, envMeanPara, env_mean_trait, nfolds = 5,
  methods = methods_no_pca,
  prediction_types = prediction_types
)

# 9.2 第二步：基因型5折交叉验证
cat("\n=== 第二步：基因型5折交叉验证 ===\n")
geno5_fold_results <- run_geno5_fold_cv_simple(
  phenotype_data_final, GD_all, GM_all, envMeanPara, env_mean_trait, nfolds = 5,
  methods = methods_no_pca,
  prediction_types = prediction_types
)

# 9.3 第三步：环境留一交叉验证
cat("\n=== 第三步：环境留一交叉验证 ===\n")
env_loo_results <- run_env_loo_cv_simple(
  phenotype_data_final, GD_all, GM_all, envMeanPara, env_mean_trait,
  methods = methods_no_pca,
  prediction_types = prediction_types
)

# ===========================================
# 10. 汇总统计和结果分析
# ===========================================

cat("\n=== 汇总统计和结果分析 ===\n")

generate_summary_statistics <- function(results_list, output_dir) {
  cat("生成汇总统计...\n")

  summary_data <- data.frame()

  for (cv_name in names(results_list)) {
    cv_data <- results_list[[cv_name]]

    if (!is.null(cv_data) && nrow(cv_data) > 0) {
      for (method in unique(cv_data$方法)) {
        for (pred_type in unique(cv_data$预测类型)) {
          method_data <- cv_data[cv_data$方法 == method & cv_data$预测类型 == pred_type, , drop = FALSE]
          if (nrow(method_data) > 0) {
            actual <- method_data$实际值
            predicted <- method_data$预测值
            valid_indices <- complete.cases(actual, predicted)

            if (sum(valid_indices) > 5) {
              actual_valid <- actual[valid_indices]
              predicted_valid <- predicted[valid_indices]

              cor_val <- suppressWarnings(cor(actual_valid, predicted_valid, use = "complete.obs"))
              rmse_val <- sqrt(mean((actual_valid - predicted_valid)^2, na.rm = TRUE))
              bias_val <- mean(predicted_valid - actual_valid, na.rm = TRUE)
              mae_val <- mean(abs(predicted_valid - actual_valid), na.rm = TRUE)

              method_parts <- strsplit(method, "_")[[1]]
              framework <- method_parts[1]
              method_type <- paste(method_parts[-1], collapse = "_")

              summary_data <- rbind(
                summary_data,
                data.frame(
                  交叉验证类型 = cv_name,
                  方法 = method,
                  框架 = framework,
                  方法类型 = method_type,
                  预测类型 = pred_type,
                  使用PCA = "否",
                  记录数 = nrow(method_data),
                  有效记录数 = sum(valid_indices),
                  平均相关性 = round(cor_val, 4),
                  平均RMSE = round(rmse_val, 4),
                  平均MAE = round(mae_val, 4),
                  平均偏差 = round(bias_val, 4),
                  备注 = "TKW性状分析",
                  stringsAsFactors = FALSE
                )
              )
            }
          }
        }
      }
    }
  }

  if (nrow(summary_data) > 0) {
    summary_file <- file.path(output_dir, "TKW交叉验证汇总统计.csv")
    write.csv(summary_data, summary_file, row.names = FALSE)
    cat("保存汇总统计:", nrow(summary_data), "条记录\n")
    cat("保存到:", summary_file, "\n")

    cat("\n汇总统计结果:\n")
    cat(strrep("=", 80), "\n")

    for (framework in c("GAPIT", "BGLR")) {
      cat("\n", framework, "框架结果:\n")
      cat(strrep("-", 60), "\n")
      framework_stats <- summary_data[summary_data$框架 == framework, , drop = FALSE]

      if (nrow(framework_stats) > 0) {
        for (method_type in c("MAS", "GBLUP", "MAS+GBLUP")) {
          # 关键修正：精确匹配，不再使用 grepl
          type_stats <- framework_stats[framework_stats$方法类型 == method_type, , drop = FALSE]

          if (nrow(type_stats) > 0) {
            cat("  ", method_type, "方法:\n")

            for (pred_type in c("截距", "截距+斜率")) {
              pred_stats <- type_stats[type_stats$预测类型 == pred_type, , drop = FALSE]

              if (nrow(pred_stats) > 0) {
                avg_cor <- mean(pred_stats$平均相关性, na.rm = TRUE)
                avg_rmse <- mean(pred_stats$平均RMSE, na.rm = TRUE)

                cat(sprintf("    %-10s: 平均相关性 = %.4f, 平均RMSE = %.4f\n",
                            pred_type, avg_cor, avg_rmse))

                for (cv_type in unique(pred_stats$交叉验证类型)) {
                  cv_row <- pred_stats[pred_stats$交叉验证类型 == cv_type, , drop = FALSE]
                  if (nrow(cv_row) > 0) {
                    cat(sprintf("      %-30s: 相关性 = %.4f, RMSE = %.4f\n",
                                cv_type, cv_row$平均相关性[1], cv_row$平均RMSE[1]))
                  }
                }
              }
            }
          }
        }
      }
    }

    cat("\n最佳预测结果（按方法类型和预测类型）:\n")
    cat(strrep("-", 60), "\n")

    for (framework in c("GAPIT", "BGLR")) {
      for (method_type in c("MAS", "GBLUP", "MAS+GBLUP")) {
        for (pred_type in c("截距", "截距+斜率")) {
          # 关键修正：精确匹配，不再使用 grepl
          type_stats <- summary_data[
            summary_data$框架 == framework &
              summary_data$方法类型 == method_type &
              summary_data$预测类型 == pred_type, , drop = FALSE
          ]

          if (nrow(type_stats) > 0) {
            best_idx <- which.max(type_stats$平均相关性)
            best_row <- type_stats[best_idx, , drop = FALSE]
            cat(sprintf("%-8s %-12s %-10s: 最佳相关性 = %.4f (%s, 记录数=%d)\n",
                        framework, method_type, pred_type,
                        best_row$平均相关性, best_row$交叉验证类型, best_row$记录数))
          }
        }
      }
    }

    if (nrow(summary_data) > 0) {
      overall_best_idx <- which.max(summary_data$平均相关性)
      overall_best <- summary_data[overall_best_idx, , drop = FALSE]
      cat("\n总体最佳预测:\n")
      cat(strrep("-", 40), "\n")
      cat(sprintf("  方法: %s\n", overall_best$方法))
      cat(sprintf("  框架: %s\n", overall_best$框架))
      cat(sprintf("  方法类型: %s\n", overall_best$方法类型))
      cat(sprintf("  预测类型: %s\n", overall_best$预测类型))
      cat(sprintf("  使用PCA: %s\n", overall_best$使用PCA))
      cat(sprintf("  交叉验证: %s\n", overall_best$交叉验证类型))
      cat(sprintf("  相关性: %.4f\n", overall_best$平均相关性))
      cat(sprintf("  RMSE: %.4f\n", overall_best$平均RMSE))
      cat(sprintf("  MAE: %.4f\n", overall_best$平均MAE))
      cat(sprintf("  记录数: %d\n", overall_best$记录数))
    }
  } else {
    cat("无有效的汇总统计数据\n")
  }

  return(summary_data)
}

results_list <- list()
if (exists("geno5_envloo_results") && !is.null(geno5_envloo_results) && nrow(geno5_envloo_results) > 0) {
  results_list[["环境留一基因型5折组合交叉验证"]] <- geno5_envloo_results
}
if (exists("geno5_fold_results") && !is.null(geno5_fold_results) && nrow(geno5_fold_results) > 0) {
  results_list[["基因型5折交叉验证"]] <- geno5_fold_results
}
if (exists("env_loo_results") && !is.null(env_loo_results) && nrow(env_loo_results) > 0) {
  results_list[["环境留一交叉验证"]] <- env_loo_results
}

if (length(results_list) > 0) {
  summary_stats <- generate_summary_statistics(results_list, gs_output_dir)

  cat("\n生成可视化比较图...\n")
  tryCatch({
    vis_data <- summary_stats

    if (nrow(vis_data) > 0) {
      library(ggplot2)

      p1 <- ggplot(vis_data, aes(x = 方法类型, y = 平均相关性, fill = 框架)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(预测类型 ~ 交叉验证类型) +
        labs(
          title = "TKW性状不同方法和预测类型的相关性比较",
          x = "方法类型", y = "平均相关性"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(file.path(gs_output_dir, "TKW相关性比较图.png"), p1, width = 16, height = 12, dpi = 300)
      cat("  相关性比较图已保存\n")

      p2 <- ggplot(vis_data, aes(x = 方法类型, y = 平均RMSE, fill = 框架)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(预测类型 ~ 交叉验证类型) +
        labs(
          title = "TKW性状不同方法和预测类型的RMSE比较",
          x = "方法类型", y = "平均RMSE"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(file.path(gs_output_dir, "TKW_RMSE比较图.png"), p2, width = 16, height = 12, dpi = 300)
      cat("  RMSE比较图已保存\n")

      method_comparison <- aggregate(
        cbind(平均相关性, 平均RMSE, 平均MAE) ~ 框架 + 方法类型 + 预测类型,
        data = vis_data, mean
      )
      write.csv(
        method_comparison,
        file.path(gs_output_dir, "TKW方法性能比较表.csv"),
        row.names = FALSE
      )
      cat("  方法性能比较表已保存\n")

      best_methods <- data.frame()
      for (cv_type in unique(vis_data$交叉验证类型)) {
        cv_data <- vis_data[vis_data$交叉验证类型 == cv_type, , drop = FALSE]
        if (nrow(cv_data) > 0) {
          best_idx <- which.max(cv_data$平均相关性)
          best_methods <- rbind(best_methods, cv_data[best_idx, , drop = FALSE])
        }
      }

      write.csv(
        best_methods,
        file.path(gs_output_dir, "TKW最佳方法推荐表.csv"),
        row.names = FALSE
      )
      cat("  最佳方法推荐表已保存\n")
    }
  }, error = function(e) {
    cat("  可视化生成失败:", e$message, "\n")
  })
} else {
  cat("没有有效的结果数据可汇总\n")
  summary_stats <- data.frame()
}

cat("\n=== 详细方法比较分析 ===\n")

if (exists("summary_stats") && nrow(summary_stats) > 0) {
  cat("\n按框架和预测类型分组比较:\n")
  cat(strrep("=", 60), "\n")

  framework_comparison <- aggregate(cbind(平均相关性, 平均RMSE) ~ 框架,
                                    data = summary_stats, mean)
  cat("\n整体框架比较:\n")
  for (i in seq_len(nrow(framework_comparison))) {
    cat(sprintf("  %-10s: 平均相关性 = %.4f, 平均RMSE = %.4f\n",
                framework_comparison$框架[i],
                framework_comparison$平均相关性[i],
                framework_comparison$平均RMSE[i]))
  }

  pred_type_comparison <- aggregate(cbind(平均相关性, 平均RMSE) ~ 预测类型,
                                    data = summary_stats, mean)
  cat("\n预测类型效果比较:\n")
  for (i in seq_len(nrow(pred_type_comparison))) {
    cat(sprintf("  %-10s: 平均相关性 = %.4f, 平均RMSE = %.4f\n",
                pred_type_comparison$预测类型[i],
                pred_type_comparison$平均相关性[i],
                pred_type_comparison$平均RMSE[i]))
  }

  method_type_comparison <- aggregate(cbind(平均相关性, 平均RMSE) ~ 方法类型,
                                      data = summary_stats, mean)
  cat("\n方法类型比较:\n")
  for (i in seq_len(nrow(method_type_comparison))) {
    cat(sprintf("  %-15s: 平均相关性 = %.4f, 平均RMSE = %.4f\n",
                method_type_comparison$方法类型[i],
                method_type_comparison$平均相关性[i],
                method_type_comparison$平均RMSE[i]))
  }

  cv_type_comparison <- aggregate(cbind(平均相关性, 平均RMSE) ~ 交叉验证类型,
                                  data = summary_stats, mean)
  cat("\n交叉验证类型比较:\n")
  for (i in seq_len(nrow(cv_type_comparison))) {
    cat(sprintf("  %-30s: 平均相关性 = %.4f, 平均RMSE = %.4f\n",
                cv_type_comparison$交叉验证类型[i],
                cv_type_comparison$平均相关性[i],
                cv_type_comparison$平均RMSE[i]))
  }

  cat("\n关键分析:\n")
  cat(strrep("-", 60), "\n")
  cat("1. 三种交叉验证的区别:\n")
  cat("   a) 环境留一基因型5折：测试基因型×测试环境（最严格）\n")
  cat("   b) 基因型5折：测试基因型×所有环境\n")
  cat("   c) 环境留一：所有基因型×测试环境\n")
  cat("   严格程度：a > b > c\n")
}

cat("\n")
cat(strrep("=", 80), "\n")
cat("TKW基因组选择分析完全结束!\n")
cat(strrep("=", 80), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", gs_output_dir, "\n")
cat("环境参数:", kPara_Name, "\n")
cat("执行顺序:", "1.环境留一基因型5折 → 2.基因型5折 → 3.环境留一\n")

cat("\n包含的配置:\n")
cat("----------------\n")
cat("方法 (6种):\n")
cat("  1. GAPIT_MAS - GAPIT MAS方法\n")
cat("  2. GAPIT_GBLUP - GAPIT GBLUP方法\n")
cat("  3. GAPIT_MAS+GBLUP - GAPIT MAS+GBLUP方法\n")
cat("  4. BGLR_MAS - BGLR MAS方法\n")
cat("  5. BGLR_GBLUP - BGLR GBLUP方法\n")
cat("  6. BGLR_MAS+GBLUP - BGLR MAS+GBLUP方法\n")
cat("预测方式 (2种):\n")
cat("  1. 截距 - 只预测截距项，直接作为表型预测值\n")
cat("  2. 截距+斜率 - 预测截距和斜率，使用FW模型计算表型预测值\n")
cat("总计: 6种方法 × 2种预测方式 = 12种配置\n")

cat("\n各交叉验证方法完成情况:\n")
cat("--------------------------------\n")
if (exists("geno5_envloo_results") && !is.null(geno5_envloo_results) && nrow(geno5_envloo_results) > 0) {
  cat("1. 环境留一基因型5折组合交叉验证:", nrow(geno5_envloo_results), "条记录\n")
} else {
  cat("1. 环境留一基因型5折组合交叉验证: 未完成或无结果\n")
}
if (exists("geno5_fold_results") && !is.null(geno5_fold_results) && nrow(geno5_fold_results) > 0) {
  cat("2. 基因型5折交叉验证:", nrow(geno5_fold_results), "条记录\n")
} else {
  cat("2. 基因型5折交叉验证: 未完成或无结果\n")
}
if (exists("env_loo_results") && !is.null(env_loo_results) && nrow(env_loo_results) > 0) {
  cat("3. 环境留一交叉验证:", nrow(env_loo_results), "条记录\n")
} else {
  cat("3. 环境留一交叉验证: 未完成或无结果\n")
}

cat("\n主要特点:\n")
cat("----------\n")
cat("1. 针对TKW性状的基因组选择分析\n")
cat("2. 简化版本：预测阶段无PCA，GWAS阶段保留PCA\n")
cat("3. 包含GAPIT和BGLR两种框架\n")
cat("4. 包含MAS、GBLUP、MAS+GBLUP三种方法类型\n")
cat("5. 包含两种预测方式：截距、截距+斜率\n")
cat("6. 三种交叉验证方法：环境留一、基因型5折、组合交叉验证\n")
cat("7. 完整的GWAS分析找到显著标记\n")
cat("8. LD过滤减少标记冗余\n")
cat("9. 详细的统计分析和可视化\n")
cat("10. 按照用户指定顺序执行：环境留一基因型5折 → 基因型5折 → 环境留一\n")
cat("11. 总计12种配置的比较分析\n")
cat("12. MAS+GBLUP背景标记过滤规则：删除显著标记3Mb内且LD>=0.7的背景标记\n")
cat("13. 最佳方法统计使用精确匹配，不再混淆 MAS 与 MAS+GBLUP\n")

# 关闭并行
stopCluster(cl)
cat("\n并行集群已关闭\n")
cat("分析结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("TKW基因组选择分析完全结束!\n")
cat("所有结果已保存到目录:", gs_output_dir, "\n")
cat("分析完成!\n")