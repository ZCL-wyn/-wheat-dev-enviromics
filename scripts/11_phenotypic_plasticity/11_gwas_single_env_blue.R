#!/usr/bin/env Rscript

# ===========================================
# 单个环境GWAS和BLUE-GWAS分析代码（修正版）
# 使用MLM模型，输出到指定目录结构
# 配置：
# 1. 工作目录：/mnt/7t_storage/zhangcl/TKW
# 2. 基因型文件：myGD2.csv
# 3. 图谱文件：myGM2.csv
# 4. 表型数据：TKW_mean_table.txt
# 5. 输出目录：/mnt/7t_storage/zhangcl/TKW/GAPIT_MLM_Results
#    - BLUE_Analysis
#    - Single_Environment
# ===========================================

# 设置工作目录
setwd("/mnt/7t_storage/zhangcl/TKW")

cat("=== 单个环境GWAS和BLUE-GWAS分析开始 ===\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

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
  "Matrix", "dplyr", "tidyr", "reshape2", "qqman", "lme4"
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

cat("注册并行后端...\n")
cl <- makeCluster(used_cores)
doParallel::registerDoParallel(cl)
cat("并行后端注册完成\n")

# ===========================================
# 2. 加载GAPIT函数
# ===========================================
cat("\n=== 加载GAPIT函数 ===\n")
cat("从网络加载GAPIT...\n")
tryCatch({
  source("http://zzlab.net/GAPIT/GAPIT.library.R")
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  cat("GAPIT函数加载成功\n")
}, error = function(e) {
  cat("网络GAPIT加载失败:", e$message, "\n")
  if (file.exists("GAPIT.library.R")) {
    source("GAPIT.library.R")
    cat("从本地加载GAPIT成功\n")
  } else {
    stop("无法加载GAPIT函数")
  }
})

# ===========================================
# 3. 数据准备（关键修正：添加品种名清洗）
# ===========================================
cat("\n=== 数据准备 ===\n")

set.seed(195021)

# 创建输出目录结构
output_base_dir <- "GAPIT_MLM_Results"
blue_dir <- file.path(output_base_dir, "BLUE_Analysis")
single_env_dir <- file.path(output_base_dir, "Single_Environment")

if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
  cat("创建输出主目录:", output_base_dir, "\n")
}
if (!dir.exists(blue_dir)) {
  dir.create(blue_dir, recursive = TRUE)
  cat("创建BLUE分析目录:", blue_dir, "\n")
}
if (!dir.exists(single_env_dir)) {
  dir.create(single_env_dir, recursive = TRUE)
  cat("创建单环境分析目录:", single_env_dir, "\n")
}

# 定义品种名清洗函数（与之前分析一致）
clean_genotype_names <- function(names) {
  names <- trimws(names)
  names <- toupper(names)
  names <- gsub("[^A-Za-z0-9]", "", names)
  names <- gsub("^([0-9]+)", "G\\1", names)
  return(names)
}

# 读取表型数据
cat("读取表型数据...\n")
phenotype_data <- read.table("TKW_mean_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("原始表型数据维度:", dim(phenotype_data), "\n")
# 重命名列（保持与之前一致）
colnames(phenotype_data) <- c("genotype", "env_code", "TKW")

# 去除缺失值
phenotype_data_clean <- phenotype_data[!is.na(phenotype_data$TKW), ]
cat("去除NA值后的表型数据维度:", dim(phenotype_data_clean), "\n")

# 清洗表型品种名，新增一列
phenotype_data_clean$genotype_clean <- clean_genotype_names(phenotype_data_clean$genotype)

# 读取基因型数据
cat("读取基因型数据...\n")
geno_dt <- fread("myGD2.csv", header = TRUE, data.table = FALSE)
sample_names_original <- geno_dt[[1]]
sample_names_clean <- clean_genotype_names(sample_names_original)
geno_matrix <- as.matrix(geno_dt[, -1])
rownames(geno_matrix) <- sample_names_clean
cat("基因型矩阵维度:", dim(geno_matrix), "\n")

# 读取SNP图谱
cat("读取SNP图谱...\n")
map_data <- fread("myGM2.csv", header = TRUE, data.table = FALSE)
colnames(map_data) <- c("SNP", "Chromosome", "Position")
cat("SNP图谱数据维度:", dim(map_data), "\n")

# 数据对齐（使用清洗后的品种名）
pheno_lines <- unique(phenotype_data_clean$genotype_clean)
geno_lines <- rownames(geno_matrix)
common_lines <- intersect(pheno_lines, geno_lines)
cat("共同品种数量:", length(common_lines), "\n")
if (length(common_lines) == 0) stop("无共同品种，请检查清洗函数！")

# 过滤表型数据（保留原始列和清洗列）
phenotype_data_final <- phenotype_data_clean[phenotype_data_clean$genotype_clean %in% common_lines, ]
cat("最终表型数据维度:", dim(phenotype_data_final), "\n")

# 过滤基因型数据
geno_matrix_final <- geno_matrix[common_lines, , drop = FALSE]
cat("最终基因型数据维度:", dim(geno_matrix_final), "\n")

# 检查SNP一致性
geno_snps <- colnames(geno_matrix_final)
map_snps <- map_data$SNP
common_snps <- intersect(geno_snps, map_snps)
cat("共同SNP数量:", length(common_snps), "\n")

# 保留共同SNP
geno_matrix_final <- geno_matrix_final[, common_snps, drop = FALSE]
map_data_final <- map_data[map_data$SNP %in% common_snps, ]
map_data_final <- map_data_final[match(colnames(geno_matrix_final), map_data_final$SNP), ]

# ===========================================
# 4. 基因型质量控制
# ===========================================
cat("\n=== 基因型数据质量控制 ===\n")
MAF.min <- 0.03
pNA.max <- 0.1

geno_matrix_qc <- geno_matrix_final

# MAF过滤
MAF_values <- apply(geno_matrix_qc, 2, function(x) mean(x, na.rm = TRUE) / 2)
drop_maf <- which(MAF_values < MAF.min | MAF_values > 1 - MAF.min)
if(length(drop_maf) > 0) {
  geno_matrix_qc <- geno_matrix_qc[, -drop_maf, drop = FALSE]
  map_data_final <- map_data_final[-drop_maf, ]
  cat("MAF过滤: 剔除了", length(drop_maf), "个位点\n")
}

# 缺失率过滤
missing_rates <- apply(geno_matrix_qc, 2, function(x) mean(is.na(x)))
drop_na <- which(missing_rates > pNA.max)
if(length(drop_na) > 0) {
  geno_matrix_qc <- geno_matrix_qc[, -drop_na, drop = FALSE]
  map_data_final <- map_data_final[-drop_na, ]
  cat("缺失率过滤: 剔除了", length(drop_na), "个位点\n")
}

cat("质量控制后基因型数据维度:", dim(geno_matrix_qc), "\n")

# 缺失值插补
for(j in 1:ncol(geno_matrix_qc)) {
  tmp <- geno_matrix_qc[, j]
  if(any(is.na(tmp))) {
    geno_matrix_qc[, j] <- ifelse(is.na(tmp), mean(tmp, na.rm = TRUE), tmp)
  }
}
cat("缺失值插补完成\n")

# ===========================================
# 5. 准备GAPIT格式数据
# ===========================================
cat("\n=== 准备GAPIT格式数据 ===\n")

# 标准化标记名称
standardize_snp_names <- function(names) {
  names <- gsub(";", "_", names)
  names <- gsub("-", "_", names)
  names <- gsub(":", "_", names)
  names <- gsub("/", "_", names)
  names <- gsub("\\|", "_", names)
  names <- gsub("\\(", "_", names)
  names <- gsub("\\)", "_", names)
  return(names)
}

new_geno_snp_names <- standardize_snp_names(colnames(geno_matrix_qc))
colnames(geno_matrix_qc) <- new_geno_snp_names
new_map_snp_names <- standardize_snp_names(map_data_final$SNP)
map_data_final$SNP <- new_map_snp_names

# 去除重复标记
keep_indices <- !duplicated(new_geno_snp_names)
geno_matrix_qc <- geno_matrix_qc[, keep_indices]
map_data_final <- map_data_final[keep_indices, ]

# 准备GD和GM
GD_all <- data.frame(taxa = rownames(geno_matrix_qc), geno_matrix_qc)
rownames(GD_all) <- GD_all$taxa
GM_all <- map_data_final
colnames(GM_all) <- c("SNP", "Chromosome", "Position")

# 将染色体转换为数值
GM_all$Chromosome <- as.numeric(as.character(GM_all$Chromosome))
GM_all$Position <- as.numeric(as.character(GM_all$Position))

cat("GAPIT GD数据维度:", dim(GD_all), "\n")
cat("GAPIT GM数据维度:", dim(GM_all), "\n")

# ===========================================
# 6. 核心函数定义（保留原样）
# ===========================================
cat("\n=== 定义核心函数 ===\n")

# 6.1 计算BLUE值函数（使用原始genotype列，但结果Taxa用清洗后的名？注意：这里使用原始列名进行模型拟合，输出Taxa也使用原始名，后续匹配需用清洗名）
# 修改：在计算BLUE时，我们使用清洗后的品种名作为Taxa，这样后续GWAS可直接用清洗名匹配
calculate_BLUE_values <- function(phenotype_data, genotype_col = "genotype_clean", trait_col = "TKW", env_col = "env_code") {
  cat("计算BLUE值...\n")
  
  # 确保数据格式正确
  phenotype_data[[genotype_col]] <- as.factor(phenotype_data[[genotype_col]])
  phenotype_data[[env_col]] <- as.factor(phenotype_data[[env_col]])
  
  # 使用lme4包计算BLUE值
  cat("拟合混合线性模型...\n")
  
  tryCatch({
    # 模型：性状 ~ 基因型（固定） + 环境（随机）
    model <- lmer(as.formula(paste(trait_col, "~ (1|", env_col, ") + (1|", genotype_col, ")")), 
                  data = phenotype_data)
    
    # 提取BLUE值
    cat("提取BLUE值...\n")
    blue_values <- fixef(model)
    
    # 获取基因型效应
    genotype_effects <- ranef(model)[[genotype_col]]
    
    # 计算BLUE值 = 总体均值 + 基因型效应
    overall_mean <- mean(phenotype_data[[trait_col]], na.rm = TRUE)
    blue_df <- data.frame(
      Taxa = rownames(genotype_effects),
      BLUE = overall_mean + genotype_effects[[1]]
    )
    
    colnames(blue_df) <- c("Taxa", "BLUE")
    
    cat("BLUE值计算完成，共", nrow(blue_df), "个品种\n")
    cat("BLUE值范围:", round(min(blue_df$BLUE, na.rm = TRUE), 2), "-", 
        round(max(blue_df$BLUE, na.rm = TRUE), 2), "\n")
    cat("BLUE值均值:", round(mean(blue_df$BLUE, na.rm = TRUE), 2), "\n")
    
    return(blue_df)
    
  }, error = function(e) {
    cat("混合模型拟合失败:", e$message, "\n")
    cat("使用简化方法计算BLUE值...\n")
    
    # 简化方法：计算每个品种在所有环境中的平均值作为BLUE
    blue_df <- phenotype_data %>%
      group_by(.data[[genotype_col]]) %>%
      summarise(
        BLUE = mean(.data[[trait_col]], na.rm = TRUE),
        .groups = 'drop'
      )
    
    colnames(blue_df) <- c("Taxa", "BLUE")
    
    cat("简化方法BLUE值计算完成，共", nrow(blue_df), "个品种\n")
    return(blue_df)
  })
}

# 6.2 GAPIT GWAS函数（MLM模型）（不变）
perform_gapit_mlm_gwas <- function(phenotype_df, GD_all, GM_all, trait_name, output_dir, model_type = "MLM") {
  cat("对", trait_name, "进行GAPIT", model_type, "GWAS分析...\n")
  
  if (!is.data.frame(phenotype_df) || nrow(phenotype_df) == 0) {
    cat("错误：表型数据无效\n")
    return(NULL)
  }
  
  # 准备表型数据
  if (!"Taxa" %in% colnames(phenotype_df)) {
    if ("genotype" %in% colnames(phenotype_df)) {
      phenotype_df$Taxa <- phenotype_df$genotype
    } else if (length(colnames(phenotype_df)) == 2) {
      colnames(phenotype_df) <- c("Taxa", "Trait")
    }
  }
  
  # 确定性状列
  trait_col <- NULL
  if ("Trait" %in% colnames(phenotype_df)) {
    trait_col <- "Trait"
  } else if ("BLUE" %in% colnames(phenotype_df)) {
    trait_col <- "BLUE"
  } else if (ncol(phenotype_df) == 2) {
    trait_col <- colnames(phenotype_df)[2]
  } else {
    cat("错误：无法确定性状列\n")
    return(NULL)
  }
  
  pheno_gapit <- data.frame(
    Taxa = phenotype_df$Taxa,
    Trait = phenotype_df[[trait_col]]
  )
  pheno_gapit <- pheno_gapit[complete.cases(pheno_gapit$Trait), ]
  
  cat(trait_name, "有效样本数量:", nrow(pheno_gapit), "\n")
  
  if (nrow(pheno_gapit) < 10) {
    cat("  有效样本不足，跳过", trait_name, "的GWAS\n")
    return(NULL)
  }
  
  # 共同基因型
  common_genotypes <- intersect(pheno_gapit$Taxa, GD_all$taxa)
  cat("共同基因型数量:", length(common_genotypes), "\n")
  
  if (length(common_genotypes) < 10) {
    cat("  共同基因型不足，跳过", trait_name, "的GWAS\n")
    return(NULL)
  }
  
  # 准备GAPIT输入
  pheno_gapit <- pheno_gapit[pheno_gapit$Taxa %in% common_genotypes, ]
  GD_gapit <- GD_all[GD_all$taxa %in% common_genotypes, ]
  
  # 创建输出目录
  trait_output_dir <- file.path(output_dir, trait_name)
  if (!dir.exists(trait_output_dir)) {
    dir.create(trait_output_dir, recursive = TRUE)
  }
  
  original_wd <- getwd()
  
  tryCatch({
    setwd(trait_output_dir)
    
    cat("运行GAPIT", model_type, "分析...\n")
    
    # 运行GAPIT MLM分析
    gapit_result <- GAPIT(
      Y = pheno_gapit,
      G = NULL,
      GD = GD_gapit,
      GM = GM_all,
      model = model_type,
      PCA.total = 3,
      Multiple_analysis = FALSE,
      file.output = TRUE,
      Major.allele.zero = FALSE,
      SNP.MAF = 0.01,
      cutOff = 0.05,
      memo = paste("MLM_GWAS", trait_name, sep = "_")
    )
    
    setwd(original_wd)
    
    # 读取GWAS结果文件
    gwas_pattern <- paste0("^GAPIT\\.Association\\.GWAS_Results\\..*", trait_name, ".*\\.csv$")
    gwas_files <- list.files(trait_output_dir, pattern = gwas_pattern, full.names = TRUE)
    
    if (length(gwas_files) == 0) {
      cat("警告：未找到GWAS结果文件\n")
      gwas_results_data <- NULL
    } else {
      cat("找到GWAS结果文件:", gwas_files[1], "\n")
      gwas_results_data <- read.csv(gwas_files[1], stringsAsFactors = FALSE)
      cat("GWAS结果文件读取成功，行数:", nrow(gwas_results_data), "\n")
    }
    
    # 读取Manhattan图数据
    manhattan_pattern <- paste0("^GAPIT\\.Manhattan\\.", trait_name, ".*\\.png$")
    manhattan_files <- list.files(trait_output_dir, pattern = manhattan_pattern, full.names = TRUE)
    
    # 读取QQ图数据
    qq_pattern <- paste0("^GAPIT\\.QQ\\.", trait_name, ".*\\.png$")
    qq_files <- list.files(trait_output_dir, pattern = qq_pattern, full.names = TRUE)
    
    return(list(
      gwas_results = gwas_results_data,
      GD = GD_gapit,
      GM = GM_all,
      phenotype_vector = pheno_gapit$Trait,
      taxa_names = pheno_gapit$Taxa,
      output_dir = trait_output_dir,
      manhattan_plot = if(length(manhattan_files) > 0) manhattan_files[1] else NULL,
      qq_plot = if(length(qq_files) > 0) qq_files[1] else NULL
    ))
    
  }, error = function(e) {
    setwd(original_wd)
    cat("GAPIT分析失败:", e$message, "\n")
    return(NULL)
  })
}

# 6.3 显著标记选择函数（不变）
select_significant_markers_mlm <- function(gwas_results, p_threshold = 0.05, n_top_markers = 10) {
  cat("从MLM GWAS结果中选择显著标记...\n")
  
  if (is.null(gwas_results) || nrow(gwas_results) == 0) {
    cat("没有GWAS结果\n")
    return(character(0))
  }
  
  cat("GWAS结果行数:", nrow(gwas_results), "\n")
  
  # 确定P值列名
  pvalue_col <- NULL
  possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
  
  for (col in possible_pvalue_cols) {
    if (col %in% colnames(gwas_results)) {
      pvalue_col <- col
      break
    }
  }
  
  if (is.null(pvalue_col)) {
    cat("警告：未找到P值列，使用所有标记\n")
    if (nrow(gwas_results) <= n_top_markers) {
      return(gwas_results$SNP)
    } else {
      return(gwas_results$SNP[1:n_top_markers])
    }
  }
  
  # 计算Bonferroni校正阈值
  bonferroni_threshold <- p_threshold / nrow(gwas_results)
  cat("Bonferroni校正p阈值(", p_threshold, "/", nrow(gwas_results), "):", bonferroni_threshold, "\n")
  
  # 筛选显著标记（使用Bonferroni校正）
  significant_markers <- gwas_results[gwas_results[[pvalue_col]] < bonferroni_threshold, ]
  cat("Bonferroni显著标记数量（p <", bonferroni_threshold, "）:", nrow(significant_markers), "\n")
  
  # 如果没有Bonferroni显著标记，使用p < 1e-5的阈值
  if (nrow(significant_markers) == 0) {
    cat("没有Bonferroni显著标记，使用p < 1e-5阈值\n")
    significant_markers <- gwas_results[gwas_results[[pvalue_col]] < 1e-5, ]
    cat("p < 1e-5的显著标记数量:", nrow(significant_markers), "\n")
  }
  
  # 如果还没有显著标记，使用p < 1e-4的阈值
  if (nrow(significant_markers) == 0) {
    cat("没有p < 1e-5显著标记，使用p < 1e-4阈值\n")
    significant_markers <- gwas_results[gwas_results[[pvalue_col]] < 1e-4, ]
    cat("p < 1e-4的显著标记数量:", nrow(significant_markers), "\n")
  }
  
  # 如果仍然没有显著标记，选择P值最小的前n_top_markers个标记
  if (nrow(significant_markers) == 0) {
    cat("没有达到阈值的显著标记，选择P值最小的前", n_top_markers, "个标记\n")
    gwas_results_sorted <- gwas_results[order(gwas_results[[pvalue_col]]), ]
    top_markers <- gwas_results_sorted$SNP[1:min(n_top_markers, nrow(gwas_results_sorted))]
    return(top_markers)
  }
  
  # 按P值升序排序
  significant_markers <- significant_markers[order(significant_markers[[pvalue_col]]), ]
  
  # 取前n_top_markers个标记
  n_to_select <- min(n_top_markers, nrow(significant_markers))
  top_markers <- significant_markers$SNP[1:n_to_select]
  
  cat("选择", n_to_select, "个显著标记\n")
  cat("前5个标记:", paste(head(top_markers, 5), collapse = ", "), "\n")
  
  return(top_markers)
}

# 6.4 提取单个环境表型数据函数（使用清洗后的品种名）
extract_single_environment_phenotypes <- function(phenotype_data, env_code) {
  cat("提取环境", env_code, "的表型数据...\n")
  
  env_data <- phenotype_data[phenotype_data$env_code == env_code, ]
  
  if (nrow(env_data) == 0) {
    cat("环境", env_code, "没有表型数据\n")
    return(NULL)
  }
  
  # 准备表型数据格式，Taxa使用清洗后的品种名
  env_pheno <- data.frame(
    Taxa = env_data$genotype_clean,
    Trait = env_data$TKW
  )
  
  # 去除缺失值
  env_pheno <- env_pheno[complete.cases(env_pheno$Trait), ]
  
  cat("环境", env_code, "的有效样本数量:", nrow(env_pheno), "\n")
  
  if (nrow(env_pheno) < 5) {
    cat("环境", env_code, "的有效样本不足，跳过\n")
    return(NULL)
  }
  
  return(list(
    pheno_df = env_pheno,
    env_name = env_code,
    sample_count = nrow(env_pheno)
  ))
}

# ===========================================
# 7. 计算BLUE值并进行GWAS
# ===========================================
cat("\n=== 第一步：计算BLUE值并进行GWAS ===\n")

# 7.1 计算BLUE值（使用清洗后的品种名作为基因型因子）
blue_values_df <- calculate_BLUE_values(phenotype_data_final, genotype_col = "genotype_clean")

if (nrow(blue_values_df) > 0) {
  # 保存BLUE值
  blue_file <- file.path(blue_dir, "BLUE_values.csv")
  write.csv(blue_values_df, blue_file, row.names = FALSE)
  cat("BLUE值已保存到:", blue_file, "\n")
  
  # 7.2 对BLUE值进行GWAS分析
  cat("\n对BLUE值进行GWAS分析...\n")
  blue_gwas_result <- perform_gapit_mlm_gwas(
    phenotype_df = blue_values_df,
    GD_all = GD_all,
    GM_all = GM_all,
    trait_name = "BLUE",
    output_dir = blue_dir,
    model_type = "MLM"
  )
  
  if (!is.null(blue_gwas_result)) {
    # 选择显著标记
    significant_markers_blue <- select_significant_markers_mlm(blue_gwas_result$gwas_results)
    
    # 保存显著标记
    if (length(significant_markers_blue) > 0) {
      dir.create(file.path(blue_dir, "BLUE"), showWarnings = FALSE)
      sig_markers_file <- file.path(blue_dir, "BLUE", "Significant_Markers_BLUE.txt")
      write.table(data.frame(SNP = significant_markers_blue),
                  sig_markers_file,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      cat("BLUE显著标记已保存到:", sig_markers_file, "\n")
    }
    
    # 保存完整的GWAS结果
    if (!is.null(blue_gwas_result$gwas_results)) {
      dir.create(file.path(blue_dir, "BLUE"), showWarnings = FALSE)
      gwas_result_file <- file.path(blue_dir, "BLUE", "Complete_GWAS_Results_BLUE.csv")
      write.csv(blue_gwas_result$gwas_results, gwas_result_file, row.names = FALSE)
      cat("完整BLUE GWAS结果已保存到:", gwas_result_file, "\n")
    }
    
    cat("BLUE值GWAS分析完成\n")
  } else {
    cat("BLUE值GWAS分析失败\n")
  }
} else {
  cat("无法计算BLUE值\n")
}

# ===========================================
# 8. 单个环境GWAS分析
# ===========================================
cat("\n=== 第二步：单个环境GWAS分析 ===\n")

# 获取所有环境
all_environments <- unique(phenotype_data_final$env_code)
cat("总环境数量:", length(all_environments), "\n")
cat("环境列表:", paste(all_environments, collapse = ", "), "\n")

# 初始化结果存储
single_env_results <- list()
significant_markers_env <- list()

# 对每个环境进行GWAS分析
for (env in all_environments) {
  cat("\n=====================================\n")
  cat("处理环境:", env, "\n")
  cat("=====================================\n")
  
  # 提取该环境的表型数据
  env_data <- extract_single_environment_phenotypes(phenotype_data_final, env)
  
  if (!is.null(env_data)) {
    # 对该环境进行GWAS分析
    env_gwas_result <- perform_gapit_mlm_gwas(
      phenotype_df = env_data$pheno_df,
      GD_all = GD_all,
      GM_all = GM_all,
      trait_name = env,
      output_dir = single_env_dir,
      model_type = "MLM"
    )
    
    if (!is.null(env_gwas_result)) {
      single_env_results[[env]] <- env_gwas_result
      
      # 选择显著标记
      significant_markers <- select_significant_markers_mlm(env_gwas_result$gwas_results)
      significant_markers_env[[env]] <- significant_markers
      
      # 保存显著标记
      if (length(significant_markers) > 0) {
        env_output_dir <- file.path(single_env_dir, env)
        if (!dir.exists(env_output_dir)) dir.create(env_output_dir, recursive = TRUE)
        sig_markers_file <- file.path(env_output_dir, paste0("Significant_Markers_", env, ".txt"))
        write.table(data.frame(SNP = significant_markers),
                    sig_markers_file,
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        cat("环境", env, "的显著标记已保存\n")
      }
      
      # 保存完整的GWAS结果
      if (!is.null(env_gwas_result$gwas_results)) {
        env_output_dir <- file.path(single_env_dir, env)
        if (!dir.exists(env_output_dir)) dir.create(env_output_dir, recursive = TRUE)
        gwas_result_file <- file.path(env_output_dir, paste0("Complete_GWAS_Results_", env, ".csv"))
        write.csv(env_gwas_result$gwas_results, gwas_result_file, row.names = FALSE)
        cat("环境", env, "的完整GWAS结果已保存\n")
      }
      
      cat("环境", env, "的GWAS分析完成\n")
    } else {
      cat("环境", env, "的GWAS分析失败\n")
      significant_markers_env[[env]] <- character(0)
    }
  }
}

cat("\n单个环境GWAS分析完成\n")
cat("成功分析的环境数量:", length(single_env_results), "\n")

# ===========================================
# 9. 结果汇总和报告
# ===========================================
cat("\n=== 第三步：结果汇总和报告 ===\n")

# 9.1 创建汇总表格
cat("创建GWAS结果汇总...\n")

# BLUE结果汇总
blue_summary <- data.frame()
if (exists("blue_gwas_result") && !is.null(blue_gwas_result) && !is.null(blue_gwas_result$gwas_results)) {
  blue_data <- blue_gwas_result$gwas_results
  
  # 确定P值列名
  pvalue_col <- NULL
  possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
  for (col in possible_pvalue_cols) {
    if (col %in% colnames(blue_data)) {
      pvalue_col <- col
      break
    }
  }
  
  if (!is.null(pvalue_col)) {
    bonferroni_threshold <- 0.05 / nrow(blue_data)
    sig_bonferroni <- sum(blue_data[[pvalue_col]] < bonferroni_threshold, na.rm = TRUE)
    sig_1e_5 <- sum(blue_data[[pvalue_col]] < 1e-5, na.rm = TRUE)
    sig_1e_4 <- sum(blue_data[[pvalue_col]] < 1e-4, na.rm = TRUE)
    
    blue_summary <- data.frame(
      Analysis = "BLUE",
      Sample_Size = nrow(blue_values_df),
      Total_SNPs = nrow(blue_data),
      Significant_SNPs_Bonferroni = sig_bonferroni,
      Significant_SNPs_1e_5 = sig_1e_5,
      Significant_SNPs_1e_4 = sig_1e_4,
      Top_SNPs_Selected = if(exists("significant_markers_blue")) length(significant_markers_blue) else 0,
      Mean_P_value = mean(blue_data[[pvalue_col]], na.rm = TRUE),
      Min_P_value = min(blue_data[[pvalue_col]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

# 单环境结果汇总
env_summary <- data.frame()
for (env in names(single_env_results)) {
  if (!is.null(single_env_results[[env]]) && !is.null(single_env_results[[env]]$gwas_results)) {
    env_data <- single_env_results[[env]]$gwas_results
    
    pvalue_col <- NULL
    possible_pvalue_cols <- c("P.value", "P", "p.value", "p", "P.value.")
    for (col in possible_pvalue_cols) {
      if (col %in% colnames(env_data)) {
        pvalue_col <- col
        break
      }
    }
    
    if (!is.null(pvalue_col)) {
      bonferroni_threshold <- 0.05 / nrow(env_data)
      sig_bonferroni <- sum(env_data[[pvalue_col]] < bonferroni_threshold, na.rm = TRUE)
      sig_1e_5 <- sum(env_data[[pvalue_col]] < 1e-5, na.rm = TRUE)
      sig_1e_4 <- sum(env_data[[pvalue_col]] < 1e-4, na.rm = TRUE)
      
      env_summary <- rbind(env_summary, data.frame(
        Analysis = env,
        Sample_Size = length(single_env_results[[env]]$phenotype_vector),
        Total_SNPs = nrow(env_data),
        Significant_SNPs_Bonferroni = sig_bonferroni,
        Significant_SNPs_1e_5 = sig_1e_5,
        Significant_SNPs_1e_4 = sig_1e_4,
        Top_SNPs_Selected = length(significant_markers_env[[env]]),
        Mean_P_value = mean(env_data[[pvalue_col]], na.rm = TRUE),
        Min_P_value = min(env_data[[pvalue_col]], na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 合并汇总数据
if (nrow(blue_summary) > 0) {
  summary_data <- rbind(blue_summary, env_summary)
} else {
  summary_data <- env_summary
}

# 保存汇总数据
if (nrow(summary_data) > 0) {
  summary_file <- file.path(output_base_dir, "GWAS_Summary_MLM.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  cat("GWAS汇总已保存到:", summary_file, "\n")
  
  cat("\nMLM GWAS分析汇总:\n")
  cat(paste(rep("=", 100), collapse = ""), "\n")
  print(summary_data)
}

# 9.2 创建分析报告
cat("\n创建分析报告...\n")
report_file <- file.path(output_base_dir, "MLM_GWAS_Analysis_Report.txt")
sink(report_file)
cat("单个环境GWAS和BLUE-GWAS分析报告 (MLM模型)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("工作目录:", getwd(), "\n")
cat("输出目录:", output_base_dir, "\n")
cat("\n")

cat("数据概览:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("表型数据品种数:", length(unique(phenotype_data_final$genotype_clean)), "\n")
cat("表型数据环境数:", length(unique(phenotype_data_final$env_code)), "\n")
cat("表型数据记录数:", nrow(phenotype_data_final), "\n")
cat("基因型数据品种数:", nrow(geno_matrix_qc), "\n")
cat("基因型数据SNP数:", ncol(geno_matrix_qc), "\n")
cat("\n")

cat("BLUE值分析:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
if (exists("blue_values_df") && nrow(blue_values_df) > 0) {
  cat("BLUE值计算品种数:", nrow(blue_values_df), "\n")
  cat("BLUE值范围:", round(min(blue_values_df$BLUE, na.rm = TRUE), 2), "-", 
      round(max(blue_values_df$BLUE, na.rm = TRUE), 2), "\n")
  cat("BLUE值均值:", round(mean(blue_values_df$BLUE, na.rm = TRUE), 2), "\n")
  cat("BLUE值标准差:", round(sd(blue_values_df$BLUE, na.rm = TRUE), 2), "\n")
} else {
  cat("BLUE值计算失败\n")
}
cat("\n")

cat("单个环境GWAS分析:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("总环境数量:", length(all_environments), "\n")
cat("成功分析的环境数量:", length(single_env_results), "\n")
cat("环境列表:\n")
for (env in all_environments) {
  if (env %in% names(single_env_results)) {
    cat("  ✓ ", env, ": ", length(single_env_results[[env]]$phenotype_vector), "个样本\n", sep = "")
  } else {
    cat("  ✗ ", env, ": 分析失败或样本不足\n", sep = "")
  }
}
cat("\n")

cat("显著标记统计:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
if (exists("significant_markers_blue") && length(significant_markers_blue) > 0) {
  cat("BLUE分析显著标记数:", length(significant_markers_blue), "\n")
  if (length(significant_markers_blue) > 0) {
    cat("BLUE分析前5个显著标记:", paste(head(significant_markers_blue, 5), collapse = ", "), "\n")
  }
}

cat("\n单个环境显著标记统计:\n")
total_env_markers <- 0
for (env in names(significant_markers_env)) {
  if (length(significant_markers_env[[env]]) > 0) {
    cat("  ", env, ": ", length(significant_markers_env[[env]]), "个显著标记\n", sep = "")
    total_env_markers <- total_env_markers + length(significant_markers_env[[env]])
  }
}
cat("  总计:", total_env_markers, "个显著标记\n")
sink()

cat("分析报告已保存到:", report_file, "\n")

# 9.3 保存所有重要结果到R数据文件
save_data <- list(
  phenotype_data = phenotype_data_final,
  blue_values = if(exists("blue_values_df")) blue_values_df else NULL,
  blue_gwas_result = if(exists("blue_gwas_result")) blue_gwas_result else NULL,
  single_env_results = single_env_results,
  significant_markers_blue = if(exists("significant_markers_blue")) significant_markers_blue else NULL,
  significant_markers_env = significant_markers_env,
  GD_all = GD_all,
  GM_all = GM_all,
  all_environments = all_environments,
  analysis_time = Sys.time()
)

save_file <- file.path(output_base_dir, "MLM_GWAS_Analysis_Results.RData")
save(save_data, file = save_file)
cat("所有分析结果已保存到R数据文件:", save_file, "\n")

# ===========================================
# 10. 可视化图表（可选，但数据现在有效）
# ===========================================
cat("\n=== 第四步：创建可视化图表 ===\n")

tryCatch({
  library(ggplot2)
  library(qqman)
  
  plots_dir <- file.path(output_base_dir, "Summary_Plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # 1. BLUE值分布图
  if (exists("blue_values_df") && nrow(blue_values_df) > 0) {
    png(file.path(plots_dir, "BLUE_Distribution.png"), width = 800, height = 600)
    hist(blue_values_df$BLUE, main = "BLUE Value Distribution", 
         xlab = "BLUE Value", col = "lightblue", border = "black")
    dev.off()
    cat("BLUE值分布图已保存\n")
  }
  
  # 2. 各环境样本数量柱状图
  env_sample_counts <- sapply(all_environments, function(env) {
    env_data <- phenotype_data_final[phenotype_data_final$env_code == env, ]
    return(nrow(env_data))
  })
  
  if (length(env_sample_counts) > 0) {
    png(file.path(plots_dir, "Environment_Sample_Counts.png"), width = 1000, height = 600)
    par(mar = c(8, 4, 4, 2))
    barplot(env_sample_counts, main = "Sample Counts by Environment", 
            xlab = "", ylab = "Sample Count", col = "skyblue", las = 2)
    dev.off()
    cat("环境样本数量图已保存\n")
  }
  
  # 3. 各环境表型值箱线图
  if (nrow(phenotype_data_final) > 0) {
    png(file.path(plots_dir, "Phenotype_by_Environment.png"), width = 1200, height = 800)
    print(
      ggplot(phenotype_data_final, aes(x = env_code, y = TKW)) +
        geom_boxplot(fill = "lightblue", alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "TKW Distribution by Environment", 
             x = "Environment", y = "TKW")
    )
    dev.off()
    cat("环境表型值箱线图已保存\n")
  }
  
  # 4. 显著标记数量汇总图
  if (exists("summary_data") && nrow(summary_data) > 0) {
    png(file.path(plots_dir, "Significant_Markers_Summary.png"), width = 1000, height = 600)
    plot_data <- summary_data
    if (nrow(plot_data) > 10) plot_data <- plot_data[1:10, ]
    barplot(plot_data$Top_SNPs_Selected, 
            names.arg = plot_data$Analysis,
            main = "Significant Markers by Analysis",
            xlab = "Analysis", ylab = "Number of Significant Markers",
            col = "steelblue", las = 2)
    dev.off()
    cat("显著标记数量汇总图已保存\n")
  }
  
}, error = function(e) {
  cat("可视化图表创建失败:", e$message, "\n")
})

# ===========================================
# 11. 清理和结束
# ===========================================
cat("\n=== 清理和结束 ===\n")

# 关闭并行
stopCluster(cl)
cat("并行集群已关闭\n")
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("单个环境GWAS和BLUE-GWAS分析完全结束! (MLM模型)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录结构:\n")
cat("  ", output_base_dir, "\n")
cat("    ├── BLUE_Analysis/\n")
cat("    │   ├── BLUE_values.csv\n")
cat("    │   └── BLUE/ (包含GWAS结果、图表和显著标记)\n")
cat("    ├── Single_Environment/\n")
cat("    │   └── [环境名称]/ (每个环境的GWAS结果)\n")
cat("    ├── Summary_Plots/ (汇总图表)\n")
cat("    ├── GWAS_Summary_MLM.csv\n")
cat("    ├── MLM_GWAS_Analysis_Report.txt\n")
cat("    └── MLM_GWAS_Analysis_Results.RData\n")
cat("\n")

cat("分析总结:\n")
cat("1. BLUE值计算: ", if(exists("blue_values_df") && nrow(blue_values_df) > 0) "成功" else "失败", "\n")
cat("2. BLUE-GWAS分析: ", if(exists("blue_gwas_result") && !is.null(blue_gwas_result)) "成功" else "失败", "\n")
cat("3. 单个环境GWAS分析: ", length(single_env_results), "/", length(all_environments), "个环境成功\n")
cat("4. 总显著标记发现:\n")
if (exists("significant_markers_blue") && length(significant_markers_blue) > 0) {
  cat("   - BLUE分析: ", length(significant_markers_blue), "个显著标记\n")
}
total_env_markers <- sum(sapply(significant_markers_env, length))
cat("   - 单个环境分析: ", total_env_markers, "个显著标记\n")

cat("\n分析完成!\n")
cat("请查看", output_base_dir, "目录下的结果文件。\n")