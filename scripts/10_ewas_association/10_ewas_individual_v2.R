# ============================================================== 
# 高性能 EWAS 分析系统（完整全组合版）
# 特点：
#   1) 不做任何组合筛选，运行基因型PC(最大6)×环境PC(最大6)全部组合
#   2) EC 名称反引号安全处理 
#   3) scale 全局一次（ECOV_scaled） 
#   4) 稳定并行（tryCatch + pass） 
#   5) bobyqa 优化器 
#   6) ECOV 零/近零方差过滤 
#   7) 最终以λ值最接近1作为最优组合选择
# ============================================================== 
suppressPackageStartupMessages({ 
  library(data.table) 
  library(doParallel) 
  library(foreach) 
  library(lme4) 
  library(poolr) 
})

cat("加载的库:\n") 
cat("  data.table v", as.character(packageVersion("data.table")), " - 高速数据处理\n") 
cat("  doParallel v", as.character(packageVersion("doParallel")), " - 并行计算框架\n") 
cat("  lme4 v", as.character(packageVersion("lme4")), " - 混合效应模型\n") 
cat("  poolr v", as.character(packageVersion("poolr")), " - 多重检验校正\n\n")

# ===============================
# 1. 基本设置（固定基因型/环境PC最大为6）
# ===============================
setwd("/mnt/7t_storage/zhangcl/TKW")

trait  <- "TKW"
comps  <- c("genotype", "env_code")

# 固定最大PC数：基因型6个，环境6个（不做筛选，运行全部组合）
max_geno_pcs <- 6
max_env_pcs  <- 6

# 移除筛选参数，直接运行全部组合（无top_n限制）
region <- "Across_regions"

num_cores <- parallel::detectCores() - 2
if(num_cores < 1) num_cores <- 1

cat("\n=== 分析设置（全组合模式）===\n")
cat("性状:", trait, "\n")
cat("随机效应:", paste(comps, collapse=" + "), "\n")
cat("基因型PC最大数:", max_geno_pcs, "（运行全部可能组合）\n")
cat("环境PC最大数:", max_env_pcs, "（运行全部可能组合）\n")
cat("并行核心数:", num_cores, "\n\n")

start_time_all <- Sys.time()

# 输出目录
output_dir <- file.path("4_ecov_WAS", "output", region)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================
# 2. 高效数据加载函数（完整保留，无删减）
# ===============================
cat("=== 数据加载 ===\n")

load_phenotype_data <- function() {
  pheno_file <- "TKW_mean_table.txt"
  cat("读取表型文件:", pheno_file, "\n")
  pheno <- fread(pheno_file, showProgress = TRUE)
  setDT(pheno)

  # 基因型列
  if (!"genotype" %in% names(pheno)) {
    if ("Genotype" %in% names(pheno)) setnames(pheno, "Genotype", "genotype")
    else if ("Line" %in% names(pheno)) setnames(pheno, "Line", "genotype")
    else if ("ID" %in% names(pheno)) setnames(pheno, "ID", "genotype")
    else stop("表型数据中未找到基因型列")
  }

  # 环境列
  if (!"env_code" %in% names(pheno)) {
    if ("Env" %in% names(pheno)) setnames(pheno, "Env", "env_code")
    else if ("Environment" %in% names(pheno)) setnames(pheno, "Environment", "env_code")
    else if ("Site" %in% names(pheno)) setnames(pheno, "Site", "env_code")
    else stop("表型数据中未找到环境代码列")
  }

  # 性状列
  if (!trait %in% names(pheno)) {
    if (paste0(trait,"_mean") %in% names(pheno)) setnames(pheno, paste0(trait,"_mean"), trait)
    else if ("Trait" %in% names(pheno)) setnames(pheno, "Trait", trait)
    else setnames(pheno, names(pheno)[1], trait)
  }

  pheno[, `:=`(
    year_loc = factor(as.character(env_code)),
    genotype = factor(as.character(genotype))
  )]

  # 删除 TKW_mean 列
  if ("TKW_mean" %in% names(pheno)) {
    pheno[, TKW_mean := NULL]
  }

  pheno <- pheno[!is.na(get(trait))]
  cat("有效观测数:", nrow(pheno), "\n")
  return(pheno)
}

load_environment_data <- function(pheno) {
  ECOV_path <- "ECs2.csv"  # 完整环境数据文件（不做筛选）
  cat("读取环境协变量文件:", ECOV_path, "\n")
  ECOV_df <- fread(ECOV_path, showProgress = TRUE)
  setDT(ECOV_df)

  envID <- ECOV_df[[1]]
  ECOV_df <- ECOV_df[envID != "" & !is.na(envID)]
  envID <- ECOV_df[[1]]
  ECOV_df <- ECOV_df[!duplicated(envID)]
  envID <- ECOV_df[[1]]

  # 转换为矩阵（保留全部环境变量，后续仅过滤近零方差）
  ECOV <- as.matrix(ECOV_df[, -1, with = FALSE])
  rownames(ECOV) <- envID

  # 对齐表型和环境数据
  common_envs <- intersect(levels(pheno$year_loc), rownames(ECOV))
  pheno <- pheno[year_loc %in% common_envs]
  pheno[, year_loc := factor(as.character(year_loc))]
  ECOV <- ECOV[levels(pheno$year_loc), , drop = FALSE]

  cat("对齐后环境数:", nrow(ECOV), "\n")
  cat("原始EC变量数:", ncol(ECOV), "\n")
  return(list(ECOV = ECOV, pheno = pheno))
}

load_genotype_data <- function(pheno) {
  geno_file <- "myGD2.csv"
  geno_rds  <- "GENO_optimized.rds"

  if (file.exists(geno_rds)) {
    cat("从RDS文件读取基因型数据...\n")
    geno <- readRDS(geno_rds)
  } else {
    cat("从CSV文件读取基因型数据...\n")
    geno_dt <- fread(geno_file, showProgress = TRUE)
    sample_names <- geno_dt[[1]]
    geno <- as.matrix(geno_dt[, -1, with = FALSE])
    rownames(geno) <- sample_names
    saveRDS(geno, geno_rds)
  }

  # 对齐基因型数据
  common_genotypes <- intersect(levels(pheno$genotype), rownames(geno))
  cat("共同基因型数量:", length(common_genotypes), "\n")
  if (length(common_genotypes) == 0) stop("没有找到共同的基因型！")

  pheno <- pheno[genotype %in% common_genotypes]
  pheno[, genotype := factor(as.character(genotype))]
  geno <- geno[common_genotypes, , drop = FALSE]
  cat("过滤后基因型矩阵维度:", dim(geno), "\n")

  return(list(geno = geno, pheno = pheno))
}

cat("\n1. 加载表型数据...\n")
pheno <- load_phenotype_data()

cat("\n2. 加载环境数据...\n")
env_data <- load_environment_data(pheno)
ECOV <- env_data$ECOV
pheno <- env_data$pheno

cat("\n3. 加载基因型数据...\n")
geno_data <- load_genotype_data(pheno)
geno <- geno_data$geno
pheno <- geno_data$pheno

cat("\n✅ 数据加载完成\n")
cat("  表型记录数:", nrow(pheno), "\n")
cat("  基因型数:", length(levels(pheno$genotype)), "\n")
cat("  环境数:", length(levels(pheno$year_loc)), "\n")
cat("  原始EC变量数:", ncol(ECOV), "\n")

# ===============================
# 3. PCA 计算（全局一次，保留全部PC）
# ===============================
cat("\n=== 计算主成分（全维度）===\n")

# ---- 环境 PCA（保留全部有效PC，最大6个）
W_scaled <- scale(ECOV)
W_scaled[is.na(W_scaled) | is.infinite(W_scaled)] <- 0
WWt <- tcrossprod(W_scaled) / ncol(W_scaled)
W_eigen <- eigen(WWt)

W_values  <- W_eigen$values
W_vectors <- W_eigen$vectors
rownames(W_vectors) <- rownames(WWt)

# 保留全部有效PC（排除近零特征值），最多不超过max_env_pcs(6)
W_index <- which(W_values > 1e-8)
W_PC_full <- W_vectors[, W_index, drop = FALSE]
W_PC_full <- sweep(W_PC_full, 2, sqrt(W_values[W_index]), FUN = "*")
# 截断到最大6个PC
W_PC_full <- W_PC_full[, 1:min(ncol(W_PC_full), max_env_pcs), drop = FALSE]
W_values <- W_values[1:min(length(W_values), max_env_pcs)]

# Meff & Bonferroni（用于后续显著性检验）
W_EVD <- list(values = W_values, vectors = W_vectors)
nt_obj <- poolr::meff(eigen = W_EVD$values[W_index], method = "galwey")
nt <- if (is.list(nt_obj)) nt_obj$Meff else nt_obj
bonferroni_threshold <- 0.05 / as.numeric(nt)
cat("环境有效独立检验数量(Meff):", as.numeric(nt), "\n")
cat("Bonferroni阈值:", bonferroni_threshold, "\n")
cat("环境PC实际维度:", ncol(W_PC_full), "（最大6个）\n")

# ---- 基因型 PCA（保留全部有效PC，最大6个）
X <- scale(geno, center = TRUE, scale = FALSE)
G <- tcrossprod(X) / ncol(X)
G_eigen <- eigen(G)

G_values  <- G_eigen$values
G_vectors <- G_eigen$vectors
rownames(G_vectors) <- rownames(G)

# 保留全部有效PC（排除近零特征值），最多不超过max_geno_pcs(6)
G_index <- which(G_values > 1e-8)
G_PC_full <- G_vectors[, G_index, drop = FALSE]
G_PC_full <- sweep(G_PC_full, 2, sqrt(G_values[G_index]), FUN = "*")
# 截断到最大6个PC
G_PC_full <- G_PC_full[, 1:min(ncol(G_PC_full), max_geno_pcs), drop = FALSE]
G_values <- G_values[1:min(length(G_values), max_geno_pcs)]

# 基因型方差解释比例
G_prop <- G_values / sum(G_values)
G_cumprop <- cumsum(G_prop)
W_cumprop_full <- cumsum(W_values / sum(W_values))

cat("基因型PC实际维度:", ncol(G_PC_full), "（最大6个）\n")
cat("基因型方差解释比例（前6个PC）:", round(tail(G_cumprop, 1)*100, 2), "%\n")
cat("环境方差解释比例（前6个PC）:", round(tail(W_cumprop_full, 1)*100, 2), "%\n")

# ---- 关键：ECOV_scaled 全局一次性计算（避免函数内重复 scale）
ECOV_scaled <- scale(ECOV)
ECOV_scaled[is.na(ECOV_scaled) | is.infinite(ECOV_scaled)] <- 0

# ---- ECOV过滤：仅去掉零方差/近零方差列（保留全部有效环境变量）
ec_sd <- apply(ECOV_scaled, 2, sd, na.rm = TRUE)
keep_ec <- which(ec_sd > 1e-8)
if (length(keep_ec) < ncol(ECOV_scaled)) {
  cat("ECOV过滤：移除近零方差变量数 =", ncol(ECOV_scaled) - length(keep_ec), "\n")
}
ECOV_scaled <- ECOV_scaled[, keep_ec, drop = FALSE]
cat("过滤后保留EC变量数:", ncol(ECOV_scaled), "（全部有效环境变量）\n")

# ---- lmer 控制：bobyqa 更稳（保留完整参数）
LME_CTRL <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# ===============================
# 4. 核心函数定义（完整保留，无筛选逻辑）
# ===============================

generate_all_pc_combinations <- function(n_pcs) {
  """生成1到n_pcs的全部可能PC组合（无智能筛选，全组合）"""
  if (n_pcs < 1) return(list())
  
  all_combinations <- list()
  
  # 生成所有可能的非空子集（1到n_pcs的全部组合）
  for (k in 1:n_pcs) {
    combos <- combn(1:n_pcs, k, simplify = FALSE)
    all_combinations <- c(all_combinations, combos)
  }
  
  # 额外添加空组合（无PC的情况，确保全覆盖）
  all_combinations <- c(list(integer(0)), all_combinations)
  
  cat(sprintf("生成 %d 个完整PC组合（含空组合）\n", length(all_combinations)))
  return(all_combinations)
}

calculate_geno_model_fit <- function(pheno_data, geno_pc_combo, trait_name) {
  if (length(geno_pc_combo) == 0) {
    G_PC_subset <- matrix(0, nrow = nrow(G_PC_full), ncol = 0)
  } else {
    G_PC_subset <- G_PC_full[, geno_pc_combo, drop = FALSE]
    colnames(G_PC_subset) <- paste0("GPC", geno_pc_combo)
  }

  geno_indices <- match(as.character(pheno_data$genotype), rownames(G_PC_subset))
  GPC_obs <- G_PC_subset[geno_indices, , drop = FALSE]

  dat <- data.table(
    y = pheno_data[[trait_name]],
    pheno_data[, comps, with = FALSE],
    GPC_obs
  )

  formula_random <- paste(paste0("(1|", comps, ")"), collapse = "+")
  if (ncol(GPC_obs) > 0) {
    formula_fixed <- paste(colnames(GPC_obs), collapse = "+")
    formula_str <- paste("y ~", formula_random, "+", formula_fixed)
  } else {
    formula_str <- paste("y ~", formula_random)
  }

  result <- tryCatch({
    fm <- lmer(as.formula(formula_str), data = dat, REML = FALSE, control = LME_CTRL)
    list(
      n_pcs = length(geno_pc_combo),
      pc_combo = paste(geno_pc_combo, collapse = ","),
      AIC = AIC(fm),
      BIC = BIC(fm),
      logLik = as.numeric(logLik(fm)),
      variance_explained = ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)], 0),
      success = TRUE
    )
  }, error = function(e) {
    list(
      n_pcs = length(geno_pc_combo),
      pc_combo = paste(geno_pc_combo, collapse = ","),
      AIC = NA_real_, BIC = NA_real_, logLik = NA_real_,
      variance_explained = ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)], 0),
      success = FALSE
    )
  })
  return(result)
}

calculate_lambda <- function(p_values) {
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0) return(NA_real_)
  chi_sq_observed <- qchisq(1 - p_values, df = 1, lower.tail = FALSE)
  median(chi_sq_observed) / qchisq(0.5, df = 1)
}

calculate_complete_model_lambda <- function(geno_pc_combo, env_pc_combo,
                                           pheno_data = pheno,
                                           eco_matrix_scaled = ECOV_scaled) {
  # 基因型PC子集（全组合对应）
  if (length(geno_pc_combo) == 0) {
    G_PC_subset <- matrix(0, nrow = nrow(G_PC_full), ncol = 0)
  } else {
    G_PC_subset <- G_PC_full[, geno_pc_combo, drop = FALSE]
    colnames(G_PC_subset) <- paste0("GPC", geno_pc_combo)
  }

  # 环境PC子集（全组合对应）
  if (length(env_pc_combo) == 0) {
    W_PC_subset <- matrix(0, nrow = nrow(W_PC_full), ncol = 0)
  } else {
    W_PC_subset <- W_PC_full[, env_pc_combo, drop = FALSE]
    colnames(W_PC_subset) <- paste0("WPC", env_pc_combo)
  }

  geno_indices <- match(as.character(pheno_data$genotype), rownames(G_PC_subset))
  env_indices  <- match(as.character(pheno_data$year_loc), rownames(W_PC_subset))

  GPC_obs <- G_PC_subset[geno_indices, , drop = FALSE]
  WPC_obs <- W_PC_subset[env_indices, , drop = FALSE]

  # EC 矩阵（观测级，保留全部过滤后环境变量）
  eco_obs <- eco_matrix_scaled[as.character(pheno_data$year_loc), , drop = FALSE]

  dat <- data.table(
    y = pheno_data[[trait]],
    pheno_data[, comps, with = FALSE],
    GPC_obs,
    WPC_obs
  )

  formula_random <- paste(paste0("(1|", comps, ")"), collapse = "+")
  fixed_terms <- c(colnames(GPC_obs), colnames(WPC_obs))
  if (length(fixed_terms) > 0) {
    formula_fixed <- paste(fixed_terms, collapse = "+")
    formula0 <- paste("y ~", formula_random, "+", formula_fixed)
  } else {
    formula0 <- paste("y ~", formula_random)
  }

  result <- tryCatch({
    fm0 <- lmer(as.formula(formula0), data = dat, REML = FALSE, control = LME_CTRL)
    ll0 <- logLik(fm0); df0 <- attr(ll0, "df")

    p_values <- rep(NA_real_, ncol(eco_obs))
    chisq_values <- rep(NA_real_, ncol(eco_obs))

    for (j in seq_len(ncol(eco_obs))) {
      ec_name <- colnames(eco_obs)[j]

      # EC 名称反引号安全处理（避免特殊字符报错）
      if (grepl("[^a-zA-Z0-9_]", ec_name)) ec_term <- paste0("`", ec_name, "`")
      else ec_term <- ec_name

      dat_ec <- copy(dat)
      dat_ec[, (ec_name) := eco_obs[, j]]
      formula1 <- paste0(formula0, " + ", ec_term)

      fit_result <- tryCatch({
        fm1 <- lmer(as.formula(formula1), data = dat_ec, REML = FALSE, control = LME_CTRL)
        ll1 <- logLik(fm1); df1 <- attr(ll1, "df")
        chisq <- 2 * (as.numeric(ll1) - as.numeric(ll0))
        p_val <- pchisq(chisq, df = df1 - df0, lower.tail = FALSE)
        list(chisq = chisq, p_val = p_val)
      }, error = function(e) {
        list(chisq = NA_real_, p_val = NA_real_)
      })

      chisq_values[j] <- fit_result$chisq
      p_values[j] <- fit_result$p_val
    }

    # 计算λ值（核心评价指标，最接近1为最优）
    lambda_val <- calculate_lambda(p_values)
    lambda_diff <- abs(lambda_val - 1)  # 与1的差异，越小越好

    # 方差解释比例
    geno_var_exp <- ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)], 0)
    env_var_exp  <- ifelse(length(env_pc_combo) > 0, W_cumprop_full[max(env_pc_combo)], 0)

    list(
      geno_pc_combo = paste(geno_pc_combo, collapse = ","),
      env_pc_combo  = paste(env_pc_combo, collapse = ","),
      n_geno_pcs = length(geno_pc_combo),
      n_env_pcs  = length(env_pc_combo),
      lambda = lambda_val,
      lambda_diff = lambda_diff,  # 排序核心依据
      n_sig_05 = sum(p_values < 0.05, na.rm = TRUE),
      n_sig_bonf = sum(p_values < bonferroni_threshold, na.rm = TRUE),
      median_p = median(p_values, na.rm = TRUE),
      mean_chisq = mean(chisq_values, na.rm = TRUE),
      AIC = AIC(fm0),
      BIC = BIC(fm0),
      geno_variance_explained = geno_var_exp * 100,
      env_variance_explained  = env_var_exp * 100,
      success = TRUE
    )
  }, error = function(e) {
    list(
      geno_pc_combo = paste(geno_pc_combo, collapse = ","),
      env_pc_combo  = paste(env_pc_combo, collapse = ","),
      n_geno_pcs = length(geno_pc_combo),
      n_env_pcs  = length(env_pc_combo),
      lambda = NA_real_,
      lambda_diff = NA_real_,
      n_sig_05 = NA_integer_,
      n_sig_bonf = NA_integer_,
      median_p = NA_real_,
      mean_chisq = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      geno_variance_explained = ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)] * 100, 0),
      env_variance_explained  = ifelse(length(env_pc_combo) > 0, W_cumprop_full[max(env_pc_combo)] * 100, 0),
      success = FALSE
    )
  })

  return(result)
}

# ===============================
# 5. 生成基因型+环境PC全部组合（无筛选）
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("生成全部PC组合（基因型6×环境6）\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# 生成基因型PC全部组合（最大6个）
geno_combinations <- generate_all_pc_combinations(n_pcs = ncol(G_PC_full))
cat("基因型PC全部组合数:", length(geno_combinations), "\n")

# 生成环境PC全部组合（最大6个）
env_combinations <- generate_all_pc_combinations(n_pcs = ncol(W_PC_full))
cat("环境PC全部组合数:", length(env_combinations), "\n")

# 构建全组合任务网格（无批次，无筛选，全部运行）
task_grid <- expand.grid(
  geno_idx = 1:length(geno_combinations),
  env_idx = 1:length(env_combinations),
  stringsAsFactors = FALSE
)
cat("总任务数（基因型×环境）:", nrow(task_grid), "（全部运行，无删减）\n")

# 解析组合（用于后续计算）
geno_combos_parsed <- lapply(1:length(geno_combinations), function(i) {
  geno_combinations[[i]]
})
env_combos_parsed <- lapply(1:length(env_combinations), function(i) {
  env_combinations[[i]]
})

# ===============================
# 6. 并行计算全部组合的λ值（核心步骤，无筛选）
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("并行计算全部组合的λ值（最接近1为最优）\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat(sprintf("启动并行计算（使用 %d 个核心）...\n", num_cores))
cl <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(cl)

# 导出所有必要变量到集群节点
clusterExport(cl, varlist = c(
  "pheno", "trait", "ECOV_scaled",
  "G_PC_full", "W_PC_full",
  "G_cumprop", "W_values", "W_cumprop_full",
  "comps", "bonferroni_threshold",
  "calculate_complete_model_lambda", "calculate_lambda",
  "geno_combos_parsed", "env_combos_parsed",
  "LME_CTRL"
), envir = environment())

cat("开始全组合λ值计算（耗时较长，请耐心等待）...\n")
full_combination_results <- foreach(
  task_idx = 1:nrow(task_grid),
  .packages = c("lme4", "data.table", "stats", "poolr"),
  .combine = "rbind",
  .errorhandling = "remove"  # 跳过拟合失败的组合，保留有效结果
) %dopar% {
  geno_idx <- task_grid$geno_idx[task_idx]
  env_idx  <- task_grid$env_idx[task_idx]

  geno_combo <- geno_combos_parsed[[geno_idx]]
  env_combo  <- env_combos_parsed[[env_idx]]

  res <- calculate_complete_model_lambda(
    geno_combo, env_combo,
    pheno_data = pheno,
    eco_matrix_scaled = ECOV_scaled
  )

  if (isTRUE(res$success)) {
    data.frame(
      geno_pc_combo = res$geno_pc_combo,
      env_pc_combo  = res$env_pc_combo,
      n_geno_pcs = res$n_geno_pcs,
      n_env_pcs  = res$n_env_pcs,
      lambda = res$lambda,
      lambda_diff = res$lambda_diff,  # 与1的差异，核心排序依据
      n_sig_05 = res$n_sig_05,
      n_sig_bonf = res$n_sig_bonf,
      median_p = res$median_p,
      mean_chisq = res$mean_chisq,
      AIC = res$AIC,
      BIC = res$BIC,
      geno_variance_explained = res$geno_variance_explained,
      env_variance_explained  = res$env_variance_explained,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

# 停止集群，释放内存
stopCluster(cl)
gc()

# 转换为data.table，整理结果
full_combination_dt <- as.data.table(full_combination_results)
full_combination_dt <- full_combination_dt[!is.na(lambda)]  # 过滤λ值为NA的无效结果

# 按lambda_diff排序（越小越接近1，最优组合排在最前）
setorder(full_combination_dt, lambda_diff)

# 保存全组合结果（无删减，完整保留）
full_combination_file <- file.path(output_dir, "full_combination_lambda_results.csv")
fwrite(full_combination_dt, full_combination_file)
cat("\n✅ 全组合λ值计算完成！\n")
cat("  有效组合数（拟合成功）:", nrow(full_combination_dt), "\n")
cat("  全组合结果已保存:", full_combination_file, "\n")

# 展示前10个最优组合（λ最接近1）
cat("\n⭐ 前10个最优组合（λ值最接近1，按lambda_diff排序）:\n")
print(head(full_combination_dt[, .(geno_pc_combo, env_pc_combo, lambda, lambda_diff, AIC, BIC)], 10))

# ===============================
# 7. 选择最优组合（λ最接近1）
# ===============================
if (nrow(full_combination_dt) == 0) {
  warning("⚠️ 所有组合拟合失败，自动fallback到：无PC + 环境PC=1")
  best_combination <- list(
    geno_pc_combo = "",
    env_pc_combo  = "1",
    n_geno_pcs = 0,
    n_env_pcs  = 1,
    lambda = NA,
    lambda_diff = NA,
    geno_variance_explained = 0,
    env_variance_explained  = W_cumprop_full[1] * 100
  )
  best_geno_pcs <- integer(0)
  best_env_pcs  <- 1
} else {
  # 选择lambda_diff最小的组合（最接近1）
  best_combination <- full_combination_dt[1]
  best_geno_pcs <- if (best_combination$geno_pc_combo == "") integer(0) else as.numeric(strsplit(best_combination$geno_pc_combo, ",")[[1]])
  best_env_pcs  <- if (best_combination$env_pc_combo == "") integer(0) else as.numeric(strsplit(best_combination$env_pc_combo, ",")[[1]])
}

cat("\n⭐ 最终最优组合（λ最接近1）:\n")
cat(sprintf("   基因型PC组合: %s\n", if (length(best_geno_pcs) == 0) "无" else paste(best_geno_pcs, collapse=",")))
cat(sprintf("   环境PC组合: %s\n", if (length(best_env_pcs) == 0) "无" else paste(best_env_pcs, collapse=",")))
cat(sprintf("   基因型PC数量: %d\n", length(best_geno_pcs)))
cat(sprintf("   环境PC数量: %d\n", length(best_env_pcs)))
if (!is.na(best_combination$lambda)) {
  cat(sprintf("   λ值: %.4f\n", best_combination$lambda))
  cat(sprintf("   与1的差异: %.4f\n", best_combination$lambda_diff))
} else {
  cat("   λ值: NA (fallback)\n")
}
cat(sprintf("   解释的基因型方差: %.1f%%\n", best_combination$geno_variance_explained))
cat(sprintf("   解释的环境方差: %.1f%%\n", best_combination$env_variance_explained))

# ===============================
# 8. 最终 EWAS 分析（基于最优组合，完整环境变量）
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("最终阶段：基于最优组合的完整EWAS分析（并行LRT）\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# ---- 准备最优组合的PC观测级矩阵
if (length(best_geno_pcs) > 0) {
  G_PC_subset <- G_PC_full[, best_geno_pcs, drop = FALSE]
  colnames(G_PC_subset) <- paste0("GPC", best_geno_pcs)
} else {
  G_PC_subset <- matrix(0, nrow = nrow(G_PC_full), ncol = 0)
}

if (length(best_env_pcs) > 0) {
  W_PC_subset <- W_PC_full[, best_env_pcs, drop = FALSE]
  colnames(W_PC_subset) <- paste0("WPC", best_env_pcs)
} else {
  W_PC_subset <- matrix(0, nrow = nrow(W_PC_full), ncol = 0)
}

geno_indices <- match(as.character(pheno$genotype), rownames(G_PC_subset))
env_indices  <- match(as.character(pheno$year_loc), rownames(W_PC_subset))

GPC_obs <- G_PC_subset[geno_indices, , drop = FALSE]
WPC_obs <- W_PC_subset[env_indices, , drop = FALSE]

# 环境协变量（观测级，保留全部过滤后有效变量）
W_obs <- ECOV_scaled[as.character(pheno$year_loc), , drop = FALSE]

# 构建分析数据集
dat <- data.table(
  y = pheno[[trait]],
  pheno[, comps, with = FALSE],
  GPC_obs,
  WPC_obs
)

# 构建基础模型公式
formula_random <- paste(paste0("(1|", comps, ")"), collapse = "+")
fixed_terms <- c(colnames(GPC_obs), colnames(WPC_obs))

if (length(fixed_terms) > 0) {
  formula_fixed <- paste(fixed_terms, collapse = "+")
  formula0 <- paste("y ~", formula_random, "+", formula_fixed)
} else {
  formula0 <- paste("y ~", formula_random)
}

cat("最优组合基础模型公式:\n  ", formula0, "\n\n")

# 拟合基础模型
fm0 <- lmer(as.formula(formula0), data = dat, REML = FALSE, control = LME_CTRL)
ll0 <- logLik(fm0); df0 <- attr(ll0, "df")

cat("基础模型拟合成功。\n")
cat("开始逐个环境变量做似然比检验（总数：", ncol(W_obs), "，全部运行）\n")

# ---- 并行运行EWAS似然比检验
cl <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(cl)

# 导出必要变量到集群节点
clusterExport(cl, varlist = c("dat", "formula0", "ll0", "df0", "W_obs", "trait", "LME_CTRL"), envir = environment())

ewas_results <- foreach(
  j = 1:ncol(W_obs),
  .packages = c("lme4", "data.table", "stats"),
  .combine = "rbind",
  .errorhandling = "pass"
) %dopar% {
  ec_name <- colnames(W_obs)[j]

  # EC名称特殊字符处理
  if (grepl("[^a-zA-Z0-9_]", ec_name)) ec_term <- paste0("`", ec_name, "`")
  else ec_term <- ec_name

  dat_ec <- copy(dat)
  dat_ec[, (ec_name) := W_obs[, j]]

  formula1 <- paste0(formula0, " + ", ec_term)

  out <- tryCatch({
    fm1 <- lmer(as.formula(formula1), data = dat_ec, REML = FALSE, control = LME_CTRL)
    ll1 <- logLik(fm1); df1 <- attr(ll1, "df")
    chisq_val <- 2*(as.numeric(ll1) - as.numeric(ll0))
    p_val <- pchisq(chisq_val, df = df1 - df0, lower.tail = FALSE)

    data.frame(
      trait = trait,
      ecov  = ec_name,
      Chisq = chisq_val,
      Chi_Df = df1 - df0,
      Pr_Chisq = p_val,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      trait = trait,
      ecov  = ec_name,
      Chisq = NA_real_,
      Chi_Df = NA_real_,
      Pr_Chisq = NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out
}

# 停止集群，释放内存
stopCluster(cl)
gc()

# 整理最终EWAS结果
out <- as.data.table(ewas_results)
out[, neg_log10_p := -log10(Pr_Chisq)]
pvals <- out[!is.na(Pr_Chisq), Pr_Chisq]
final_lambda_value <- calculate_lambda(pvals)

cat(sprintf("\n✅ EWAS分析完成！最终λ值: %.4f\n", final_lambda_value))
setorder(out, Pr_Chisq)

# 补充最优组合信息到结果中
out[, `:=`(
  lambda = final_lambda_value,
  geno_pc_combo = if (length(best_geno_pcs) == 0) "" else paste(best_geno_pcs, collapse=","),
  env_pc_combo  = if (length(best_env_pcs) == 0) "" else paste(best_env_pcs, collapse=","),
  n_geno_pcs = length(best_geno_pcs),
  n_env_pcs  = length(best_env_pcs),
  geno_variance_explained = best_combination$geno_variance_explained,
  env_variance_explained  = best_combination$env_variance_explained,
  bonferroni_threshold = bonferroni_threshold
)]

# ===============================
# 9. 保存完整结果 + 详细报告（无删减）
# ===============================
cat("\n=== 保存完整结果（无删减）===\n")

# 保存最终EWAS结果
result_file_csv <- file.path(output_dir, paste0("ecov_WAS_results_", trait, "_final_full_combination.csv"))
fwrite(out, result_file_csv)
cat("最终EWAS全结果已保存:", result_file_csv, "\n")

# 保存RData格式（包含所有中间结果）
result_file_rdata <- file.path(output_dir, paste0("ecov_WAS_results_", trait, "_full_combination.RData"))
save(out, fm0, full_combination_dt, best_combination, final_lambda_value,
     file = result_file_rdata)
cat("RData格式完整结果已保存:", result_file_rdata, "\n")

# 生成详细最终报告
final_report_file <- file.path(output_dir, paste0("EWAS_final_report_", trait, "_full_combination.txt"))
sink(final_report_file)

cat("高性能EWAS分析报告（全组合版，λ最接近1为最优）\n")
cat("==================================================\n")
cat("分析日期:", date(), "\n")
cat("性状:", trait, "\n")
cat("随机效应:", paste(comps, collapse=" + "), "\n")
cat("并行核心数:", num_cores, "\n")
cat("基因型PC最大数:", max_geno_pcs, "\n")
cat("环境PC最大数:", max_env_pcs, "\n")
cat("全组合任务数:", nrow(task_grid), "\n")
cat("有效拟合组合数:", nrow(full_combination_dt), "\n\n")

cat("最优组合（λ最接近1）:\n")
cat("----------------------\n")
cat("基因型PC组合:", if (length(best_geno_pcs) == 0) "无" else paste(best_geno_pcs, collapse=","), "\n")
cat("环境PC组合:",  if (length(best_env_pcs) == 0) "无" else paste(best_env_pcs, collapse=","), "\n")
cat("基因型PC数量:", length(best_geno_pcs), "\n")
cat("环境PC数量:",  length(best_env_pcs), "\n")
cat("解释的基因型方差:", round(best_combination$geno_variance_explained, 1), "%\n")
cat("解释的环境方差:",  round(best_combination$env_variance_explained, 1), "%\n")
cat("最终λ值:", round(final_lambda_value, 4), "\n")
cat("λ与1的差异:", round(best_combination$lambda_diff, 4), "\n\n")

cat("基础模型公式:\n")
cat("-------------\n")
cat(formula0, "\n\n")

# 显著性统计（全部环境变量）
sig_05   <- sum(out$Pr_Chisq < 0.05, na.rm = TRUE)
sig_01   <- sum(out$Pr_Chisq < 0.01, na.rm = TRUE)
sig_001  <- sum(out$Pr_Chisq < 0.001, na.rm = TRUE)
sig_bonf <- sum(out$Pr_Chisq < bonferroni_threshold, na.rm = TRUE)
total_ec <- ncol(W_obs)

cat("显著性结果统计（全部", total_ec, "个EC变量）:\n")
cat("--------------------------------------------\n")
cat("p < 0.05:", sig_05, "/", total_ec, "(", round(sig_05/total_ec*100, 2), "%)\n")
cat("p < 0.01:", sig_01, "/", total_ec, "(", round(sig_01/total_ec*100, 2), "%)\n")
cat("p < 0.001:", sig_001, "/", total_ec, "(", round(sig_001/total_ec*100, 2), "%)\n")
cat("p < Bonferroni:", sig_bonf, "/", total_ec, "(", round(sig_bonf/total_ec*100, 2), "%)\n\n")

# 展示前20个最显著的EC变量
if (sig_05 > 0) {
  cat("前20个最显著的ECOV变量:\n")
  cat("------------------------\n")
  for (i in 1:min(20, nrow(out))) {
    row <- out[i]
    if (!is.na(row$Pr_Chisq)) {
      cat(sprintf("%2d. %-40s p=%.2e  -log10(p)=%.2f\n",
                  i, substr(row$ecov, 1, 40), row$Pr_Chisq, row$neg_log10_p))
    }
  }
}

# 模型评估
cat("\n模型评估:\n")
cat("---------\n")
if (final_lambda_value > 1.05) {
  cat("⚠️ λ偏高(>1.05)，可能存在残留群体结构或未校正混杂因素\n")
} else if (final_lambda_value < 0.95) {
  cat("⚠️ λ偏低(<0.95)，可能存在过度校正或数据异常\n")
} else {
  cat("✅ λ在合理范围(0.95-1.05)，模型校正效果良好\n")
}

sink()
cat("最终详细报告已保存:", final_report_file, "\n")

# ===============================
# 10. 分析完成总结
# ===============================
end_time_all <- Sys.time()
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("✅ 高性能EWAS全组合分析全部完成！\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("总运行时间:", difftime(end_time_all, start_time_all, units = "hours"), "（小时）\n")
cat("输出目录:", output_dir, "\n")
cat("核心结果：\n")
cat("  1. 全组合λ值结果：full_combination_lambda_results.csv\n")
cat("  2. 最终EWAS结果：ecov_WAS_results_", trait, "_final_full_combination.csv\n")
cat("  3. 详细分析报告：EWAS_final_report_", trait, "_full_combination.txt\n")