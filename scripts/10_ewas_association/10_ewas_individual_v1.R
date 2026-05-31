# ============================================================== 
# 高性能 EWAS 分析系统（终极完整修复版） 
# Fix:
#   1) EC 名称反引号 
#   2) scale 全局一次（ECOV_scaled） 
#   3) Stage2 fallback 
#   4) 更稳并行（tryCatch + pass） 
#   5) bobyqa 优化器 
#   6) ECOV 零/近零方差过滤 
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
# 1. 基本设置
# ===============================
setwd("/mnt/7t_storage/zhangcl/TKW")

trait  <- "TKW"
comps  <- c("genotype", "env_code")

max_geno_pcs <- 5
max_env_pcs  <- 6

top_n_geno_combinations <- 50
top_n_env_combinations  <- 20

region <- "Across_regions"

num_cores <- parallel::detectCores() - 2
if(num_cores < 1) num_cores <- 1

cat("\n=== 分析设置 ===\n")
cat("性状:", trait, "\n")
cat("随机效应:", paste(comps, collapse=" + "), "\n")
cat("最大基因型PC数:", max_geno_pcs, "\n")
cat("最大环境PC数:", max_env_pcs, "\n")
cat("并行核心数:", num_cores, "\n\n")

start_time_all <- Sys.time()

# 输出目录
output_dir <- file.path("4_ecov_WAS", "output", region)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================
# 2. 高效数据加载函数
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
  ECOV_path <- "ECs2.csv"  # 更新为新的环境文件
  cat("读取环境协变量文件:", ECOV_path, "\n")
  ECOV_df <- fread(ECOV_path, showProgress = TRUE)
  setDT(ECOV_df)

  envID <- ECOV_df[[1]]
  ECOV_df <- ECOV_df[envID != "" & !is.na(envID)]
  envID <- ECOV_df[[1]]
  ECOV_df <- ECOV_df[!duplicated(envID)]
  envID <- ECOV_df[[1]]

  ECOV <- as.matrix(ECOV_df[, -1, with = FALSE])
  rownames(ECOV) <- envID

  common_envs <- intersect(levels(pheno$year_loc), rownames(ECOV))
  pheno <- pheno[year_loc %in% common_envs]
  pheno[, year_loc := factor(as.character(year_loc))]
  ECOV <- ECOV[levels(pheno$year_loc), , drop = FALSE]

  cat("对齐后环境数:", nrow(ECOV), "\n")
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
cat("  EC变量数:", ncol(ECOV), "\n")

# ===============================
# 3. PCA 计算（全局一次）
# ===============================
cat("\n=== 计算主成分 ===\n")

# ---- 环境 PCA
W_scaled <- scale(ECOV)
W_scaled[is.na(W_scaled) | is.infinite(W_scaled)] <- 0
WWt <- tcrossprod(W_scaled) / ncol(W_scaled)
W_eigen <- eigen(WWt)

W_values  <- W_eigen$values
W_vectors <- W_eigen$vectors
rownames(W_vectors) <- rownames(WWt)

W_index <- which(W_values > 1e-8)
W_PC_full <- W_vectors[, W_index, drop = FALSE]
W_PC_full <- sweep(W_PC_full, 2, sqrt(W_values[W_index]), FUN = "*")

# Meff & Bonferroni
W_EVD <- list(values = W_values, vectors = W_vectors)
nt_obj <- poolr::meff(eigen = W_EVD$values[W_index], method = "galwey")
nt <- if (is.list(nt_obj)) nt_obj$Meff else nt_obj
bonferroni_threshold <- 0.05 / as.numeric(nt)
cat("有效独立检验数量(Meff):", as.numeric(nt), "\n")
cat("Bonferroni阈值:", bonferroni_threshold, "\n\n")

# ---- 基因型 PCA
X <- scale(geno, center = TRUE, scale = FALSE)
G <- tcrossprod(X) / ncol(X)
G_eigen <- eigen(G)

G_values  <- G_eigen$values
G_vectors <- G_eigen$vectors
rownames(G_vectors) <- rownames(G)

G_index <- which(G_values > 1e-8)
G_PC_full <- G_vectors[, G_index, drop = FALSE]
G_PC_full <- sweep(G_PC_full, 2, sqrt(G_values[G_index]), FUN = "*")

G_prop <- G_values / sum(G_values)
G_cumprop <- cumsum(G_prop)

max_geno_pcs <- min(max_geno_pcs, ncol(G_PC_full))
max_env_pcs  <- min(max_env_pcs, ncol(W_PC_full))
cat("实际使用的最大基因型PC数:", max_geno_pcs, "\n")
cat("实际使用的最大环境PC数:", max_env_pcs, "\n")

# ---- 关键：ECOV_scaled 全局一次性计算（避免函数内重复 scale）
ECOV_scaled <- scale(ECOV)
ECOV_scaled[is.na(ECOV_scaled) | is.infinite(ECOV_scaled)] <- 0

# ---- ECOV过滤：去掉零方差/近零方差列（加速+减少拟合失败）
ec_sd <- apply(ECOV_scaled, 2, sd, na.rm = TRUE)
keep_ec <- which(ec_sd > 1e-8)
if (length(keep_ec) < ncol(ECOV_scaled)) {
  cat("ECOV过滤：移除近零方差变量数 =", ncol(ECOV_scaled) - length(keep_ec), "\n")
}
ECOV_scaled <- ECOV_scaled[, keep_ec, drop = FALSE]
cat("过滤后EC变量数:", ncol(ECOV_scaled), "\n")

# ---- lmer 控制：bobyqa 更稳
LME_CTRL <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)
# ===============================
# 4. 核心函数定义（修复版）
# ===============================

generate_smart_pc_combinations <- function(n_pcs, max_to_include = NULL, values = NULL) {
  if (is.null(max_to_include)) max_to_include <- n_pcs
  max_to_include <- min(max_to_include, n_pcs)

  cat(sprintf("智能生成PC组合 (n=%d, max_k=%d)\n", n_pcs, max_to_include))
  all_combinations <- list()

  # 连续组合
  for (k in 1:max_to_include) all_combinations[[length(all_combinations)+1]] <- 1:k

  # 基于解释方差阈值
  if (!is.null(values)) {
    prop <- values[1:n_pcs] / sum(values[1:n_pcs])
    cumprop <- cumsum(prop)
    thresholds <- seq(0.1, 0.95, by = 0.05)
    for (threshold in thresholds) {
      k <- which(cumprop >= threshold)[1]
      if (!is.na(k) && k <= max_to_include) {
        combo <- 1:k
        combo_str <- paste(combo, collapse = ",")
        if (!combo_str %in% sapply(all_combinations, function(x) paste(x, collapse=","))) {
          all_combinations[[length(all_combinations)+1]] <- combo
        }
      }
    }
  }

  # 前5个PC任意组合
  if (max_to_include >= 3) {
    important_pcs <- 1:min(5, max_to_include)
    for (k in 2:length(important_pcs)) {
      combos <- combn(important_pcs, k, simplify = FALSE)
      for (combo in combos) {
        combo_str <- paste(combo, collapse = ",")
        if (!combo_str %in% sapply(all_combinations, function(x) paste(x, collapse=","))) {
          all_combinations[[length(all_combinations)+1]] <- combo
        }
      }
    }
  }

  # 特征值断点
  if (max_to_include >= 3 && !is.null(values)) {
    ratios <- values[1:(max_to_include-1)] / values[2:max_to_include]
    cutoff <- which(ratios > 2)
    for (cut in cutoff) {
      combo <- 1:cut
      combo_str <- paste(combo, collapse = ",")
      if (!combo_str %in% sapply(all_combinations, function(x) paste(x, collapse=","))) {
        all_combinations[[length(all_combinations)+1]] <- combo
      }
    }
  }

  cat(sprintf("  生成 %d 个独特的PC组合\n", length(all_combinations)))
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
  # 基因型PC
  if (length(geno_pc_combo) == 0) {
    G_PC_subset <- matrix(0, nrow = nrow(G_PC_full), ncol = 0)
  } else {
    G_PC_subset <- G_PC_full[, geno_pc_combo, drop = FALSE]
    colnames(G_PC_subset) <- paste0("GPC", geno_pc_combo)
  }

  # 环境PC
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

  # EC 矩阵（观测级）
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

      # EC 名称反引号安全处理
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

    lambda_val <- calculate_lambda(p_values)

    n_sig_05  <- sum(p_values < 0.05, na.rm = TRUE)
    n_sig_bonf <- sum(p_values < bonferroni_threshold, na.rm = TRUE)

    geno_var_exp <- ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)], 0)
    W_cumprop_full <- cumsum(W_values / sum(W_values))
    env_var_exp  <- ifelse(length(env_pc_combo) > 0, W_cumprop_full[max(env_pc_combo)], 0)

    list(
      geno_pc_combo = paste(geno_pc_combo, collapse = ","),
      env_pc_combo  = paste(env_pc_combo, collapse = ","),
      n_geno_pcs = length(geno_pc_combo),
      n_env_pcs  = length(env_pc_combo),
      lambda = lambda_val,
      lambda_diff = abs(lambda_val - 1),
      n_sig_05 = n_sig_05,
      n_sig_bonf = n_sig_bonf,
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
      lambda = NA_real_, lambda_diff = NA_real_,
      n_sig_05 = NA_integer_, n_sig_bonf = NA_integer_,
      median_p = NA_real_, mean_chisq = NA_real_,
      AIC = NA_real_, BIC = NA_real_,
      geno_variance_explained = ifelse(length(geno_pc_combo) > 0, G_cumprop[max(geno_pc_combo)] * 100, 0),
      env_variance_explained  = ifelse(length(env_pc_combo) > 0, cumsum(W_values/sum(W_values))[max(env_pc_combo)] * 100, 0),
      success = FALSE
    )
  })

  return(result)
}

# ===============================
# 5. 第一阶段：筛选基因型PC组合
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("第一阶段：基于AIC/BIC筛选基因型PC组合 (并行)\n")
cat(paste(rep("=", 80), collapse=""), "\n")

geno_combinations <- generate_smart_pc_combinations(
  n_pcs = max_geno_pcs,
  max_to_include = max_geno_pcs,
  values = G_values
)
cat("总基因型PC组合数:", length(geno_combinations), "\n")

cat(sprintf("启动并行计算 (使用 %d 个核心)...\n", num_cores))
cl <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(cl)

clusterExport(cl, varlist = c("pheno", "trait", "G_PC_full", "G_cumprop", "comps",
                              "calculate_geno_model_fit", "LME_CTRL"),
              envir = environment())

cat("并行计算AIC/BIC...\n")
stage1_results <- foreach(
  i = 1:length(geno_combinations),
  .packages = c("lme4", "data.table", "stats"),
  .combine = "rbind",
  .errorhandling = "remove"
) %dopar% {
  combo <- geno_combinations[[i]]
  res <- calculate_geno_model_fit(pheno, combo, trait)
  if (isTRUE(res$success)) {
    data.frame(
      n_pcs = res$n_pcs,
      pc_combo = res$pc_combo,
      AIC = res$AIC,
      BIC = res$BIC,
      logLik = res$logLik,
      variance_explained = res$variance_explained,
      stringsAsFactors = FALSE
    )
  } else NULL
}

stopCluster(cl); gc()

stage1_dt <- as.data.table(stage1_results)
stage1_dt <- stage1_dt[!is.na(AIC)]
setorder(stage1_dt, AIC)
top_geno_combinations <- head(stage1_dt, top_n_geno_combinations)

cat("\n✅ 第一阶段完成！有效组合数:", nrow(stage1_dt), "\n")
cat("前10个最优基因型PC组合（按AIC）:\n")
print(head(top_geno_combinations[, .(pc_combo, n_pcs, AIC, BIC, variance_explained)], 10))

stage1_file <- file.path(output_dir, "stage1_geno_pc_selection.csv")
fwrite(stage1_dt, stage1_file)
cat("第一阶段结果已保存:", stage1_file, "\n")

# ===============================
# 6. 第二阶段：环境PC组合 + λ优化
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("第二阶段：生成环境PC组合 & 计算λ值 (并行)\n")
cat(paste(rep("=", 80), collapse=""), "\n")

env_combinations <- generate_smart_pc_combinations(
  n_pcs = max_env_pcs,
  max_to_include = max_env_pcs,
  values = W_values
)
cat("总环境PC组合数:", length(env_combinations), "\n")

geno_indices <- 1:nrow(top_geno_combinations)
env_indices  <- 1:length(env_combinations)

task_grid <- expand.grid(
  geno_idx = geno_indices,
  env_idx = env_indices,
  stringsAsFactors = FALSE
)
cat("需要计算的组合总数:", nrow(task_grid), "\n")

geno_combos_parsed <- lapply(geno_indices, function(i) {
  as.numeric(strsplit(top_geno_combinations$pc_combo[i], ",")[[1]])
})

batch_size <- 500
total_batches <- ceiling(nrow(task_grid) / batch_size)

all_stage2_results <- list()
batch_counter <- 1

for (batch in 1:total_batches) {
  cat(sprintf("\n批次 %d/%d (%.1f%%)...\n", batch, total_batches, batch/total_batches*100))

  start_idx <- (batch - 1) * batch_size + 1
  end_idx   <- min(batch * batch_size, nrow(task_grid))
  batch_tasks <- task_grid[start_idx:end_idx, ]

  cl <- makeCluster(num_cores, type = "PSOCK")
  registerDoParallel(cl)

  clusterExport(cl, varlist = c(
    "pheno","trait","ECOV_scaled",
    "G_PC_full","W_PC_full",
    "G_cumprop","W_values",
    "comps","bonferroni_threshold",
    "calculate_complete_model_lambda",
    "calculate_lambda",
    "geno_combos_parsed","env_combinations",
    "LME_CTRL"
  ), envir = environment())

  batch_results <- foreach(
    task_idx = 1:nrow(batch_tasks),
    .packages = c("lme4","data.table","stats","poolr"),
    .combine = "rbind",
    .errorhandling = "remove"
  ) %dopar% {
    geno_idx <- batch_tasks$geno_idx[task_idx]
    env_idx  <- batch_tasks$env_idx[task_idx]

    geno_combo <- geno_combos_parsed[[geno_idx]]
    env_combo  <- env_combinations[[env_idx]]

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
        lambda_diff = res$lambda_diff,
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
    } else NULL
  }

  stopCluster(cl); gc()

  if (!is.null(batch_results) && nrow(batch_results) > 0) {
    batch_dt <- as.data.table(batch_results)
    batch_file <- file.path(output_dir, sprintf("stage2_batch_%03d.csv", batch))
    fwrite(batch_dt, batch_file)

    all_stage2_results[[batch_counter]] <- batch_dt
    batch_counter <- batch_counter + 1

    cat(sprintf("  批次完成，有效结果: %d\n", nrow(batch_dt)))
  } else {
    cat("  批次无有效结果\n")
  }
}

cat("\n合并所有批次结果...\n")
if (length(all_stage2_results) > 0) {
  stage2_dt <- rbindlist(all_stage2_results, use.names = TRUE, fill = TRUE)
} else {
  stage2_dt <- data.table()
}

stage2_dt <- stage2_dt[!is.na(lambda)]
setorder(stage2_dt, lambda_diff)

stage2_file <- file.path(output_dir, "stage2_lambda_optimization.csv")
fwrite(stage2_dt, stage2_file)
cat("第二阶段结果已保存:", stage2_file, "\n")

# ===============================
# 7. 最优组合选择（含 fallback）
# ===============================
if (nrow(stage2_dt) == 0) {
  warning("⚠️ 第二阶段 lambda 全为 NA（可能仍存在模型拟合问题），自动 fallback 到：AIC最优基因型组合 + 环境PC=1")
  best_combination <- list(
    geno_pc_combo = top_geno_combinations$pc_combo[1],
    env_pc_combo  = "1",
    n_geno_pcs = length(strsplit(top_geno_combinations$pc_combo[1], ",")[[1]]),
    n_env_pcs  = 1,
    lambda = NA,
    lambda_diff = NA,
    geno_variance_explained = top_geno_combinations$variance_explained[1] * 100,
    env_variance_explained  = cumsum(W_values/sum(W_values))[1] * 100
  )
  best_geno_pcs <- as.numeric(strsplit(best_combination$geno_pc_combo, ",")[[1]])
  best_env_pcs  <- 1
} else {
  best_combination <- stage2_dt[1]
  best_geno_pcs <- as.numeric(strsplit(best_combination$geno_pc_combo, ",")[[1]])
  best_env_pcs  <- as.numeric(strsplit(best_combination$env_pc_combo, ",")[[1]])
}

cat("\n⭐ 最终选择的最佳组合:\n")
cat(sprintf("   基因型PC组合: %s\n", paste(best_geno_pcs, collapse=",")))
cat(sprintf("   环境PC组合: %s\n", paste(best_env_pcs, collapse=",")))
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
# 8. 最终 EWAS 分析
# ===============================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("最终阶段：完整EWAS分析（并行 LRT）\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# ---- 准备 PC 观测级矩阵
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

# 环境协变量（观测级）
W_obs <- ECOV_scaled[as.character(pheno$year_loc), , drop = FALSE]

dat <- data.table(
  y = pheno[[trait]],
  pheno[, comps, with = FALSE],
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

cat("完整模型公式:\n  ", formula0, "\n\n")

fm0 <- lmer(as.formula(formula0), data = dat, REML = FALSE, control = LME_CTRL)
ll0 <- logLik(fm0); df0 <- attr(ll0, "df")

cat("基础模型拟合成功。\n")
cat("开始逐个环境变量做似然比检验，总数：", ncol(W_obs), "\n")

# ---- 并行 EWAS
cl <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(cl)

clusterExport(cl, varlist = c("dat","formula0","ll0","df0","W_obs","trait","LME_CTRL"), envir = environment())

ewas_results <- foreach(
  j = 1:ncol(W_obs),
  .packages = c("lme4","data.table","stats"),
  .combine = "rbind",
  .errorhandling = "pass"
) %dopar% {
  ec_name <- colnames(W_obs)[j]

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

stopCluster(cl); gc()

out <- as.data.table(ewas_results)
out[, neg_log10_p := -log10(Pr_Chisq)]
pvals <- out[!is.na(Pr_Chisq), Pr_Chisq]
final_lambda_value <- calculate_lambda(pvals)

cat(sprintf("\n✅ EWAS分析完成！最终λ值: %.4f\n", final_lambda_value))
setorder(out, Pr_Chisq)

out[, `:=`(
  lambda = final_lambda_value,
  geno_pc_combo = paste(best_geno_pcs, collapse=","),
  env_pc_combo  = paste(best_env_pcs, collapse=","),
  n_geno_pcs = length(best_geno_pcs),
  n_env_pcs  = length(best_env_pcs),
  geno_variance_explained = best_combination$geno_variance_explained,
  env_variance_explained  = best_combination$env_variance_explained,
  bonferroni_threshold = bonferroni_threshold
)]

# ===============================
# 9. 保存最终结果 + 报告
# ===============================
cat("\n=== 保存最终结果 ===\n")

result_file_csv <- file.path(output_dir, paste0("ecov_WAS_results_", trait, "_final.csv"))
fwrite(out, result_file_csv)
cat("最终EWAS结果已保存:", result_file_csv, "\n")

result_file_rdata <- file.path(output_dir, paste0("ecov_WAS_results_", trait, ".RData"))
save(out, fm0, stage1_dt, stage2_dt, best_combination, final_lambda_value,
     file = result_file_rdata)
cat("RData格式结果已保存:", result_file_rdata, "\n")

final_report_file <- file.path(output_dir, paste0("EWAS_final_report_", trait, ".txt"))
sink(final_report_file)

cat("高性能EWAS分析报告\n")
cat("==================\n")
cat("分析日期:", date(), "\n")
cat("性状:", trait, "\n")
cat("随机效应:", paste(comps, collapse=" + "), "\n")
cat("并行核心数:", num_cores, "\n")
cat("基因型最大PC数:", max_geno_pcs, "\n")
cat("环境最大PC数:", max_env_pcs, "\n\n")

cat("最佳PC组合:\n")
cat("-----------\n")
cat("基因型PC组合:", paste(best_geno_pcs, collapse=","), "\n")
cat("环境PC组合:",  paste(best_env_pcs, collapse=","), "\n")
cat("基因型PC数量:", length(best_geno_pcs), "\n")
cat("环境PC数量:",  length(best_env_pcs), "\n")
cat("解释的基因型方差:", round(best_combination$geno_variance_explained, 1), "%\n")
cat("解释的环境方差:",  round(best_combination$env_variance_explained, 1), "%\n")
cat("最终λ值:", round(final_lambda_value, 4), "\n\n")

cat("模型公式:\n")
cat("---------\n")
cat(formula0, "\n\n")

sig_05   <- sum(out$Pr_Chisq < 0.05, na.rm = TRUE)
sig_01   <- sum(out$Pr_Chisq < 0.01, na.rm = TRUE)
sig_001  <- sum(out$Pr_Chisq < 0.001, na.rm = TRUE)
sig_bonf <- sum(out$Pr_Chisq < bonferroni_threshold, na.rm = TRUE)

cat("显著性结果统计:\n")
cat("---------------\n")
cat("总EC变量数:", ncol(W_obs), "\n")
cat("p < 0.05:", sig_05, "/", ncol(W_obs), "(", round(sig_05/ncol(W_obs)*100, 2), "%)\n")
cat("p < 0.01:", sig_01, "/", ncol(W_obs), "(", round(sig_01/ncol(W_obs)*100, 2), "%)\n")
cat("p < 0.001:", sig_001, "/", ncol(W_obs), "(", round(sig_001/ncol(W_obs)*100, 2), "%)\n")
cat("p < Bonferroni:", sig_bonf, "/", ncol(W_obs), "(", round(sig_bonf/ncol(W_obs)*100, 2), "%)\n\n")

if (sig_05 > 0) {
  cat("前20个最显著的ECOV:\n")
  for (i in 1:min(20, nrow(out))) {
    row <- out[i]
    if (!is.na(row$Pr_Chisq)) {
      cat(sprintf("%2d. %-40s p=%.2e  -log10(p)=%.2f\n",
                  i, substr(row$ecov, 1, 40), row$Pr_Chisq, row$neg_log10_p))
    }
  }
}

cat("\n模型评估:\n")
cat("---------\n")
if (final_lambda_value > 1.05) {
  cat("⚠️ λ偏高(>1.05)，可能残留结构\n")
} else if (final_lambda_value < 0.95) {
  cat("⚠️ λ偏低(<0.95)，可能过度校正\n")
} else {
  cat("✅ λ在合理范围(0.95-1.05)\n")
}

sink()
cat("最终报告已保存:", final_report_file, "\n")

end_time_all <- Sys.time()
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("✅ 高性能EWAS分析全部完成！\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("总运行时间:", end_time_all - start_time_all, "\n")
cat("输出目录:", output_dir, "\n")
