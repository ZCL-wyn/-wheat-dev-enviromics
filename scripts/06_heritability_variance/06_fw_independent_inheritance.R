#!/usr/bin/env Rscript
# ==============================================================================
# 小麦千粒重（TKW）FW参数计算 + 二元GBLUP分析（无Bootstrap；含LRT+p敏感性梯度）
#
# 目标性状：
#   α = Intcp_para_adj  （以DTR环境指数为自变量的FW截距，并校正到全局均值kPara）
#   β = Slope_para      （以DTR环境指数为自变量的FW斜率）
#
# 分析内容：
#   1) 读取并清洗表型数据（TKW_mean_table.txt）
#   2) 读取环境因子（EC8.csv）并提取目标环境指数（DTR&mean&PostFlowering1_14）
#   3) 若FW参数文件不存在，则计算并保存7个FW参数（每个基因型一套）
#   4) 读取基因型矩阵（myGD2.csv），MAF过滤 + 缺失插补
#   5) 构建基因组关系矩阵（GRM, G）
#   6) 单性状遗传力估计（BGLR-RKHS）: h2(alpha), h2(beta)
#   7) 二元GBLUP（sommer::mmer）主模型：
#        full：遗传协方差自由估计（unsm(2)）
#        null：遗传协方差固定为0（diag(2)）
#      输出：r_g, Wald CI, LRT(AIC反推), p(χ²1), p(mixture)
#   8) 敏感性分析（每个条件都计算 full vs null 的 LRT & p）：
#        - 标记密度：10K/20K/30K/full
#        - PC校正：0/1/3/5/7 PCs
#   9) 保存所有结果到 output_dir
#
# 备注：
#   - 去除参数自助法（bootstrap），因为：主模型 Wald + LRT + 系统敏感性分析已足够支撑稳健性。
#   - LRT 使用 AIC 等价反推（当 logLik 不可用时的实用方案）：
#       AIC = -2logLik + 2k
#       LRT = 2(logLik_full - logLik_null) = -(AIC_full - AIC_null) + 2*(k_full - k_null)
#     对本检验 k_diff = 1（full多一个遗传协方差参数）。
# ==============================================================================

# ------------------------------ 0. 起始信息 ----------------------------------
cat("========== FW参数计算 + 二元GBLUP分析（无Bootstrap；含LRT+p敏感性）==========\n")
start_time <- Sys.time()
cat("开始时间：", format(start_time), "\n\n")

# ------------------------------ 1. 工作目录 ----------------------------------
# !!! 根据实际情况修改
setwd("/mnt/7t_storage/zhangcl/TKW")
cat("工作目录：", getwd(), "\n\n")

# ------------------------------ 2. 加载依赖包 --------------------------------
required_packages <- c(
  "data.table", "dplyr", "tidyr", "reshape2", "parallel",
  "BGLR", "sommer", "ggplot2", "coda", "Matrix", "stringr"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
cat("所有包加载完成\n\n")

# ------------------------------ 3. 输出目录 ----------------------------------
output_dir <- "S4.1_Bivariate_GBLUP_Results_NoBootstrap"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("输出目录：", output_dir, "\n\n")

# ------------------------------ 4. 分析参数 ----------------------------------
# BGLR参数
bglr_nIter  <- 20000
bglr_burnIn <- 5000
bglr_thin   <- 10

# MAF阈值
maf_threshold <- 0.03

# 环境指数（必须存在于EC8.csv列名中）
target_factor <- "DTR&mean&PostFlowering1_14"

# 敏感性分析梯度
marker_grid <- c(10000, 20000, 30000)   # + full
pc_grid     <- c(0, 1, 3, 5, 7)

cat("参数设置：\n")
cat("  BGLR: nIter=", bglr_nIter, " burnIn=", bglr_burnIn, " thin=", bglr_thin, "\n")
cat("  MAF阈值：", maf_threshold, "\n")
cat("  环境指数：", target_factor, "\n")
cat("  标记密度梯度：", paste(marker_grid, collapse = ", "), " + full\n")
cat("  PC梯度：", paste(pc_grid, collapse = ", "), "\n\n")

# ==============================================================================
# 5. 工具函数
# ==============================================================================

# 5.1 基因型命名清洗：去空格/大写/去非字母数字/数字开头加G前缀
clean_genotype_names <- function(names) {
  names <- trimws(names)
  names <- toupper(names)
  names <- gsub("[^A-Za-z0-9]", "", names)
  names <- gsub("^([0-9]+)", "G\\1", names)
  return(names)
}

# 5.2 在表型文件中自动识别列名（更稳健）
detect_pheno_cols <- function(df) {
  lower_names <- tolower(colnames(df))
  pick_col <- function(candidates) {
    idx <- which(lower_names %in% candidates)
    if (length(idx) >= 1) return(colnames(df)[idx[1]])
    NA_character_
  }

  geno_col <- pick_col(c("genotype","geno","line","variety","cultivar","gid","id"))
  env_col  <- pick_col(c("env_code","env","environment","site","trial","loc","location"))
  tkw_col  <- pick_col(c("tkw","tkw_mean","tkwmean","trait","value","phenotype"))

  if (is.na(geno_col)) {
    idx <- grep("geno|line|var|cult", lower_names)[1]
    if (!is.na(idx)) geno_col <- colnames(df)[idx]
  }
  if (is.na(env_col)) {
    idx <- grep("env|site|loc|trial", lower_names)[1]
    if (!is.na(idx)) env_col <- colnames(df)[idx]
  }
  if (is.na(tkw_col)) {
    idx <- grep("tkw", lower_names)[1]
    if (!is.na(idx)) tkw_col <- colnames(df)[idx]
  }

  if (any(is.na(c(geno_col, env_col, tkw_col)))) {
    stop("无法自动识别 TKW_mean_table.txt 的 genotype/env_code/TKW 列，请检查列名。")
  }
  list(geno_col = geno_col, env_col = env_col, tkw_col = tkw_col)
}

# 5.3 计算FW七参数（每个基因型）
calculate_FW_7params <- function(pheno_data, env_factors, env_means) {
  # 输入：
  #   pheno_data: data.frame(genotype_clean, env_code, TKW)
  #   env_factors: data.frame(env_code, kPara)
  #   env_means: data.frame(env_code, meanY)
  #
  # 输出：
  #   data.frame 每行一个基因型，包含7个FW参数

  ph <- pheno_data %>%
    dplyr::transmute(
      line_code = genotype_clean,
      env_code = as.character(env_code),
      PH = as.numeric(TKW)
    )

  line_codes <- unique(ph$line_code)
  all_envs <- unique(ph$env_code)

  mean_kPara <- mean(env_factors$kPara[env_factors$env_code %in% all_envs], na.rm = TRUE)
  mean_meanY <- mean(env_means$meanY[env_means$env_code %in% all_envs], na.rm = TRUE)

  fw_7params <- data.frame()

  for (line in line_codes) {
    line_data <- subset(ph, line_code == line)
    if (nrow(line_data) < 3) next

    merged <- merge(line_data, env_factors, by = "env_code")
    merged <- merge(merged, env_means, by = "env_code")
    clean <- merged[complete.cases(merged$PH, merged$meanY, merged$kPara), ]
    if (nrow(clean) < 2) next

    # 两条回归：PH~meanY 和 PH~kPara
    tmp <- tryCatch({
      lm_mean <- lm(PH ~ meanY, data = clean)
      Intcp_mean <- predict(lm_mean, data.frame(meanY = mean_meanY))
      Slope_mean <- coef(lm_mean)[2]
      R2_mean <- summary(lm_mean)$r.squared

      lm_para <- lm(PH ~ kPara, data = clean)
      Intcp_para <- coef(lm_para)[1]
      Slope_para <- coef(lm_para)[2]
      Intcp_para_adj <- predict(lm_para, data.frame(kPara = mean_kPara))
      R2_para <- summary(lm_para)$r.squared

      data.frame(
        genotype = line,
        Intcp_mean = as.numeric(Intcp_mean),
        Slope_mean = as.numeric(Slope_mean),
        R2_mean = as.numeric(R2_mean),
        Intcp_para_adj = as.numeric(Intcp_para_adj),
        Intcp_para = as.numeric(Intcp_para),
        Slope_para = as.numeric(Slope_para),
        R2_para = as.numeric(R2_para),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(tmp)) fw_7params <- rbind(fw_7params, tmp)
  }

  fw_7params
}

# 5.4 计算MAF
calc_maf <- function(x) {
  p <- mean(x, na.rm = TRUE) / 2
  if (is.na(p)) return(NA_real_)
  ifelse(p > 0.5, 1 - p, p)
}

# 5.5 缺失插补（列均值）
impute_mean_by_col <- function(M) {
  for (j in 1:ncol(M)) {
    if (anyNA(M[, j])) {
      M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
    }
  }
  M
}

# 5.6 统一函数：拟合二元GBLUP并输出检验统计量
fit_bivar_and_tests <- function(dat_in, G_in, fixed_form) {
  fit_full <- try(
    sommer::mmer(
      fixed = fixed_form,
      random = ~ sommer::vsr(gid, Gu = G_in, Gtc = sommer::unsm(2)),
      rcov   = ~ sommer::vsr(units, Gtc = sommer::unsm(2)),
      data = dat_in,
      verbose = FALSE
    ),
    silent = TRUE
  )

  fit_null <- try(
    sommer::mmer(
      fixed = fixed_form,
      random = ~ sommer::vsr(gid, Gu = G_in, Gtc = diag(2)),
      rcov   = ~ sommer::vsr(units, Gtc = sommer::unsm(2)),
      data = dat_in,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(fit_full, "try-error") || inherits(fit_null, "try-error")) {
    return(list(
      ok = FALSE,
      out = data.frame(
        rg = NA, wald_lower = NA, wald_upper = NA,
        AIC_full = NA, AIC_null = NA, BIC_full = NA, BIC_null = NA,
        LRT = NA, p_chisq = NA, p_mix = NA,
        stringsAsFactors = FALSE
      )
    ))
  }

  vc <- summary(fit_full)$varcomp
  Vg11 <- as.numeric(vc[1, "VarComp"])
  Vg12 <- as.numeric(vc[2, "VarComp"])
  Vg22 <- as.numeric(vc[3, "VarComp"])
  rg_hat <- Vg12 / sqrt(Vg11 * Vg22)

  # Wald CI (Fisher Z)
  n_loc <- nrow(dat_in)
  z_loc <- atanh(rg_hat)
  se_z_loc <- 1 / sqrt(n_loc - 3)
  ci_z_loc <- c(z_loc - 1.96 * se_z_loc, z_loc + 1.96 * se_z_loc)
  ci_rg_loc <- tanh(ci_z_loc)

  # LRT via AIC (k_diff=1)
  k_diff <- 1
  AIC_full <- as.numeric(fit_full$AIC)
  AIC_null <- as.numeric(fit_null$AIC)
  BIC_full <- as.numeric(fit_full$BIC)
  BIC_null <- as.numeric(fit_null$BIC)

  LRT_loc <- -(AIC_full - AIC_null) + 2 * k_diff
  p_chisq_loc <- pchisq(LRT_loc, df = 1, lower.tail = FALSE)
  p_mix_loc <- 0.5 * p_chisq_loc

  list(
    ok = TRUE,
    out = data.frame(
      rg = rg_hat,
      wald_lower = ci_rg_loc[1],
      wald_upper = ci_rg_loc[2],
      AIC_full = AIC_full, AIC_null = AIC_null,
      BIC_full = BIC_full, BIC_null = BIC_null,
      LRT = LRT_loc, p_chisq = p_chisq_loc, p_mix = p_mix_loc,
      stringsAsFactors = FALSE
    ),
    fit_full = fit_full,
    fit_null = fit_null
  )
}

# ==============================================================================
# 6. 第一部分：读取数据 + 计算/读取FW参数
# ==============================================================================

cat("=== 第一部分：读取数据并获取FW参数 ===\n")

# 6.1 读取表型
pheno_file <- "TKW_mean_table.txt"
if (!file.exists(pheno_file)) stop("缺少表型文件：", pheno_file)

ph_raw <- read.table(pheno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cols <- detect_pheno_cols(ph_raw)

pheno <- ph_raw %>%
  dplyr::transmute(
    genotype = trimws(as.character(.data[[cols$geno_col]])),
    env_code  = trimws(as.character(.data[[cols$env_col]])),
    TKW       = as.numeric(.data[[cols$tkw_col]])
  ) %>%
  dplyr::filter(!is.na(genotype), !is.na(env_code), !is.na(TKW))

cat("  表型数据维度：", nrow(pheno), " x ", ncol(pheno), "\n")
cat("  表型列名：", paste(colnames(pheno), collapse = ", "), "\n")

# 6.2 读取环境数据（EC8）
ec8_file <- "EC8.csv"
if (!file.exists(ec8_file)) stop("缺少环境文件：", ec8_file)

ec8 <- read.csv(ec8_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(ec8)[1] <- "env_code"

if (!(target_factor %in% colnames(ec8))) {
  stop("EC8.csv 中找不到目标环境指数列：", target_factor)
}

env_factor <- ec8[, c("env_code", target_factor)]
colnames(env_factor)[2] <- "kPara"
env_factor$env_code <- trimws(as.character(env_factor$env_code))
env_factor$kPara <- as.numeric(env_factor$kPara)

cat("  环境数据维度：", nrow(ec8), " x ", ncol(ec8), "\n")
cat("  使用环境指数：", target_factor, "\n")

# 6.3 计算各环境表型均值（meanY）
env_mean <- pheno %>%
  dplyr::group_by(env_code) %>%
  dplyr::summarise(meanY = mean(TKW, na.rm = TRUE), .groups = "drop")

cat("  环境均值表维度：", nrow(env_mean), " x ", ncol(env_mean), "\n\n")

# 6.4 读取基因型文件（用于匹配品种名，也用于后续建GRM）
geno_file <- "myGD2.csv"
if (!file.exists(geno_file)) stop("缺少基因型文件：", geno_file)

geno_dt0 <- data.table::fread(geno_file, data.table = FALSE)
if (ncol(geno_dt0) < 2) stop("myGD2.csv 列数不足，至少应包含样本列 + SNP列。")

geno_gids_raw <- geno_dt0[[1]]
geno_gids <- clean_genotype_names(geno_gids_raw)

cat("  基因型原始维度：", nrow(geno_dt0), " x ", ncol(geno_dt0), "\n")
cat("  基因型样本数：", length(geno_gids), "\n\n")

# 6.5 清洗表型品种名并匹配到基因型
pheno$genotype_clean <- clean_genotype_names(pheno$genotype)
common_lines_phenoGeno <- intersect(unique(pheno$genotype_clean), unique(geno_gids))
if (length(common_lines_phenoGeno) == 0) stop("表型与基因型无共同品种，终止。")

pheno2 <- pheno %>%
  dplyr::filter(genotype_clean %in% common_lines_phenoGeno)

cat("  表型与基因型共同品种数：", length(common_lines_phenoGeno), "\n")
cat("  过滤后表型记录数：", nrow(pheno2), "\n\n")

# 6.6 FW参数文件：若存在直接读；否则计算并保存
fw_dir <- "FW_7Params_Results"
fw_file <- file.path(fw_dir, "TKW_FW_7Core_Parameters.csv")

if (file.exists(fw_file)) {
  fw_params <- read.csv(fw_file, stringsAsFactors = FALSE)
  cat("  从文件读取FW参数：", fw_file, "\n")
  cat("  FW参数维度：", nrow(fw_params), " x ", ncol(fw_params), "\n\n")
} else {
  cat("  未发现FW参数文件，开始计算7参数...\n")

  # 只保留 env_factor 和 env_mean 有定义的环境
  common_envs <- intersect(intersect(unique(pheno2$env_code), unique(env_factor$env_code)), unique(env_mean$env_code))
  pheno_fw <- pheno2 %>% dplyr::filter(env_code %in% common_envs)

  env_factor_fw <- env_factor %>% dplyr::filter(env_code %in% common_envs)
  env_mean_fw   <- env_mean   %>% dplyr::filter(env_code %in% common_envs)

  cat("  FW计算使用环境数：", length(common_envs), "\n")
  cat("  FW计算使用表型记录数：", nrow(pheno_fw), "\n")

  fw_params <- calculate_FW_7params(pheno_fw, env_factor_fw, env_mean_fw)
  cat("  成功计算FW参数的品种数：", nrow(fw_params), "\n")

  if (!dir.exists(fw_dir)) dir.create(fw_dir, recursive = TRUE)
  write.csv(fw_params, fw_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("  FW参数已保存：", fw_file, "\n\n")
}

cat("  FW参数列名：", paste(colnames(fw_params), collapse = ", "), "\n\n")

# ==============================================================================
# 7. 第二部分：提取α/β并构建GRM
# ==============================================================================

cat("=== 第二部分：提取α/β并构建GRM ===\n")

# 7.1 提取 α / β
need_cols <- c("genotype", "Intcp_para_adj", "Slope_para")
if (!all(need_cols %in% colnames(fw_params))) {
  stop("FW参数文件缺少必要列：", paste(setdiff(need_cols, colnames(fw_params)), collapse = ", "))
}

fw_params$alpha <- as.numeric(fw_params$Intcp_para_adj)
fw_params$beta  <- as.numeric(fw_params$Slope_para)
fw_params$gid_clean <- clean_genotype_names(fw_params$genotype)

alpha_values <- fw_params$alpha
names(alpha_values) <- fw_params$gid_clean
beta_values <- fw_params$beta
names(beta_values) <- fw_params$gid_clean

alpha_values <- alpha_values[is.finite(alpha_values)]
beta_values  <- beta_values[is.finite(beta_values)]

# 7.2 构建基因型矩阵
geno_matrix <- as.matrix(geno_dt0[, -1])
rownames(geno_matrix) <- geno_gids

# 7.3 三方共同品种
common_gids <- intersect(intersect(names(alpha_values), names(beta_values)), rownames(geno_matrix))
cat("  α/β/基因型共同品种数：", length(common_gids), "\n")
if (length(common_gids) < 10) stop("共同品种太少（<10），终止。")

alpha_common <- alpha_values[common_gids]
beta_common  <- beta_values[common_gids]
geno_common  <- geno_matrix[common_gids, , drop = FALSE]

cat("  共同品种示例：", paste(head(common_gids, 5), collapse = ", "), "\n")
cat("  基因型矩阵维度（共同品种）：", nrow(geno_common), " x ", ncol(geno_common), "\n\n")

# 7.4 MAF过滤
maf <- apply(geno_common, 2, calc_maf)
keep_maf <- which(is.finite(maf) & maf >= maf_threshold)

geno_qc <- geno_common[, keep_maf, drop = FALSE]
cat("  MAF过滤后标记数：", ncol(geno_qc), "\n")
if (ncol(geno_qc) < 1000) cat("  警告：过滤后标记数较少，请确认maf_threshold是否过高。\n")

# 7.5 缺失插补
geno_qc <- impute_mean_by_col(geno_qc)

# 7.6 构建GRM
G <- sommer::A.mat(geno_qc)
G <- (G + t(G)) / 2 + diag(1e-6, nrow(G))
rownames(G) <- common_gids
colnames(G) <- common_gids
cat("  GRM维度：", nrow(G), " x ", ncol(G), "\n\n")

# ==============================================================================
# 8. 第三部分：单性状遗传力估计（BGLR-RKHS）
# ==============================================================================

cat("=== 第三部分：单性状遗传力估计（BGLR-RKHS）===\n")

y_alpha <- scale(alpha_common, center = TRUE, scale = FALSE)[, 1]
y_beta  <- scale(beta_common,  center = TRUE, scale = FALSE)[, 1]

ETA <- list(geno = list(K = G, model = "RKHS"))

set.seed(123)
fit_alpha <- BGLR::BGLR(
  y = y_alpha,
  ETA = ETA,
  nIter = bglr_nIter,
  burnIn = bglr_burnIn,
  thin = bglr_thin,
  verbose = FALSE
)
h2_alpha <- fit_alpha$ETA$geno$varU / (fit_alpha$ETA$geno$varU + fit_alpha$varE)

set.seed(456)
fit_beta <- BGLR::BGLR(
  y = y_beta,
  ETA = ETA,
  nIter = bglr_nIter,
  burnIn = bglr_burnIn,
  thin = bglr_thin,
  verbose = FALSE
)
h2_beta <- fit_beta$ETA$geno$varU / (fit_beta$ETA$geno$varU + fit_beta$varE)

cat("  h2(alpha) =", round(h2_alpha, 4), "\n")
cat("  h2(beta)  =", round(h2_beta,  4), "\n\n")

# ==============================================================================
# 9. 第四部分：二元GBLUP主模型（full vs null）
# ==============================================================================

cat("=== 第四部分：二元GBLUP主模型拟合（full vs null）===\n")

dat <- data.frame(
  gid = factor(common_gids, levels = common_gids),
  alpha = as.numeric(alpha_common),
  beta  = as.numeric(beta_common),
  stringsAsFactors = FALSE
)

fixed_form_main <- cbind(alpha, beta) ~ 1

main_fit <- fit_bivar_and_tests(dat_in = dat, G_in = G, fixed_form = fixed_form_main)
if (!main_fit$ok) stop("主模型拟合失败，请检查数据/GRM/是否数值不稳定。")

main_out <- main_fit$out
rg_hat <- main_out$rg
ci_rg  <- c(main_out$wald_lower, main_out$wald_upper)
LRT_main <- main_out$LRT
p_chisq_main <- main_out$p_chisq
p_mix_main <- main_out$p_mix

cat("  遗传相关 r_g =", round(rg_hat, 4), "\n")
cat("  Wald 95% CI: [", round(ci_rg[1], 4), ", ", round(ci_rg[2], 4), "]\n")
cat("  LRT(AIC反推) =", round(LRT_main, 4), "\n")
cat("  p(χ²1) =", format(p_chisq_main, scientific = TRUE, digits = 3), "\n")
cat("  p(mixture) =", format(p_mix_main, scientific = TRUE, digits = 3), "\n")
cat("  AIC full =", round(main_out$AIC_full, 3), " | null =", round(main_out$AIC_null, 3), "\n")
cat("  BIC full =", round(main_out$BIC_full, 3), " | null =", round(main_out$BIC_null, 3), "\n\n")

# ==============================================================================
# 10. 第五部分：敏感性分析（每个条件输出 rg + WaldCI + LRT + p）
# ==============================================================================

cat("=== 第五部分：敏感性分析（标记密度 & PC梯度；均输出LRT与p）===\n")

# ------------------------------
# 10.1 标记密度梯度：10K / 20K / 30K / full
# ------------------------------
cat("\n--- 10.1 标记密度梯度（10K/20K/30K/full）---\n")

set.seed(2026)
n_markers_all <- ncol(geno_qc)

pick_k <- function(k, n) {
  if (k >= n) return(1:n)
  sample(n, k)
}

idx_full <- 1:n_markers_all
idx_10k  <- pick_k(10000, n_markers_all)
idx_20k  <- pick_k(20000, n_markers_all)
idx_30k  <- pick_k(30000, n_markers_all)

density_list <- list(
  full = idx_full,
  k30  = idx_30k,
  k20  = idx_20k,
  k10  = idx_10k
)

sensitivity_marker_LRT <- data.frame(
  density = names(density_list),
  n_markers = sapply(density_list, length),
  rg = NA, wald_lower = NA, wald_upper = NA,
  LRT = NA, p_chisq = NA, p_mix = NA,
  AIC_full = NA, AIC_null = NA, BIC_full = NA, BIC_null = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(density_list)) {
  lab <- names(density_list)[i]
  idx <- density_list[[i]]

  M_sub <- geno_qc[, idx, drop = FALSE]
  G_sub <- sommer::A.mat(M_sub)
  G_sub <- (G_sub + t(G_sub)) / 2 + diag(1e-6, nrow(G_sub))
  rownames(G_sub) <- common_gids
  colnames(G_sub) <- common_gids

  res <- fit_bivar_and_tests(dat_in = dat, G_in = G_sub, fixed_form = fixed_form_main)

  sensitivity_marker_LRT[i, c("rg","wald_lower","wald_upper",
                              "LRT","p_chisq","p_mix",
                              "AIC_full","AIC_null","BIC_full","BIC_null")] <- res$out[1, c(
                                "rg","wald_lower","wald_upper",
                                "LRT","p_chisq","p_mix",
                                "AIC_full","AIC_null","BIC_full","BIC_null"
                              )]

  cat("  ", lab,
      " markers=", sensitivity_marker_LRT$n_markers[i],
      " rg=", round(sensitivity_marker_LRT$rg[i], 4),
      " LRT=", round(sensitivity_marker_LRT$LRT[i], 3),
      " p_mix=", format(sensitivity_marker_LRT$p_mix[i], scientific = TRUE, digits = 3), "\n")
}

print(sensitivity_marker_LRT)

# ------------------------------
# 10.2 PC梯度：0 / 1 / 3 / 5 / 7 PCs（固定效应加入PC）
# ------------------------------
cat("\n--- 10.2 PC梯度（0/1/3/5/7 PCs）---\n")

eig <- eigen(G)
# 至少取到 max(pc_grid) 个PC
max_pc_need <- max(pc_grid)
PCs <- eig$vectors[, 1:max(10, max_pc_need), drop = FALSE]
colnames(PCs) <- paste0("PC", 1:ncol(PCs))

pc_df <- data.frame(gid = common_gids, PCs, stringsAsFactors = FALSE)
dat_pc <- merge(dat, pc_df, by = "gid")

sensitivity_pc_LRT <- data.frame(
  nPC = pc_grid,
  rg = NA, wald_lower = NA, wald_upper = NA,
  LRT = NA, p_chisq = NA, p_mix = NA,
  AIC_full = NA, AIC_null = NA, BIC_full = NA, BIC_null = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(pc_grid)) {
  npc <- pc_grid[i]
  if (npc == 0) {
    fixed_form <- cbind(alpha, beta) ~ 1
  } else {
    pc_terms <- paste(paste0("PC", 1:npc), collapse = " + ")
    fixed_form <- as.formula(paste("cbind(alpha, beta) ~ 1 +", pc_terms))
  }

  res <- fit_bivar_and_tests(dat_in = dat_pc, G_in = G, fixed_form = fixed_form)

  sensitivity_pc_LRT[i, c("rg","wald_lower","wald_upper",
                          "LRT","p_chisq","p_mix",
                          "AIC_full","AIC_null","BIC_full","BIC_null")] <- res$out[1, c(
                            "rg","wald_lower","wald_upper",
                            "LRT","p_chisq","p_mix",
                            "AIC_full","AIC_null","BIC_full","BIC_null"
                          )]

  cat("  PCs=", npc,
      " rg=", round(sensitivity_pc_LRT$rg[i], 4),
      " LRT=", round(sensitivity_pc_LRT$LRT[i], 3),
      " p_mix=", format(sensitivity_pc_LRT$p_mix[i], scientific = TRUE, digits = 3), "\n")
}

print(sensitivity_pc_LRT)

# ==============================================================================
# 11. 第六部分：保存结果
# ==============================================================================

cat("\n=== 第六部分：保存结果 ===\n")

main_results <- data.frame(
  metric = c(
    "h2_alpha", "h2_beta",
    "rg_hat", "wald_lower", "wald_upper",
    "LRT_main", "p_chisq_main", "p_mix_main",
    "AIC_full", "AIC_null", "BIC_full", "BIC_null",
    "n_common_gids", "n_markers_after_maf"
  ),
  value = c(
    h2_alpha, h2_beta,
    rg_hat, ci_rg[1], ci_rg[2],
    LRT_main, p_chisq_main, p_mix_main,
    main_out$AIC_full, main_out$AIC_null, main_out$BIC_full, main_out$BIC_null,
    length(common_gids), ncol(geno_qc)
  ),
  stringsAsFactors = FALSE
)

write.csv(main_results, file.path(output_dir, "main_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(sensitivity_marker_LRT, file.path(output_dir, "sensitivity_marker_LRT.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(sensitivity_pc_LRT, file.path(output_dir, "sensitivity_pc_LRT.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 保存关键对象（可复现）
saveRDS(list(
  fw_params = fw_params,
  common_gids = common_gids,
  alpha_common = alpha_common,
  beta_common = beta_common,
  maf_threshold = maf_threshold,
  geno_qc_dim = dim(geno_qc),
  G = G,
  main_results = main_results,
  sensitivity_marker_LRT = sensitivity_marker_LRT,
  sensitivity_pc_LRT = sensitivity_pc_LRT
), file.path(output_dir, "analysis_objects.rds"))

cat("  结果已保存至：", output_dir, "\n")

# ==============================================================================
# 12. 结束
# ==============================================================================
end_time <- Sys.time()
cat("\n========== 全部完成 ==========\n")
cat("结束时间：", format(end_time), "\n")
cat("总耗时：", format(difftime(end_time, start_time, units = "mins")), "分钟\n")

# ==============================================================================
# 将分析结果格式化为论文 Table 2 发表版格式
# ==============================================================================

cat("========== 生成 Table 2 发表版格式 ==========\n")

# 设置结果目录（与主分析脚本一致）
output_dir <- "S4.1_Bivariate_GBLUP_Results_NoBootstrap"

# ------------------------------
# 1. 读取结果文件
# ------------------------------

main_results <- read.csv(file.path(output_dir, "main_results.csv"),
                         stringsAsFactors = FALSE)

sens_marker <- read.csv(file.path(output_dir, "sensitivity_marker_LRT.csv"),
                        stringsAsFactors = FALSE)

sens_pc <- read.csv(file.path(output_dir, "sensitivity_pc_LRT.csv"),
                    stringsAsFactors = FALSE)

# ------------------------------
# 2. 主模型行
# ------------------------------

rg_main <- main_results$value[main_results$metric == "rg_hat"]
wald_lower_main <- main_results$value[main_results$metric == "wald_lower"]
wald_upper_main <- main_results$value[main_results$metric == "wald_upper"]
LRT_main <- main_results$value[main_results$metric == "LRT_main"]
p_mix_main <- main_results$value[main_results$metric == "p_mix_main"]

table_main <- data.frame(
  Section = "Main model (no correction)",
  Condition = "",
  rg_hat = round(rg_main, 3),
  Wald_CI = paste0("[", round(wald_lower_main, 3), ", ",
                    round(wald_upper_main, 3), "]"),
  LRT = round(LRT_main, 3),
  p_mixture = format(p_mix_main, scientific = TRUE, digits = 3),
  stringsAsFactors = FALSE
)

# ------------------------------
# 3. 标记密度部分
# ------------------------------

marker_section <- sens_marker

# 美化密度标签
marker_section$Condition <- ifelse(
  marker_section$density == "full", "Full",
  ifelse(marker_section$density == "k30", "30K",
  ifelse(marker_section$density == "k20", "20K",
  ifelse(marker_section$density == "k10", "10K",
         marker_section$density)))
)

table_marker <- data.frame(
  Section = "Marker density",
  Condition = marker_section$Condition,
  rg_hat = round(marker_section$rg, 3),
  Wald_CI = paste0("[",
                   round(marker_section$wald_lower, 3), ", ",
                   round(marker_section$wald_upper, 3), "]"),
  LRT = round(marker_section$LRT, 3),
  p_mixture = format(marker_section$p_mix, scientific = TRUE, digits = 3),
  stringsAsFactors = FALSE
)

# ------------------------------
# 4. PC校正部分
# ------------------------------

pc_section <- sens_pc

table_pc <- data.frame(
  Section = "Principal component correction",
  Condition = paste0(pc_section$nPC, " PC"),
  rg_hat = round(pc_section$rg, 3),
  Wald_CI = paste0("[",
                   round(pc_section$wald_lower, 3), ", ",
                   round(pc_section$wald_upper, 3), "]"),
  LRT = round(pc_section$LRT, 3),
  p_mixture = format(pc_section$p_mix, scientific = TRUE, digits = 3),
  stringsAsFactors = FALSE
)

# ------------------------------
# 5. 合并成最终Table 2
# ------------------------------

Table2 <- rbind(
  table_main,
  table_marker,
  table_pc
)

# ------------------------------
# 6. 保存发表版表格
# ------------------------------

write.csv(Table2,
          file.path(output_dir, "Table2_GeneticCorrelation_Formatted.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")

cat("Table 2 已生成：",
    file.path(output_dir, "Table2_GeneticCorrelation_Formatted.csv"), "\n")


