# ==============================================================
# MLM (Mixed Linear Model) GWAS with crossed random intercepts + PCs fixed
# Phenotype: genotype + env_code + trait (multiple env records per genotype)
# Genotype matrix: rows=genotype, cols=SNP (0/1/2), first col=ID
# MAP file columns: SNP, Chromosome, Position
#
# Goal:
#   1) Compute genotype PCs once (from GRM)
#   2) Use LD-pruned SNPs (pure R) to compute lambda for each k (PC number)
#   3) Choose best k by score = |lambda-1| + penalty*k
#   4) Final full scan on ALL SNPs using Score test (variance components fixed from null MLM)
#
# Notes:
#   - LD pruning uses the Prune() function from SFSI repo + the LD_prune() you provided
#   - No PLINK needed
# ==============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(lme4)
})

# -----------------------
# 0) Settings (edit as needed)
# -----------------------
setwd("/mnt/7t_storage/zhangcl/TKW")

trait <- "TKW"
pheno_file <- "TKW_mean_table.txt"   # must contain genotype/env_code/trait (names can be mapped)
geno_file  <- "myGD2.csv"            # first col = genotype ID, remaining columns = SNPs (0/1/2)
geno_rds   <- "GENO_optimized.rds"   # cache for genotype matrix

# MAP columns are exactly: SNP,Chromosome,Position
map_file <- "myGD2.csv"          # <- your map file (CSV/TSV both ok with fread)

out_dir <- file.path("MLM_PCA_lambda_LDprune_FULL")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# PC screening range: k = 0..max_pcs_to_test (will be clipped by rank)
max_pcs_to_test <- 2

# SNP QC on genotype-level matrix
maf_min <- 0.01
miss_max <- 0.05
impute_missing <- TRUE

# LD prune params (for lambda screening only)
ld_r2_threshold <- 0.2       # recommended 0.1~0.2; 0.95 is too loose
window_size <- 2000          # SNPs per window (adjust)
mc_cores_ld <- 1L            # Windows use 1; Linux can use >1

# Score scan block size
block_size <- 3000

# PC selection criterion
penalty_k <- 0.001           # score = |lambda-1| + penalty*k ; set 0 for no penalty

# lmer control (stability)
LME_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# -----------------------
# 1) Load Prune.R (SFSI) and define LD_prune (your exact function)
# -----------------------
source("https://raw.githubusercontent.com/MarcooLopez/SFSI/main/R/Prune.R")

LD_prune <- function(X, MAP, threshold=0.95, n.windows=NULL, window.size=NULL,
                     d.max=NULL, mc.cores=1L, verbose=TRUE)
{
    colnames(MAP) <- toupper(colnames(MAP))
    CHRs <- sort(unique(MAP$CHR))

    if(is.null(colnames(X))){
      colnames(X) <- paste0("S",MAP$CHR,"_",MAP$POS)
    }
    index <- NULL
    for(j in 1:ncol(MAP)){
      if(all(MAP[,j]==colnames(X))) index <- j
    }
    if(is.null(index)){
      MAP <- data.frame(NAME=colnames(X),MAP)
    }else{
      colnames(MAP)[index] <- "NAME"
    }

    if(!is.null(n.windows) & !is.null(window.size)){
      window.size <- NULL
      warning("Input 'window.size' is ignored when also 'n.windows' is provided",
              .immediate=TRUE)
    }

    if("MAF" %in% colnames(MAP)){
      message("SNP with the highest MAF will be selected from each group sharing an R2>threshold")
    }else{
      message("First-appeared SNP will be selected from each group sharing an R2>threshold")
    }

    if(any(is.na(X))){
      message("Correlation is obtained from complete pairwise observations")
    }

    for(chr in CHRs){
      index <- which(MAP$CHR==chr)
      tmp <- sprintf('%.2f',diff(range(MAP$POS[index]))/1E6)
      cat(" - Chromosome ",chr,": size=",tmp," Mb. ",length(index)," SNPs\n",sep="")
    }

    if(!is.null(d.max)) d.max <- as.integer(d.max)

    compApply <- function(i)
    {
      chr <- CHRs[i]
      indexCHR <- which(MAP$CHR==chr)

      if(is.null(n.windows) & is.null(window.size)){
        n.windows0 <- 1L
        window.size0 <- length(indexCHR)
      }else{
        if(is.null(n.windows)){
          n.windows0 <- ceiling(length(indexCHR)/window.size)
          window.size0 <- table(rep(1:n.windows0, each=window.size)[1:length(indexCHR)])

          if(window.size0[n.windows0] < window.size*0.30){ # If the last window is very small
            window.size0[n.windows0-1] <- window.size0[n.windows0-1] + window.size0[n.windows0]
            window.size0 <- window.size0[1:(n.windows0-1)]
            n.windows0 <- n.windows0 - 1
          }

        }else{
          if(is.null(window.size)){
            n.windows0 <- n.windows
            window.size0 <- table(rep(1:n.windows,ceiling(length(indexCHR)/n.windows))[1:length(indexCHR)])
          }
        }
      }
      stopifnot(length(window.size0) == n.windows0)
      stopifnot(sum(window.size0) == length(indexCHR))
      window.size.cum <- c(0, cumsum(window.size0))

      id0 <- c()
      out0  <- vector('list', length(threshold))
      nConn0 <- vector('list', length(threshold))
      #names(out0) <- names(nConn0) <- threshold
      for(k in 1:n.windows0)
      {
        tmp <- c(window.size.cum[k] + 1, window.size.cum[k+1])
        index0 <- indexCHR[seq(tmp[1], tmp[2])]
        if(any(is.na(index0))){ stop("Marker matrix could not be subset at window ",k) }

        if(any(is.na(X[,index0]))){
          R0 <- cor(X[,index0], use="pairwise.complete.obs")^2
        }else{
          R0 <- (crossprod(scale(X[,index0]))/(nrow(X) - 1))^2
        }

        if(!is.null(d.max)){
          # D0 <- as.matrix(dist(MAP[index0,"POS",drop=F], method='manhattan'))
          D0 <- sapply(MAP$POS[index0],function(x)abs(x-MAP$POS[index0]))
        }else {
          D0 <- NULL
        }

        if("MAF" %in% colnames(MAP)){
          MAF0 <- MAP$MAF[index0]
        }else{
          MAF0 <- NULL
        }
        for(tr in seq_along(threshold)){
          ans <- Prune(R0, threshold=threshold[tr], D=D0, d.max=d.max, MAF=MAF0, verbose=FALSE)
          ##tt <- match(ans$prune.in,colnames(R))
          ##A <- (R[tt,tt] > threshold) & (D[tt,tt] <= d.max); diag(A) <- FALSE
          ##range(colSums(A))
          out0[[tr]] <- rbind(out0[[tr]], data.frame(CHR=chr,window=k,NAME=ans$prune.in))
          nConn0[[tr]] <- c(nConn0[[tr]], ans$nConn)

          message("    CHR=",sprintf('%2d',chr),". Window ",ifelse(n.windows0==1,NA,k)," [",
                  sprintf('%3d',tmp[1]),",",tmp[2],"]: ",sprintf('%5d',length(ans$prune.in)),
                  " SNPs selected with R2=",threshold[tr])
        }
        id0 <- c(id0, MAP[index0,"NAME"])
      }
      for(tr in seq_along(threshold)){
        out0[[tr]] <- out0[[tr]][order(match(out0[[tr]]$NAME, MAP$NAME[indexCHR])),]
        rownames(out0[[tr]]) <- NULL
        if(n.windows0 == 1L){
          out0[[tr]]$window <- NA
        }
      }

      return(list(pruneIn=out0, id=id0, nConn=nConn0))
    }

    if(mc.cores == 1L){
       res = lapply(X=seq_along(CHRs),FUN=compApply)
    }else{
       res = parallel::mclapply(X=seq_along(CHRs),FUN=compApply,mc.cores=mc.cores)
    }
    if(any(unlist(lapply(res, function(x)x$id)) != MAP$NAME)){
      stop("Something went wrong during the prunning procedure failed")
    }

    pruneIn <- lapply(seq_along(threshold), function(tr){
      do.call(rbind, lapply(res, function(x)x$pruneIn[[tr]]))
    })

    nConn <- lapply(seq_along(threshold), function(tr){
      unlist(lapply(res, function(x)x$nConn[[tr]]))
    })
    names(pruneIn) <- names(nConn) <- threshold

    tmp <- unlist(lapply(pruneIn, nrow))
    for(i in 1:length(tmp)){
     message(" - ",sprintf('%5d',tmp[i])," of ",ncol(X)," SNPs selected across the ",
             length(CHRs)," CHRs for R2=",names(tmp)[i])
    }

    return(list(pruneIn=pruneIn, nConn=nConn))
}

# -----------------------
# 2) Helper functions (lambda + score)
# -----------------------
calc_lambda <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  chisq <- qchisq(1 - p, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

score_block <- function(W, solveV, X, ViX, XtViX_inv, Py) {
  ViW <- solveV(W)                            # nObs x m
  U   <- as.numeric(crossprod(W, Py))         # m
  XtViW <- crossprod(X, ViW)                  # p x m
  wViW <- as.numeric(colSums(W * ViW))        # diag(W'ViW)
  mid <- XtViX_inv %*% XtViW                  # p x m
  quad <- as.numeric(colSums(XtViW * mid))    # diag(XtViW' XtViX_inv XtViW)
  I <- wViW - quad
  I[!is.finite(I) | I <= 0] <- NA_real_
  chisq <- (U * U) / I
  chisq[!is.finite(chisq)] <- NA_real_
  p <- pchisq(chisq, df = 1, lower.tail = FALSE)
  list(chisq = chisq, p = p)
}

# -----------------------
# 3) Load phenotype (genotype/env_code/trait)
# -----------------------
cat("=== Load phenotype ===\n")
pheno <- fread(pheno_file)
setDT(pheno)

# map genotype column
if (!"genotype" %in% names(pheno)) {
  if ("Genotype" %in% names(pheno)) setnames(pheno, "Genotype", "genotype")
  else if ("Line" %in% names(pheno)) setnames(pheno, "Line", "genotype")
  else if ("ID" %in% names(pheno)) setnames(pheno, "ID", "genotype")
  else stop("表型中未找到 genotype 列（或可映射列）")
}
# map env_code column
if (!"env_code" %in% names(pheno)) {
  if ("Env" %in% names(pheno)) setnames(pheno, "Env", "env_code")
  else if ("Environment" %in% names(pheno)) setnames(pheno, "Environment", "env_code")
  else if ("Site" %in% names(pheno)) setnames(pheno, "Site", "env_code")
  else stop("表型中未找到 env_code 列（或可映射列）")
}
# map trait column
if (!trait %in% names(pheno)) {
  if (paste0(trait, "_mean") %in% names(pheno)) setnames(pheno, paste0(trait, "_mean"), trait)
  else stop(paste0("表型中未找到性状列: ", trait))
}

pheno <- pheno[!is.na(get(trait))]
pheno[, genotype := factor(as.character(genotype))]
pheno[, env_code := factor(as.character(env_code))]

cat("pheno n(obs) =", nrow(pheno), "\n")
cat("genotypes =", nlevels(pheno$genotype), " envs =", nlevels(pheno$env_code), "\n")

# -----------------------
# 4) Load genotype matrix (rows=genotype)
# -----------------------
cat("\n=== Load genotype matrix ===\n")
if (file.exists(geno_rds)) {
  cat("Read RDS:", geno_rds, "\n")
  geno <- readRDS(geno_rds)
} else {
  cat("Read CSV:", geno_file, "\n")
  gd <- fread(geno_file)
  ids <- gd[[1]]
  geno <- as.matrix(gd[, -1, with = FALSE])
  rownames(geno) <- ids
  saveRDS(geno, geno_rds)
  cat("Saved RDS:", geno_rds, "\n")
}
storage.mode(geno) <- "double"
cat("geno dim =", paste(dim(geno), collapse=" x "), "\n")

# align genotype IDs between phenotype and genotype matrix
common <- intersect(levels(pheno$genotype), rownames(geno))
if (length(common) < 5) stop("共同基因型太少，请检查 phenotype 的 genotype 与基因型矩阵行名是否一致")

pheno <- pheno[genotype %in% common]
pheno[, genotype := factor(as.character(genotype), levels = common)]
geno <- geno[common, , drop = FALSE]
cat("aligned genotypes =", nrow(geno), " SNPs =", ncol(geno), "\n")

# -----------------------
# 5) SNP QC on genotype-level matrix (0/1/2): missing + MAF + impute
# -----------------------
cat("\n=== SNP QC (missing + MAF) ===\n")
miss_rate <- colMeans(!is.finite(geno))

af <- rep(NA_real_, ncol(geno))
for (j in seq_len(ncol(geno))) {
  x <- geno[, j]
  ok <- is.finite(x)
  if (sum(ok) > 0) af[j] <- mean(x[ok]) / 2
}
maf <- pmin(af, 1 - af)

keep <- which(miss_rate <= miss_max & is.finite(maf) & maf >= maf_min)
cat("keep SNPs:", length(keep), "/", ncol(geno), "\n")
if (length(keep) < 100) stop("QC后SNP过少，放宽 maf_min 或 miss_max")

geno <- geno[, keep, drop = FALSE]
maf  <- maf[keep]
miss_rate <- miss_rate[keep]
snp_names <- colnames(geno)

if (impute_missing) {
  na_cnt <- sum(!is.finite(geno))
  if (na_cnt > 0) {
    cat("Impute missing by SNP mean, NA count =", na_cnt, "\n")
    for (j in seq_len(ncol(geno))) {
      x <- geno[, j]
      bad <- !is.finite(x)
      if (any(bad)) {
        mu <- mean(x[!bad])
        x[bad] <- mu
        geno[, j] <- x
      }
    }
  }
}

cat("post-QC geno dim =", paste(dim(geno), collapse=" x "), "\n")

# save QC summary
qc_dt <- data.table(
  snp = colnames(geno),
  maf = maf,
  miss_rate = miss_rate
)
fwrite(qc_dt, file.path(out_dir, paste0("SNP_QC_", trait, ".csv")))

# -----------------------
# 6) Load MAP (SNP,Chromosome,Position) and align to geno columns
# -----------------------
cat("\n=== Load MAP (SNP,Chromosome,Position) ===\n")
MAP0 <- fread(map_file)
setDT(MAP0)

# Normalize column names exactly as user said
# Expect: SNP, Chromosome, Position (case-insensitive)
nms <- names(MAP0)
nms_upper <- toupper(nms)
names(MAP0) <- nms_upper

if (!all(c("SNP", "CHROMOSOME", "POSITION") %in% names(MAP0))) {
  stop("MAP 文件必须包含三列：SNP, Chromosome, Position（大小写不敏感）")
}

MAP <- MAP0[, .(
  NAME = as.character(SNP),
  CHR  = as.integer(Chromosome),
  POS  = as.integer(Position)
)]

# Align MAP rows to geno SNP order using SNP name matching
idxm <- match(colnames(geno), MAP$NAME)
if (any(is.na(idxm))) {
  missing_snps <- colnames(geno)[is.na(idxm)]
  stop(paste0(
    "MAP$SNP 与基因型矩阵列名不一致，以下 SNP 在 MAP 中找不到（展示前20个）：\n",
    paste(head(missing_snps, 20), collapse = ", ")
  ))
}
MAP <- MAP[idxm]
MAP[, MAF := maf]   # add MAF for tie-break in pruning

cat("MAP aligned rows =", nrow(MAP), "\n")

# -----------------------
# 7) LD prune SNPs (pure R) for lambda screening
# -----------------------
cat("\n=== LD pruning for lambda screening (pure R) ===\n")
pr <- LD_prune(
  X = geno,
  MAP = MAP[, .(CHR, POS, NAME, MAF)],
  threshold = ld_r2_threshold,
  window.size = window_size,
  mc.cores = mc_cores_ld,
  verbose = TRUE
)

# Only one threshold used -> take first element
pruned_names <- pr$pruneIn[[1]]$NAME
snp_idx_pruned <- match(pruned_names, colnames(geno))
snp_idx_pruned <- snp_idx_pruned[!is.na(snp_idx_pruned)]

cat("LD-pruned SNPs:", length(snp_idx_pruned), "/", ncol(geno), "\n")
fwrite(
  data.table(snp = colnames(geno)[snp_idx_pruned]),
  file.path(out_dir, paste0("LD_pruned_SNPs_r2_", ld_r2_threshold, ".csv"))
)

if (length(snp_idx_pruned) < 500) {
  warning("LD-pruned SNP数量过少，lambda筛k可能不稳定。可放宽阈值(0.2->0.3)或增大window_size。")
}

# -----------------------
# 8) Compute genotype PCs once (GRM eigen)
# -----------------------
cat("\n=== Compute genotype PCs (once) ===\n")
max_k <- min(max_pcs_to_test, nrow(geno) - 2)

Xc <- scale(geno, center = TRUE, scale = FALSE)
Xc[!is.finite(Xc)] <- 0

GRM <- tcrossprod(Xc) / ncol(Xc)   # nG x nG
eig <- eigen(GRM, symmetric = TRUE)

vals <- eig$values
vecs <- eig$vectors
idx <- which(vals > 1e-8)
vals <- vals[idx]
vecs <- vecs[, idx, drop = FALSE]
max_k <- min(max_k, ncol(vecs))

PC_full <- vecs[, 1:max_k, drop = FALSE]
PC_full <- sweep(PC_full, 2, sqrt(vals[1:max_k]), `*`)
colnames(PC_full) <- paste0("GPC", 1:max_k)
rownames(PC_full) <- rownames(geno)

cat("PC computed:", ncol(PC_full), "\n")

# Save PCs
pc_out <- as.data.table(PC_full)
pc_out[, genotype := rownames(PC_full)]
setcolorder(pc_out, c("genotype", paste0("GPC", 1:ncol(PC_full))))
fwrite(pc_out, file.path(out_dir, paste0("Genotype_PCs_", trait, ".csv")))

# -----------------------
# 9) Fit null MLM + score scan (core function)
# -----------------------
fit_null_and_scan <- function(k, snp_idx_use) {
  # 9.1 observation-level PCs
  if (k == 0) {
    pc_obs <- NULL
  } else {
    pc_obs <- PC_full[as.character(pheno$genotype), 1:k, drop = FALSE]
    colnames(pc_obs) <- paste0("GPC", 1:k)
  }

  dat <- data.table(
    y = pheno[[trait]],
    genotype = pheno$genotype,
    env_code = pheno$env_code
  )
  if (k > 0) dat <- cbind(dat, as.data.table(pc_obs))

  fixed_part <- if (k == 0) "1" else paste(colnames(pc_obs), collapse = " + ")
  fml <- as.formula(paste0("y ~ ", fixed_part, " + (1|genotype) + (1|env_code)"))

  # 9.2 fit null MLM once
  fm0 <- lmer(fml, data = dat, REML = FALSE, control = LME_CTRL)

  # 9.3 variance components
  vc <- VarCorr(fm0)
  s2g <- as.numeric(vc[["genotype"]])^2
  s2e <- as.numeric(vc[["env_code"]])^2
  s2r <- sigma(fm0)^2

  # 9.4 build V (sparse)
  Zg <- sparse.model.matrix(~0 + genotype, data = dat)
  Ze <- sparse.model.matrix(~0 + env_code, data = dat)

  nObs <- nrow(dat)
  V <- s2r * Diagonal(nObs) + s2g * (Zg %*% t(Zg)) + s2e * (Ze %*% t(Ze))

  cholV <- Cholesky(V, LDL = FALSE, super = TRUE)
  solveV <- function(B) as.matrix(solve(cholV, B))

  # 9.5 fixed design matrix X and y from lme4
  X <- getME(fm0, "X")
  y <- getME(fm0, "y")

  # 9.6 precompute for P = V^-1 - V^-1 X (X'V^-1X)^-1 X' V^-1
  ViX <- solveV(X)
  XtViX <- crossprod(X, ViX)
  XtViX_inv <- solve(XtViX)

  Viy <- solveV(y)
  beta_hat <- XtViX_inv %*% crossprod(X, Viy)
  Py <- Viy - ViX %*% beta_hat  # P*y

  # 9.7 prepare SNP matrix: genotype-level -> obs-level expansion
  g_idx <- match(as.character(dat$genotype), rownames(geno))
  if (any(is.na(g_idx))) stop("Internal error: genotype index NA during SNP expansion")

  G_use <- geno[, snp_idx_use, drop = FALSE]
  # center SNP at genotype-level mean
  G_use <- scale(G_use, center = TRUE, scale = FALSE)

  m <- ncol(G_use)
  pvals <- rep(NA_real_, m)
  chisq <- rep(NA_real_, m)

  idx0 <- 1
  while (idx0 <= m) {
    j2 <- min(m, idx0 + block_size - 1)
    Wg_blk <- G_use[, idx0:j2, drop = FALSE]      # nG x mblk
    W_obs  <- Wg_blk[g_idx, , drop = FALSE]       # nObs x mblk

    sb <- score_block(W_obs, solveV, X, ViX, XtViX_inv, Py)
    pvals[idx0:j2] <- sb$p
    chisq[idx0:j2] <- sb$chisq

    idx0 <- j2 + 1
  }

  lam <- calc_lambda(pvals)

  list(
    k = k,
    lambda = lam,
    lambda_diff = abs(lam - 1),
    s2g = s2g, s2e = s2e, s2r = s2r,
    pvals = pvals,
    chisq = chisq
  )
}

# -----------------------
# 10) Screen k by lambda using LD-pruned SNPs
# -----------------------
cat("\n=== Screen k by lambda using LD-pruned SNPs ===\n")
screen_dt <- data.table(k = 0:max_k, lambda = NA_real_, lambda_diff = NA_real_, score = NA_real_)

for (k in 0:max_k) {
  cat(sprintf(" -> k=%d\n", k))
  rr <- fit_null_and_scan(k, snp_idx_pruned)
  screen_dt[k == rr$k, `:=`(
    lambda = rr$lambda,
    lambda_diff = rr$lambda_diff,
    score = rr$lambda_diff + penalty_k * k
  )]
}

setorder(screen_dt, score)
best_k <- screen_dt$k[1]
best_lambda_pruned <- screen_dt$lambda[1]

cat("\nBEST k =", best_k, " (lambda on pruned SNPs =", round(best_lambda_pruned, 4), ")\n")

screen_file <- file.path(out_dir, paste0("SCREEN_lambda_LDprune_", trait, ".csv"))
fwrite(screen_dt, screen_file)
cat("Saved:", screen_file, "\n")

# -----------------------
# 11) Final full scan with best_k on ALL SNPs
# -----------------------
cat("\n=== Final full scan with best_k on ALL SNPs ===\n")
final_res <- fit_null_and_scan(best_k, seq_len(ncol(geno)))

final_lambda <- final_res$lambda
cat("Final lambda (ALL SNPs) =", round(final_lambda, 4), "\n")

out <- data.table(
  trait = trait,
  snp = colnames(geno),
  chisq = final_res$chisq,
  p = final_res$pvals,
  neg_log10_p = -log10(final_res$pvals),
  k = best_k,
  lambda = final_lambda
)
setorder(out, p)

final_file <- file.path(out_dir, paste0("GWAS_MLM_SCORE_", trait, "_bestk", best_k, ".csv"))
fwrite(out, final_file)
cat("Saved:", final_file, "\n")

# -----------------------
# 12) Save report
# -----------------------
rep_file <- file.path(out_dir, paste0("REPORT_", trait, ".txt"))
sink(rep_file)

cat("MLM GWAS Score Test Report (FULL)\n")
cat("================================\n")
cat("Trait:", trait, "\n")
cat("Phenotype file:", pheno_file, "\n")
cat("Genotype file:", geno_file, "\n")
cat("Map file:", map_file, "\n\n")

cat("Data summary:\n")
cat("  Obs n:", nrow(pheno), "\n")
cat("  Genotypes:", nlevels(pheno$genotype), "\n")
cat("  Envs:", nlevels(pheno$env_code), "\n")
cat("  SNPs after QC:", ncol(geno), "\n\n")

cat("Model:\n")
cat("  y ~ GPC1..GPCk + (1|genotype) + (1|env_code)\n\n")

cat("SNP QC:\n")
cat("  maf_min:", maf_min, "\n")
cat("  miss_max:", miss_max, "\n")
cat("  impute_missing:", impute_missing, "\n\n")

cat("LD pruning (lambda screening):\n")
cat("  r2 threshold:", ld_r2_threshold, "\n")
cat("  window_size:", window_size, "\n")
cat("  mc_cores_ld:", mc_cores_ld, "\n")
cat("  pruned SNP count:", length(snp_idx_pruned), "\n\n")

cat("PC screening:\n")
cat("  max_pcs_to_test:", max_pcs_to_test, "\n")
cat("  effective max_k:", max_k, "\n")
cat("  penalty_k:", penalty_k, "\n")
cat("  best_k:", best_k, "\n")
cat("  lambda(pruned) at best_k:", round(best_lambda_pruned, 4), "\n")
cat("  lambda(ALL) at best_k:", round(final_lambda, 4), "\n\n")

cat("Variance components (null @ best k):\n")
cat("  sigma_g^2:", final_res$s2g, "\n")
cat("  sigma_e^2:", final_res$s2e, "\n")
cat("  sigma_r^2:", final_res$s2r, "\n\n")

cat("Top 20 hits:\n")
print(head(out[, .(snp, chisq, p, neg_log10_p)], 20))

sink()

cat("Saved:", rep_file, "\n")
cat("\nDONE.\n")
