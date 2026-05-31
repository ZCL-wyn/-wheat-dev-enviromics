# ==============================================================
# C-Path GRM-K MLM GWAS (HAND-WRITTEN) with:
#  - REML variance estimation via numerical optimization (nlminb)
#  - Score-test GWAS (block scan)
#  - Lambda-based PC-set screening (ARBITRARY PC combinations allowed)
#  - Print per PC-set: lambda, logLik(REML), AIC, BIC, nll
#  - Supports:
#      * RAW mode: phenotype has genotype+env_code+trait, env as random (if >1 env)
#      * BLUE mode: compute BLUE across env (env fixed), then GWAS (no env random)
#
# Inputs:
#  phenotype file: contains genotype, env_code, trait
#  genotype file: rows=genotype, cols=SNP (0/1/2), first col=ID
#  map file: SNP,Chromosome,Position (optional)
# ==============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# -----------------------
# 0) Settings
# -----------------------
setwd("/mnt/7t_storage/zhangcl/TKW")

trait <- "TKW"
pheno_file <- "TKW_mean_table.txt"
geno_file  <- "myGD2.csv"
geno_rds   <- "GENO_optimized.rds"
map_file   <- "myGM2.csv"   # optional check

out_dir <- file.path("GWAS_GRM_K_REML_SCORE_ULTIMATE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# choose analysis mode: "RAW" or "BLUE" or "BOTH"
run_mode <- "BOTH"  # <- "RAW" / "BLUE" / "BOTH"

# PC range
max_pcs_to_use <- 20

# PC set screening controls
max_pc_set_size <- 10
max_pc_sets_to_test <- 200
penalty_k <- 0.001

# SNP QC
do_qc <- FALSE
maf_min <- 0.01
miss_max <- 0.05
impute_missing <- TRUE

# score scan block size
block_size <- 2000

# optimization controls
opt_maxit <- 100
opt_rel_tol <- 1e-6

# numerical stability
jitter <- 1e-8

set.seed(1)

# -----------------------
# 1) Utilities
# -----------------------
calc_lambda <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  chisq <- qchisq(1 - p, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

chol_safe <- function(V) {
  for (i in 0:6) {
    eps <- jitter * (10^i)
    V2 <- V
    diag(V2) <- diag(V2) + eps
    ok <- TRUE
    R <- tryCatch(chol(V2), error = function(e) { ok <<- FALSE; NULL })
    if (ok) return(list(R=R, eps=eps))
  }
  stop("Cholesky failed even after adding jitter. V might be indefinite.")
}

solve_chol <- function(R, B) {
  backsolve(R, forwardsolve(t(R), B))
}

# REML negative log-likelihood
nll_reml <- function(theta, y, X, AKAt, Ee, use_env, n, p) {
  s2g <- exp(theta[1])
  if (use_env) {
    s2e <- exp(theta[2])
    s2r <- exp(theta[3])
    V <- s2g * AKAt + s2e * Ee
    diag(V) <- diag(V) + s2r
  } else {
    s2r <- exp(theta[2])
    V <- s2g * AKAt
    diag(V) <- diag(V) + s2r
  }

  ch <- chol_safe(V)
  R <- ch$R
  logdetV <- 2 * sum(log(diag(R)))

  ViX <- solve_chol(R, X)
  XtViX <- crossprod(X, ViX)

  ch2 <- tryCatch(chol(XtViX), error = function(e) NULL)
  if (is.null(ch2)) return(1e100)
  logdetXtViX <- 2 * sum(log(diag(ch2)))

  Viy <- solve_chol(R, y)
  beta_hat <- solve(XtViX, crossprod(X, Viy))
  Py <- Viy - ViX %*% beta_hat
  yPy <- as.numeric(crossprod(y, Py))

  nll <- 0.5 * (logdetV + logdetXtViX + yPy + (n - p) * log(2*pi))
  if (!is.finite(nll)) nll <- 1e100
  nll
}

score_block <- function(W, solveV, X, ViX, XtViX_inv, Py) {
  ViW <- solveV(W)
  U <- as.numeric(crossprod(W, Py))
  XtViW <- crossprod(X, ViW)

  wViW <- as.numeric(colSums(W * ViW))
  mid <- XtViX_inv %*% XtViW
  quad <- as.numeric(colSums(XtViW * mid))

  I <- wViW - quad
  I[!is.finite(I) | I <= 0] <- NA_real_

  chisq <- (U * U) / I
  chisq[!is.finite(chisq)] <- NA_real_
  p <- pchisq(chisq, df = 1, lower.tail = FALSE)
  list(chisq=chisq, p=p)
}

# ----------------------------------------
# 2) PC-set generator
# ----------------------------------------
generate_pc_sets <- function(max_pc, max_set_size=10, max_sets=200) {
  sets <- list(integer(0))

  for (i in 1:max_pc) sets[[length(sets)+1]] <- i
  for (k in 2:min(max_set_size, max_pc)) sets[[length(sets)+1]] <- 1:k

  if (max_pc >= 8) {
    patterns <- list(c(2,5,8), c(1,3,5), c(2,4,6,8), c(3,6,9), c(5,10,15), c(2,5,10,15))
    for (s in patterns) if (all(s <= max_pc)) sets[[length(sets)+1]] <- s
  }

  key <- function(v) paste(v, collapse=",")
  seen <- unique(sapply(sets, key))

  while (length(sets) < max_sets) {
    kmax <- min(max_set_size, max_pc)
    k <- sample(1:kmax, 1, prob = rev(seq_len(kmax)))
    s <- sort(sample.int(max_pc, k))
    ks <- key(s)
    if (!ks %in% seen) {
      sets[[length(sets)+1]] <- s
      seen <- c(seen, ks)
    }
  }
  sets
}

# ----------------------------------------
# 3) BLUE calculator
# ----------------------------------------
compute_BLUE <- function(pheno_dt, trait_col) {
  if (nlevels(pheno_dt$env_code) <= 1) {
    blue <- pheno_dt[, .(BLUE = mean(get(trait_col), na.rm=TRUE)), by=.(genotype)]
    return(blue)
  }

  df <- as.data.frame(pheno_dt[, .(y=get(trait_col), genotype=genotype, env_code=env_code)])
  fit <- lm(y ~ genotype + env_code, data=df)

  env_levels <- levels(pheno_dt$env_code)
  geno_levels <- levels(pheno_dt$genotype)

  w_env <- prop.table(table(pheno_dt$env_code))
  w_env <- w_env[env_levels]

  grid <- expand.grid(genotype=geno_levels, env_code=env_levels, stringsAsFactors=FALSE)
  pred <- predict(fit, newdata=grid)

  pred_mat <- matrix(pred, nrow=length(geno_levels), byrow=TRUE,
                     dimnames=list(geno_levels, env_levels))

  blue_val <- as.numeric(pred_mat %*% as.numeric(w_env))
  data.table(genotype=factor(geno_levels, levels=geno_levels), BLUE=blue_val)
}

# ----------------------------------------
# helper: print metrics nicely
# ----------------------------------------
print_fit_metrics <- function(pc_set_str, rr, nObs, p, use_env) {
  # rr must include: lambda, logLik, AIC, BIC, nll, s2g,s2e,s2r, conv
  cat(sprintf(
    "    [PCset={%s}] nObs=%d p=%d use_env=%s | lambda=%.4f | logLik(REML)=%.3f | AIC=%.3f | BIC=%.3f | nll=%.3f | conv=%d | s2g=%.6f s2e=%.6f s2r=%.6f\n",
    pc_set_str, nObs, p, ifelse(use_env, "TRUE", "FALSE"),
    rr$lambda, rr$logLik, rr$AIC, rr$BIC, rr$nll, rr$conv,
    rr$s2g, ifelse(is.na(rr$s2e), NaN, rr$s2e), rr$s2r
  ))
}

# -----------------------
# 4) Load phenotype
# -----------------------
cat("=== Load phenotype ===\n")
pheno <- fread(pheno_file)
setDT(pheno)

if (!"genotype" %in% names(pheno)) {
  if ("Genotype" %in% names(pheno)) setnames(pheno, "Genotype", "genotype")
  else if ("Line" %in% names(pheno)) setnames(pheno, "Line", "genotype")
  else if ("ID" %in% names(pheno)) setnames(pheno, "ID", "genotype")
  else stop("表型中未找到 genotype 列（或可映射列）")
}
if (!"env_code" %in% names(pheno)) {
  if ("Env" %in% names(pheno)) setnames(pheno, "Env", "env_code")
  else if ("Environment" %in% names(pheno)) setnames(pheno, "Environment", "env_code")
  else if ("Site" %in% names(pheno)) setnames(pheno, "Site", "env_code")
  else stop("表型中未找到 env_code 列（或可映射列）")
}
if (!trait %in% names(pheno)) {
  if (paste0(trait, "_mean") %in% names(pheno)) setnames(pheno, paste0(trait, "_mean"), trait)
  else stop(paste0("表型中未找到性状列: ", trait))
}

pheno <- pheno[!is.na(get(trait))]
pheno[, genotype := factor(as.character(genotype))]
pheno[, env_code := factor(as.character(env_code))]

cat("Obs n =", nrow(pheno), "\n")
cat("Genotypes =", nlevels(pheno$genotype), " Envs =", nlevels(pheno$env_code), "\n")

# -----------------------
# 5) Load genotype matrix
# -----------------------
cat("\n=== Load genotype matrix ===\n")
if (file.exists(geno_rds)) {
  cat("Read RDS:", geno_rds, "\n")
  geno <- readRDS(geno_rds)
} else {
  cat("Read CSV:", geno_file, "\n")
  gd <- fread(geno_file)
  ids <- gd[[1]]
  geno <- as.matrix(gd[, -1, with=FALSE])
  rownames(geno) <- ids
  saveRDS(geno, geno_rds)
  cat("Saved RDS:", geno_rds, "\n")
}
storage.mode(geno) <- "double"
cat("geno dim =", paste(dim(geno), collapse=" x "), "\n")

# align IDs
common <- intersect(levels(pheno$genotype), rownames(geno))
if (length(common) < 5) stop("共同基因型太少，请检查 genotype ID 是否一致")
pheno <- pheno[genotype %in% common]
pheno[, genotype := factor(as.character(genotype), levels=common)]
geno <- geno[common, , drop=FALSE]
cat("Aligned genotypes =", nrow(geno), " SNPs =", ncol(geno), "\n")

# optional MAP check
if (file.exists(map_file)) {
  cat("\n=== MAP check (optional) ===\n")
  MAP0 <- fread(map_file)
  setDT(MAP0)
  names(MAP0) <- toupper(names(MAP0))
  if (all(c("SNP","CHROMOSOME","POSITION") %in% names(MAP0))) {
    miss_in_map <- setdiff(colnames(geno), as.character(MAP0$SNP))
    if (length(miss_in_map) > 0) {
      cat("Warning: some geno SNP not found in MAP. (show 10):\n")
      cat(paste(head(miss_in_map, 10), collapse=", "), "\n")
    } else {
      cat("MAP covers all geno SNP names.\n")
    }
  }
}

# -----------------------
# 6) Optional QC + impute
# -----------------------
if (do_qc) {
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
  geno <- geno[, keep, drop=FALSE]
  cat("post-QC SNPs =", ncol(geno), "\n")
}

if (impute_missing) {
  na_cnt <- sum(!is.finite(geno))
  if (na_cnt > 0) {
    cat("\nImpute missing by SNP mean, NA count =", na_cnt, "\n")
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

# -----------------------
# 7) Build GRM K and PCs (max 20)
# -----------------------
cat("\n=== Build GRM K and PCs (max 20) ===\n")
nG <- nrow(geno)
mS <- ncol(geno)

Xc <- scale(geno, center=TRUE, scale=FALSE)
Xc[!is.finite(Xc)] <- 0
K <- tcrossprod(Xc) / mS
K <- (K + t(K)) / 2

eig <- eigen(K, symmetric=TRUE)
vals <- eig$values
vecs <- eig$vectors

idx_pos <- which(vals > 1e-10)
if (length(idx_pos) < 2) stop("K 的正特征值太少，无法做 PCs。")

vals <- vals[idx_pos]
vecs <- vecs[, idx_pos, drop=FALSE]

max_pc <- min(max_pcs_to_use, ncol(vecs))
cat("PC available:", ncol(vecs), " use max_pc =", max_pc, "\n")

PC_full <- vecs[, 1:max_pc, drop=FALSE]
PC_full <- sweep(PC_full, 2, sqrt(vals[1:max_pc]), `*`)
colnames(PC_full) <- paste0("GPC", 1:max_pc)
rownames(PC_full) <- rownames(geno)

fwrite(
  cbind(data.table(genotype=rownames(PC_full)), as.data.table(PC_full)),
  file.path(out_dir, paste0("Genotype_PCs_fromK_", trait, ".csv"))
)

# -----------------------
# 8) Core runner for one dataset (RAW or BLUE)
# -----------------------
run_one <- function(mode_label, pheno_dt) {
  cat("\n================================================\n")
  cat("RUN MODE:", mode_label, "\n")
  cat("================================================\n")

  use_env <- (mode_label == "RAW") && (nlevels(pheno_dt$env_code) > 1)

  nObs <- nrow(pheno_dt)
  y <- as.matrix(pheno_dt[[trait]])
  storage.mode(y) <- "double"

  pheno_dt[, genotype := factor(as.character(genotype), levels=rownames(geno))]
  if (any(is.na(pheno_dt$genotype))) stop("Some phenotype genotypes not found in geno rownames after alignment.")

  A <- sparse.model.matrix(~0 + genotype, data=pheno_dt)
  AK <- as.matrix(A %*% K)
  AKAt <- AK %*% t(as.matrix(A))

  if (use_env) {
    Ze <- sparse.model.matrix(~0 + env_code, data=pheno_dt)
    Ee <- as.matrix(Ze %*% t(Ze))
  } else {
    Ee <- NULL
  }

  pc_sets <- generate_pc_sets(max_pc=max_pc, max_set_size=max_pc_set_size, max_sets=max_pc_sets_to_test)

  fit_set_and_scan <- function(pc_set, snp_idx_use) {
    # X = intercept + selected PCs expanded to obs
    if (length(pc_set) == 0) {
      X <- matrix(1, nObs, 1)
      colnames(X) <- "(Intercept)"
    } else {
      pc_obs <- PC_full[as.character(pheno_dt$genotype), pc_set, drop=FALSE]
      X <- cbind(1, pc_obs)
      colnames(X)[1] <- "(Intercept)"
      colnames(X)[-1] <- paste0("GPC", pc_set)
    }
    p_fix <- ncol(X)

    # Fit REML
    yvar <- var(y)
    if (use_env) {
      theta0 <- log(c(yvar*0.3, yvar*0.3, yvar*0.4) + 1e-6)
      obj <- function(th) nll_reml(th, y, X, AKAt, Ee, TRUE, nObs, p_fix)
      opt <- nlminb(start=theta0, objective=obj,
                    control=list(iter.max=opt_maxit, eval.max=opt_maxit*5, rel.tol=opt_rel_tol))
      s2g <- exp(opt$par[1]); s2e <- exp(opt$par[2]); s2r <- exp(opt$par[3])
      V <- s2g * AKAt + s2e * Ee
      diag(V) <- diag(V) + s2r
      k_param <- p_fix + 3  # fixed effects + variance params
    } else {
      theta0 <- log(c(yvar*0.6, yvar*0.4) + 1e-6)
      obj <- function(th) nll_reml(th, y, X, AKAt, NULL, FALSE, nObs, p_fix)
      opt <- nlminb(start=theta0, objective=obj,
                    control=list(iter.max=opt_maxit, eval.max=opt_maxit*5, rel.tol=opt_rel_tol))
      s2g <- exp(opt$par[1]); s2e <- NA_real_; s2r <- exp(opt$par[2])
      V <- s2g * AKAt
      diag(V) <- diag(V) + s2r
      k_param <- p_fix + 2
    }

    ch <- chol_safe(V)
    R <- ch$R
    solveV <- function(B) solve_chol(R, B)

    # Projection pieces
    ViX <- solveV(X)
    XtViX <- crossprod(X, ViX)
    XtViX_inv <- tryCatch(solve(XtViX), error=function(e) NULL)
    if (is.null(XtViX_inv)) {
      return(list(ok=FALSE,
                  lambda=NA_real_, lambda_diff=NA_real_,
                  logLik=NA_real_, AIC=NA_real_, BIC=NA_real_, nll=as.numeric(opt$objective),
                  pvals=NULL, chisq=NULL,
                  s2g=s2g, s2e=s2e, s2r=s2r, conv=as.integer(opt$convergence),
                  opt=opt, p_fix=p_fix))
    }

    Viy <- solveV(y)
    beta_hat <- XtViX_inv %*% crossprod(X, Viy)
    Py <- Viy - ViX %*% beta_hat

    # Score scan
    G_use <- geno[, snp_idx_use, drop=FALSE]
    G_use <- scale(G_use, center=TRUE, scale=FALSE)
    g_idx <- match(as.character(pheno_dt$genotype), rownames(geno))

    m <- ncol(G_use)
    pvals <- rep(NA_real_, m)
    chisq <- rep(NA_real_, m)

    idx0 <- 1
    while (idx0 <= m) {
      j2 <- min(m, idx0 + block_size - 1)
      Wg_blk <- G_use[, idx0:j2, drop=FALSE]
      W_obs <- Wg_blk[g_idx, , drop=FALSE]
      sb <- score_block(W_obs, solveV, X, ViX, XtViX_inv, Py)
      pvals[idx0:j2] <- sb$p
      chisq[idx0:j2] <- sb$chisq
      idx0 <- j2 + 1
    }

    lam <- calc_lambda(pvals)

    # logLik/AIC/BIC based on REML objective
    nll <- as.numeric(opt$objective)
    logLik <- -nll
    AIC <- 2 * k_param - 2 * logLik
    BIC <- log(nObs) * k_param - 2 * logLik

    list(ok=TRUE,
         lambda=lam, lambda_diff=abs(lam-1),
         logLik=logLik, AIC=AIC, BIC=BIC, nll=nll,
         pvals=pvals, chisq=chisq,
         s2g=s2g, s2e=s2e, s2r=s2r,
         conv=as.integer(opt$convergence),
         opt=opt, p_fix=p_fix)
  }

  snp_idx_all <- seq_len(ncol(geno))

  screen <- data.table(
    pc_set = character(0),
    n_pc = integer(0),
    lambda = numeric(0),
    lambda_diff = numeric(0),
    score = numeric(0),
    logLik = numeric(0),
    AIC = numeric(0),
    BIC = numeric(0),
    nll = numeric(0),
    s2g = numeric(0),
    s2e = numeric(0),
    s2r = numeric(0),
    conv = integer(0),
    p_fix = integer(0)
  )

  cat("PC sets to test =", length(pc_sets), "\n")

  for (i in seq_along(pc_sets)) {
    ps <- pc_sets[[i]]
    ps_str <- if (length(ps)==0) "" else paste(ps, collapse=",")
    cat(sprintf(" -> PC set %d/%d : {%s}\n", i, length(pc_sets), ps_str))

    rr <- fit_set_and_scan(ps, snp_idx_all)

    if (isTRUE(rr$ok) && is.finite(rr$lambda)) {
      screen <- rbind(
        screen,
        data.table(
          pc_set = ps_str,
          n_pc = length(ps),
          lambda = rr$lambda,
          lambda_diff = rr$lambda_diff,
          score = rr$lambda_diff + penalty_k * length(ps),
          logLik = rr$logLik,
          AIC = rr$AIC,
          BIC = rr$BIC,
          nll = rr$nll,
          s2g = rr$s2g,
          s2e = rr$s2e,
          s2r = rr$s2r,
          conv = rr$conv,
          p_fix = rr$p_fix
        )
      )

      # <<< 关键：每次都打印 >>>
      print_fit_metrics(pc_set_str = ps_str, rr = rr, nObs = nObs, p = rr$p_fix, use_env = use_env)

    } else {
      cat("    [FAILED] XtViX singular or optimization failed.\n")
    }
  }

  if (nrow(screen) == 0) stop("All PC-set fits failed; check data / increase jitter / reduce PC-set size.")

  setorder(screen, score)
  best <- screen[1]
  best_set <- if (best$pc_set == "") integer(0) else as.integer(strsplit(best$pc_set, ",")[[1]])

  cat("\nBEST PC set:",
      ifelse(best$pc_set=="","{ }", paste0("{", best$pc_set, "}")),
      sprintf("| score=%.6f lambda=%.4f AIC=%.2f BIC=%.2f logLik=%.2f\n", best$score, best$lambda, best$AIC, best$BIC, best$logLik)
  )

  # Save screen
  screen_file <- file.path(out_dir, paste0("SCREEN_PCsets_metrics_", trait, "_", mode_label, ".csv"))
  fwrite(screen, screen_file)

  # Final full scan with best set
  cat("\n=== Final scan with best PC set ===\n")
  final <- fit_set_and_scan(best_set, snp_idx_all)

  out <- data.table(
    trait = trait,
    snp = colnames(geno),
    chisq = final$chisq,
    p = final$pvals,
    neg_log10_p = -log10(final$pvals),
    pc_set = best$pc_set,
    n_pc = length(best_set),
    lambda = final$lambda,
    mode = mode_label,
    logLik = final$logLik,
    AIC = final$AIC,
    BIC = final$BIC,
    nll = final$nll,
    s2g = final$s2g,
    s2e = final$s2e,
    s2r = final$s2r,
    conv = final$conv
  )
  setorder(out, p)

  out_file <- file.path(out_dir, paste0("GWAS_GRM_K_SCORE_", trait, "_", mode_label, "_PCset_", gsub(",", "_", best$pc_set), ".csv"))
  fwrite(out, out_file)

  # Report
  rep_file <- file.path(out_dir, paste0("REPORT_", trait, "_", mode_label, ".txt"))
  sink(rep_file)
  cat("GWAS GRM-K MLM (C-Path) + REML(nlminb) + Score test\n")
  cat("===================================================\n")
  cat("MODE:", mode_label, "\n\n")
  cat("Trait:", trait, "\n")
  cat("Obs n:", nObs, "\n")
  cat("Genotypes:", nlevels(pheno_dt$genotype), "\n")
  cat("Envs:", nlevels(pheno_dt$env_code), "\n")
  cat("SNPs:", ncol(geno), "\n\n")

  cat("Random components:\n")
  cat("  Genetic: s2g * A K A'\n")
  if (use_env) cat("  Env:     s2e * Ze Ze'\n") else cat("  Env:     (not used)\n")
  cat("  Resid:   s2r * I\n\n")

  cat("Best PC set:", ifelse(best$pc_set=="","(none)", best$pc_set), "\n")
  cat("nPC:", length(best_set), "\n")
  cat(sprintf("lambda: %.4f | logLik(REML): %.3f | AIC: %.3f | BIC: %.3f | nll: %.3f | conv: %d\n\n",
              final$lambda, final$logLik, final$AIC, final$BIC, final$nll, final$conv))

  cat("Variance components:\n")
  cat("  s2g:", final$s2g, "\n")
  cat("  s2e:", final$s2e, "\n")
  cat("  s2r:", final$s2r, "\n\n")

  cat("Top 20 hits:\n")
  print(head(out[, .(snp, chisq, p, neg_log10_p)], 20))
  sink()

  cat("Saved screen:", screen_file, "\n")
  cat("Saved GWAS:", out_file, "\n")
  cat("Saved report:", rep_file, "\n")

  list(screen=screen, out=out, best_set=best_set)
}

# -----------------------
# 9) Run modes
# -----------------------
results <- list()

if (run_mode %in% c("RAW","BOTH")) {
  results$RAW <- run_one("RAW", copy(pheno))
}

if (run_mode %in% c("BLUE","BOTH")) {
  cat("\n=== Compute BLUE phenotype ===\n")
  blue_dt <- compute_BLUE(copy(pheno), trait)

  ph_blue <- data.table(
    genotype = factor(as.character(blue_dt$genotype), levels=rownames(geno)),
    env_code = factor(rep("BLUE", nrow(blue_dt))),
    TKW = blue_dt$BLUE
  )
  results$BLUE <- run_one("BLUE", ph_blue)
}

cat("\nALL DONE.\n")
