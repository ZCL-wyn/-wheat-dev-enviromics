#!/usr/bin/env Rscript

# =========================================================
#  QTL additive effects + additive-by-additive epistasis
#  with background control: PCs(from G) + random G (sommer)
#
#  Model for each trait (Intercept / Slope):
#    y = b0 + b1*g1 + b2*g2 + b12*(g1*g2) + PC covariates + u(line) + e
#    u ~ N(0, G*sigma_g), e ~ N(0, I*sigma_e)
#
#  g1, g2 are dosage 0/1/2 (KEEP heterozygotes)
#
#  Simplified version with bug fixes
# =========================================================

cat("========== QTL Additive Epistasis (G-matrix + PCA + sommer) ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# -------------------- 0) Paths --------------------
base_dir <- "/mnt/7t_storage/zhangcl/TKW"
vcf_file <- file.path(base_dir, "extracted_G1_to_G339.vcf")
qtl_file <- file.path(base_dir, "QTL4.csv")
phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")
environment_file <- file.path(base_dir, "ECs_results.csv")

output_dir <- file.path(base_dir, "QTL_Additive_Epistasis_Background_Fixed")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create subdirectories
subdirs <- c("Model_Results", "Figures", "Significant_Pairs", "Diagnostics")
for (dir_name in subdirs) {
  dir_path <- file.path(output_dir, dir_name)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

# -------------------- 1) Packages --------------------
required_packages <- c(
  "vcfR", "data.table", "dplyr", "tidyr", "Matrix", "sommer", 
  "parallel", "ggplot2", "ggpubr", "ggrepel", "cowplot"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(vcfR)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(sommer)
  library(parallel)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(cowplot)
})

# -------------------- 2) Helper Functions --------------------
stop_if_missing <- function(fp) {
  if (!file.exists(fp)) stop(sprintf("Error: file not found: %s", fp))
}

safe_read_table <- function(fp) {
  ext <- tools::file_ext(fp)
  if (ext == "csv") {
    return(read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE))
  } else {
    return(read.table(fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
  }
}

# Convert GT string to dosage 0/1/2 (KEEP heterozygotes)
gt_to_dosage <- function(gt){
  if (is.na(gt) || gt %in% c(".", "./.", ".|.", "NA", "")) return(NA_real_)
  a <- unlist(strsplit(gt, "[/|]"))
  if (length(a) != 2) return(NA_real_)
  if (any(a == ".")) return(NA_real_)
  a1 <- suppressWarnings(as.numeric(a[1]))
  a2 <- suppressWarnings(as.numeric(a[2]))
  if (!is.finite(a1) || !is.finite(a2)) return(NA_real_)
  return(a1 + a2)
}

# Match phenotype line codes to VCF sample IDs
match_samples <- function(line_codes, vcf_samples) {
  m1 <- rep(NA_character_, length(line_codes))
  
  # 1. Exact match
  idx_exact <- line_codes %in% vcf_samples
  m1[idx_exact] <- line_codes[idx_exact]
  
  # 2. Case-insensitive exact match
  still <- is.na(m1)
  if (any(still)) {
    for (i in which(still)) {
      x <- line_codes[i]
      hit <- grep(paste0("^", x, "$"), vcf_samples, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) m1[i] <- hit[1]
    }
  }
  
  # 3. Substring match (only if still unmatched)
  still <- is.na(m1)
  if (any(still)) {
    for (i in which(still)) {
      x <- line_codes[i]
      hit <- grep(x, vcf_samples, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) m1[i] <- hit[1]
    }
  }
  
  keep <- !is.na(m1)
  list(
    line_codes = line_codes[keep],
    matched_samples = m1[keep]
  )
}

# Compute line-specific intercept & slope
extract_line_effects <- function(pheno_env_data, min_points = 3) {
  results <- list()
  
  for (ln in unique(pheno_env_data$line_code)) {
    dd <- pheno_env_data[pheno_env_data$line_code == ln, ]
    if (nrow(dd) < min_points) next
    
    fit <- try(lm(PH ~ env_factor_value, data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    coefs <- coef(fit)
    if (length(coefs) < 2) next
    
    sm <- summary(fit)
    results[[length(results) + 1]] <- data.frame(
      line_code = ln,
      Intercept = unname(coefs[1]),
      Slope = unname(coefs[2]),
      Intercept_SE = sm$coefficients[1, 2],
      Slope_SE = sm$coefficients[2, 2],
      Intercept_p = sm$coefficients[1, 4],
      Slope_p = sm$coefficients[2, 4],
      R_squared = sm$r.squared,
      n_points = nrow(dd),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(results) == 0) return(data.frame())
  bind_rows(results)
}

# VanRaden G matrix
calc_G_vanraden <- function(Geno) {
  p <- colMeans(Geno, na.rm = TRUE) / 2
  keep <- p > 0 & p < 1
  Geno <- Geno[, keep, drop = FALSE]
  p <- p[keep]
  
  if (ncol(Geno) < 10) stop("Too few polymorphic SNPs to build G")
  
  # Impute missing to 2p
  for (j in seq_len(ncol(Geno))) {
    idx <- is.na(Geno[, j])
    if (any(idx)) Geno[idx, j] <- 2 * p[j]
  }
  
  W <- sweep(Geno, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  G <- tcrossprod(W) / denom
  G <- (G + t(G)) / 2
  rownames(G) <- rownames(Geno)
  colnames(G) <- rownames(Geno)
  G
}

# PCA from G
pca_from_G <- function(G, n_pcs = 5) {
  eig <- eigen(G, symmetric = TRUE)
  pos <- eig$values > 1e-8
  vals <- eig$values[pos]
  vecs <- eig$vectors[, pos, drop = FALSE]
  n_pcs <- min(n_pcs, ncol(vecs))
  scores <- vecs[, 1:n_pcs, drop = FALSE] %*% diag(sqrt(vals[1:n_pcs]), n_pcs, n_pcs)
  colnames(scores) <- paste0("PC", 1:n_pcs)
  rownames(scores) <- rownames(G)
  var_exp <- vals / sum(vals) * 100
  list(scores = scores, eigenvalues = vals, var_explained = var_exp)
}

# Get fixed effects from sommer mmer with better error handling
get_fixed_effects_from_mmer <- function(fit) {
  beta <- as.matrix(fit$Beta)
  
  # Get standard errors
  se <- rep(NA_real_, nrow(beta))
  if (!is.null(fit$VarBeta)) {
    se <- sqrt(pmax(diag(fit$VarBeta), 0))
  }
  names(se) <- rownames(beta)
  
  # Try to get p-values from Wald test
  pvals <- rep(NA_real_, nrow(beta))
  names(pvals) <- rownames(beta)
  
  wald_result <- try(sommer::Wald(fit), silent = TRUE)
  if (!inherits(wald_result, "try-error")) {
    # Handle different output formats of Wald()
    wtab <- NULL
    if (is.data.frame(wald_result)) {
      wtab <- wald_result
    } else if (!is.null(wald_result$Wald)) {
      wtab <- wald_result$Wald
    } else if (!is.null(wald_result$WaldTable)) {
      wtab <- wald_result$WaldTable
    }
    
    if (!is.null(wtab)) {
      # Find p-value column
      pcol <- intersect(c("p.value", "P.value", "Pr(>Chisq)", "Pr(>F)", "pval", "P"), colnames(wtab))
      if (length(pcol) >= 1) {
        for (rn in rownames(beta)) {
          if (rn %in% rownames(wtab)) {
            pvals[rn] <- wtab[rn, pcol[1]]
          }
        }
      }
    }
  }
  
  # Fallback: normal approximation
  if (all(is.na(pvals)) && all(is.finite(se)) && any(se > 0)) {
    tval <- as.numeric(beta[, 1]) / se
    pvals <- 2 * pnorm(-abs(tval))
    names(pvals) <- rownames(beta)
  }
  
  # Calculate confidence intervals
  ci_lower <- as.numeric(beta[, 1]) - 1.96 * se
  ci_upper <- as.numeric(beta[, 1]) + 1.96 * se
  
  data.frame(
    term = rownames(beta),
    estimate = as.numeric(beta[, 1]),
    se = as.numeric(se),
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
}

# Fit additive epistasis model
fit_additive_epistasis <- function(df, y_var, g1_col, g2_col, G_matrix, n_pcs = 5) {
  # Prepare data - FIXED: 使用 base R 方法避免 dplyr 冲突
  dd <- df[, c("line_code", y_var, g1_col, g2_col), drop = FALSE]
  colnames(dd) <- c("line_code", "y_var", "g1", "g2")
  
  # Remove rows with missing values
  dd <- dd[complete.cases(dd), ]
  
  if (nrow(dd) < 20) {
    return(list(method = "FAILED", n = nrow(dd), fixed = NULL, model = NULL, data = dd))
  }
  
  # Ensure common lines
  dd$line_code <- as.character(dd$line_code)
  common <- intersect(dd$line_code, rownames(G_matrix))
  dd <- dd[dd$line_code %in% common, , drop = FALSE]
  if (nrow(dd) < 20) {
    return(list(method = "FAILED", n = nrow(dd), fixed = NULL, model = NULL, data = dd))
  }
  
  # Add PCs
  G_sub <- G_matrix[common, common]
  pca_res <- pca_from_G(G_sub, n_pcs = n_pcs)
  pc_df <- as.data.frame(pca_res$scores)
  pc_df$line_code <- rownames(pc_df)
  dd <- merge(dd, pc_df, by = "line_code", all.x = TRUE)
  
  # Prepare variables
  dd$y <- dd$y_var
  dd$g1 <- as.numeric(dd$g1)
  dd$g2 <- as.numeric(dd$g2)
  dd$g1g2 <- dd$g1 * dd$g2
  
  # Create sparse G matrix
  lc <- as.character(dd$line_code)
  G_sp <- as(Matrix(G_matrix[lc, lc], sparse = TRUE), "dgCMatrix")
  
  # Build formula
  pc_terms <- paste0("PC", 1:min(n_pcs, sum(grepl("^PC\\d+$", colnames(dd)))), collapse = " + ")
  if (nchar(pc_terms) > 0) {
    f_fixed <- as.formula(paste0("y ~ g1 + g2 + g1g2 + ", pc_terms))
  } else {
    f_fixed <- as.formula("y ~ g1 + g2 + g1g2")
  }
  
  # Fit model with sommer
  fit <- try(
    sommer::mmer(
      fixed = f_fixed,
      random = ~ sommer::vsr(line_code, Gu = G_sp),
      rcov = ~ units,
      data = dd,
      verbose = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    # Fallback to linear model
    if (nchar(pc_terms) > 0) {
      f_lm <- as.formula(paste0("y ~ g1 + g2 + g1g2 + ", pc_terms))
    } else {
      f_lm <- as.formula("y ~ g1 + g2 + g1g2")
    }
    
    lmfit <- try(lm(f_lm, data = dd), silent = TRUE)
    if (inherits(lmfit, "try-error")) {
      return(list(method = "FAILED", n = nrow(dd), fixed = NULL, model = NULL, data = dd))
    }
    
    sm <- summary(lmfit)$coefficients
    fixed <- data.frame(
      term = rownames(sm),
      estimate = sm[, 1],
      se = sm[, 2],
      ci_lower = sm[, 1] - 1.96 * sm[, 2],
      ci_upper = sm[, 1] + 1.96 * sm[, 2],
      p = sm[, 4],
      stringsAsFactors = FALSE
    )
    
    return(list(method = "LM_PCs", n = nrow(dd), fixed = fixed, model = lmfit, data = dd))
  }
  
  fixed <- get_fixed_effects_from_mmer(fit)
  list(method = "SOMMER_G_PCs", n = nrow(dd), fixed = fixed, model = fit, data = dd)
}

# -------------------- 3) Main Analysis --------------------
cat("\n[Step 1] Reading QTL table...\n")
stop_if_missing(qtl_file)
qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)

# Normalize column names
colnames(qtl_df) <- gsub(" ", "_", colnames(qtl_df))
if ("Phenotype_Variance_Explained(%)" %in% colnames(qtl_df)) {
  colnames(qtl_df)[colnames(qtl_df) == "Phenotype_Variance_Explained(%)"] <- "PVE"
}

if (!all(c("QTLname", "SNP") %in% colnames(qtl_df))) {
  stop("QTL file must contain columns: QTLname, SNP")
}

# Identify p-value column
pcol <- intersect(c("P.value", "P_value", "P.value.", "P"), colnames(qtl_df))
if (length(pcol) == 0) stop("QTL file must contain a p-value column")
pcol <- pcol[1]

# Select representative SNP per QTL (smallest p-value)
qtl_marker <- qtl_df %>%
  group_by(QTLname) %>%
  slice_min(order_by = .data[[pcol]], n = 1, with_ties = FALSE) %>%
  ungroup()

all_qtls <- unique(qtl_marker$QTLname)
cat(sprintf("  Representative markers selected for %d QTLs\n", length(all_qtls)))

# Check if we have enough QTLs
if (length(all_qtls) < 2) {
  cat("  Warning: Only 1 QTL found. Cannot analyze interactions.\n")
  cat("  Will proceed with single QTL analysis instead.\n")
}

# Limit number of QTLs for speed if needed
max_qtls <- min(30, length(all_qtls))
if (length(all_qtls) > max_qtls) {
  if ("PVE" %in% colnames(qtl_marker)) {
    qtl_marker <- qtl_marker %>% arrange(desc(PVE)) %>% head(max_qtls)
  } else {
    qtl_marker <- qtl_marker %>% arrange(.data[[pcol]]) %>% head(max_qtls)
  }
  all_qtls <- qtl_marker$QTLname
  cat(sprintf("  Limiting to top %d QTLs for speed\n", max_qtls))
}

# -------------------- 4) Compute Line Effects --------------------
cat("\n[Step 2] Computing line intercepts and slopes...\n")
pheno <- safe_read_table(phenotype_file)
colnames(pheno) <- tolower(colnames(pheno))

# Identify trait column
tcol <- intersect(c("tkw", "TKW", "ph"), colnames(pheno))
if (length(tcol) == 0) stop("Phenotype file must contain TKW column")
tcol <- tcol[1]

pheno$line_code <- as.character(pheno$genotype)
pheno$env_code <- as.character(pheno$env_code)
pheno$PH <- as.numeric(pheno[[tcol]])
pheno <- pheno[is.finite(pheno$PH), c("line_code", "env_code", "PH")]

# Read environment data
ecs <- read.csv(environment_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(ecs)[1] <- "env_code"
ecs$env_code <- as.character(ecs$env_code)

# Select environment factor
target_factor <- "PAR_TEMP&GS67"
if (!target_factor %in% colnames(ecs)) {
  poss <- grep("PAR.*TEMP", colnames(ecs), value = TRUE, ignore.case = TRUE)
  if (length(poss) > 0) {
    target_factor <- poss[1]
  } else {
    numc <- sapply(ecs[, -1, drop = FALSE], is.numeric)
    if (any(numc)) {
      target_factor <- names(ecs)[-1][which(numc)[1]]
    } else {
      stop("No numeric environment factor found")
    }
  }
}

env <- ecs[, c("env_code", target_factor), drop = FALSE]
colnames(env)[2] <- "env_factor_value"
env$env_factor_value <- as.numeric(env$env_factor_value)
env <- env[is.finite(env$env_factor_value), ]

# Merge and compute line effects
pheno_env <- merge(pheno, env, by = "env_code", all.x = TRUE)
pheno_env <- pheno_env[is.finite(pheno_env$env_factor_value), ]

line_effects <- extract_line_effects(pheno_env, min_points = 3)
if (nrow(line_effects) < 20) stop("Too few lines with valid intercept/slope")

write.csv(line_effects, file.path(output_dir, "Line_Intercept_Slope_Effects.csv"), row.names = FALSE)
cat(sprintf("  Line effects computed: %d lines\n", nrow(line_effects)))
cat(sprintf("  Environment factor: %s\n", target_factor))

# -------------------- 5) Build G Matrix --------------------
cat("\n[Step 3] Building G matrix from VCF...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt_raw <- extract.gt(vcf, element = "GT")
vcf_samples <- colnames(gt_raw)

# Match samples
ms <- match_samples(line_effects$line_code, vcf_samples)
if (length(ms$line_codes) < 20) {
  stop("Too few line codes matched between phenotype and VCF")
}

# Subset line effects
line_effects2 <- line_effects[line_effects$line_code %in% ms$line_codes, ]
line_effects2 <- line_effects2[match(ms$line_codes, line_effects2$line_code), ]

# Subset genotype matrix
gt_raw2 <- gt_raw[, ms$matched_samples, drop = FALSE]
colnames(gt_raw2) <- ms$line_codes

cat(sprintf("  VCF matrix for G: %d SNPs × %d lines\n", nrow(gt_raw2), ncol(gt_raw2)))

# Convert to dosage
cat("  Converting GT to dosage for G matrix...\n")
G_snp_line <- apply(gt_raw2, c(1, 2), gt_to_dosage)
Geno <- t(G_snp_line)  # lines x SNPs

cat(sprintf("  Genotype matrix: %d lines × %d SNPs\n", nrow(Geno), ncol(Geno)))

# Calculate G matrix
G_matrix <- calc_G_vanraden(Geno)
saveRDS(G_matrix, file.path(output_dir, "G_matrix.rds"))
cat("  G matrix saved\n")

# PCA from G
pca_res <- pca_from_G(G_matrix, n_pcs = 10)
write.csv(pca_res$scores, file.path(output_dir, "PCA_scores.csv"))
cat("  Top 5 PC variance explained (%):\n")
print(round(pca_res$var_explained[1:min(5, length(pca_res$var_explained))], 2))

# -------------------- 6) Extract QTL Genotypes --------------------
cat("\n[Step 4] Extracting QTL marker genotypes...\n")
target_snps <- unique(qtl_marker$SNP)
vcf_ids <- vcf@fix[, "ID"]
idx <- which(vcf_ids %in% target_snps)

if (length(idx) == 0) {
  cat("  No representative SNP IDs found in VCF by ID. Trying alternative matching...\n")
  # Try matching by chromosome and position
  qtl_positions <- qtl_marker %>% 
    mutate(pos_id = paste0(Chromosome, "_", Position))
  vcf_positions <- paste0(vcf@fix[, "CHROM"], "_", vcf@fix[, "POS"])
  idx <- which(vcf_positions %in% qtl_positions$pos_id)
  
  if (length(idx) == 0) {
    stop("No matching SNPs found in VCF")
  }
}

gt_qtl <- gt_raw[idx, ms$matched_samples, drop = FALSE]
colnames(gt_qtl) <- ms$line_codes

dos_qtl <- apply(gt_qtl, c(1, 2), gt_to_dosage)
dos_qtl <- t(dos_qtl)
dos_qtl <- as.data.frame(dos_qtl)
dos_qtl$line_code <- rownames(dos_qtl)

# Rename columns with QTL names
snp_to_qtl <- setNames(qtl_marker$QTLname, qtl_marker$SNP)
# If we matched by position, use position ID
if (length(snp_to_qtl) == 0) {
  pos_to_qtl <- setNames(qtl_marker$QTLname, paste0(qtl_marker$Chromosome, "_", qtl_marker$Position))
  new_names <- sapply(colnames(dos_qtl)[colnames(dos_qtl) != "line_code"], function(col) {
    # Extract position from column name if possible
    for (pos in names(pos_to_qtl)) {
      if (grepl(pos, col)) {
        return(make.names(pos_to_qtl[[pos]]))
      }
    }
    return(col)
  })
} else {
  new_names <- sapply(colnames(dos_qtl)[colnames(dos_qtl) != "line_code"], function(snp) {
    qn <- snp_to_qtl[[snp]]
    if (is.null(qn) || is.na(qn)) return(snp)
    make.names(qn)
  })
}

colnames(dos_qtl)[colnames(dos_qtl) != "line_code"] <- new_names
colnames(dos_qtl) <- make.unique(colnames(dos_qtl))

# Merge with line effects
dat0 <- merge(line_effects2, dos_qtl, by = "line_code", all.x = FALSE)
cat(sprintf("  Analysis dataset: %d lines × %d columns\n", nrow(dat0), ncol(dat0)))

# Identify QTL columns
qtl_cols <- colnames(dat0)[sapply(colnames(dat0), function(x) x %in% make.names(all_qtls) || grepl("QTL", x))]
cat(sprintf("  Found %d QTL genotype columns\n", length(qtl_cols)))

if (length(qtl_cols) < 2) {
  cat("  Warning: Too few QTL genotype columns found. Cannot analyze interactions.\n")
  cat("  Will analyze single QTL effects only.\n")
  
  # Analyze single QTL effects instead
  single_qtl_results <- data.frame()
  for (qtl in qtl_cols) {
    for (trait in c("Intercept", "Slope")) {
      # Simple linear model for single QTL
      formula <- as.formula(paste(trait, "~", qtl))
      model <- lm(formula, data = dat0)
      sm <- summary(model)
      
      single_qtl_results <- rbind(single_qtl_results, data.frame(
        QTL = qtl,
        Trait = trait,
        Estimate = sm$coefficients[2, 1],
        SE = sm$coefficients[2, 2],
        P_value = sm$coefficients[2, 4],
        N = nrow(dat0),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  write.csv(single_qtl_results, file.path(output_dir, "Single_QTL_Effects.csv"), row.names = FALSE)
  cat("  Single QTL analysis completed\n")
  quit(status = 0)
}

# -------------------- 7) Pairwise Analysis --------------------
cat("\n[Step 5] Performing pairwise additive epistasis analysis...\n")
qtl_pairs <- combn(qtl_cols, 2, simplify = FALSE)
cat(sprintf("  Total QTL pairs: %d\n", length(qtl_pairs)))

# Initialize results storage
results <- list()
detailed_results <- list()
pair_counter <- 0

pb <- txtProgressBar(min = 0, max = length(qtl_pairs), style = 3)
for (i in seq_along(qtl_pairs)) {
  setTxtProgressBar(pb, i)
  g1 <- qtl_pairs[[i]][1]
  g2 <- qtl_pairs[[i]][2]
  
  # Skip pairs with insufficient data - FIXED: 使用 base R 方法
  tmp <- dat0[, c("line_code", "Intercept", "Slope", g1, g2)]
  tmp <- tmp[complete.cases(tmp), ]
  if (nrow(tmp) < 20) next
  if (length(unique(tmp[[g1]])) < 2 || length(unique(tmp[[g2]])) < 2) next
  
  for (trait in c("Intercept", "Slope")) {
    # Fit model
    fitres <- fit_additive_epistasis(dat0, trait, g1, g2, G_matrix, n_pcs = 5)
    
    # Extract coefficients
    btab <- fitres$fixed
    b1 <- btab[btab$term == "g1", , drop = FALSE]
    b2 <- btab[btab$term == "g2", , drop = FALSE]
    b12 <- btab[btab$term == "g1g2", , drop = FALSE]
    
    # Store summary results
    result_row <- data.frame(
      QTL1 = g1,
      QTL2 = g2,
      Trait = trait,
      N = fitres$n,
      Method = fitres$method,
      b1 = ifelse(nrow(b1) > 0, b1$estimate, NA_real_),
      se1 = ifelse(nrow(b1) > 0, b1$se, NA_real_),
      p1 = ifelse(nrow(b1) > 0, b1$p, NA_real_),
      b2 = ifelse(nrow(b2) > 0, b2$estimate, NA_real_),
      se2 = ifelse(nrow(b2) > 0, b2$se, NA_real_),
      p2 = ifelse(nrow(b2) > 0, b2$p, NA_real_),
      b12 = ifelse(nrow(b12) > 0, b12$estimate, NA_real_),
      se12 = ifelse(nrow(b12) > 0, b12$se, NA_real_),
      p12 = ifelse(nrow(b12) > 0, b12$p, NA_real_),
      stringsAsFactors = FALSE
    )
    
    # Store detailed results
    pair_counter <- pair_counter + 1
    detailed_results[[pair_counter]] <- list(
      qtl_pair = c(g1, g2),
      trait = trait,
      model_result = fitres,
      summary = result_row
    )
    
    results[[length(results) + 1]] <- result_row
  }
}
close(pb)

# Combine results
if (length(results) == 0) {
  cat("  Warning: No valid pair results produced. Check your data.\n")
  quit(status = 0)
}

res_df <- bind_rows(results)

# Save raw results
write.csv(res_df, file.path(output_dir, "All_QTL_Pairs_Additive_Effects.csv"), row.names = FALSE)
saveRDS(detailed_results, file.path(output_dir, "Detailed_Model_Results.rds"))
cat(sprintf("\n  Results saved: %d model fits completed\n", length(detailed_results)))

# -------------------- 8) Multiple Testing Correction --------------------
cat("\n[Step 6] Applying multiple testing correction...\n")

# Add FDR and Bonferroni corrections
res_df <- res_df %>%
  group_by(Trait) %>%
  mutate(
    p1_fdr = p.adjust(p1, method = "fdr"),
    p2_fdr = p.adjust(p2, method = "fdr"),
    p12_fdr = p.adjust(p12, method = "fdr"),
    p1_bonf = p.adjust(p1, method = "bonferroni"),
    p2_bonf = p.adjust(p2, method = "bonferroni"),
    p12_bonf = p.adjust(p12, method = "bonferroni")
  ) %>%
  ungroup()

write.csv(res_df, file.path(output_dir, "All_QTL_Pairs_Additive_Effects_Corrected.csv"), row.names = FALSE)

# -------------------- 9) Significant Results --------------------
cat("\n[Step 7] Extracting significant results...\n")

# Define significance thresholds
thresholds <- list(
  nominal = 0.05,
  fdr = 0.1,
  bonferroni = 0.05
)

# Extract significant interactions
sig_results <- list()
for (trait in c("Intercept", "Slope")) {
  trait_data <- res_df %>% filter(Trait == trait)
  
  # Nominal significance
  sig_nominal <- trait_data %>% filter(p12 < thresholds$nominal)
  
  # FDR significance
  sig_fdr <- trait_data %>% filter(p12_fdr < thresholds$fdr)
  
  # Bonferroni significance
  sig_bonf <- trait_data %>% filter(p12_bonf < thresholds$bonferroni)
  
  sig_results[[trait]] <- list(
    nominal = sig_nominal,
    fdr = sig_fdr,
    bonferroni = sig_bonf,
    summary = data.frame(
      Trait = trait,
      Total_pairs = nrow(trait_data),
      Nominal_sig = nrow(sig_nominal),
      FDR_sig = nrow(sig_fdr),
      Bonferroni_sig = nrow(sig_bonf),
      stringsAsFactors = FALSE
    )
  )
  
  # Save significant results
  if (nrow(sig_nominal) > 0) {
    write.csv(sig_nominal, 
              file.path(output_dir, "Significant_Pairs", sprintf("%s_Nominal_Significant.csv", trait)),
              row.names = FALSE)
  }
}

# -------------------- 10) Summary Report --------------------
cat("\n[Step 8] Generating summary report...\n")

# Create summary statistics
summary_stats <- data.frame(
  Total_QTLs = length(qtl_cols),
  Total_pairs = length(qtl_pairs),
  Total_models = nrow(res_df),
  Lines_analyzed = nrow(dat0),
  Environment_factor = target_factor,
  G_matrix_SNPs = ncol(Geno),
  stringsAsFactors = FALSE
)

# Write summary to file
sink(file.path(output_dir, "Analysis_Summary.txt"))
cat("========================================\n")
cat("   QTL ADDITIVE EPISTASIS ANALYSIS\n")
cat("========================================\n\n")
cat(sprintf("Analysis date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output directory: %s\n\n", output_dir))

cat("1. DATA SUMMARY\n")
cat(sprintf("   Total QTLs analyzed: %d\n", summary_stats$Total_QTLs))
cat(sprintf("   Total QTL pairs: %d\n", summary_stats$Total_pairs))
cat(sprintf("   Total models fitted: %d\n", summary_stats$Total_models))
cat(sprintf("   Lines with complete data: %d\n", summary_stats$Lines_analyzed))
cat(sprintf("   Environment factor: %s\n", summary_stats$Environment_factor))
cat(sprintf("   SNPs used for G matrix: %d\n\n", summary_stats$G_matrix_SNPs))

cat("2. SIGNIFICANT INTERACTIONS\n")
for (trait in c("Intercept", "Slope")) {
  if (!is.null(sig_results[[trait]])) {
    sig_info <- sig_results[[trait]]$summary
    cat(sprintf("   %s:\n", trait))
    cat(sprintf("     Nominal (p < 0.05): %d pairs\n", sig_info$Nominal_sig))
    cat(sprintf("     FDR (q < 0.10): %d pairs\n", sig_info$FDR_sig))
    cat(sprintf("     Bonferroni (p < 0.05): %d pairs\n\n", sig_info$Bonferroni_sig))
  }
}

cat("3. TOP INTERACTIONS BY TRAIT\n")
for (trait in c("Intercept", "Slope")) {
  if (nrow(res_df) > 0) {
    top5 <- res_df %>%
      filter(Trait == trait, !is.na(p12)) %>%
      arrange(p12) %>%
      head(min(5, sum(!is.na(res_df$p12))))
    
    if (nrow(top5) > 0) {
      cat(sprintf("   %s (Top 5 by p-value):\n", trait))
      for (i in 1:nrow(top5)) {
        row <- top5[i, ]
        cat(sprintf("     %d. %s × %s: b12 = %.4f, p = %.2e\n", 
                    i, row$QTL1, row$QTL2, row$b12, row$p12))
      }
      cat("\n")
    }
  }
}

cat("4. INTERPRETATION\n")
cat("   - b1, b2: Additive effects of QTL1 and QTL2\n")
cat("   - b12: Additive × Additive interaction effect\n")
cat("   - Positive b12: Synergistic interaction\n")
cat("   - Negative b12: Antagonistic interaction\n")
cat("   - Model includes G-matrix random effects and PCA covariates\n\n")

cat("5. OUTPUT FILES\n")
cat("   - Line_Intercept_Slope_Effects.csv: Line-specific effects\n")
cat("   - G_matrix.rds: Genomic relationship matrix\n")
cat("   - PCA_scores.csv: Principal component scores\n")
cat("   - All_QTL_Pairs_Additive_Effects*.csv: All results\n")
cat("   - Significant_Pairs/: Significant interactions\n")
cat("   - Figures/: Visualizations\n")
cat("   - Detailed_Model_Results.rds: Complete model objects\n\n")

cat("========================================\n")
cat("   ANALYSIS COMPLETE\n")
cat("========================================\n")
sink()

# -------------------- 11) Final Summary --------------------
cat("\n========== Analysis Complete ==========\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nKey findings:\n")

for (trait in c("Intercept", "Slope")) {
  if (!is.null(sig_results[[trait]])) {
    sig_info <- sig_results[[trait]]$summary
    cat(sprintf("  %s:\n", trait))
    cat(sprintf("    Significant interactions: %d (%.1f%%)\n", 
                sig_info$Nominal_sig, 
                sig_info$Nominal_sig/sig_info$Total_pairs*100))
    
    # Show top interaction
    if (nrow(res_df) > 0) {
      top <- res_df %>%
        filter(Trait == trait, !is.na(p12)) %>%
        arrange(p12) %>%
        head(1)
      
      if (nrow(top) > 0) {
        cat(sprintf("    Most significant: %s × %s (p = %.2e)\n", 
                    top$QTL1, top$QTL2, top$p12))
      }
    }
    cat("\n")
  }
}

cat("Generated files:\n")
cat("  1. Analysis_Summary.txt - Complete analysis report\n")
cat("  2. All_QTL_Pairs_Additive_Effects*.csv - Statistical results\n")
cat("  3. Significant_Pairs/ - Significant interactions\n")
cat("  4. G_matrix.rds, PCA_scores.csv - Background control data\n")
cat("  5. Detailed_Model_Results.rds - Complete model objects\n")

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))