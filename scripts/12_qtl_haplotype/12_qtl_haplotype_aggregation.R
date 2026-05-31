# ============================================================
# Line vs Haplotype Fitted Lines (Top6 Color + Others Gray + Haplotype-based Multiple Comparison)
# 最终修复版：兼容无roman()函数环境 + 鲁棒性优化 + 完整功能
# ============================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(grid)
  library(car)       # For ANOVA
  library(multcomp)  # For Tukey HSD
  library(cowplot)   # For legend extraction
})

# ==================== 1. Paths ====================
base_dir <- "/mnt/7t_storage/zhangcl/TKW/"
output_dir <- file.path(base_dir, "Haplotype_Environment_Analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

vcf_file <- file.path(base_dir, "extracted_G1_to_G339.vcf")
qtl_file <- file.path(base_dir, "QTL4.csv")
phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")
environment_file <- file.path(base_dir, "ECs_results.csv")

target_factor <- "PAR_TEMP&GS67"

# Top6 legend placement (normalized 0~1)
legend_x <- 0.98
legend_y <- 0.98
legend_width <- 0.3
legend_height <- 0.3

# ==================== 2. Read QTL markers ====================
cat("Reading QTL markers...\n")
qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)

if ("SNP" %in% colnames(qtl_df)) {
  qtl_snps <- unique(qtl_df$SNP)
} else if ("Marker" %in% colnames(qtl_df)) {
  qtl_snps <- unique(qtl_df$Marker)
} else {
  snp_col <- grep("SNP|snp|Marker|marker|ID|id", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)
  if (length(snp_col) > 0) {
    qtl_snps <- unique(qtl_df[[snp_col[1]]])
  } else {
    stop("Cannot find SNP/Marker column in QTL.csv")
  }
}
qtl_snps <- qtl_snps[!is.na(qtl_snps) & qtl_snps != ""]
cat(sprintf("Found %d QTL markers\n", length(qtl_snps)))

# ==================== 3. Read VCF and extract target SNPs ====================
cat("Reading VCF...\n")
vcf_data <- read.vcfR(vcf_file, verbose = FALSE)

vcf_snp_ids <- vcf_data@fix[, "ID"]
matched_indices <- which(vcf_snp_ids %in% qtl_snps)

if (length(matched_indices) == 0) {
  cat("No match by ID. Trying CHROM+POS matching...\n")
  chr_col <- grep("CHR|Chr|chr|chromosome|Chromosome", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)
  pos_col <- grep("POS|Pos|pos|position|Position", colnames(qtl_df),
                  value = TRUE, ignore.case = TRUE)

  if (length(chr_col) > 0 && length(pos_col) > 0) {
    qtl_pos_ids <- paste(qtl_df[[chr_col[1]]], qtl_df[[pos_col[1]]], sep = "_")
    vcf_pos_ids <- paste(vcf_data@fix[, "CHROM"], vcf_data@fix[, "POS"], sep = "_")
    matched_indices <- which(vcf_pos_ids %in% qtl_pos_ids)
  }
}

if (length(matched_indices) == 0) stop("No matched SNP markers found in VCF.")
cat(sprintf("Matched %d markers\n", length(matched_indices)))

extracted_vcf <- vcf_data[matched_indices, ]

# ==================== 4. Extract genotypes and convert to 0/1/2 ====================
gt_matrix <- extract.gt(extracted_vcf, element = "GT")  # rows=SNP, cols=Samples

convert_gt_to_numeric <- function(gt_mat) {
  out <- matrix(NA_real_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  for (i in 1:nrow(gt_mat)) {
    for (j in 1:ncol(gt_mat)) {
      gt_str <- gt_mat[i, j]
      if (is.na(gt_str) || gt_str %in% c("./.", ".|.", ".")) {
        out[i, j] <- NA_real_
      } else {
        alleles <- strsplit(gt_str, "[|/]")[[1]]
        if (length(alleles) == 2) {
          a1 <- suppressWarnings(as.numeric(alleles[1]))
          a2 <- suppressWarnings(as.numeric(alleles[2]))
          if (is.finite(a1) && is.finite(a2)) out[i, j] <- a1 + a2
        }
      }
    }
  }
  rownames(out) <- rownames(gt_mat)
  colnames(out) <- colnames(gt_mat)
  out
}

gt_numeric <- convert_gt_to_numeric(gt_matrix)

# ==================== 5. Build haplotypes (one row per sample) ====================
cat("Building haplotypes...\n")
haplotype_df <- as.data.frame(t(gt_numeric), stringsAsFactors = FALSE) # samples × SNPs
haplotype_df$Sample <- rownames(haplotype_df)

snp_cols <- setdiff(colnames(haplotype_df), "Sample")
haplotype_df$Haplotype <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
  paste(ifelse(is.na(x), "N", x), collapse = "-")
})

# 修复：替换dplyr::select为基础R列选择（兼容低版本）
haplotype_df <- haplotype_df[, c("Sample", "Haplotype", snp_cols)]

# 统计单倍型频率
haplotype_counts <- as.data.frame(table(haplotype_df$Haplotype), stringsAsFactors = FALSE)
colnames(haplotype_counts) <- c("Haplotype", "Frequency")
haplotype_counts <- haplotype_counts[order(haplotype_counts$Frequency, decreasing = TRUE), ]

top6_haps <- head(haplotype_counts$Haplotype, 6)

# 修复：替换roman()函数，手动生成HapI~HapVI命名（兼容所有R环境）
top6_haps_named <- paste0("Hap", c("I", "II", "III", "IV", "V", "VI")[1:length(top6_haps)])
names(top6_haps_named) <- top6_haps
cat(sprintf("Detected %d haplotypes | Top 6: %s\n", 
            nrow(haplotype_counts), paste(top6_haps_named, collapse = ", ")))

# ==================== 6. Read phenotype and environment factor ====================
cat("Reading phenotype and environment data...\n")

phenotype_data <- read.table(
  phenotype_file, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
) 
# 修复：替换dplyr管道为基础R操作（兼容低版本）
phenotype_data$genotype <- as.character(phenotype_data$genotype)
phenotype_data$env_code <- as.character(phenotype_data$env_code)
phenotype_data$TKW <- as.numeric(phenotype_data$TKW)
colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
colnames(phenotype_data)[colnames(phenotype_data) == "TKW"] <- "PH"
phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]

ecs_data <- read.csv(environment_file, header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
colnames(ecs_data)[1] <- "env_code"  # 替换rename
ecs_data$env_code <- as.character(ecs_data$env_code)

if (!target_factor %in% colnames(ecs_data)) {
  cat("WARNING: target factor not found. Trying fuzzy match...\n")
  possible <- grep("PAR.*TEMP", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
  if (length(possible) > 0) {
    target_factor <- possible[1]
    cat("Using:", target_factor, "\n")
  } else {
    stop("No matching environment factor found.")
  }
}

# 修复：替换dplyr::select/rename为基础R
env_factor_data <- ecs_data[, c("env_code", target_factor)]
colnames(env_factor_data)[colnames(env_factor_data) == target_factor] <- "env_factor_value"
env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)

# ==================== 7. Merge PH + haplotype + env factor ====================
cat("Merging...\n")

# 修复：替换dplyr::inner_join为merge（基础R）
pheno_haplo <- merge(
  phenotype_data, 
  haplotype_df[, c("Sample", "Haplotype")],  # 基础R列选择
  by.x = "line_code", by.y = "Sample", all.x = FALSE, all.y = FALSE
)

pheno_haplo_env <- merge(
  pheno_haplo, 
  env_factor_data, 
  by = "env_code", all.x = FALSE, all.y = FALSE
)
pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$env_factor_value), ]

# 为观测数据添加单倍型命名（HapI~VI）
pheno_haplo_env$Haplotype_Named <- ifelse(
  pheno_haplo_env$Haplotype %in% top6_haps,
  top6_haps_named[match(pheno_haplo_env$Haplotype, names(top6_haps_named))],
  "Other"
)

cat(sprintf("Merged rows: %d\n", nrow(pheno_haplo_env)))
cat(sprintf("Lines: %d | Haplotypes: %d\n",
            length(unique(pheno_haplo_env$line_code)),
            length(unique(pheno_haplo_env$Haplotype))))

# ==================== 8. Fit per-line models (dashed) ====================
cat("Fitting per-line models...\n")

# 修复：替换dplyr分组统计为基础R
line_env_count <- aggregate(env_code ~ line_code, data = pheno_haplo_env, FUN = function(x) length(unique(x)))
colnames(line_env_count)[2] <- "n_env"
valid_lines <- line_env_count$line_code[line_env_count$n_env >= 3]

line_models <- data.frame()
for (line in valid_lines) {
  line_data <- pheno_haplo_env[pheno_haplo_env$line_code == line, ]
  if (nrow(line_data) >= 3) {
    m <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
    if (!inherits(m, "try-error")) {
      hap <- unique(line_data$Haplotype)[1]
      hap_named <- ifelse(hap %in% top6_haps,
                          top6_haps_named[match(hap, names(top6_haps_named))],
                          "Other")
      line_models <- rbind(line_models, data.frame(
        line_code = line,
        Haplotype = hap,
        Haplotype_Named = hap_named,
        Slope = as.numeric(coef(m)[2]),
        Intercept = as.numeric(coef(m)[1]),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 鲁棒性优化：处理无有效模型的情况
if (nrow(line_models) == 0) {
  cat("WARNING: No valid line models found. Reducing env requirement to 2...\n")
  valid_lines <- line_env_count$line_code[line_env_count$n_env >= 2]
  for (line in valid_lines) {
    line_data <- pheno_haplo_env[pheno_haplo_env$line_code == line, ]
    if (nrow(line_data) >= 2) {
      m <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
      if (!inherits(m, "try-error")) {
        hap <- unique(line_data$Haplotype)[1]
        hap_named <- ifelse(hap %in% top6_haps,
                            top6_haps_named[match(hap, names(top6_haps_named))],
                            "Other")
        line_models <- rbind(line_models, data.frame(
          line_code = line,
          Haplotype = hap,
          Haplotype_Named = hap_named,
          Slope = as.numeric(coef(m)[2]),
          Intercept = as.numeric(coef(m)[1]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

cat(sprintf("Line models: %d\n", nrow(line_models)))

# ==================== 9. Fit per-haplotype mean models (solid) ====================
cat("Fitting haplotype mean models...\n")

haplotype_models <- data.frame()
for (hap in unique(pheno_haplo_env$Haplotype)) {
  dat <- pheno_haplo_env[pheno_haplo_env$Haplotype == hap, ]
  n_lines <- length(unique(dat$line_code))
  n_obs <- nrow(dat)

  # 鲁棒性优化：降低模型拟合门槛
  if (n_lines >= 2 && n_obs >= 4) {
    m <- try(lm(PH ~ env_factor_value, data = dat), silent = TRUE)
    if (!inherits(m, "try-error")) {
      hap_named <- ifelse(hap %in% top6_haps,
                          top6_haps_named[match(hap, names(top6_haps_named))],
                          "Other")
      haplotype_models <- rbind(haplotype_models, data.frame(
        Haplotype = hap,
        Haplotype_Named = hap_named,
        Slope = as.numeric(coef(m)[2]),
        Intercept = as.numeric(coef(m)[1]),
        Is_Top6 = ifelse(hap %in% top6_haps, "Top6", "Other"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("Haplotype models: %d (Top6: %d)\n", 
            nrow(haplotype_models), sum(haplotype_models$Is_Top6 == "Top6")))

write.csv(line_models, file.path(output_dir, "Line_Models_by_Haplotype.csv"), row.names = FALSE)
write.csv(haplotype_models, file.path(output_dir, "Haplotype_Models.csv"), row.names = FALSE)

# ==================== 10. Significance Test: ANOVA + Tukey HSD (Haplotype as Group) ====================
cat("Testing line slope/intercept significance (ANOVA + Tukey HSD, grouped by haplotype)...\n")

# 鲁棒性优化：仅当有足够数据时执行ANOVA
if (nrow(line_models) >= 5 && length(unique(line_models$Haplotype)) >= 2) {
  # 修复1：基础R列选择（替换dplyr::select）
  line_params <- line_models[, c("Haplotype", "Slope", "Intercept")]
  # 修复2：Haplotype转为因子（ANOVA/Tukey必需）
  line_params$Haplotype <- as.factor(line_params$Haplotype)

  # -------------------- Slope: ANOVA + Tukey HSD --------------------
  slope_anova <- aov(Slope ~ Haplotype, data = line_params)
  slope_tukey <- glht(slope_anova, linfct = mcp(Haplotype = "Tukey"))
  slope_tukey_res <- summary(slope_tukey)

  # 提取Tukey结果（基础R方式）
  slope_tukey_df <- data.frame(
    Comparison = names(slope_tukey_res$test$coefficients),
    Slope_Difference = as.numeric(slope_tukey_res$test$coefficients),
    P_Value = as.numeric(slope_tukey_res$test$pvalues),
    stringsAsFactors = FALSE
  )
  slope_tukey_df$Significance <- ifelse(slope_tukey_df$P_Value < 0.05, "Significant", "Not Significant")

  # -------------------- Intercept: ANOVA + Tukey HSD --------------------
  intercept_anova <- aov(Intercept ~ Haplotype, data = line_params)
  intercept_tukey <- glht(intercept_anova, linfct = mcp(Haplotype = "Tukey"))
  intercept_tukey_res <- summary(intercept_tukey)

  # 提取Tukey结果（基础R方式）
  intercept_tukey_df <- data.frame(
    Comparison = names(intercept_tukey_res$test$coefficients),
    Intercept_Difference = as.numeric(intercept_tukey_res$test$coefficients),
    P_Value = as.numeric(intercept_tukey_res$test$pvalues),
    stringsAsFactors = FALSE
  )
  intercept_tukey_df$Significance <- ifelse(intercept_tukey_df$P_Value < 0.05, "Significant", "Not Significant")

  # 保存多重比较结果
  write.csv(slope_tukey_df, file.path(output_dir, "Line_Slope_Tukey_HSD_Results.csv"), row.names = FALSE)
  write.csv(intercept_tukey_df, file.path(output_dir, "Line_Intercept_Tukey_HSD_Results.csv"), row.names = FALSE)
} else {
  cat("WARNING: Insufficient data for ANOVA/Tukey test. Skipping...\n")
  # 创建空结果文件避免后续报错
  write.csv(data.frame(), file.path(output_dir, "Line_Slope_Tukey_HSD_Results.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(output_dir, "Line_Intercept_Tukey_HSD_Results.csv"), row.names = FALSE)
}

cat("Line slope/intercept Tukey HSD results saved.\n")

# ==================== 11. Color Map (Top6 Colored, Others Gray) ====================
cat("Creating color map (Top6 colored, others gray)...\n")

# Top6颜色（使用命名后的HapI~VI）
top6_cols <- brewer.pal(min(6, length(top6_haps)), "Set2")
names(top6_cols) <- top6_haps_named

# 所有单倍型的颜色映射（原始名→颜色）
all_haps <- unique(c(as.character(line_models$Haplotype), as.character(haplotype_models$Haplotype)))
hap_colors <- ifelse(
  all_haps %in% top6_haps,
  top6_cols[match(top6_haps_named[match(all_haps, names(top6_haps_named))], names(top6_cols))],
  "gray80"  # 其他单倍型用灰色
)
names(hap_colors) <- all_haps

# 命名后的单倍型颜色映射（用于图例）
hap_colors_named <- c(top6_cols, "Other" = "gray80")

# ==================== 12. Prediction Data ====================
env_range <- range(pheno_haplo_env$env_factor_value, na.rm = TRUE)
env_seq <- seq(env_range[1], env_range[2], length.out = 120)

# 株系预测（所有单倍型，Top6彩色/其他灰色）
line_pred_list <- list()
if (nrow(line_models) > 0) {
  for (i in 1:nrow(line_models)) {
    one <- line_models[i, ]
    line_pred_list[[i]] <- data.frame(
      line_code = one$line_code,
      Haplotype = one$Haplotype,
      Haplotype_Named = one$Haplotype_Named,
      env_factor_value = env_seq,
      Predicted_PH = one$Intercept + one$Slope * env_seq,
      Is_Top6 = ifelse(one$Haplotype %in% top6_haps, "Top6", "Other"),
      stringsAsFactors = FALSE
    )
  }
}
line_pred <- if (length(line_pred_list) > 0) do.call(rbind, line_pred_list) else data.frame()

# 单倍型平均预测（Top6彩色/其他灰色）
hap_pred_list <- list()
if (nrow(haplotype_models) > 0) {
  for (i in 1:nrow(haplotype_models)) {
    one <- haplotype_models[i, ]
    hap_pred_list[[i]] <- data.frame(
      Haplotype = one$Haplotype,
      Haplotype_Named = one$Haplotype_Named,
      env_factor_value = env_seq,
      Predicted_PH = one$Intercept + one$Slope * env_seq,
      Is_Top6 = one$Is_Top6,
      stringsAsFactors = FALSE
    )
  }
}
hap_pred <- if (length(hap_pred_list) > 0) do.call(rbind, hap_pred_list) else data.frame()

# ==================== 13. Top6 Haplotype Legend Plot ====================
cat("Building Top6 haplotype legend...\n")
p_legend <- ggplot() +
  geom_line(
    data = data.frame(
      Haplotype_Named = top6_haps_named,
      x = 0, y = 0
    ),
    aes(x = x, y = y, color = Haplotype_Named, group = Haplotype_Named)
  ) +
  # 修改图例标题为"Haplotypes"，使用命名后的HapI~VI
  scale_color_manual(values = top6_cols, name = "Haplotypes") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

legend_grob <- get_legend(p_legend)

# ==================== 14. MAIN Plot (Top6 Colored, Others Gray) ====================
cat("Drawing main plot...\n")

p_main <- ggplot() +
  # 1. 添加株系实际观测值（散点）- 核心展示内容
  geom_point(
    data = pheno_haplo_env,
    aes(x = env_factor_value, y = PH, color = Haplotype, shape = env_code),
    alpha = 0.6, size = 2, stroke = 0.5
  ) +
  # 2. 株系拟合线（仅当有数据时添加）
  {if (nrow(line_pred) > 0) 
    geom_line(
      data = line_pred,
      aes(x = env_factor_value, y = Predicted_PH, group = line_code, color = Haplotype),
      linetype = "dashed", alpha = 0.35, linewidth = 0.35
    )
  } +
  # 3. 单倍型平均拟合线（仅当有数据时添加）
  {if (nrow(hap_pred) > 0)
    geom_line(
      data = hap_pred,
      aes(x = env_factor_value, y = Predicted_PH, group = Haplotype, color = Haplotype),
      linetype = "solid", alpha = 0.95, linewidth = 1.6
    )
  } +
  # 颜色映射（Top6彩色，其他灰色）
  scale_color_manual(values = hap_colors, name = "Haplotypes") +
  # 添加环境形状图例（区分不同环境）
  scale_shape_discrete(name = "Environment") +
  labs(
    title = "TKW Fitted Lines: Line vs Haplotype (Top6 Colored)",
    x = target_factor,
    y = "TKW (Thousand Kernel Weight)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(size = 11, color = "black"),
    legend.position = "none",  # 隐藏主图例，仅显示自定义Top6图例
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
    plot.margin = margin(14, 14, 14, 14),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ==================== 15. Save: Main Plot + Top6 Legend Overlay ====================
main_grob <- ggplotGrob(p_main)

png_out <- file.path(output_dir, "Combined_Main_Top6_Colored_Others_Gray.png")
pdf_out <- file.path(output_dir, "Combined_Main_Top6_Colored_Others_Gray.pdf")

# Save PNG
png(png_out, width = 18, height = 12, units = "in", res = 300)
grid.newpage()
grid.draw(main_grob)
# Add Top6 legend at top-right
pushViewport(viewport(
  x = legend_x,
  y = legend_y,
  width = unit(legend_width, "npc"),
  height = unit(legend_height, "npc"),
  just = c("right", "top")
))
grid.draw(legend_grob)
upViewport()
dev.off()

# Save PDF
pdf(pdf_out, width = 18, height = 12)
grid.newpage()
grid.draw(main_grob)
pushViewport(viewport(
  x = legend_x,
  y = legend_y,
  width = unit(legend_width, "npc"),
  height = unit(legend_height, "npc"),
  just = c("right", "top")
))
grid.draw(legend_grob)
upViewport()
dev.off()

# ==================== 16. Final Output ====================
cat("\n========== Done ==========\n")
cat("Outputs:\n")
cat("  - ", png_out, "\n", sep = "")
cat("  - ", pdf_out, "\n", sep = "")
cat("  - ", file.path(output_dir, "Haplotype_Models.csv"), "\n", sep = "")
cat("  - ", file.path(output_dir, "Line_Slope_Tukey_HSD_Results.csv"), "\n", sep = "")
cat("  - ", file.path(output_dir, "Line_Intercept_Tukey_HSD_Results.csv"), "\n", sep = "")
cat(sprintf("Env factor: %s | Range: %.3f to %.3f\n", target_factor, env_range[1], env_range[2]))
cat(sprintf("Saved to: %s\n", output_dir))