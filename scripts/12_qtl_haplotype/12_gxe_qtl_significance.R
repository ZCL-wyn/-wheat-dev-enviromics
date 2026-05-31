# ============================================================
# 图：-log10(P)_slope vs -log10(P)_intercept
#  - X轴：-log10(P) intercept
#  - Y轴：-log10(P) slope
#  - 灰色对角线 y=x
#  - 从 /mnt/7t_storage/zhangcl/TKW/QTL.csv 读取QTL对应SNP
#    -> 这些SNP点：红色 + 更大
#  - 输出PNG/PDF
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ==================== 1. 路径 ====================
base_dir <- "/mnt/7t_storage/zhangcl/TKW/"
output_dir <- file.path(base_dir, "QTL_plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

slope_gwas_file <- file.path(
  base_dir,
  "FW_GWAS_Results/GWAS_Slope_para/GAPIT.Association.GWAS_Results.GWAS_Slope_para.BLINK.Trait(NYC).csv"
)
intcp_gwas_file <- file.path(
  base_dir,
  "FW_GWAS_Results/GWAS_Intcp_para_adj/GAPIT.Association.GWAS_Results.GWAS_Intcp_para_adj.BLINK.Trait(NYC).csv"
)

qtl_file <- file.path(base_dir, "QTL4.csv")

png_file <- file.path(output_dir, "Slope_vs_Intercept_-log10P_QTL_highlight.png")
pdf_file <- file.path(output_dir, "Slope_vs_Intercept_-log10P_QTL_highlight.pdf")

# ==================== 2. 读取数据 ====================
cat("Reading slope GWAS:\n", slope_gwas_file, "\n")
slope_df <- read_csv(slope_gwas_file, show_col_types = FALSE)

cat("Reading intercept GWAS:\n", intcp_gwas_file, "\n")
intcp_df <- read_csv(intcp_gwas_file, show_col_types = FALSE)

cat("Reading QTL list:\n", qtl_file, "\n")
qtl_df <- read_csv(qtl_file, show_col_types = FALSE)

# ==================== 3. 自动识别P列 ====================
p_col_candidates <- c("P.value", "P.Value", "Pvalue", "P", "p.value", "p_value", "pvalue")

find_p_col <- function(df) {
  hit <- p_col_candidates[p_col_candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  hit2 <- grep("^p(\\.|_|\\s)?value$|^p$|pvalue", colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hit2) > 0) return(hit2[1])
  stop("找不到P值列。请检查GWAS结果列名是否包含：P.value / P.Value / P 等。")
}

p_slope_col <- find_p_col(slope_df)
p_intcp_col <- find_p_col(intcp_df)

# ==================== 4. 校验SNP列 ====================
if (!("SNP" %in% colnames(slope_df)) || !("SNP" %in% colnames(intcp_df))) {
  stop("两个GWAS文件都必须包含 SNP 列用于合并。")
}

# QTL.csv 里也需要 SNP 列（若不是SNP列名，请在这里改）
if (!("SNP" %in% colnames(qtl_df))) {
  # 尝试自动找一下类似列名
  snp_hit <- grep("^snp$|marker|rs|^id$", colnames(qtl_df), ignore.case = TRUE, value = TRUE)
  if (length(snp_hit) > 0) {
    message(sprintf("QTL.csv 未找到列 'SNP'，使用 '%s' 代替", snp_hit[1]))
    qtl_df <- qtl_df %>% rename(SNP = all_of(snp_hit[1]))
  } else {
    stop("QTL.csv 中未找到 SNP 列（或类似列名）。请确认QTL.csv里SNP列叫什么。")
  }
}

# ==================== 5. 计算 -log10(P) 并合并 ====================
slope_df2 <- slope_df %>%
  transmute(
    SNP = SNP,
    slope_p = .data[[p_slope_col]],
    slope_logp = -log10(pmax(.data[[p_slope_col]], 1e-300))
  )

intcp_df2 <- intcp_df %>%
  transmute(
    SNP = SNP,
    intcp_p = .data[[p_intcp_col]],
    intcp_logp = -log10(pmax(.data[[p_intcp_col]], 1e-300))
  )

plot_df <- slope_df2 %>%
  inner_join(intcp_df2, by = "SNP") %>%
  filter(is.finite(slope_logp), is.finite(intcp_logp))

cat(sprintf("Merged SNP count: %d\n", nrow(plot_df)))

# ==================== 6. 标记QTL对应的点（红色变大） ====================
qtl_snps <- qtl_df %>%
  filter(!is.na(SNP)) %>%
  distinct(SNP) %>%
  pull(SNP)

plot_df <- plot_df %>%
  mutate(is_qtl = SNP %in% qtl_snps)

cat(sprintf("QTL SNPs in QTL.csv: %d\n", length(qtl_snps)))
cat(sprintf("QTL SNPs matched in GWAS merge: %d\n", sum(plot_df$is_qtl, na.rm = TRUE)))

# ==================== 7. 坐标范围（方形好看） ====================
max_lim <- max(plot_df$slope_logp, plot_df$intcp_logp, na.rm = TRUE)
max_lim <- ceiling(max_lim)
max_lim <- max(max_lim, 10)

# ==================== 8. 作图 ====================
# 先画非QTL点，再叠加QTL点（确保红点在上层）
p <- ggplot() +
  # 非QTL点（黑色）
  geom_point(
    data = plot_df %>% filter(!is_qtl),
    aes(x = intcp_logp, y = slope_logp),
    size = 1.1, alpha = 0.85, color = "black"
  ) +
  # QTL点（红色、变大）
  geom_point(
    data = plot_df %>% filter(is_qtl),
    aes(x = intcp_logp, y = slope_logp),
    size = 2.8, alpha = 0.95, color = "red"
  ) +
  # y=x 对角线
  geom_abline(intercept = 0, slope = 1, linewidth = 0.8, color = "grey70") +
  coord_equal(xlim = c(0, max_lim), ylim = c(0, max_lim), expand = FALSE) +
  labs(
    x = expression(-log[10](italic(P))~intercept),
    y = expression(-log[10](italic(P))~slope)
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 1.0, colour = "black"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(size = 12, colour = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==================== 9. 保存 ====================
ggsave(png_file, p, width = 5.2, height = 4.6, dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 5.2, height = 4.6, bg = "white")

cat("\nSaved:\n")
cat("  - ", png_file, "\n", sep = "")
cat("  - ", pdf_file, "\n", sep = "")
