# ==================== 加载所有必需包 ====================
suppressPackageStartupMessages({
  # 基础数据处理
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(grid)
  library(gridExtra)
  
  # 绘图核心包
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  
  # 特殊分析包
  library(vcfR)
  library(RColorBrewer)
  library(car)
  library(multcomp)
  
  # scales包用于alpha函数
  library(scales)
  
  # PDF转换包
  if (!requireNamespace("pdftools", quietly = TRUE)) {
    cat("安装pdftools包用于PDF转图片...\n")
    install.packages("pdftools", quiet = TRUE)
  }
  if (!requireNamespace("magick", quietly = TRUE)) {
    cat("安装magick包用于图片处理...\n")
    install.packages("magick", quiet = TRUE)
  }
  library(pdftools)
  library(magick)
  
  # 色盲模拟包
  if (!requireNamespace("dichromat", quietly = TRUE)) {
    cat("安装dichromat包用于色盲模拟...\n")
    install.packages("dichromat", quiet = TRUE)
  }
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    cat("安装colorspace包用于颜色分析...\n")
    install.packages("colorspace", quiet = TRUE)
  }
  library(dichromat)
  library(colorspace)
})

# ==================== 全局参数设置 ====================
base_dir <- "/mnt/7t_storage/zhangcl/TKW"
setwd(base_dir)
cat("当前工作目录:", getwd(), "\n")

# 创建输出目录结构
output_main_dir <- file.path(base_dir, "Final_Figure_Output")
if (!dir.exists(output_main_dir)) dir.create(output_main_dir, recursive = TRUE)

cat("创建输出目录结构...\n")
sub_dirs <- c(
  "Intermediate_Results",
  "Manhattan_Plots",
  "QTL_Plots",
  "Haplotype_Analysis",
  "Correlation_Analysis",
  "Final_Figures"
)

for (dir_name in sub_dirs) {
  dir_path <- file.path(output_main_dir, dir_name)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("创建目录:", dir_path, "\n")
  }
}

# 输出路径
final_pdf <- file.path(output_main_dir, "Final_Figures", "Final_Combined_Figure.pdf")
final_png <- file.path(output_main_dir, "Final_Figures", "Final_Combined_Figure.png")
final_tiff <- file.path(output_main_dir, "Final_Figures", "Final_Combined_Figure.tiff")

# ==================== 环境列表和形状映射 ====================
env_dirs <- c("2024NY", "2024PY", "2024YL", "2024ZMD",
              "2025YLW", "2025YLZ", "2025ZMDW", "2025ZMDZ")
env_shape_values <- setNames(c(1, 2, 3, 4, 5, 6, 7, 8), env_dirs)

# ==================== 统一参数设置 ====================
label_params <- list(
  label_size = 8.5,
  label_fontface = "bold",
  label_hjust = 1.25,
  label_vjust_top = 1.15,
  label_x = -Inf,
  label_y = Inf,
  ab_label_y = 0.95,
  cde_label_y = 0.95,
  qtl_label_fontface = "bold.italic"
)

legend_params <- list(
  title_size = 11.5,
  text_size  = 10,
  title_face = "bold",
  text_face  = "plain"
)

axis_title_size <- 14
axis_text_size <- 12
axis_line_width <- 0.8
axis_tick_width <- 0.6

# d/e图X轴刻度 - 使用用户提供的 DTR&mean&PostFlowering1_14 数值
x_ticks_de <- c(13.64333333, 15.70304348, 14.65909091, 14.16857143,
                13.35590909, 13.41952381, 14.58863636, 14.73227273)
x_tick_labels_de <- format(round(x_ticks_de, 1), nsmall = 1)  # 保留1位小数

panel_margin_de <- margin(14, 40, 14, 30)
panel_margin_c <- panel_margin_de

# ==================== 极致统一：panel label + legend + axis label 的工具函数 ====================
panel_label_layer <- function(letter) {
  annotate(
    "text",
    x = label_params$label_x,
    y = label_params$label_y,
    label = letter,
    size = label_params$label_size,
    fontface = label_params$label_fontface,
    hjust = label_params$label_hjust,
    vjust = label_params$label_vjust_top,
    color = "black"
  )
}

theme_legend_std <- function(position = "right", direction = "vertical") {
  theme(
    legend.position = position,
    legend.direction = direction,
    legend.title = element_text(face = legend_params$title_face, size = legend_params$title_size),
    legend.text  = element_text(face = legend_params$text_face,  size = legend_params$text_size),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
}

theme_axis_std <- function(base_size = 13, plot_margin = margin(14, 14, 14, 30)) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line  = element_line(color = "black", linewidth = axis_line_width),
      axis.ticks = element_line(color = "black", linewidth = axis_tick_width),
      axis.title = element_text(face = "bold", size = axis_title_size),
      axis.text  = element_text(size = axis_text_size, color = "black"),
      plot.margin = plot_margin,
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

# ==================== 辅助函数：保存中间结果 ====================
save_intermediate_data <- function(data, filename, description = "") {
  if (!missing(description) && description != "") {
    cat("保存中间结果:", description, "\n")
  }
  file_path <- file.path(output_main_dir, "Intermediate_Results", filename)
  if (is.data.frame(data)) {
    write.csv(data, file_path, row.names = FALSE)
    cat("✓ 保存到:", file_path, "\n")
  } else if (is.list(data)) {
    saveRDS(data, file = sub("\\.csv$", ".rds", file_path))
    cat("✓ 保存列表为RDS格式:", sub("\\.csv$", ".rds", file_path), "\n")
  } else {
    cat("✗ 无法保存的数据类型:", class(data), "\n")
  }
}

# ==================== 辅助函数：PDF转图片 ====================
convert_pdf_to_images <- function(pdf_path, output_dir, base_name, formats = c("png", "tiff")) {
  cat("\n转换PDF到其他图片格式...\n")
  if (!file.exists(pdf_path)) {
    cat("✗ PDF文件不存在:", pdf_path, "\n")
    return(FALSE)
  }
  pdf_img <- image_read_pdf(pdf_path, density = 300)
  results <- list()
  for (fmt in formats) {
    output_path <- file.path(output_dir, paste0(base_name, ".", fmt))
    tryCatch({
      if (fmt == "png") {
        image_write(pdf_img, output_path, format = "png", density = 300, quality = 100)
        cat("✓ PNG保存完成:", output_path, " (300 DPI)\n")
        results[[fmt]] <- output_path
      } else if (fmt == "tiff") {
        image_write(pdf_img, output_path, format = "tiff", density = 300, compression = "lzw")
        cat("✓ TIFF保存完成:", output_path, " (300 DPI, LZW压缩)\n")
        results[[fmt]] <- output_path
      } else if (fmt == "jpeg") {
        image_write(pdf_img, output_path, format = "jpeg", density = 300, quality = 95)
        cat("✓ JPEG保存完成:", output_path, " (300 DPI, 质量95%)\n")
        results[[fmt]] <- output_path
      }
    }, error = function(e) {
      cat("✗ 转换", fmt, "失败:", e$message, "\n")
    })
  }
  info_path <- file.path(output_dir, paste0(base_name, "_image_info.txt"))
  writeLines(c(
    paste("PDF源文件:", pdf_path),
    paste("生成时间:", Sys.time()),
    paste("图像宽度:", image_info(pdf_img)$width, "像素"),
    paste("图像高度:", image_info(pdf_img)$height, "像素"),
    paste("分辨率: 300 DPI"),
    paste("生成的文件:", paste(names(results), collapse = ", "))
  ), info_path)
  return(results)
}

# ==================== 辅助函数：检查配色对色盲友好性 ====================
check_colorblind_friendly <- function(color_vector, plot_type = "Manhattan") {
  cat("\n检查配色对色盲友好性...\n")
  colors <- unname(color_vector)
  cat("模拟色盲视觉效果：\n")
  deut_colors <- dichromat(colors, type = "deutan")
  cat("  - 红绿色盲（绿色盲）下：\n")
  color_comparison <- data.frame(
    Original = colors,
    Deuteranopia = deut_colors,
    stringsAsFactors = FALSE
  )
  print(color_comparison)
  prot_colors <- dichromat(colors, type = "protan")
  cat("  - 红色盲下：\n")
  color_comparison$Protanopia <- prot_colors
  print(color_comparison[, c("Original", "Protanopia")])
  trit_colors <- dichromat(colors, type = "tritan")
  cat("  - 蓝黄色盲下：\n")
  color_comparison$Tritanopia <- trit_colors
  print(color_comparison[, c("Original", "Tritanopia")])
  cat("\n颜色差异度分析：\n")
  rgb_matrix <- col2rgb(colors)
  color_distances <- as.matrix(dist(t(rgb_matrix)))
  diag(color_distances) <- NA
  min_dist <- min(color_distances, na.rm = TRUE)
  mean_dist <- mean(color_distances, na.rm = TRUE)
  cat(sprintf("  最小颜色距离：%.1f (建议 > 50)\n", min_dist))
  cat(sprintf("  平均颜色距离：%.1f\n", mean_dist))
  if (min_dist < 50) {
    cat("  ⚠️  警告：有些颜色在色盲视角下可能难以区分\n")
  } else {
    cat("  ✅ 颜色差异度良好\n")
  }
  cat("\n推荐的红绿色盲友好配色（使用分类颜色）：\n")
  if (plot_type == "Manhattan") {
    recommended <- list(
      "Set3" = c(
        brewer.pal(8, "Set3")[1:7],
        brewer.pal(8, "Set2")[1:7],
        brewer.pal(8, "Set1")[1:7],
        "#999999"
      ),
      "Paired" = c(
        brewer.pal(12, "Paired")[c(1,3,5,7,9,11,2,4,6,8,10,12,1,3,5,7,9,11,2,4,6)],
        "#999999"
      ),
      "Dark2" = c(
        brewer.pal(8, "Dark2")[1:7],
        brewer.pal(8, "Accent")[1:7],
        brewer.pal(8, "Set1")[1:7],
        "#999999"
      )
    )
    for (scheme_name in names(recommended)) {
      cat(sprintf("  %s: %d种颜色\n", scheme_name, length(recommended[[scheme_name]])))
    }
  }
  return(list(
    original_colors = colors,
    deuteranopia = deut_colors,
    protanopia = prot_colors,
    tritanopia = trit_colors,
    min_distance = min_dist,
    mean_distance = mean_dist,
    color_comparison = color_comparison
  ))
}

# ==================== 通用函数：从VCF构建单倍型数据 ====================
build_haplotype_data <- function(qtl_file, vcf_file) {
  cat("读取QTL标记信息...\n")
  qtl_df <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)

  snp_col <- NULL
  if ("SNP" %in% colnames(qtl_df)) {
    snp_col <- "SNP"
  } else if ("Marker" %in% colnames(qtl_df)) {
    snp_col <- "Marker"
  } else {
    snp_col <- grep("SNP|snp|Marker|marker|ID|id", colnames(qtl_df),
                    value = TRUE, ignore.case = TRUE)[1]
  }
  if (is.null(snp_col)) stop("无法在QTL.csv中找到SNP/Marker列")

  qtl_snps <- unique(qtl_df[[snp_col]])
  qtl_snps <- qtl_snps[!is.na(qtl_snps)]
  qtl_snps <- qtl_snps[qtl_snps != ""]
  if (length(qtl_snps) == 0) {
    stop("QTL.csv中没有找到有效的QTL SNPs")
  }
  cat(sprintf("找到 %d 个QTL标记\n", length(qtl_snps)))

  cat("读取VCF文件...\n")
  vcf_data <- read.vcfR(vcf_file, verbose = FALSE)

  vcf_snp_ids <- vcf_data@fix[, "ID"]
  matched_indices <- which(vcf_snp_ids %in% qtl_snps)

  if (length(matched_indices) == 0) {
    cat("通过ID没有匹配到。尝试通过CHROM+POS匹配...\n")
    chr_col <- grep("CHR|Chr|chr|chromosome|Chromosome", colnames(qtl_df),
                    value = TRUE, ignore.case = TRUE)[1]
    pos_col <- grep("POS|Pos|pos|position|Position", colnames(qtl_df),
                    value = TRUE, ignore.case = TRUE)[1]
    if (!is.null(chr_col) && !is.null(pos_col)) {
      qtl_pos_ids <- paste(qtl_df[[chr_col]], qtl_df[[pos_col]], sep = "_")
      vcf_pos_ids <- paste(vcf_data@fix[, "CHROM"], vcf_data@fix[, "POS"], sep = "_")
      matched_indices <- which(vcf_pos_ids %in% qtl_pos_ids)
    }
  }
  if (length(matched_indices) == 0) stop("在VCF中没有找到匹配的SNP标记。")
  cat(sprintf("匹配到 %d 个标记\n", length(matched_indices)))

  extracted_vcf <- vcf_data[matched_indices, ]

  matched_markers <- data.frame(
    VCF_ID = vcf_data@fix[matched_indices, "ID"],
    CHROM = vcf_data@fix[matched_indices, "CHROM"],
    POS = vcf_data@fix[matched_indices, "POS"],
    Matched_from_QTL = "Yes",
    stringsAsFactors = FALSE
  )
  save_intermediate_data(matched_markers, "Matched_QTL_Markers_VCF.csv", "VCF中匹配的QTL标记")

  gt_matrix <- extract.gt(extracted_vcf, element = "GT")

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
    return(out)
  }

  gt_numeric <- convert_gt_to_numeric(gt_matrix)

  genotype_summary <- data.frame(
    Total_samples = ncol(gt_numeric),
    Total_markers = nrow(gt_numeric),
    Missing_rate = sum(is.na(gt_numeric)) / (nrow(gt_numeric) * ncol(gt_numeric)),
    Homozygous_0 = sum(gt_numeric == 0, na.rm = TRUE),
    Heterozygous_1 = sum(gt_numeric == 1, na.rm = TRUE),
    Homozygous_2 = sum(gt_numeric == 2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(genotype_summary, "Genotype_Matrix_Summary.csv", "基因型矩阵摘要")

  cat("构建单倍型...\n")
  haplotype_df <- as.data.frame(t(gt_numeric), stringsAsFactors = FALSE)
  haplotype_df$Sample <- rownames(haplotype_df)

  snp_cols <- setdiff(colnames(haplotype_df), "Sample")

  haplotype_df$has_heterozygous <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
    any(x == 1, na.rm = TRUE)
  })
  haplotype_df_original <- haplotype_df
  haplotype_df <- haplotype_df[!haplotype_df$has_heterozygous, ]
  cat(sprintf("过滤杂合子后剩余样本数: %d (原始: %d)\n", nrow(haplotype_df), nrow(haplotype_df_original)))

  haplotype_df$Haplotype <- apply(haplotype_df[, snp_cols, drop = FALSE], 1, function(x) {
    paste(ifelse(is.na(x), "N", x), collapse = "-")
  })

  cat("保存品种名-基因型-单倍型关系...\n")
  genotype_haplotype_data <- haplotype_df[, c("Sample", snp_cols, "Haplotype")]
  for (snp in snp_cols) {
    genotype_haplotype_data[[paste0(snp, "_genotype")]] <- sapply(genotype_haplotype_data[[snp]], function(x) {
      if (is.na(x)) return("Missing")
      if (x == 0) return("Homozygous_ref")
      if (x == 1) return("Heterozygous")
      if (x == 2) return("Homozygous_alt")
      return("Other")
    })
  }
  colnames(genotype_haplotype_data) <- c("Sample_ID", snp_cols, "Haplotype",
                                         paste0(snp_cols, "_genotype_desc"))

  haplotype_output_dir <- file.path(output_main_dir, "Haplotype_Analysis")
  if (!dir.exists(haplotype_output_dir)) dir.create(haplotype_output_dir, recursive = TRUE)

  genotype_haplotype_file <- file.path(haplotype_output_dir, "Sample_Genotype_Haplotype_Relationship.csv")
  write.csv(genotype_haplotype_data, genotype_haplotype_file, row.names = FALSE)
  cat("✓ 品种名-基因型-单倍型关系保存:", genotype_haplotype_file, "\n")
  save_intermediate_data(genotype_haplotype_data, "Sample_Genotype_Haplotype_Relationship.csv",
                         "品种名-基因型-单倍型关系")

  haplotype_counts <- as.data.frame(table(haplotype_df$Haplotype), stringsAsFactors = FALSE)
  colnames(haplotype_counts) <- c("Haplotype", "Frequency")
  haplotype_counts <- haplotype_counts[order(haplotype_counts$Frequency, decreasing = TRUE), ]
  haplotype_counts$Proportion <- haplotype_counts$Frequency / sum(haplotype_counts$Frequency)
  haplotype_counts$Cumulative_Proportion <- cumsum(haplotype_counts$Proportion)

  all_haps <- haplotype_counts$Haplotype
  n_haps <- length(all_haps)

  roman_nums <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
  if (n_haps <= length(roman_nums)) {
    haps_named <- paste0("Hap", roman_nums[1:n_haps])
  } else {
    haps_named <- paste0("Hap", 1:n_haps)
  }
  names(haps_named) <- all_haps
  cat(sprintf("检测到 %d 个单倍型 | 全部展示: %s\n",
              n_haps, paste(haps_named, collapse = ", ")))

  haplotype_counts$Haplotype_Named <- haps_named[match(haplotype_counts$Haplotype, names(haps_named))]
  save_intermediate_data(haplotype_counts, "Haplotype_Frequency_Table.csv", "单倍型频率表")

  # 清理样本名：去除空格、引号，转换为大写
  haplotype_df$Sample <- toupper(trimws(gsub('"', '', haplotype_df$Sample)))

  return(list(
    haplotype_df = haplotype_df,
    all_haps = all_haps,
    haps_named = haps_named
  ))
}

# ==================== 1. 定义图a/b (曼哈顿图) 绘制函数 ====================
plot_ab <- function() {
  # ============ 参数设置 ============
  slope_gwas_file <- "FW_GWAS5_Results/GWAS_All_Traits/GWAS_Slope_para/GAPIT.Association.GWAS_Results.GWAS_Slope_para.MLM.Trait(NYC).csv"
  slope_qtl_file  <- "QTLslope.csv"
  intcp_gwas_file <- "FW_GWAS5_Results/GWAS_All_Traits/GWAS_Intcp_para_adj/GAPIT.Association.GWAS_Results.GWAS_Intcp_para_adj.MLM.Trait(NYC).csv"
  intcp_qtl_file  <- "QTLIntcp.csv"
  p_thresh <- 0.0001
  log_thresh <- -log10(p_thresh)

  # ============ 染色体映射 ============
  chr_map <- data.frame(
    Chr = 1:22,
    ChrName = c(
      "Chr1A", "Chr1B", "Chr1D",
      "Chr2A", "Chr2B", "Chr2D",
      "Chr3A", "Chr3B", "Chr3D",
      "Chr4A", "Chr4B", "Chr4D",
      "Chr5A", "Chr5B", "Chr5D",
      "Chr6A", "Chr6B", "Chr6D",
      "Chr7A", "Chr7B", "Chr7D",
      "ChrUN"
    ),
    stringsAsFactors = FALSE
  )

  # ============ 简化颜色方案（3种颜色+灰色） - 色盲友好版本 ============
  set2_colors <- brewer.pal(3, "Set2")
  color_mapping <- c(
    "Chr1A" = set2_colors[1], "Chr2A" = set2_colors[1], "Chr3A" = set2_colors[1],
    "Chr4A" = set2_colors[1], "Chr5A" = set2_colors[1], "Chr6A" = set2_colors[1],
    "Chr7A" = set2_colors[1],
    "Chr1B" = set2_colors[2], "Chr2B" = set2_colors[2], "Chr3B" = set2_colors[2],
    "Chr4B" = set2_colors[2], "Chr5B" = set2_colors[2], "Chr6B" = set2_colors[2],
    "Chr7B" = set2_colors[2],
    "Chr1D" = set2_colors[3], "Chr2D" = set2_colors[3], "Chr3D" = set2_colors[3],
    "Chr4D" = set2_colors[3], "Chr5D" = set2_colors[3], "Chr6D" = set2_colors[3],
    "Chr7D" = set2_colors[3],
    "ChrUN" = "#9E9E9E"
  )

  color_groups <- data.frame(
    Chromosome = names(color_mapping),
    Color = color_mapping,
    Group = ifelse(grepl("A$", names(color_mapping)), "A组",
                   ifelse(grepl("B$", names(color_mapping)), "B组",
                          ifelse(grepl("D$", names(color_mapping)), "D组", "UN"))),
    stringsAsFactors = FALSE
  )

  cat("\n检查曼哈顿图配色色盲友好性...\n")
  cat("使用Set2调色板的前3种颜色 + 灰色，对色盲友好\n")
  cat("颜色分配：\n")
  cat("  A组染色体（1A-7A）：", set2_colors[1], "\n")
  cat("  B组染色体（1B-7B）：", set2_colors[2], "\n")
  cat("  D组染色体（1D-7D）：", set2_colors[3], "\n")
  cat("  ChrUN：#9E9E9E（灰色）\n")

  color_check_result <- check_colorblind_friendly(color_mapping, plot_type = "Manhattan")
  save_intermediate_data(
    color_check_result$color_comparison,
    "Manhattan_Colorblind_Check_Results.csv",
    "曼哈顿图色盲友好性检查结果"
  )
  save_intermediate_data(color_groups, "Simplified_Color_Scheme.csv", "简化颜色方案")

  # ============ 定义函数：处理单个GWAS数据 ============
  process_gwas_data <- function(gwas_file, qtl_file, plot_title) {
    cat("处理数据:", plot_title, "\n")
    gwas <- fread(gwas_file, data.table = FALSE)
    qtl  <- fread(qtl_file, data.table = FALSE)
    colnames(gwas) <- gsub(" ", "", colnames(gwas))
    colnames(qtl)  <- gsub(" ", "", colnames(qtl))

    required_gwas_cols <- c("SNP", "Chr", "Pos", "P.value")
    required_qtl_cols <- c("QTLname", "SNP")

    if (!"P.value" %in% colnames(gwas)) {
      cat(plot_title, "文件缺少P.value列，尝试计算P值...\n")
      if ("tValue" %in% colnames(gwas) && "DF" %in% colnames(gwas)) {
        gwas$P.value <- 2 * pt(-abs(gwas$tValue), df = gwas$DF)
        cat(sprintf("计算了P值：范围 = [%.2e, %.2e]\n",
                    min(gwas$P.value, na.rm = TRUE),
                    max(gwas$P.value, na.rm = TRUE)))
      } else {
        stop(paste("GWAS文件缺少必要列，需要:", paste(required_gwas_cols, collapse = ", ")))
      }
    }

    if (!all(required_qtl_cols %in% colnames(qtl))) {
      stop(paste("QTL文件缺少必要列，需要:", paste(required_qtl_cols, collapse = ", ")))
    }

    gwas <- gwas[gwas$Chr %in% 1:22, ]
    gwas$ChrName <- chr_map$ChrName[match(gwas$Chr, chr_map$Chr)]
    gwas$Pos <- as.numeric(gwas$Pos)
    gwas$P.value <- as.numeric(gwas$P.value)

    gwas$ChrGroup <- substr(gwas$ChrName, nchar(gwas$ChrName), nchar(gwas$ChrName))
    gwas$ChrGroup <- ifelse(gwas$ChrName == "ChrUN", "UN", gwas$ChrGroup)

    gwas$logP <- -log10(gwas$P.value)
    gwas <- gwas[order(gwas$Chr, gwas$Pos), ]

    gwas$chr_index <- ave(1:nrow(gwas), gwas$Chr, FUN = seq_along)
    gwas$chr_max <- ave(gwas$Pos, gwas$Chr, FUN = function(x) max(x, na.rm = TRUE))
    gwas$chr_min <- ave(gwas$Pos, gwas$Chr, FUN = function(x) min(x, na.rm = TRUE))

    chr_lengths <- tapply(gwas$Pos, gwas$Chr, function(x) max(x, na.rm = TRUE))
    chr_offsets_df <- data.frame(
      Chr = as.numeric(names(chr_lengths)),
      chr_len = as.numeric(chr_lengths),
      stringsAsFactors = FALSE
    )
    chr_offsets_df <- chr_offsets_df[order(chr_offsets_df$Chr), ]
    chr_offsets_df$offset <- c(0, cumsum(chr_offsets_df$chr_len[-nrow(chr_offsets_df)]))
    chr_offsets_df$cumulative <- chr_offsets_df$offset + chr_offsets_df$chr_len / 2

    gwas$offset <- chr_offsets_df$offset[match(gwas$Chr, chr_offsets_df$Chr)]
    gwas$CumPos <- gwas$Pos + gwas$offset

    chr_centers <- aggregate(CumPos ~ Chr + ChrName, data = gwas, FUN = median, na.rm = TRUE)
    chr_centers <- chr_centers[order(chr_centers$Chr), ]
    colnames(chr_centers)[3] <- "center"

    ann <- merge(qtl[, c("QTLname", "SNP")], gwas, by = "SNP")
    ann <- unique(ann)
    if (nrow(ann) > 0) {
      ann$point_color <- color_mapping[ann$ChrName]
    }

    # 计算Genomic inflation factor (lambda) 但不用于绘图
    p_for_lambda <- pmax(pmin(gwas$P.value, 1 - 1e-15), 1e-300)
    chisq_vals <- qchisq(1 - p_for_lambda, df = 1)
    lambda_gc <- median(chisq_vals, na.rm = TRUE) / qchisq(0.5, df = 1)
    cat(sprintf("%s的基因组膨胀因子λ = %.3f (范围: %.3f~%.3f)\n",
                plot_title, lambda_gc, min(chisq_vals, na.rm = TRUE) / qchisq(0.5, df = 1),
                max(chisq_vals, na.rm = TRUE) / qchisq(0.5, df = 1)))

    gwas_info <- data.frame(
      Parameter = plot_title,
      Total_SNPs = nrow(gwas),
      Max_logP = max(gwas$logP, na.rm = TRUE),
      Min_Pvalue = min(gwas$P.value, na.rm = TRUE),
      Mean_logP = mean(gwas$logP, na.rm = TRUE),
      Significant_SNPs = sum(gwas$P.value < p_thresh, na.rm = TRUE),
      QTL_markers = nrow(ann),
      Lambda = lambda_gc,
      Lambda_formatted = sprintf("%.3f", lambda_gc),
      stringsAsFactors = FALSE
    )
    save_intermediate_data(gwas_info,
                           paste0("GWAS_Summary_", gsub(" ", "_", plot_title), ".csv"),
                           paste0(plot_title, " GWAS汇总统计"))

    return(list(
      gwas = gwas,
      ann = ann,
      chr_centers = chr_centers,
      chr_offsets = chr_offsets_df,
      plot_title = plot_title,
      max_logP = max(gwas$logP, na.rm = TRUE),
      max_CumPos = max(gwas$CumPos, na.rm = TRUE),
      lambda = lambda_gc,
      lambda_formatted = sprintf("%.3f", lambda_gc)
    ))
  }

  # ============ 处理两个数据集 ============
  intercept_data <- process_gwas_data(
    gwas_file = intcp_gwas_file,
    qtl_file = intcp_qtl_file,
    plot_title = "Intercept"
  )
  slope_data <- process_gwas_data(
    gwas_file = slope_gwas_file,
    qtl_file = slope_qtl_file,
    plot_title = "Slope"
  )

  # ============ 计算统一的y轴上限 ============
  raw_max_logP <- max(intercept_data$max_logP, slope_data$max_logP, na.rm = TRUE)
  max_logP_both <- ceiling(raw_max_logP) + 1
  y_break_by <- ifelse(max_logP_both <= 10, 1, 2)
  y_breaks <- seq(0, max_logP_both, by = y_break_by)
  max_CumPos_both <- max(intercept_data$max_CumPos, slope_data$max_CumPos, na.rm = TRUE)

  # ============ 定义函数：创建单个曼哈顿图（已移除lambda标注） ============
  create_manhattan_plot <- function(data_list, subplot_label = "a", show_title = TRUE) {
    gwas <- data_list$gwas
    ann <- data_list$ann
    chr_centers <- data_list$chr_centers
    plot_title <- data_list$plot_title

    p <- ggplot(gwas, aes(x = CumPos, y = logP)) +
      geom_point(aes(color = ChrGroup), size = 0.8, alpha = 0.75) +
      scale_color_manual(
        values = c("A" = set2_colors[1], "B" = set2_colors[2], "D" = set2_colors[3], "UN" = "#9E9E9E"),
        guide = "none"
      ) +
      geom_hline(yintercept = log_thresh,
                 linetype = "dashed",
                 color = "#E63946",
                 linewidth = 0.6) +
      annotate("text",
               x = max_CumPos_both * 0.98,
               y = log_thresh + 0.2,
               label = paste0("P = ", p_thresh),
               color = "#E63946",
               size = 3,
               hjust = 1,
               vjust = -0.2) +
      scale_x_continuous(
        name = NULL,
        breaks = chr_centers$center,
        labels = chr_centers$ChrName,
        expand = expansion(mult = c(0.005, 0.005)),
        limits = c(0, max_CumPos_both)
      ) +
      scale_y_continuous(
        name = expression(-log[10](italic(P))),
        limits = c(0, max_logP_both),
        breaks = y_breaks,
        expand = expansion(mult = c(0, 0))
      ) +
      panel_label_layer(subplot_label) +
      theme_axis_std(base_size = 13, plot_margin = margin(5, 8, 8, 30)) +
      theme(
        axis.ticks.length = unit(1.5, "mm"),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = axis_text_size,
          color = "black"
        ),
        axis.title.y = element_text(margin = margin(r = 5)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank()
      ) +
      coord_cartesian(xlim = c(0, max_CumPos_both), clip = "off")

    if (show_title) {
      p <- p +
        annotate("text",
                 x = max_CumPos_both * 0.01,
                 y = max_logP_both * 0.85,
                 label = plot_title,
                 size = 4.5,
                 fontface = "bold",
                 hjust = 0,
                 vjust = 0.5)
    }

    if (nrow(ann) > 0) {
      p <- p +
        geom_point(
          data = ann,
          aes(x = CumPos, y = logP),
          shape = 21,
          fill = color_mapping[ann$ChrName],
          color = "black",
          stroke = 0.35,
          size = 3.5,
          alpha = 0.95
        ) +
        geom_text_repel(
          data = ann,
          aes(x = CumPos, y = logP, label = QTLname),
          size = 3.2,
          fontface = label_params$qtl_label_fontface,
          color = "black",
          segment.color = "gray50",
          segment.size = 0.3,
          segment.alpha = 0.6,
          box.padding = 0.3,
          point.padding = 0.2,
          min.segment.length = 0.1,
          force = 0.8,
          nudge_y = 0.1,
          max.overlaps = Inf,
          family = "sans"
        )
    }
    return(p)
  }

  intercept_plot <- create_manhattan_plot(
    data_list = intercept_data,
    subplot_label = "a",
    show_title = TRUE
  )
  slope_plot <- create_manhattan_plot(
    data_list = slope_data,
    subplot_label = "b",
    show_title = TRUE
  )

  combined_plot <- intercept_plot + slope_plot + plot_layout(ncol = 2, widths = c(1, 1))

  # 单独保存曼哈顿图
  manhattan_output_dir <- file.path(output_main_dir, "Manhattan_Plots")
  if (!dir.exists(manhattan_output_dir)) dir.create(manhattan_output_dir, recursive = TRUE)

  intercept_pdf <- file.path(manhattan_output_dir, "Manhattan_Intercept.pdf")
  ggsave(intercept_pdf, intercept_plot, width = 12, height = 6, units = "in", bg = "white")
  cat("✓ Intercept曼哈顿图保存:", intercept_pdf, "\n")

  slope_pdf <- file.path(manhattan_output_dir, "Manhattan_Slope.pdf")
  ggsave(slope_pdf, slope_plot, width = 12, height = 6, units = "in", bg = "white")
  cat("✓ Slope曼哈顿图保存:", slope_pdf, "\n")

  convert_pdf_to_images(intercept_pdf, manhattan_output_dir, "Manhattan_Intercept", c("png", "tiff"))
  convert_pdf_to_images(slope_pdf, manhattan_output_dir, "Manhattan_Slope", c("png", "tiff"))

  # ============ 生成Supplementary QQ Plots (S1) ============
  cat("\n生成Supplementary QQ Plots...\n")
  plot_qq_with_lambda <- function(p_values, lambda, title = "QQ plot") {
    p <- pmax(pmin(p_values, 1 - 1e-15), 1e-300)
    p <- p[is.finite(p) & !is.na(p)]
    p <- sort(p)
    n <- length(p)
    exp_p <- (1:n) / (n + 1)
    df <- data.frame(
      exp = -log10(exp_p),
      obs = -log10(p)
    )
    max_lim <- ceiling(max(df$exp, df$obs, na.rm = TRUE)) + 1
    lambda_label_qq <- paste0("lambda = ", round(lambda, 3))

    p_plot <- ggplot(df, aes(x = exp, y = obs)) +
      geom_point(size = 0.9, alpha = 0.8, color = "#2171B5") +
      geom_abline(intercept = 0, slope = 1, linewidth = 0.7, color = "grey60") +
      coord_equal(xlim = c(0, max_lim), ylim = c(0, max_lim), expand = FALSE) +
      labs(x = expression(Expected~~-log[10](italic(P))),
           y = expression(Observed~~-log[10](italic(P)))) +
      theme_axis_std(base_size = 12, plot_margin = margin(8, 8, 8, 8)) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11)
      ) +
      ggtitle(title) +
      annotate(
        "text",
        x = max_lim * 0.05,
        y = max_lim * 0.95,
        label = lambda_label_qq,
        hjust = 0, vjust = 1,
        size = 3.2,
        color = "grey30"
      )
    return(p_plot)
  }

  qq_intcp <- plot_qq_with_lambda(
    p_values = intercept_data$gwas$P.value,
    lambda = intercept_data$lambda,
    title = "Intercept"
  ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = "S1a",
      size = 8.5,
      fontface = "bold",
      hjust = 1.25,
      vjust = 1.15,
      color = "black"
    )

  qq_slope <- plot_qq_with_lambda(
    p_values = slope_data$gwas$P.value,
    lambda = slope_data$lambda,
    title = "Slope"
  ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = "S1b",
      size = 8.5,
      fontface = "bold",
      hjust = 1.25,
      vjust = 1.15,
      color = "black"
    )

  qq_supp <- qq_intcp + qq_slope + plot_layout(ncol = 2, widths = c(1, 1))

  supp_dir <- file.path(output_main_dir, "Final_Figures")
  if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

  supp_pdf <- file.path(supp_dir, "Supplementary_Figure_S1_QQplots.pdf")
  ggsave(supp_pdf, qq_supp, width = 10, height = 5, units = "in", bg = "white")
  cat("✓ Supplementary QQ plots保存:", supp_pdf, "\n")
  convert_pdf_to_images(supp_pdf, supp_dir, "Supplementary_Figure_S1_QQplots", c("png", "tiff"))

  lambda_stats <- data.frame(
    Parameter = c("Intercept", "Slope"),
    Lambda = c(intercept_data$lambda, slope_data$lambda),
    Lambda_formatted = c(sprintf("%.3f", intercept_data$lambda),
                         sprintf("%.3f", slope_data$lambda)),
    Max_logP = c(intercept_data$max_logP, slope_data$max_logP),
    SNPs_count = c(nrow(intercept_data$gwas), nrow(slope_data$gwas)),
    P_threshold = p_thresh,
    stringsAsFactors = FALSE
  )
  save_intermediate_data(lambda_stats, "Genomic_Inflation_Factors.csv", "基因组膨胀因子统计")

  color_info <- data.frame(
    Chromosome = names(color_mapping),
    Color_Hex = color_mapping,
    Color_Group = c(rep("A组 (Set2颜色1)", 7), rep("B组 (Set2颜色2)", 7), rep("D组 (Set2颜色3)", 7), "UN (灰色)"),
    Colorblind_Safe = "是（Set2调色板，色盲友好）",
    stringsAsFactors = FALSE
  )
  save_intermediate_data(color_info, "Manhattan_Color_Scheme.csv", "曼哈顿图配色方案")

  cat("\n生成色盲模拟预览图数据...\n")
  color_df <- data.frame(
    Color_Name = names(color_mapping),
    Original = color_mapping,
    Deuteranopia = dichromat(color_mapping, "deutan"),
    Protanopia = dichromat(color_mapping, "protan"),
    Tritanopia = dichromat(color_mapping, "tritan")
  )
  save_intermediate_data(color_df, "Colorblind_Simulation_Preview.csv", "色盲模拟预览数据")

  return(combined_plot)
}

# ==================== 2. 定义图c (Slope vs Intercept) 绘制函数 ====================
plot_c <- function() {
  slope_gwas_file <- file.path(
    base_dir,
    "FW_GWAS5_Results/GWAS_All_Traits/GWAS_Slope_para/GAPIT.Association.GWAS_Results.GWAS_Slope_para.MLM.Trait(NYC).csv"
  )
  intcp_gwas_file <- file.path(
    base_dir,
    "FW_GWAS5_Results/GWAS_All_Traits/GWAS_Intcp_para_adj/GAPIT.Association.GWAS_Results.GWAS_Intcp_para_adj.MLM.Trait(NYC).csv"
  )
  qtl_file <- file.path(base_dir, "QTL.csv")

  cat("读取slope GWAS数据:\n", slope_gwas_file, "\n")
  slope_df <- read_csv(slope_gwas_file, show_col_types = FALSE)
  cat("读取intercept GWAS数据:\n", intcp_gwas_file, "\n")
  intcp_df <- read_csv(intcp_gwas_file, show_col_types = FALSE)
  cat("读取QTL列表:\n", qtl_file, "\n")
  qtl_df <- read_csv(qtl_file, show_col_types = FALSE)

  colnames(slope_df) <- gsub(" ", "", colnames(slope_df))
  colnames(intcp_df) <- gsub(" ", "", colnames(intcp_df))
  colnames(qtl_df) <- gsub(" ", "", colnames(qtl_df))

  find_p_col <- function(df) {
    p_col_candidates <- c("P.value", "P.Value", "Pvalue", "P", "p.value", "p_value", "pvalue")
    hit <- p_col_candidates[p_col_candidates %in% colnames(df)]
    if (length(hit) > 0) return(hit[1])
    hit2 <- grep("^p(\\.|_|\\s)?value$|^p$|pvalue", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hit2) > 0) return(hit2[1])
    stop("找不到P值列。请检查GWAS结果列名是否包含：P.value / P.Value / P 等。")
  }

  if (!"P.value" %in% colnames(slope_df)) {
    cat("Slope GWAS文件缺少P.value列，尝试计算P值...\n")
    if ("tValue" %in% colnames(slope_df) && "DF" %in% colnames(slope_df)) {
      slope_df$P.value <- 2 * pt(-abs(slope_df$tValue), df = slope_df$DF)
      cat(sprintf("计算了Slope P值：范围 = [%.2e, %.2e]\n",
                  min(slope_df$P.value, na.rm = TRUE),
                  max(slope_df$P.value, na.rm = TRUE)))
    } else {
      p_slope_col <- find_p_col(slope_df)
      if (p_slope_col %in% colnames(slope_df)) {
        slope_df$P.value <- slope_df[[p_slope_col]]
      } else {
        stop("Slope GWAS文件缺少P值列且无法计算。")
      }
    }
  }

  if (!"P.value" %in% colnames(intcp_df)) {
    p_intcp_col <- find_p_col(intcp_df)
    if (p_intcp_col %in% colnames(intcp_df)) {
      intcp_df$P.value <- intcp_df[[p_intcp_col]]
      cat("使用Intercept P值列:", p_intcp_col, "\n")
    } else {
      stop("Intercept GWAS文件缺少P值列。")
    }
  }

  if (!("SNP" %in% colnames(slope_df)) || !("SNP" %in% colnames(intcp_df))) {
    stop("两个GWAS文件都必须包含 SNP 列用于合并。")
  }

  if (!("SNP" %in% colnames(qtl_df))) {
    snp_hit <- grep("^snp$|marker|rs|^id$", colnames(qtl_df), ignore.case = TRUE, value = TRUE)
    if (length(snp_hit) > 0) {
      message(sprintf("QTL.csv 未找到列 'SNP'，使用 '%s' 代替", snp_hit[1]))
      colnames(qtl_df)[colnames(qtl_df) == snp_hit[1]] <- "SNP"
    } else {
      stop("QTL.csv 中未找到 SNP 列（或类似列名）。")
    }
  }

  slope_df2 <- data.frame(
    SNP = slope_df$SNP,
    slope_p = slope_df$P.value,
    slope_logp = -log10(pmax(slope_df$P.value, 1e-300))
  )
  intcp_df2 <- data.frame(
    SNP = intcp_df$SNP,
    intcp_p = intcp_df$P.value,
    intcp_logp = -log10(pmax(intcp_df$P.value, 1e-300))
  )

  plot_df <- merge(slope_df2, intcp_df2, by = "SNP")
  plot_df <- plot_df[is.finite(plot_df$slope_logp) & is.finite(plot_df$intcp_logp), ]
  cat(sprintf("合并后的SNP数量: %d\n", nrow(plot_df)))

  qtl_snps <- unique(qtl_df$SNP[!is.na(qtl_df$SNP)])
  plot_df$is_qtl <- plot_df$SNP %in% qtl_snps
  cat(sprintf("QTL.csv中的QTL SNPs: %d\n", length(qtl_snps)))
  cat(sprintf("GWAS合并数据中匹配的QTL SNPs: %d\n", sum(plot_df$is_qtl, na.rm = TRUE)))

  cat("\n计算斜率和截距-log10(P)值的相关系数和p值...\n")
  cor_test <- NULL
  cor_test_spearman <- NULL
  correlation_results <- NULL
  summary_stats <- NULL

  if (nrow(plot_df) >= 3) {
    cor_test <- tryCatch({
      cor.test(plot_df$intcp_logp, plot_df$slope_logp,
               method = "pearson", use = "complete.obs")
    }, error = function(e) {
      cat("Pearson相关系数计算失败:", e$message, "\n")
      return(NULL)
    })

    cor_test_spearman <- tryCatch({
      cor.test(plot_df$intcp_logp, plot_df$slope_logp,
               method = "spearman", use = "complete.obs", exact = FALSE)
    }, error = function(e) {
      cat("Spearman相关系数计算失败:", e$message, "\n")
      return(NULL)
    })

    if (!is.null(cor_test)) {
      correlation_results <- data.frame(
        Analysis = c("Pearson", "Spearman"),
        Correlation = c(if (!is.null(cor_test)) cor_test$estimate else NA,
                        if (!is.null(cor_test_spearman)) cor_test_spearman$estimate else NA),
        P_value = c(if (!is.null(cor_test)) cor_test$p.value else NA,
                    if (!is.null(cor_test_spearman)) cor_test_spearman$p.value else NA),
        Statistic = c(if (!is.null(cor_test)) cor_test$statistic else NA,
                      if (!is.null(cor_test_spearman)) cor_test_spearman$statistic else NA),
        DF = c(if (!is.null(cor_test)) cor_test$parameter else NA, NA),
        CI_lower = c(if (!is.null(cor_test)) cor_test$conf.int[1] else NA, NA),
        CI_upper = c(if (!is.null(cor_test)) cor_test$conf.int[2] else NA, NA),
        Method = c(if (!is.null(cor_test)) cor_test$method else NA,
                   if (!is.null(cor_test_spearman)) cor_test_spearman$method else NA),
        Alternative = c(if (!is.null(cor_test)) cor_test$alternative else NA,
                        if (!is.null(cor_test_spearman)) cor_test_spearman$alternative else NA),
        N = c(if (!is.null(cor_test)) cor_test$parameter + 2 else NA,
              length(na.omit(plot_df$intcp_logp))),
        stringsAsFactors = FALSE
      )

      summary_stats <- data.frame(
        Statistic = c("Mean_intcp_logp", "Mean_slope_logp", "SD_intcp_logp", "SD_slope_logp",
                      "Min_intcp_logp", "Max_intcp_logp", "Min_slope_logp", "Max_slope_logp",
                      "N_total", "N_QTL"),
        Value = c(mean(plot_df$intcp_logp, na.rm = TRUE),
                  mean(plot_df$slope_logp, na.rm = TRUE),
                  sd(plot_df$intcp_logp, na.rm = TRUE),
                  sd(plot_df$slope_logp, na.rm = TRUE),
                  min(plot_df$intcp_logp, na.rm = TRUE),
                  max(plot_df$intcp_logp, na.rm = TRUE),
                  min(plot_df$slope_logp, na.rm = TRUE),
                  max(plot_df$slope_logp, na.rm = TRUE),
                  nrow(plot_df),
                  sum(plot_df$is_qtl, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )

      save_intermediate_data(correlation_results,
                             "Slope_Intercept_Correlation_Analysis.csv",
                             "斜率和截距-log10(P)值相关系数分析")
      save_intermediate_data(summary_stats,
                             "Slope_Intercept_Summary_Statistics.csv",
                             "斜率和截距-log10(P)值汇总统计")
      cor_data <- plot_df[, c("SNP", "intcp_logp", "slope_logp", "is_qtl")]
      save_intermediate_data(cor_data,
                             "Slope_Intercept_Correlation_Data.csv",
                             "斜率和截距-log10(P)值相关分析数据")

      cat("\n斜率和截距-log10(P)值相关分析结果:\n")
      cat("========================================\n")
      if (!is.null(cor_test)) {
        cat(sprintf("Pearson相关系数: r = %.4f, p = %.2e\n",
                    cor_test$estimate, cor_test$p.value))
      }
      if (!is.null(cor_test_spearman)) {
        cat(sprintf("Spearman相关系数: ρ = %.4f, p = %.2e\n",
                    cor_test_spearman$estimate, cor_test_spearman$p.value))
      }
      cat(sprintf("样本数量: %d (其中QTL标记: %d)\n", nrow(plot_df), sum(plot_df$is_qtl)))
      cat(sprintf("截距-log10(P)均值: %.3f ± %.3f\n",
                  mean(plot_df$intcp_logp, na.rm = TRUE),
                  sd(plot_df$intcp_logp, na.rm = TRUE)))
      cat(sprintf("斜率-log10(P)均值: %.3f ± %.3f\n",
                  mean(plot_df$slope_logp, na.rm = TRUE),
                  sd(plot_df$slope_logp, na.rm = TRUE)))

      if (!is.null(cor_test)) {
        cor_label <- sprintf("r = %.3f", cor_test$estimate)
      } else {
        cor_label <- "Correlation not calculated"
      }
    }
  } else {
    cat("警告：数据点不足，无法计算相关系数\n")
    cor_label <- "Insufficient data for correlation"
  }

  max_lim <- max(plot_df$slope_logp, plot_df$intcp_logp, na.rm = TRUE)
  max_lim <- ceiling(max_lim)
  max_lim <- max(max_lim, 10)

  p <- ggplot() +
    geom_point(
      data = plot_df[!plot_df$is_qtl, ],
      aes(x = intcp_logp, y = slope_logp),
      size = 1.1, alpha = 0.85, color = "gray60"
    ) +
    geom_point(
      data = plot_df[plot_df$is_qtl, ],
      aes(x = intcp_logp, y = slope_logp),
      size = 2.8, alpha = 0.95, color = "#D55E00"
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.8, color = "grey70") +
    coord_equal(xlim = c(0, max_lim), ylim = c(0, max_lim), expand = FALSE) +
    labs(
      x = expression(-log[10](italic(P)) ~ intercept),
      y = expression(-log[10](italic(P)) ~ slope)
    ) +
    panel_label_layer("c") +
    theme_axis_std(base_size = 13, plot_margin = panel_margin_c) +
    coord_cartesian(xlim = c(0, max_lim), ylim = c(0, max_lim), clip = "off")

  if (exists("cor_label")) {
    p <- p +
      annotate("text",
               x = max_lim * 0.05,
               y = max_lim * 0.95,
               label = cor_label,
               hjust = 0, vjust = 1,
               size = 4, color = "blue",
               fontface = "italic")
  }

  plot_c_dir <- file.path(output_main_dir, "QTL_Plots")
  if (!dir.exists(plot_c_dir)) dir.create(plot_c_dir, recursive = TRUE)

  plot_c_pdf <- file.path(plot_c_dir, "Slope_vs_Intercept_Correlation.pdf")
  ggsave(plot_c_pdf, p, width = 8, height = 8, units = "in", bg = "white")
  cat("✓ 图c保存:", plot_c_pdf, "\n")

  convert_pdf_to_images(plot_c_pdf, plot_c_dir, "Slope_vs_Intercept_Correlation", c("png", "tiff"))

  return(p)
}

# ==================== 3. 定义图d (QTL Effect by DTR&mean&PostFlowering1_14) 绘制函数 ====================
plot_d <- function() {
  slope_qtl_path <- file.path(base_dir, "QTL.csv")
  base_env_path  <- file.path(base_dir, "GAPIT_MLM_Results/Single_Environment")
  ecs_path       <- file.path(base_dir, "EC8.csv")

  cat("读取QTL数据...\n")
  slope_qtl <- read_csv(slope_qtl_path, show_col_types = FALSE)

  cat("读取环境协变量数据...\n")
  ecs_data <- read_csv(ecs_path, show_col_types = FALSE)
  if ("env" %in% colnames(ecs_data)) colnames(ecs_data)[colnames(ecs_data) == "env"] <- "env_code"

  pve_col <- "Phenotype_Variance_Explained(%)"
  if (!(pve_col %in% colnames(slope_qtl))) {
    hit <- grep("Phenotype_Variance_Explained", colnames(slope_qtl), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) {
      message(sprintf("未找到列 '%s'，使用 '%s' 代替", pve_col, hit[1]))
      pve_col <- hit[1]
    } else {
      stop("QTL.csv 中未找到 Phenotype_Variance_Explained(%) 列")
    }
  }

  slope_qtl_valid <- slope_qtl[!is.na(slope_qtl$QTLname), c("QTLname", pve_col)]
  slope_qtl_valid <- slope_qtl_valid[!is.na(slope_qtl_valid[[pve_col]]), ]

  qtl_pve <- data.frame()
  for (qtl in unique(slope_qtl_valid$QTLname)) {
    qtl_data <- slope_qtl_valid[slope_qtl_valid$QTLname == qtl, ]
    if (nrow(qtl_data) > 0) {
      pve_val <- qtl_data[[pve_col]][1]
      pve_fmt <- sprintf("%.2f%%", as.numeric(pve_val))
      qtl_label <- paste0(qtl, " (", pve_fmt, ")")
      qtl_pve <- rbind(qtl_pve, data.frame(
        QTLname = qtl,
        PVE = pve_val,
        PVE_fmt = pve_fmt,
        QTL_label = qtl_label,
        stringsAsFactors = FALSE
      ))
    }
  }
  save_intermediate_data(qtl_pve, "QTL_PVE_Summary.csv", "QTL表型方差解释率汇总")

  all_effects <- unique(slope_qtl[, c("QTLname", "SNP", "Chromosome", "Position")])
  all_effects <- all_effects[!duplicated(all_effects$SNP), ]
  cat(sprintf("基础QTL数据: %d 个SNP\n", nrow(all_effects)))

  cat("\n正在从各个环境文件中提取效应值...\n")
  for (env in env_dirs) {
    file_path1 <- file.path(base_env_path, env,
                            paste0("GAPIT.Association.GWAS_Results.MLM_GWAS_", env, ".MLM.Trait(NYC).csv"))
    file_path2 <- file.path(base_env_path, env,
                            "GAPIT.Association.GWAS_Results.TKW.MLM.Trait(NYC).csv")

    if (file.exists(file_path1)) {
      file_path <- file_path1
    } else if (file.exists(file_path2)) {
      file_path <- file_path2
    } else {
      warning(paste("环境", env, "文件不存在，跳过"))
      all_effects[[env]] <- NA
      next
    }

    cat(sprintf("读取环境 %s 文件: %s\n", env, basename(file_path)))
    env_data <- read_csv(file_path, show_col_types = FALSE)

    colnames(env_data) <- gsub(" ", "", colnames(env_data))
    colnames(env_data) <- gsub("&", ".", colnames(env_data))

    if ("Chromosome" %in% colnames(env_data)) colnames(env_data)[colnames(env_data) == "Chromosome"] <- "Chr"
    if ("Position" %in% colnames(env_data)) colnames(env_data)[colnames(env_data) == "Position"] <- "Pos"
    if ("effect" %in% colnames(env_data)) colnames(env_data)[colnames(env_data) == "effect"] <- "Effect"

    if ("Effect" %in% colnames(env_data) && "SNP" %in% colnames(env_data)) {
      env_effects <- env_data[, c("SNP", "Effect")]
      env_effects <- env_effects[!duplicated(env_effects$SNP), ]
      colnames(env_effects)[2] <- env

      all_effects <- merge(all_effects, env_effects, by = "SNP", all.x = TRUE)

      matched <- sum(!is.na(all_effects[[env]]))
      total   <- nrow(all_effects)
      cat(sprintf("环境 %s: 匹配 %d/%d (%.1f%%)\n", env, matched, total, 100 * matched / total))
    } else {
      warning(paste("环境", env, "文件缺少Effect或SNP列"))
      all_effects[[env]] <- NA
    }
  }
  save_intermediate_data(all_effects, "QTL_Effect_Values_All_Environments.csv", "所有环境的QTL效应值")

  long_format <- data.frame()
  for (i in 1:nrow(all_effects)) {
    row_data <- all_effects[i, ]
    qtl_name <- row_data$QTLname
    snp <- row_data$SNP
    chr <- row_data$Chromosome
    pos <- row_data$Position

    for (env in env_dirs) {
      effect_val <- row_data[[env]]
      if (!is.na(effect_val)) {
        long_format <- rbind(long_format, data.frame(
          QTLname = qtl_name,
          SNP = snp,
          Chromosome = chr,
          Position = pos,
          env_code = env,
          Effect = effect_val,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  x_var_name <- "DTR&mean&PostFlowering1_14"
  if (!(x_var_name %in% colnames(ecs_data))) {
    hit <- grep("DTR.*mean.*PostFlowering", colnames(ecs_data), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) {
      message(sprintf("未找到列 '%s'，使用 '%s' 代替", x_var_name, hit[1]))
      x_var_name <- hit[1]
    } else {
      stop("EC8.csv 中未找到 DTR&mean&PostFlowering1_14 列")
    }
  }

  env_xvar <- ecs_data[, c("env_code", x_var_name)]
  colnames(env_xvar)[2] <- "Xvar"
  long_format <- merge(long_format, env_xvar, by = "env_code", all.x = TRUE)

  data_clean <- long_format[!is.na(long_format$Effect) & !is.na(long_format$Xvar), ]
  save_intermediate_data(data_clean, "QTL_Effect_Clean_Data.csv", "清理后的QTL效应数据")

  qtl_summary <- data.frame()
  for (qtl in unique(data_clean$QTLname)) {
    qtl_data <- data_clean[data_clean$QTLname == qtl, ]
    for (env in unique(qtl_data$env_code)) {
      env_data <- qtl_data[qtl_data$env_code == env, ]
      mean_effect <- mean(env_data$Effect, na.rm = TRUE)
      n_snps <- nrow(env_data)
      xvar <- unique(env_data$Xvar)

      qtl_summary <- rbind(qtl_summary, data.frame(
        QTLname = qtl,
        env_code = env,
        Xvar = xvar,
        mean_effect = mean_effect,
        n_snps = n_snps,
        stringsAsFactors = FALSE
      ))
    }
  }
  save_intermediate_data(qtl_summary, "QTL_Summary_by_Environment.csv", "按环境汇总的QTL效应")

  available_qtl <- unique(qtl_summary$QTLname)
  selected_qtl <- available_qtl
  cat(sprintf("检测到 %d 个QTL，全部展示\n", length(selected_qtl)))

  plot_df <- qtl_summary[qtl_summary$QTLname %in% selected_qtl, ]
  plot_df$QTLname <- factor(plot_df$QTLname, levels = selected_qtl)

  plot_df <- merge(plot_df, qtl_pve[, c("QTLname", "QTL_label")], by = "QTLname", all.x = TRUE)
  plot_df$QTL_label <- ifelse(is.na(plot_df$QTL_label), as.character(plot_df$QTLname), plot_df$QTL_label)

  env_order <- unique(plot_df[, c("env_code", "Xvar")])
  env_order <- env_order[order(env_order$Xvar), ]

  lm_data <- data.frame()
  lm_summary <- data.frame()

  for (qtl in selected_qtl) {
    qtl_data <- plot_df[plot_df$QTLname == qtl, ]
    if (nrow(qtl_data) >= 2) {
      lm_model <- lm(mean_effect ~ Xvar, data = qtl_data)
      lm_sum <- summary(lm_model)

      lm_stats <- data.frame(
        QTLname = qtl,
        Intercept = coef(lm_model)[1],
        Slope = coef(lm_model)[2],
        R_squared = lm_sum$r.squared,
        Adjusted_R_squared = lm_sum$adj.r.squared,
        F_statistic = lm_sum$fstatistic[1],
        P_value = pf(lm_sum$fstatistic[1], lm_sum$fstatistic[2], lm_sum$fstatistic[3], lower.tail = FALSE),
        N = nrow(qtl_data),
        stringsAsFactors = FALSE
      )
      lm_summary <- rbind(lm_summary, lm_stats)

      x_range <- range(env_order$Xvar, na.rm = TRUE)
      x_seq <- seq(x_range[1], x_range[2], length.out = 250)
      y_pred <- predict(lm_model, newdata = data.frame(Xvar = x_seq))
      lm_data <- rbind(lm_data, data.frame(
        QTLname = qtl,
        Xvar = x_seq,
        mean_effect = y_pred,
        stringsAsFactors = FALSE
      ))
    }
  }
  if (nrow(lm_summary) > 0) {
    save_intermediate_data(lm_summary, "QTL_Regression_Analysis.csv", "QTL效应回归分析")
  }

  n_qtl <- length(selected_qtl)
  if (n_qtl <= 8) {
    colors <- brewer.pal(max(3, n_qtl), "Set2")[1:n_qtl]
  } else {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_qtl)
  }
  names(colors) <- selected_qtl

  color_map <- data.frame(
    QTLname = names(colors),
    Color = colors,
    stringsAsFactors = FALSE
  )
  save_intermediate_data(color_map, "QTL_Color_Mapping.csv", "QTL颜色映射")

  label_positions <- data.frame()
  for (qtl in selected_qtl) {
    qtl_lm <- lm_data[lm_data$QTLname == qtl, ]
    if (nrow(qtl_lm) > 0) {
      max_x_row <- qtl_lm[qtl_lm$Xvar == max(qtl_lm$Xvar), ]
      label_positions <- rbind(label_positions, max_x_row)
    }
  }
  label_positions <- merge(label_positions, plot_df[, c("QTLname", "QTL_label")], by = "QTLname")
  label_positions <- unique(label_positions)

  if (exists("lm_summary") && nrow(lm_summary) > 0) {
    lm_summary_sub <- lm_summary[, c("QTLname", "Slope", "P_value")]
    label_positions <- merge(label_positions, lm_summary_sub, by = "QTLname", all.x = TRUE)
    label_positions$Slope_formatted <- sprintf("%.3f", label_positions$Slope)
    label_positions$sig_stars <- cut(label_positions$P_value,
                                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                                     labels = c("***", "**", "*", ""))
    label_positions$QTL_label_new <- ifelse(
      !is.na(label_positions$Slope_formatted),
      paste0(label_positions$QTL_label, "\nSlope = ", label_positions$Slope_formatted, label_positions$sig_stars),
      label_positions$QTL_label
    )
  } else {
    label_positions$QTL_label_new <- label_positions$QTL_label
  }

  x_range <- range(env_order$Xvar, na.rm = TRUE)
  label_positions$label_x <- label_positions$Xvar + (x_range[2] - x_range[1]) * 0.15

  x_min <- min(env_order$Xvar, na.rm = TRUE)
  x_max <- max(env_order$Xvar, na.rm = TRUE)
  x_pad_right <- (x_max - x_min) * 0.45

  y_range <- range(plot_df$mean_effect, na.rm = TRUE)
  y_pad <- (y_range[2] - y_range[1]) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1

  # 使用统一的X轴刻度（用户提供的数值）
  custom_theme <- theme_axis_std(base_size = 12, plot_margin = panel_margin_de) +
    theme(
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12)),
      axis.text.x  = element_text(hjust = 0.5, vjust = 0.5, margin = margin(t = 5)),
      axis.text.y  = element_text(colour = "black")
    ) +
    theme_legend_std(position = c(0.98, 0.98), direction = "vertical") +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1)
    )

  p <- ggplot() +
    geom_point(
      data = plot_df,
      aes(x = Xvar, y = mean_effect, color = QTLname, shape = env_code),
      size = 2.6, alpha = 0.85, stroke = 0.7
    ) +
    geom_line(
      data = lm_data,
      aes(x = Xvar, y = mean_effect, color = QTLname, group = QTLname),
      linewidth = 1.0,
      alpha = 0.95
    ) +
    geom_text_repel(
      data = label_positions,
      aes(x = label_x, y = mean_effect, label = QTL_label_new, color = QTLname),
      hjust = 0,
      vjust = 0.5,
      size = 3.5,
      fontface = label_params$qtl_label_fontface,
      direction = "y",
      box.padding = 0.6,
      point.padding = 0.5,
      segment.size = 0.3,
      min.segment.length = 0,
      segment.color = "grey50",
      seed = 123,
      max.iter = 5000,
      family = "sans"
    ) +
    scale_x_continuous(
      name = expression("DTR&mean&PostFlowering1_14 ("*degree*C*")"),  # 修改为摄氏度
      breaks = x_ticks_de,
      labels = x_tick_labels_de,
      expand = expansion(mult = c(0.00, 0.00)),
      limits = c(x_min, x_max + x_pad_right)
    ) +
    labs(y = "QTL Effect (g)") +
    scale_color_manual(values = colors, guide = "none") +
    scale_shape_manual(values = env_shape_values, breaks = env_dirs, name = "Environments") +
    guides(shape = guide_legend(override.aes = list(color = "black", alpha = 1, size = 3))) +
    custom_theme +
    panel_label_layer("d") +
    coord_cartesian(
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      clip = "off"
    )

  plot_d_dir <- file.path(output_main_dir, "QTL_Plots")
  if (!dir.exists(plot_d_dir)) dir.create(plot_d_dir, recursive = TRUE)

  plot_d_pdf <- file.path(plot_d_dir, "QTL_Effect_by_DTR_mean_PostFlowering1_14.pdf")
  ggsave(plot_d_pdf, p, width = 12, height = 8, units = "in", bg = "white")
  cat("✓ 图d保存:", plot_d_pdf, "\n")

  convert_pdf_to_images(plot_d_pdf, plot_d_dir, "QTL_Effect_by_DTR_mean_PostFlowering1_14", c("png", "tiff"))

  return(p)
}

# ==================== 4. 定义图e (Haplotype Fitted Lines) 绘制函数（使用新基因型文件）====================
plot_e <- function() {
  qtl_file <- file.path(base_dir, "QTL.csv")
  vcf_file <- file.path(base_dir, "5Kfiltered_alt_geno_imputed_renamed.vcf.gz")  # 新基因型文件
  phenotype_file <- file.path(base_dir, "TKW_mean_table.txt")  # 原始表型文件
  environment_file <- file.path(base_dir, "EC8.csv")
  target_factor <- "DTR&mean&PostFlowering1_14"

  # 构建单倍型数据
  haplo_data <- build_haplotype_data(qtl_file, vcf_file)
  haplotype_df <- haplo_data$haplotype_df
  all_haps <- haplo_data$all_haps
  haps_named <- haplo_data$haps_named

  cat("读取表型和环境数据...\n")
  phenotype_data <- read.table(
    phenotype_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, na.strings = c("", "NA")
  )
  phenotype_data$genotype <- as.character(phenotype_data$genotype)
  phenotype_data$env_code <- as.character(phenotype_data$env_code)
  phenotype_data$TKW <- as.numeric(phenotype_data$TKW)
  colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
  colnames(phenotype_data)[colnames(phenotype_data) == "TKW"] <- "PH"
  phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]

  # 清理表型品种名：去除两端空格、引号，转换为大写
  phenotype_data$line_code <- toupper(trimws(gsub('"', '', phenotype_data$line_code)))

  phenotype_summary <- data.frame(
    Total_lines = length(unique(phenotype_data$line_code)),
    Total_environments = length(unique(phenotype_data$env_code)),
    Total_observations = nrow(phenotype_data),
    Mean_PH = mean(phenotype_data$PH, na.rm = TRUE),
    SD_PH = sd(phenotype_data$PH, na.rm = TRUE),
    Min_PH = min(phenotype_data$PH, na.rm = TRUE),
    Max_PH = max(phenotype_data$PH, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(phenotype_summary, "Phenotype_Data_Summary.csv", "表型数据摘要")

  ecs_data <- read.csv(environment_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ecs_data)[1] <- "env_code"
  ecs_data$env_code <- as.character(ecs_data$env_code)

  if (!target_factor %in% colnames(ecs_data)) {
    cat("警告: 目标环境因子未找到。尝试模糊匹配...\n")
    possible <- grep("DTR.*mean.*PostFlowering", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
    if (length(possible) > 0) {
      target_factor <- possible[1]
      cat("使用:", target_factor, "\n")
    } else {
      stop("没有找到匹配的环境因子。")
    }
  }

  env_factor_data <- ecs_data[, c("env_code", target_factor)]
  colnames(env_factor_data)[2] <- "env_factor_value"
  env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)

  env_factor_summary <- data.frame(
    Environment_factor = target_factor,
    N_environments = nrow(env_factor_data),
    Mean_value = mean(env_factor_data$env_factor_value, na.rm = TRUE),
    SD_value = sd(env_factor_data$env_factor_value, na.rm = TRUE),
    Min_value = min(env_factor_data$env_factor_value, na.rm = TRUE),
    Max_value = max(env_factor_data$env_factor_value, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(env_factor_summary, "Environment_Factor_Summary.csv", "环境因子摘要")

  cat("合并数据...\n")
  # 左连接保留所有表型数据
  pheno_haplo <- merge(
    phenotype_data,
    haplotype_df[, c("Sample", "Haplotype")],
    by.x = "line_code", by.y = "Sample",
    all.x = TRUE, all.y = FALSE
  )

  # 合并环境数据（内连接，环境必须存在）
  pheno_haplo_env <- merge(
    pheno_haplo,
    env_factor_data,
    by = "env_code",
    all.x = FALSE, all.y = FALSE
  )
  pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$env_factor_value), ]

  # 过滤掉单倍型为NA的行
  pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$Haplotype), ]

  pheno_haplo_env$Haplotype_Named <- ifelse(
    pheno_haplo_env$Haplotype %in% all_haps,
    haps_named[match(pheno_haplo_env$Haplotype, names(haps_named))],
    "Other"
  )

  cat(sprintf("合并后的数据行数: %d\n", nrow(pheno_haplo_env)))
  cat(sprintf("品种数: %d | 单倍型数: %d\n",
              length(unique(pheno_haplo_env$line_code)),
              length(unique(pheno_haplo_env$Haplotype))))
  save_intermediate_data(pheno_haplo_env, "Phenotype_Haplotype_Environment_Merged.csv",
                         "表型-单倍型-环境合并数据")

  cat("拟合株系模型...\n")
  line_env_count <- aggregate(env_code ~ line_code, data = pheno_haplo_env, FUN = function(x) length(unique(x)))
  colnames(line_env_count)[2] <- "n_env"
  valid_lines <- line_env_count$line_code[line_env_count$n_env >= 2]

  line_models <- data.frame()
  for (line in valid_lines) {
    line_data <- pheno_haplo_env[pheno_haplo_env$line_code == line, ]
    if (nrow(line_data) >= 2) {
      m <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
      if (!inherits(m, "try-error")) {
        fit_summary <- summary(m)
        intercept_fit <- fit_summary$coefficients["(Intercept)", "Estimate"]
        slope_fit <- fit_summary$coefficients["env_factor_value", "Estimate"]
        r_squared <- fit_summary$r.squared
        p_value <- pf(fit_summary$fstatistic[1], fit_summary$fstatistic[2],
                      fit_summary$fstatistic[3], lower.tail = FALSE)

        hap <- unique(line_data$Haplotype)[1]
        hap_named <- ifelse(hap %in% all_haps,
                            haps_named[match(hap, names(haps_named))],
                            "Other")
        line_models <- rbind(line_models, data.frame(
          line_code = line,
          Haplotype = hap,
          Haplotype_Named = hap_named,
          Slope = slope_fit,
          Intercept = intercept_fit,
          R_squared = r_squared,
          P_value = p_value,
          N_points = nrow(line_data),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  cat(sprintf("株系模型数量: %d\n", nrow(line_models)))

  cat("拟合单倍型平均模型...\n")
  haplotype_models <- data.frame()
  for (hap in all_haps) {
    dat <- pheno_haplo_env[pheno_haplo_env$Haplotype == hap, ]
    n_lines <- length(unique(dat$line_code))
    n_obs <- nrow(dat)

    if (n_lines >= 2 && n_obs >= 4) {
      m <- try(lm(PH ~ env_factor_value, data = dat), silent = TRUE)
      if (!inherits(m, "try-error")) {
        fit_summary <- summary(m)
        intercept_fit <- fit_summary$coefficients["(Intercept)", "Estimate"]
        slope_fit <- fit_summary$coefficients["env_factor_value", "Estimate"]
        r_squared <- fit_summary$r.squared
        p_value <- pf(fit_summary$fstatistic[1], fit_summary$fstatistic[2],
                      fit_summary$fstatistic[3], lower.tail = FALSE)

        hap_named <- haps_named[match(hap, names(haps_named))]
        haplotype_models <- rbind(haplotype_models, data.frame(
          Haplotype = hap,
          Haplotype_Named = hap_named,
          Slope = slope_fit,
          Intercept = intercept_fit,
          R_squared = r_squared,
          P_value = p_value,
          N_lines = n_lines,
          N_observations = n_obs,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  cat(sprintf("单倍型模型数量: %d\n", nrow(haplotype_models)))

  haplotype_output_dir <- file.path(output_main_dir, "Haplotype_Analysis")
  write.csv(line_models, file.path(haplotype_output_dir, "Line_Models_by_Haplotype.csv"), row.names = FALSE)
  write.csv(haplotype_models, file.path(haplotype_output_dir, "Haplotype_Models.csv"), row.names = FALSE)
  save_intermediate_data(line_models, "Line_Models_by_Haplotype.csv", "按单倍型的株系模型")
  save_intermediate_data(haplotype_models, "Haplotype_Mean_Models.csv", "单倍型平均模型")

  cat("检验株系斜率和截距的显著性(ANOVA + Tukey HSD)...\n")
  if (nrow(line_models) >= 5 && length(unique(line_models$Haplotype)) >= 2) {
    line_params <- line_models[, c("Haplotype", "Slope", "Intercept")]
    line_params$Haplotype <- as.factor(line_params$Haplotype)

    slope_anova <- aov(Slope ~ Haplotype, data = line_params)
    slope_tukey <- glht(slope_anova, linfct = mcp(Haplotype = "Tukey"))
    slope_tukey_res <- summary(slope_tukey)
    slope_tukey_df <- data.frame(
      Comparison = names(slope_tukey_res$test$coefficients),
      Slope_Difference = as.numeric(slope_tukey_res$test$coefficients),
      Std_Error = as.numeric(slope_tukey_res$test$sigma),
      t_Value = as.numeric(slope_tukey_res$test$tstat),
      P_Value = as.numeric(slope_tukey_res$test$pvalues),
      stringsAsFactors = FALSE
    )
    slope_tukey_df$Significance <- ifelse(slope_tukey_df$P_Value < 0.05, "Significant", "Not Significant")
    slope_tukey_df$Significance_Level <- cut(slope_tukey_df$P_Value,
                                              breaks = c(0, 0.001, 0.01, 0.05, 1),
                                              labels = c("***", "**", "*", "NS"))

    intercept_anova <- aov(Intercept ~ Haplotype, data = line_params)
    intercept_tukey <- glht(intercept_anova, linfct = mcp(Haplotype = "Tukey"))
    intercept_tukey_res <- summary(intercept_tukey)
    intercept_tukey_df <- data.frame(
      Comparison = names(intercept_tukey_res$test$coefficients),
      Intercept_Difference = as.numeric(intercept_tukey_res$test$coefficients),
      Std_Error = as.numeric(intercept_tukey_res$test$sigma),
      t_Value = as.numeric(intercept_tukey_res$test$tstat),
      P_Value = as.numeric(intercept_tukey_res$test$pvalues),
      stringsAsFactors = FALSE
    )
    intercept_tukey_df$Significance <- ifelse(intercept_tukey_df$P_Value < 0.05, "Significant", "Not Significant")
    intercept_tukey_df$Significance_Level <- cut(intercept_tukey_df$P_Value,
                                                  breaks = c(0, 0.001, 0.01, 0.05, 1),
                                                  labels = c("***", "**", "*", "NS"))

    write.csv(slope_tukey_df, file.path(haplotype_output_dir, "Line_Slope_Tukey_HSD_Results.csv"), row.names = FALSE)
    write.csv(intercept_tukey_df, file.path(haplotype_output_dir, "Line_Intercept_Tukey_HSD_Results.csv"), row.names = FALSE)
    save_intermediate_data(slope_tukey_df, "Line_Slope_Tukey_HSD_Results.csv", "株系斜率Tukey HSD分析")
    save_intermediate_data(intercept_tukey_df, "Line_Intercept_Tukey_HSD_Results.csv", "株系截距Tukey HSD分析")

    cat("\n单倍型间差异显著性摘要:\n")
    cat("==========================\n")
    cat(sprintf("斜率差异显著对比数: %d/%d\n",
                sum(slope_tukey_df$P_Value < 0.05), nrow(slope_tukey_df)))
    cat(sprintf("截距差异显著对比数: %d/%d\n",
                sum(intercept_tukey_df$P_Value < 0.05), nrow(intercept_tukey_df)))
  } else {
    cat("警告: 数据不足，无法进行ANOVA/Tukey检验。跳过...\n")
    write.csv(data.frame(), file.path(haplotype_output_dir, "Line_Slope_Tukey_HSD_Results.csv"), row.names = FALSE)
    write.csv(data.frame(), file.path(haplotype_output_dir, "Line_Intercept_Tukey_HSD_Results.csv"), row.names = FALSE)
  }

  cat("创建颜色映射...\n")
  n_colors <- length(all_haps)
  if (n_colors <= 8) {
    hap_cols <- brewer.pal(max(3, n_colors), "Set2")[1:n_colors]
  } else if (n_colors <= 12) {
    hap_cols <- brewer.pal(n_colors, "Paired")[1:n_colors]
  } else {
    hap_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  }
  names(hap_cols) <- haps_named

  haplotype_color_map <- data.frame(
    Haplotype = all_haps,
    Haplotype_Named = haps_named,
    Color = hap_cols[haps_named],
    stringsAsFactors = FALSE
  )
  save_intermediate_data(haplotype_color_map, "Haplotype_Color_Mapping.csv", "单倍型颜色映射")

  haplotype_labels <- character(length(haps_named))
  names(haplotype_labels) <- haps_named
  for (hap_name in haps_named) {
    idx <- which(haplotype_models$Haplotype_Named == hap_name)
    if (length(idx) > 0) {
      intercept_fit <- haplotype_models$Intercept[idx]
      slope_fit <- haplotype_models$Slope[idx]
      haplotype_labels[hap_name] <- sprintf("%s (%.2f, %.2f)", hap_name, intercept_fit, slope_fit)
    } else {
      haplotype_labels[hap_name] <- hap_name
    }
  }

  # 环境范围根据实际数据计算
  env_range <- range(env_factor_data$env_factor_value, na.rm = TRUE)
  env_seq <- seq(env_range[1], env_range[2], length.out = 120)

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
        stringsAsFactors = FALSE
      )
    }
  }
  line_pred <- if (length(line_pred_list) > 0) do.call(rbind, line_pred_list) else data.frame()

  hap_pred_list <- list()
  if (nrow(haplotype_models) > 0) {
    for (i in 1:nrow(haplotype_models)) {
      one <- haplotype_models[i, ]
      hap_pred_list[[i]] <- data.frame(
        Haplotype = one$Haplotype,
        Haplotype_Named = one$Haplotype_Named,
        env_factor_value = env_seq,
        Predicted_PH = one$Intercept + one$Slope * env_seq,
        stringsAsFactors = FALSE
      )
    }
  }
  hap_pred <- if (length(hap_pred_list) > 0) do.call(rbind, hap_pred_list) else data.frame()

  if (nrow(line_pred) > 0) {
    save_intermediate_data(line_pred, "Line_Predictions.csv", "株系表型预测值")
  }
  if (nrow(hap_pred) > 0) {
    save_intermediate_data(hap_pred, "Haplotype_Predictions.csv", "单倍型表型预测值")
  }

  cat("绘制主图...\n")
  y_range <- range(pheno_haplo_env$PH, na.rm = TRUE)
  y_pad <- diff(y_range) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1

  p_final <- ggplot() +
    geom_point(
      data = pheno_haplo_env,
      aes(x = env_factor_value, y = PH, color = Haplotype_Named),
      alpha = 0.6, size = 2, stroke = 0.5
    ) +
    {if (nrow(line_pred) > 0)
      geom_line(
        data = line_pred,
        aes(x = env_factor_value, y = Predicted_PH, group = line_code, color = Haplotype_Named),
        linetype = "dashed", alpha = 0.35, linewidth = 0.35
      )
    } +
    {if (nrow(hap_pred) > 0)
      geom_line(
        data = hap_pred,
        aes(x = env_factor_value, y = Predicted_PH, group = Haplotype_Named, color = Haplotype_Named),
        linetype = "solid", alpha = 0.95, linewidth = 1.6
      )
    } +
    scale_color_manual(values = hap_cols, labels = haplotype_labels, name = "Haplotypes") +
    scale_x_continuous(
      name = expression("DTR&mean&PostFlowering1_14 ("*degree*C*")"),  # 修改为摄氏度
      limits = env_range,
      breaks = x_ticks_de,
      labels = x_tick_labels_de,
      expand = expansion(mult = c(0.00, 0.00))
    ) +
    labs(y = "Thousand Kernel Weight (g)") +
    panel_label_layer("e") +
    theme_axis_std(base_size = 13, plot_margin = panel_margin_de) +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.75), color = NA),
      legend.key = element_blank(),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      panel.border = element_blank()
    ) +
    coord_cartesian(
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      clip = "off"
    )

  plot_e_dir <- file.path(output_main_dir, "Haplotype_Analysis")
  if (!dir.exists(plot_e_dir)) dir.create(plot_e_dir, recursive = TRUE)

  plot_e_pdf <- file.path(plot_e_dir, "Haplotype_Fitted_Lines.pdf")
  ggsave(plot_e_pdf, p_final, width = 12, height = 8, units = "in", bg = "white")
  cat("✓ 图e保存:", plot_e_pdf, "\n")

  convert_pdf_to_images(plot_e_pdf, plot_e_dir, "Haplotype_Fitted_Lines", c("png", "tiff"))

  return(p_final)
}

# ==================== 新增函数：对第二个表型 testTKW.txt 进行单倍型分析并生成T检验图（含株系过滤）====================
plot_extra_t_test <- function() {
  qtl_file <- file.path(base_dir, "QTL.csv")
  vcf_file <- file.path(base_dir, "5Kfiltered_alt_geno_imputed_renamed.vcf.gz")
  phenotype_file <- file.path(base_dir, "testTKW.txt")  # 第二个表型文件
  environment_file <- file.path(base_dir, "EC8.csv")
  target_factor <- "DTR&mean&PostFlowering1_14"

  # 构建单倍型数据（与图e相同）
  haplo_data <- build_haplotype_data(qtl_file, vcf_file)
  haplotype_df <- haplo_data$haplotype_df
  all_haps <- haplo_data$all_haps
  haps_named <- haplo_data$haps_named

  cat("\n=== 处理第二个表型 testTKW.txt ===\n")
  cat("读取表型和环境数据...\n")
  phenotype_data <- read.table(
    phenotype_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, na.strings = c("", "NA")
  )
  # 标准化列名：假设文件包含 genotype, env_code, TKW
  colnames(phenotype_data) <- tolower(colnames(phenotype_data))
  if (!("genotype" %in% colnames(phenotype_data)) || !("env_code" %in% colnames(phenotype_data)) || !("tkw" %in% colnames(phenotype_data))) {
    stop("testTKW.txt 必须包含列: genotype, env_code, TKW")
  }
  phenotype_data$genotype <- as.character(phenotype_data$genotype)
  phenotype_data$env_code <- as.character(phenotype_data$env_code)
  phenotype_data$tkw <- as.numeric(phenotype_data$tkw)
  colnames(phenotype_data)[colnames(phenotype_data) == "genotype"] <- "line_code"
  colnames(phenotype_data)[colnames(phenotype_data) == "tkw"] <- "PH"
  phenotype_data <- phenotype_data[!is.na(phenotype_data$PH), ]

  # 清理品种名
  phenotype_data$line_code <- toupper(trimws(gsub('"', '', phenotype_data$line_code)))

  phenotype_summary <- data.frame(
    Total_lines = length(unique(phenotype_data$line_code)),
    Total_environments = length(unique(phenotype_data$env_code)),
    Total_observations = nrow(phenotype_data),
    Mean_PH = mean(phenotype_data$PH, na.rm = TRUE),
    SD_PH = sd(phenotype_data$PH, na.rm = TRUE),
    Min_PH = min(phenotype_data$PH, na.rm = TRUE),
    Max_PH = max(phenotype_data$PH, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(phenotype_summary, "Phenotype_Data_Summary_testTKW.csv", "表型数据摘要 (testTKW)")

  # 读取环境数据（与之前相同）
  ecs_data <- read.csv(environment_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ecs_data)[1] <- "env_code"
  ecs_data$env_code <- as.character(ecs_data$env_code)

  if (!target_factor %in% colnames(ecs_data)) {
    possible <- grep("DTR.*mean.*PostFlowering", colnames(ecs_data), value = TRUE, ignore.case = TRUE)
    if (length(possible) > 0) {
      target_factor <- possible[1]
      cat("使用环境因子:", target_factor, "\n")
    } else {
      stop("没有找到匹配的环境因子。")
    }
  }

  env_factor_data <- ecs_data[, c("env_code", target_factor)]
  colnames(env_factor_data)[2] <- "env_factor_value"
  env_factor_data$env_factor_value <- as.numeric(env_factor_data$env_factor_value)

  env_factor_summary <- data.frame(
    Environment_factor = target_factor,
    N_environments = nrow(env_factor_data),
    Mean_value = mean(env_factor_data$env_factor_value, na.rm = TRUE),
    SD_value = sd(env_factor_data$env_factor_value, na.rm = TRUE),
    Min_value = min(env_factor_data$env_factor_value, na.rm = TRUE),
    Max_value = max(env_factor_data$env_factor_value, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(env_factor_summary, "Environment_Factor_Summary_testTKW.csv", "环境因子摘要 (testTKW)")

  cat("合并数据...\n")
  pheno_haplo <- merge(
    phenotype_data,
    haplotype_df[, c("Sample", "Haplotype")],
    by.x = "line_code", by.y = "Sample",
    all.x = TRUE, all.y = FALSE
  )

  pheno_haplo_env <- merge(
    pheno_haplo,
    env_factor_data,
    by = "env_code",
    all.x = FALSE, all.y = FALSE
  )
  pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$env_factor_value), ]
  pheno_haplo_env <- pheno_haplo_env[!is.na(pheno_haplo_env$Haplotype), ]

  # ==================== 关键步骤：过滤环境数据不完整的株系 ====================
  # 使用方案 B：保留环境数等于最大值的株系（适用于“少数株系缺一个环境”的情况）
  line_env_n <- pheno_haplo_env %>%
    dplyr::distinct(line_code, env_code) %>%
    dplyr::count(line_code, name = "n_env")
  
  max_env <- max(line_env_n$n_env)
  keep_lines <- line_env_n %>%
    dplyr::filter(n_env == max_env) %>%
    dplyr::pull(line_code)
  
  pheno_haplo_env <- pheno_haplo_env %>% dplyr::filter(line_code %in% keep_lines)
  
  cat(sprintf("过滤后：保留株系 %d 个（要求 n_env == %d）\n",
              length(unique(pheno_haplo_env$line_code)), max_env))
  
  # 方案 A（按全部环境集合过滤）如果需要，可取消注释并注释掉上面的方案 B
  # env_set <- sort(unique(pheno_haplo_env$env_code))
  # line_env_n <- pheno_haplo_env %>%
  #   dplyr::distinct(line_code, env_code) %>%
  #   dplyr::count(line_code, name = "n_env")
  # keep_lines <- line_env_n %>%
  #   dplyr::filter(n_env == length(env_set)) %>%
  #   dplyr::pull(line_code)
  # pheno_haplo_env <- pheno_haplo_env %>% dplyr::filter(line_code %in% keep_lines)
  # cat(sprintf("过滤后：保留株系 %d 个（要求覆盖 %d 个环境）\n",
  #             length(unique(pheno_haplo_env$line_code)), length(env_set)))
  # ========================================================================

  pheno_haplo_env$Haplotype_Named <- ifelse(
    pheno_haplo_env$Haplotype %in% all_haps,
    haps_named[match(pheno_haplo_env$Haplotype, names(haps_named))],
    "Other"
  )

  cat(sprintf("合并后数据行数: %d\n", nrow(pheno_haplo_env)))
  cat(sprintf("品种数: %d | 单倍型数: %d\n",
              length(unique(pheno_haplo_env$line_code)),
              length(unique(pheno_haplo_env$Haplotype))))
  save_intermediate_data(pheno_haplo_env, "Phenotype_Haplotype_Environment_Merged_testTKW.csv",
                         "表型-单倍型-环境合并数据 (testTKW)")

  cat("拟合株系模型...\n")
  line_env_count <- aggregate(env_code ~ line_code, data = pheno_haplo_env, FUN = function(x) length(unique(x)))
  colnames(line_env_count)[2] <- "n_env"
  valid_lines <- line_env_count$line_code[line_env_count$n_env >= 2]

  line_models <- data.frame()
  for (line in valid_lines) {
    line_data <- pheno_haplo_env[pheno_haplo_env$line_code == line, ]
    if (nrow(line_data) >= 2) {
      m <- try(lm(PH ~ env_factor_value, data = line_data), silent = TRUE)
      if (!inherits(m, "try-error")) {
        fit_summary <- summary(m)
        intercept_fit <- fit_summary$coefficients["(Intercept)", "Estimate"]
        slope_fit <- fit_summary$coefficients["env_factor_value", "Estimate"]
        r_squared <- fit_summary$r.squared
        p_value <- pf(fit_summary$fstatistic[1], fit_summary$fstatistic[2],
                      fit_summary$fstatistic[3], lower.tail = FALSE)

        hap <- unique(line_data$Haplotype)[1]
        hap_named <- ifelse(hap %in% all_haps,
                            haps_named[match(hap, names(haps_named))],
                            "Other")
        line_models <- rbind(line_models, data.frame(
          line_code = line,
          Haplotype = hap,
          Haplotype_Named = hap_named,
          Slope = slope_fit,
          Intercept = intercept_fit,
          R_squared = r_squared,
          P_value = p_value,
          N_points = nrow(line_data),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  cat(sprintf("株系模型数量: %d\n", nrow(line_models)))

  # 如果模型数量不足，无法进行T检验，则返回空图
  if (nrow(line_models) < 5 || length(unique(line_models$Haplotype)) < 2) {
    cat("警告：数据不足，无法进行T检验，生成空白图。\n")
    p_empty <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for T-test", size = 6) +
      theme_void()
    ggsave(file.path(output_main_dir, "Haplotype_Analysis", "TTest_Intercept_Slope_testTKW.pdf"), p_empty, width = 8, height = 6)
    return(p_empty)
  }

  # 准备用于ANOVA的数据
  line_params <- line_models[, c("Haplotype_Named", "Slope", "Intercept")]
  line_params$Haplotype_Named <- as.factor(line_params$Haplotype_Named)

  # 关键修复：过滤掉样本量少于2的单倍型组
  hap_counts <- table(line_params$Haplotype_Named)
  valid_hap <- names(hap_counts[hap_counts >= 2])
  if (length(valid_hap) < 2) {
    cat("警告：有效的单倍型组少于2个，无法进行T检验，生成空白图。\n")
    p_empty <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient haplotypes for T-test", size = 6) +
      theme_void()
    ggsave(file.path(output_main_dir, "Haplotype_Analysis", "TTest_Intercept_Slope_testTKW.pdf"), p_empty, width = 8, height = 6)
    return(p_empty)
  }
  line_params <- line_params[line_params$Haplotype_Named %in% valid_hap, ]
  line_params$Haplotype_Named <- droplevels(line_params$Haplotype_Named)

  cat("进行截距和斜率的单倍型间T检验...\n")

  # 使用Tukey HSD进行多重比较
  slope_anova <- aov(Slope ~ Haplotype_Named, data = line_params)
  slope_tukey <- glht(slope_anova, linfct = mcp(Haplotype_Named = "Tukey"))
  slope_tukey_res <- summary(slope_tukey)
  slope_sig <- data.frame(
    Comparison = names(slope_tukey_res$test$coefficients),
    p_value = as.numeric(slope_tukey_res$test$pvalues),
    stringsAsFactors = FALSE
  )

  intercept_anova <- aov(Intercept ~ Haplotype_Named, data = line_params)
  intercept_tukey <- glht(intercept_anova, linfct = mcp(Haplotype_Named = "Tukey"))
  intercept_tukey_res <- summary(intercept_tukey)
  intercept_sig <- data.frame(
    Comparison = names(intercept_tukey_res$test$coefficients),
    p_value = as.numeric(intercept_tukey_res$test$pvalues),
    stringsAsFactors = FALSE
  )

  # 生成显著性字母
  library(multcompView)
  slope_cld <- multcompLetters4(slope_anova, slope_tukey)
  intercept_cld <- multcompLetters4(intercept_anova, intercept_tukey)

  slope_letters <- data.frame(
    Haplotype_Named = names(slope_cld$Haplotype_Named$Letters),
    letter = slope_cld$Haplotype_Named$Letters,
    stringsAsFactors = FALSE
  )
  intercept_letters <- data.frame(
    Haplotype_Named = names(intercept_cld$Haplotype_Named$Letters),
    letter = intercept_cld$Haplotype_Named$Letters,
    stringsAsFactors = FALSE
  )

  # 计算每个单倍型的均值用于放置字母位置
  slope_means <- aggregate(Slope ~ Haplotype_Named, data = line_params, FUN = mean)
  intercept_means <- aggregate(Intercept ~ Haplotype_Named, data = line_params, FUN = mean)
  slope_letters <- merge(slope_letters, slope_means, by = "Haplotype_Named")
  intercept_letters <- merge(intercept_letters, intercept_means, by = "Haplotype_Named")

  # 颜色映射（与图e保持一致）
  n_colors <- length(unique(line_params$Haplotype_Named))
  if (n_colors <= 8) {
    hap_cols <- brewer.pal(max(3, n_colors), "Set2")[1:n_colors]
  } else if (n_colors <= 12) {
    hap_cols <- brewer.pal(n_colors, "Paired")[1:n_colors]
  } else {
    hap_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  }
  names(hap_cols) <- levels(line_params$Haplotype_Named)

  # 绘制截距箱线图
  p_intercept <- ggplot(line_params, aes(x = Haplotype_Named, y = Intercept, fill = Haplotype_Named)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_manual(values = hap_cols, guide = "none") +
    labs(x = "Haplotype", y = "Intercept") +
    theme_axis_std(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = intercept_letters, aes(x = Haplotype_Named, y = Intercept + 0.1 * max(Intercept), label = letter),
              size = 5, color = "black")

  # 绘制斜率箱线图
  p_slope <- ggplot(line_params, aes(x = Haplotype_Named, y = Slope, fill = Haplotype_Named)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_manual(values = hap_cols, guide = "none") +
    labs(x = "Haplotype", y = "Slope") +
    theme_axis_std(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = slope_letters, aes(x = Haplotype_Named, y = Slope + 0.1 * max(Slope), label = letter),
              size = 5, color = "black")

  # 组合两个图
  p_combined <- p_intercept + p_slope + plot_annotation(title = "T-test of Intercept and Slope among Haplotypes (testTKW)")

  # 保存
  ttest_dir <- file.path(output_main_dir, "Haplotype_Analysis")
  ttest_pdf <- file.path(ttest_dir, "TTest_Intercept_Slope_testTKW.pdf")
  ggsave(ttest_pdf, p_combined, width = 12, height = 6, units = "in", bg = "white")
  cat("✓ T检验图保存:", ttest_pdf, "\n")

  # 同时保存显著性表格
  write.csv(slope_sig, file.path(ttest_dir, "Slope_Tukey_HSD_testTKW.csv"), row.names = FALSE)
  write.csv(intercept_sig, file.path(ttest_dir, "Intercept_Tukey_HSD_testTKW.csv"), row.names = FALSE)

  return(p_combined)
}

# ==================== 5. 组合所有图表 ====================
cat("\n")
cat(paste0(rep("=", 80), collapse = ""))
cat("\n")
cat("开始绘制所有子图...\n")
cat(paste0(rep("=", 80), collapse = ""))
cat("\n\n")

start_time <- Sys.time()

tryCatch({
  p_ab <- plot_ab()
  cat("✓ 图a/b绘制成功\n")
}, error = function(e) {
  cat("✗ 图a/b绘制失败:", e$message, "\n")
  p_ab <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "图a/b绘制失败", size = 10, color = "red") +
    theme_void()
})

tryCatch({
  p_c <- plot_c()
  cat("✓ 图c绘制成功\n")
}, error = function(e) {
  cat("✗ 图c绘制失败:", e$message, "\n")
  p_c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "图c绘制失败", size = 10, color = "red") +
    theme_void()
})

tryCatch({
  p_d <- plot_d()
  cat("✓ 图d绘制成功\n")
}, error = function(e) {
  cat("✗ 图d绘制失败:", e$message, "\n")
  p_d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "图d绘制失败", size = 10, color = "red") +
    theme_void()
})

tryCatch({
  p_e <- plot_e()
  cat("✓ 图e绘制成功\n")
}, error = function(e) {
  cat("✗ 图e绘制失败:", e$message, "\n")
  p_e <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "图e绘制失败", size = 10, color = "red") +
    theme_void()
})

# 使用cowplot强制底部基线对齐
cat("\n使用cowplot进行基线强制对齐...\n")

bottom_row <- cowplot::plot_grid(
  p_c, p_d, p_e,
  nrow = 1,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1, 1)
)

final_layout <- cowplot::plot_grid(
  p_ab,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.65, 1),
  align = "v"
)

# ==================== 6. 保存最终组合图 ====================
cat("\n")
cat(paste0(rep("=", 80), collapse = ""))
cat("\n")
cat("保存最终组合图...\n")
cat(paste0(rep("=", 80), collapse = ""))
cat("\n")

final_figures_dir <- file.path(output_main_dir, "Final_Figures")
if (!dir.exists(final_figures_dir)) dir.create(final_figures_dir, recursive = TRUE)

cat("\n第一步：保存PDF文件...\n")
tryCatch({
  ggsave(final_pdf, final_layout, width = 24, height = 18, units = "in", bg = "white")
  cat("✅ PDF保存完成:", final_pdf, "\n")
}, error = function(e) {
  cat("❌ PDF保存失败:", e$message, "\n")
})

# ==================== 7. 转换PDF到其他格式 ====================
cat("\n第二步：转换PDF到PNG和TIFF格式...\n")
if (file.exists(final_pdf)) {
  image_files <- convert_pdf_to_images(final_pdf, final_figures_dir,
                                        "Final_Combined_Figure", c("png", "tiff", "jpeg"))

  summary_info <- data.frame(
    File_Type = c("PDF", "PNG", "TIFF", "JPEG"),
    File_Path = c(final_pdf,
                  ifelse(!is.null(image_files$png), image_files$png, "Not generated"),
                  ifelse(!is.null(image_files$tiff), image_files$tiff, "Not generated"),
                  ifelse(!is.null(image_files$jpeg), image_files$jpeg, "Not generated")),
    Status = c(ifelse(file.exists(final_pdf), "Generated", "Failed"),
               ifelse(!is.null(image_files$png) && file.exists(image_files$png), "Generated", "Failed"),
               ifelse(!is.null(image_files$tiff) && file.exists(image_files$tiff), "Generated", "Failed"),
               ifelse(!is.null(image_files$jpeg) && file.exists(image_files$jpeg), "Generated", "Failed")),
    stringsAsFactors = FALSE
  )
  save_intermediate_data(summary_info, "Final_Figure_File_Summary.csv", "最终图片文件汇总")
} else {
  cat("❌ PDF文件不存在，无法转换\n")
}

# ==================== 8. 绘制额外的T检验图（针对第二个表型）====================
cat("\n第三步：对第二个表型 testTKW.txt 进行T检验分析并绘图...\n")
tryCatch({
  p_ttest <- plot_extra_t_test()
  cat("✓ T检验图绘制成功\n")
}, error = function(e) {
  cat("✗ T检验图绘制失败:", e$message, "\n")
})

# ==================== 9. 创建分析总结报告 ====================
cat("\n第四步：创建分析总结报告...\n")
end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

summary_report <- c(
  paste0(rep("=", 80), collapse = ""),
  "TKW GWAS和QTL分析总结报告（色盲友好版本）",
  paste0(rep("=", 80), collapse = ""),
  paste("生成时间:", Sys.time()),
  paste("运行时长:", round(as.numeric(run_time), 2), "分钟"),
  "",
  "✅ λ符号显示修复：",
  "1. 使用最简单的显示方法：lambda = 1.026",
  "2. 避免使用Unicode字符或表达式解析导致的兼容性问题",
  "3. 曼哈顿图和QQ图中的λ符号现在应该能正确显示为'lambda'",
  "",
  "✅ 色盲友好修改内容：",
  "1. 曼哈顿图配色：",
  "   - A组染色体（1A-7A）：Set2颜色1",
  "   - B组染色体（1B-7B）：Set2颜色2",
  "   - D组染色体（1D-7D）：Set2颜色3",
  "   - ChrUN：灰色 #9E9E9E",
  "   - 只使用3种颜色+灰色，避免红绿色盲混淆",
  "2. 图c（斜率和截距相关性图）：",
  "   - 非QTL点：灰色（gray60）",
  "   - QTL点：橙色（#D55E00，对红绿色盲友好）",
  "   - 相关系数标注简化为'r = 0.017'（无p值）",
  "3. 图d（QTL效应随环境变化）：",
  "   - Environment图例标题改为Environments（复数）",
  "   - 使用Set2分类颜色，避免渐变色",
  "   - X轴单位改为摄氏度",
  "4. 图e（单倍型拟合线）：",
  "   - 去除Environment图例，只保留Haplotypes图例",
  "   - 图例位置与d图保持一致（右上角）",
  "   - 使用Set2分类颜色，避免渐变色",
  "   - X轴单位改为摄氏度",
  "5. 新增分析：使用第二个表型 testTKW.txt 进行单倍型分析，并对截距和斜率进行T检验，生成箱线图（含显著性字母）",
  "6. 基因型文件更新为 5Kfiltered_alt_geno_imputed_renamed.vcf.gz",
  "7. 所有颜色方案均考虑色盲读者，使用高对比度分类颜色",
  "",
  "1. 输出目录结构:",
  paste("  主目录:", output_main_dir),
  paste("  - Intermediate_Results: 中间分析结果"),
  paste("  - Manhattan_Plots: 曼哈顿图"),
  paste("  - QTL_Plots: QTL相关图表"),
  paste("  - Haplotype_Analysis: 单倍型分析结果（包含T检验图）"),
  paste("  - Correlation_Analysis: 相关分析结果"),
  paste("  - Final_Figures: 最终组合图"),
  "",
  "2. 主要分析内容:",
  "   a) Intercept和Slope的曼哈顿图（使用Set2分类颜色，λ符号显示为'lambda'）",
  "   b) 斜率和截距-log10(P)值的相关分析（灰/橙配色）",
  "   c) QTL效应随环境因子(DTR&mean&PostFlowering1_14)的变化",
  "   d) 单倍型在不同环境下的表型响应（原始表型 TKW_mean_table.txt）",
  "   e) 对第二个表型 testTKW.txt 的单倍型截距和斜率进行T检验（箱线图+显著性字母）",
  "",
  "3. 色盲友好视觉设计:",
  "   - 所有图表使用分类颜色，避免渐变色",
  "   - 高对比度颜色方案，适合红绿色盲读者",
  "   - 图例位置统一（d/e图均在右上角）",
  "   - 使用cowplot::plot_grid强制对齐所有图表",
  "",
  "4. 主要输出文件:",
  paste("  最终组合图PDF:", final_pdf),
  paste("  最终组合图PNG:", ifelse(file.exists(final_png), final_png, "未生成")),
  paste("  最终组合图TIFF:", ifelse(file.exists(final_tiff), final_tiff, "未生成")),
  paste("  Supplementary QQ图S1:", file.path(output_main_dir, "Final_Figures", "Supplementary_Figure_S1_QQplots.pdf")),
  paste("  品种名-基因型-单倍型关系:", file.path(output_main_dir, "Haplotype_Analysis", "Sample_Genotype_Haplotype_Relationship.csv")),
  paste("  T检验图 (testTKW):", file.path(output_main_dir, "Haplotype_Analysis", "TTest_Intercept_Slope_testTKW.pdf")),
  "",
  paste0(rep("=", 80), collapse = "")
)

report_path <- file.path(output_main_dir, "Analysis_Summary_Report.txt")
writeLines(summary_report, report_path)
cat("✅ 分析总结报告保存:", report_path, "\n")

# ==================== 10. 最终输出总结 ====================
cat("\n")
cat(paste0(rep("*", 80), collapse = ""))
cat("\n")
cat("✅ 所有分析完成！λ符号已修复（显示为'lambda'），图表已优化为色盲友好版本！\n")
cat("✅ 已针对第二个表型 testTKW.txt 生成T检验图。\n")
cat(paste0(rep("*", 80), collapse = ""))
cat("\n\n")

cat("修复总结:\n")
cat("============\n")
cat("1. λ符号显示修复：使用最简单的显示方法 'lambda = 1.026'，避免所有兼容性问题\n")
cat("2. 色盲友好配色：所有图表去除渐变色，改用分类颜色（Set2调色板）\n")
cat("3. 曼哈顿图颜色简化：\n")
cat("   - A组染色体：Set2颜色1\n")
cat("   - B组染色体：Set2颜色2\n")
cat("   - D组染色体：Set2颜色3\n")
cat("   - ChrUN：灰色 #9E9E9E\n")
cat("4. 图c优化：\n")
cat("   - 非QTL点：灰色（gray60）\n")
cat("   - QTL点：橙色（#D55E00，对红绿色盲友好）\n")
cat("   - 相关系数标注：简化为'r = 0.017'\n")
cat("5. 图d优化：Environment图例标题改为Environments（复数），X轴单位改为摄氏度\n")
cat("6. 图e优化：去除Environment图例，只保留Haplotypes图例，X轴单位改为摄氏度，图例位置与d图一致\n")
cat("7. 新增T检验图：对第二个表型 testTKW.txt 的单倍型截距和斜率进行显著性分析，生成箱线图并标注显著性字母\n")
cat("8. 基因型文件更新为 5Kfiltered_alt_geno_imputed_renamed.vcf.gz\n")
cat("9. 输出目录: ", output_main_dir, "\n")
cat("10. 主要文件:\n")
cat("   - 最终组合图PDF: ", final_pdf, "\n")
if (file.exists(final_png)) cat("   - 最终组合图PNG: ", final_png, "\n")
if (file.exists(final_tiff)) cat("   - 最终组合图TIFF: ", final_tiff, "\n")
cat("   - 分析总结报告: ", report_path, "\n")
cat("   - T检验图: ", file.path(output_main_dir, "Haplotype_Analysis", "TTest_Intercept_Slope_testTKW.pdf"), "\n")
cat("\n运行时间: ", round(as.numeric(run_time), 2), "分钟\n")
cat(paste0(rep("*", 80), collapse = ""))
cat("\n")

cat("\n输出目录结构:\n")
cat("------------\n")
print(list.dirs(output_main_dir, recursive = TRUE, full.names = FALSE))

cat("\n✅ 脚本执行完毕！所有图表已生成，包含新增T检验图。\n")