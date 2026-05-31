#!/usr/bin/env Rscript

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
  
  # 特殊分析包（图e需要）
  library(vcfR)
  library(RColorBrewer)
  library(car)
  library(multcomp)
})

# ==================== 全局参数设置 ====================
base_dir <- "/mnt/7t_storage/zhangcl/TKW"
setwd(base_dir)
cat("当前工作目录:", getwd(), "\n")

# 输出最终组合图路径
final_png <- "Final_Combined_Figure.png"
final_pdf <- "Final_Combined_Figure.pdf"

# 全局标签样式参数（核心修改：增大偏移量）
label_params <- list(
  ab_label_y = 0.95,        # ab标签统一y轴高度（相对坐标）
  cde_label_y = 0.95,       # cde标签统一y轴高度（相对坐标）
  label_size = 5,
  label_fontface = "bold",
  label_hjust = 1.8,        # 核心修改：标签向外偏移更多（原1.2→1.8）
  qtl_label_fontface = "bold.italic"  # QTL名称：斜体+加粗
)

# ==================== 1. 定义图a/b (曼哈顿图) 绘制函数 ====================
plot_ab <- function() {
  # ============ 参数设置 ============
  # 斜率(Slope) GWAS文件
  slope_gwas_file <- "FW_GWAS_Results/GWAS_Slope_para/GAPIT.Association.GWAS_Results.GWAS_Slope_para.BLINK.Trait(NYC).csv"
  slope_qtl_file  <- "SlopeQTL.csv"
  
  # 截距(Intcp) GWAS文件
  intcp_gwas_file <- "FW_GWAS_Results/GWAS_Intcp_para_adj/GAPIT.Association.GWAS_Results.GWAS_Intcp_para_adj.BLINK.Trait(NYC).csv"
  intcp_qtl_file  <- "IntcpQTL.csv"
  
  # 显著性阈值
  p_thresh <- 0.0001
  log_thresh <- -log10(p_thresh)
  
  # ============ 染色体映射 ============
  chr_map <- data.frame(
    Chr = 1:22,
    ChrName = c(
      "Chr1A","Chr1B","Chr1D",
      "Chr2A","Chr2B","Chr2D",
      "Chr3A","Chr3B","Chr3D",
      "Chr4A","Chr4B","Chr4D",
      "Chr5A","Chr5B","Chr5D",
      "Chr6A","Chr6B","Chr6D",
      "Chr7A","Chr7B","Chr7D",
      "ChrUN"
    ),
    stringsAsFactors = FALSE
  )
  
  # ============ 定义颜色方案 ============
  # A组染色体 - 蓝色系
  blue_palette <- colorRampPalette(c("#4A7BFF", "#6C9EFF", "#8EC1FF"))(7)
  # B组染色体 - 绿色系
  green_palette <- colorRampPalette(c("#36A852", "#58CA74", "#7AEC96"))(7)
  # D组染色体 - 橙色系
  orange_palette <- colorRampPalette(c("#FF8C42", "#FFA566", "#FFBE8A"))(7)
  
  # 创建颜色映射
  color_mapping <- c(
    # A组
    "Chr1A" = blue_palette[1], "Chr2A" = blue_palette[2], "Chr3A" = blue_palette[3],
    "Chr4A" = blue_palette[4], "Chr5A" = blue_palette[5], "Chr6A" = blue_palette[6],
    "Chr7A" = blue_palette[7],
    
    # B组
    "Chr1B" = green_palette[1], "Chr2B" = green_palette[2], "Chr3B" = green_palette[3],
    "Chr4B" = green_palette[4], "Chr5B" = green_palette[5], "Chr6B" = green_palette[6],
    "Chr7B" = green_palette[7],
    
    # D组
    "Chr1D" = orange_palette[1], "Chr2D" = orange_palette[2], "Chr3D" = orange_palette[3],
    "Chr4D" = orange_palette[4], "Chr5D" = orange_palette[5], "Chr6D" = orange_palette[6],
    "Chr7D" = orange_palette[7