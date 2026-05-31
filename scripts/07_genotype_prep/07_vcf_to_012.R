# 加载必要的包，如果未安装则自动安装
if (!require(vcfR)) install.packages("vcfR")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tibble)) install.packages("tibble")
library(vcfR)
library(dplyr)
library(tibble)

# 读取 VCF 文件，替换为实际的 VCF 文件路径
vcf_file <- "C:\\Users\\Lenovo\\Desktop\\小麦千粒重文章\\extracted_G1_to_G339.vcf"  
cat("正在读取VCF文件:", vcf_file, "\n")
vcf <- read.vcfR(vcf_file)

# 显示VCF文件基本信息
cat("\nVCF文件基本信息:\n")
cat("SNP数量:", nrow(vcf@fix), "\n")
cat("样本数量:", ncol(vcf@gt) - 1, "\n")
cat("染色体数量:", length(unique(vcf@fix[, "CHROM"])), "\n")

# 提取基因型数据（GT 字段）并转换为数值矩阵
cat("\n正在提取基因型数据...\n")
gt <- extract.gt(vcf, element = "GT")

# 将基因型转换为数值（0,1,2）
cat("正在转换基因型为数值格式...\n")
myGD <- apply(gt, c(1, 2), function(x) {
  if (all(is.na(x))) return(NA)
  x <- gsub("/", "", x)  # 处理斜杠分隔
  x <- gsub("\\|", "", x)  # 处理竖线分隔
  x <- gsub("-", "", x)  # 处理连字符分隔
  
  # 转换基因型为数值
  if (x %in% c("00", "0|0","0/0", "0-0")) return(0)
  if (x %in% c("01", "10", "0|1", "0/1","1|0","1/0", "0-1", "1-0")) return(1)
  if (x %in% c("11", "1|1","1/1", "1-1")) return(2)
  return(NA)
})

cat("基因型转换完成!\n")
cat("缺失值比例:", round(sum(is.na(myGD)) / (nrow(myGD) * ncol(myGD)) * 100, 2), "%\n")

# 转置矩阵，使行 = 样本，列 = SNP
myGD <- t(myGD)

# 设置行名为样本名，列名为 SNP ID
rownames(myGD) <- colnames(gt)  # 样本名来自 VCF 的列名（即样本 ID）
colnames(myGD) <- vcf@fix[, "ID"]  # SNP 名来自 VCF 的 ID 字段

# 将矩阵转换为数据框，添加 taxa 列作为第一列
cat("\n正在创建GD数据框...\n")
myGD_df <- as.data.frame(myGD)

# 使用tibble包添加taxa列作为第一列
myGD_df <- tibble::rownames_to_column(myGD_df, var = "taxa")

# 显示GD数据框信息
cat("GD数据框维度:", dim(myGD_df), "\n")
cat("GD数据框列数:", ncol(myGD_df), " (taxa +", ncol(myGD_df)-1, "个SNP)\n")
cat("GD数据框行数:", nrow(myGD_df), " (", nrow(myGD_df), "个样本)\n")

# 提取 SNP 元数据（SNP 名称、染色体、物理位置）
cat("\n正在创建GM数据框...\n")
myGM_df <- data.frame(
  SNP = vcf@fix[, "ID"],
  Chromosome = vcf@fix[, "CHROM"],
  Position = as.numeric(vcf@fix[, "POS"]),
  stringsAsFactors = FALSE
)

# 过滤 SNP 名称为空的行
myGM_df <- myGM_df[!is.na(myGM_df$SNP), ]
myGM_df <- myGM_df[myGM_df$SNP != ".", ]

# 显示GM数据框信息
cat("GM数据框维度:", dim(myGM_df), "\n")
cat("染色体分布:\n")
print(table(myGM_df$Chromosome))

# 验证输出
cat("\n=== 验证输出 ===\n")
cat("myGD前10行前5列:\n")
print(myGD_df[1:10, 1:5])
cat("\nmyGM前10行:\n")
print(myGM_df[1:10, ])

# 检查数据质量
cat("\n=== 数据质量检查 ===\n")
cat("1. GD缺失值统计:\n")
missing_count <- apply(myGD_df[, -1], 2, function(x) sum(is.na(x)))
cat("   SNP缺失值范围:", min(missing_count), "-", max(missing_count), "\n")
cat("   平均每个SNP缺失样本数:", round(mean(missing_count), 2), "\n")

cat("\n2. 基因型频率统计:\n")
# 只统计数值型的基因型（排除taxa列）
geno_values <- as.numeric(as.matrix(myGD_df[, -1]))
geno_freq <- table(geno_values, useNA = "always")
names(geno_freq)[is.na(names(geno_freq))] <- "Missing"
print(geno_freq)

# 保存为CSV文件
cat("\n=== 保存文件 ===\n")

# 保存GD文件
gd_filename <- "myGD2.csv"
write.csv(myGD_df, gd_filename, row.names = FALSE, na = "NA")
cat("已保存GD文件:", gd_filename, "\n")
cat("文件大小:", file.info(gd_filename)$size / 1024 / 1024, "MB\n")

# 保存GM文件
gm_filename <- "myGM2.csv"
write.csv(myGM_df, gm_filename, row.names = FALSE, na = "NA")
cat("已保存GM文件:", gm_filename, "\n")
cat("文件大小:", file.info(gm_filename)$size / 1024 / 1024, "MB\n")

# 创建数据摘要报告
cat("\n=== 数据转换完成 ===\n")
cat("输入文件:", vcf_file, "\n")
cat("输出文件:\n")
cat("  - myGD2.csv: 基因型数据文件 (", nrow(myGD_df), "个样本 × ", ncol(myGD_df)-1, "个SNP)\n")
cat("  - myGM2.csv: SNP元数据文件 (", nrow(myGM_df), "个SNP信息)\n")
cat("\n数据统计:\n")
cat("1. 样本数量:", nrow(myGD_df), "\n")
cat("2. SNP数量:", ncol(myGD_df)-1, "\n")
cat("3. 染色体数量:", length(unique(myGM_df$Chromosome)), "\n")
cat("4. 总基因型数:", nrow(myGD_df) * (ncol(myGD_df)-1), "\n")
cat("5. 缺失值比例:", round(sum(is.na(myGD_df[, -1])) / (nrow(myGD_df) * (ncol(myGD_df)-1)) * 100, 2), "%\n")
cat("6. 基因型频率:\n")
for (i in 0:2) {
  count <- sum(geno_values == i, na.rm = TRUE)
  prop <- count / sum(!is.na(geno_values)) * 100
  cat("   基因型", i, ":", count, "(", round(prop, 2), "%)\n")
}

# 可选：保存为文本文件（如果需要）
cat("\n=== 可选：保存为文本文件 ===\n")
save_txt <- FALSE  # 设为TRUE如果需要同时保存为文本文件
if (save_txt) {
  write.table(myGD_df, "myGD2.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(myGM_df, "myGM2.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("已保存为文本文件: myGD2.txt 和 myGM2.txt\n")
}

# 清理内存
cat("\n正在清理内存...\n")
rm(vcf, gt, myGD)
gc()

cat("\n=== 所有操作完成! ===\n")
cat("生成的CSV文件可以直接导入GAPIT进行GWAS分析。\n")
cat("文件格式说明:\n")
cat("1. myGD2.csv: 第一列'taxa'为样本名称，后续各列为SNP基因型(0,1,2)\n")
cat("2. myGM2.csv: 包含SNP、Chromosome、Position三列\n")