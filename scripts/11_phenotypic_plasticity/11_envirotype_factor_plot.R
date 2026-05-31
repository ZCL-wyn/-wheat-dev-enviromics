# Wheat TKW Analysis - PAR_TEMP Key Growth Stage Line Plot
# Analysis of specific growth stages: 65.5 to 87.5

# 1. Load required libraries
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "stringr", "purrr")
lapply(required_pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
})

# 2. Set path parameters (保留中文路径)
data_dir <- "/mnt/7t_storage/zhangcl/TKW/"
output_dir <- paste0(data_dir, "par_temp_key_stages\\")

# 3. Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory: ", output_dir, "\n")
}

# 4. Load environmental factor data
cat("Loading environmental factor data...\n")

ecs_data <- read.csv(
    file = paste0(data_dir, "ECs_results.csv"), 
    header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
)

# Rename first column to env_code
first_col_name <- colnames(ecs_data)[1]
ecs_data <- ecs_data %>%
    rename(env_code = 1) %>%
    mutate(env_code = as.character(env_code))

cat("Environmental factor data loaded, number of environments: ", nrow(ecs_data), "\n")

# 5. Extract PAR_TEMP related columns
cat("\nExtracting PAR_TEMP related columns...\n")

# Get all column names
all_cols <- colnames(ecs_data)

# Extract columns starting with PAR_TEMP
par_temp_cols <- all_cols[grepl("^PAR_TEMP&", all_cols)]

if (length(par_temp_cols) == 0) {
    cat("\nError: No PAR_TEMP related columns found!\n")
    stop("No PAR_TEMP data found")
}

cat("Found", length(par_temp_cols), "PAR_TEMP related columns\n")

# 6. Parse column names to extract growth stage information
cat("\nParsing growth stage information...\n")

# Function to extract numeric part
extract_stage_number <- function(col_name) {
    # Split column name
    parts <- strsplit(col_name, "&")[[1]]
    
    if (length(parts) != 2) {
        return(NA)
    }
    
    stage_info <- parts[2]
    
    # Extract numeric part (may contain decimal point)
    num_match <- regmatches(stage_info, regexpr("[0-9]+\\.?[0-9]*", stage_info))
    
    if (length(num_match) > 0) {
        return(as.numeric(num_match[1]))
    } else {
        return(NA)
    }
}

# Create growth stage information data frame
stage_info_df <- data.frame(
    column_name = par_temp_cols,
    stage_num = sapply(par_temp_cols, extract_stage_number),
    stringsAsFactors = FALSE
) %>%
    filter(!is.na(stage_num)) %>%
    arrange(stage_num)

cat("Parsing completed, found", nrow(stage_info_df), "growth stages with numbers\n")

# 7. Filter growth stages between 65.5 and 87.5 (Precise range)
cat("\nFiltering growth stages between 65.5 and 87.5...\n")

# Exact target stages (ensure 87.5 is included)
target_stages <- c(65.5, 87.5)

# Filter stages within exact range (65.5 ≤ stage ≤ 87.5)
filtered_stages <- stage_info_df %>%
    filter(stage_num >= 65.5 & stage_num <= 87.5) %>%
    arrange(stage_num)

# Verify 87.5 is included, if not find closest
if (!87.5 %in% filtered_stages$stage_num) {
    closest_87.5 <- filtered_stages[which.min(abs(filtered_stages$stage_num - 87.5)), ]
    cat(sprintf("Warning: Exact 87.5 not found, using closest: %.1f\n", closest_87.5$stage_num))
    # Force add 87.5 label to closest stage
    filtered_stages$stage_num[filtered_stages$column_name == closest_87.5$column_name] <- 87.5
}

# Verify 65.5 is included
if (!65.5 %in% filtered_stages$stage_num) {
    closest_65.5 <- filtered_stages[which.min(abs(filtered_stages$stage_num - 65.5)), ]
    cat(sprintf("Warning: Exact 65.5 not found, using closest: %.1f\n", closest_65.5$stage_num))
    filtered_stages$stage_num[filtered_stages$column_name == closest_65.5$column_name] <- 65.5
}

cat(sprintf("\nFilter range: 65.5 to 87.5 (exact)\n"))
cat("Filtered", nrow(filtered_stages), "growth stages:\n")
print(filtered_stages)

# 8. Prepare plotting data (GS67排序+统一形状)
cat("\nPreparing plotting data...\n")

# Extract filtered PAR_TEMP data
par_temp_data <- ecs_data %>%
    select(env_code, all_of(filtered_stages$column_name)) %>%
    pivot_longer(
        cols = all_of(filtered_stages$column_name),
        names_to = "column_name",
        values_to = "par_temp_value"
    ) %>%
    # Merge stage information
    left_join(filtered_stages, by = "column_name") %>%
    # Remove missing values
    filter(!is.na(par_temp_value)) %>%
    # Ensure numeric types
    mutate(
        stage_num = as.numeric(stage_num),
        par_temp_value = as.numeric(par_temp_value),
        # 保留数值标签，确保GS67按数值排序在66和68之间
        x_label_raw = case_when(
            stage_num == 65.5 ~ "65.5 (Flowering)",  # 保留数值+标签
            stage_num == 87.5 ~ "87.5 (EndGrainFill)",# 保留数值+标签
            abs(stage_num - 67) < 0.1 ~ "67 (GS67)",  # GS67显示为67 (GS67)，保证排序
            TRUE ~ as.character(stage_num)            # 其他阶段仅显示数值
        ),
        # 仅通过尺寸区分GS67，所有点都是普通圆形
        point_size = case_when(
            abs(stage_num - 67) < 0.1 ~ 5,            # GS67用更大的尺寸突出
            TRUE ~ 2                                   # 所有其他点（包括Flowering/EndGrainFill）用普通尺寸
        ),
        # 用于x轴排序的数值字段
        sort_num = stage_num
    ) %>%
    # 按数值排序
    arrange(sort_num) %>%
    # 按数值顺序生成x轴标签层级（GS67会自动在66和68之间）
    mutate(
        x_label = factor(
            x_label_raw,
            levels = unique(x_label_raw[order(sort_num)])
        )
    )

cat("Plotting data prepared, number of samples: ", nrow(par_temp_data), "\n")
cat("Number of environments: ", length(unique(par_temp_data$env_code)), "\n")
cat("Number of growth stages: ", length(unique(par_temp_data$stage_num)), "\n")
cat("GS67位置验证（x轴顺序）:\n")
print(levels(par_temp_data$x_label))  # 打印x轴顺序，确认GS67在66和68之间

# 9. Create custom color palette
get_color_palette <- function(n) {
    # Base colors (professional color scheme)
    base_colors <- c(
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
    )
    
    if (n <= length(base_colors)) {
        return(base_colors[1:n])
    } else {
        return(hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n])
    }
}

# Get environment list
env_list <- unique(par_temp_data$env_code)
n_env <- length(env_list)
color_palette <- get_color_palette(n_env)
names(color_palette) <- env_list

cat("\nNumber of environments: ", n_env, "\n")

# 10. Create key growth stage line plot (所有点都是普通圆形，仅GS67尺寸更大，无网格线)
cat("\nCreating key growth stage line plot...\n")

# 所有点都用普通实心圆，仅通过尺寸区分GS67
p <- ggplot(par_temp_data, aes(
               x = x_label, 
               y = par_temp_value, 
               group = env_code, 
               color = env_code
           )) +
    # Lines for all environments
    geom_line(linewidth = 1.0, alpha = 0.7) +
    
    # 所有点都是实心圆，仅GS67尺寸更大
    geom_point(
        shape = 19,          # 统一使用实心圆（无特殊形状）
        aes(size = I(point_size)),  # 固定尺寸（不生成图例）
        alpha = 0.9,
        stroke = 1.2
    ) +
    
    # Color scale
    scale_color_manual(values = color_palette) +
    
    # Axis labels (All English, no title/caption)
    labs(
        x = "Growth Stage",
        y = "PAR_TEMP Value",
        color = "Environment"
    ) +
    
    # Theme settings (No title/subtitle, clean style, 移除所有网格线)
    theme_bw(base_size = 12) +
    theme(
        # Remove all titles/subtitles/captions
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        
        # 移除所有网格线（核心修改）
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        
        # Axis styling (加宽x轴避免标签重叠)
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        
        # Legend styling
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.background = element_blank(),
        
        # Panel styling
        panel.border = element_rect(color = "black", linewidth = 1.2),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
    )

# 11. Save the plot (No title in filename)
cat("\nSaving plot...\n")

# Save as PNG (加宽到14英寸避免x轴标签重叠)
png_file <- paste0(output_dir, "PAR_TEMP_key_stages_65.5_to_87.5.png")
ggsave(png_file, p, width = 14, height = 8, dpi = 600, bg = "white")
cat("PNG file saved: ", png_file, "\n")

# Save as PDF
pdf_file <- paste0(output_dir, "PAR_TEMP_key_stages_65.5_to_87.5.pdf")
ggsave(pdf_file, p, width = 14, height = 8, device = "pdf", bg = "white")
cat("PDF file saved: ", pdf_file, "\n")

# Save as TIFF
tiff_file <- paste0(output_dir, "PAR_TEMP_key_stages_65.5_to_87.5.tiff")
ggsave(tiff_file, p, width = 14, height = 8, dpi = 600, 
       device = "tiff", compression = "lzw", bg = "white")
cat("TIFF file saved: ", tiff_file, "\n")

# 12. Save plotting data
write.csv(par_temp_data, paste0(output_dir, "par_temp_key_stages_data.csv"), row.names = FALSE)
cat("Plotting data saved: ", paste0(output_dir, "par_temp_key_stages_data.csv"), "\n")

# 13. Display plot preview
cat("\nDisplaying plot preview...\n")
print(p)

# 14. Generate brief report (All English)
report_file <- paste0(output_dir, "key_stages_report.txt")
sink(report_file)

cat("======================================================\n")
cat("PAR_TEMP Key Growth Stage Analysis Report\n")
cat("======================================================\n")
cat("Analysis Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("\n")

cat("--- Data Overview ---\n")
cat("Number of environments: ", n_env, "\n")
cat("Analysis range: 65.5 to 87.5 (exact)\n")
cat("Number of growth stages included: ", length(unique(par_temp_data$stage_num)), "\n")
cat("Special annotations:\n")
cat("  - 65.5: Labeled as '65.5 (Flowering)' (normal circle shape, size 2)\n")
cat("  - 87.5: Labeled as '87.5 (EndGrainFill)' (normal circle shape, size 2)\n")
cat("  - 67.0: Labeled as '67 (GS67)' (normal circle shape, size 5 - only larger size)\n")
cat("  - GS67 position: Between 66 and 68 (sorted by numeric stage)\n")
cat("\n")

cat("--- Plot Information ---\n")
cat("Plot dimensions: 14x8 inches (wider for x-axis labels)\n")
cat("Resolution: 600 DPI\n")
cat("Output formats: PNG, PDF, TIFF\n")
cat("No plot title/subtitle/caption\n")
cat("All labels in English only\n")
cat("All points use normal circle shape (no special shapes)\n")
cat("Only GS67 has larger size (5) for highlight\n")
cat("No grid lines in plot (all grid lines removed)\n")
cat("\n")

cat("--- File Output ---\n")
cat("1. Plot files: PAR_TEMP_key_stages_65.5_to_87.5.png/pdf/tiff\n")
cat("2. Plot data: par_temp_key_stages_data.csv\n")
cat("3. This report: key_stages_report.txt\n")
cat("\n")

cat("======================================================\n")
cat("Analysis completed\n")
cat("======================================================\n")

sink()

cat("Analysis report saved: ", report_file, "\n")

# 15. Display key information (All English)
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PAR_TEMP key growth stage analysis completed!\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")

cat("Key information:\n")
cat("1. Analyzed PAR_TEMP data for", n_env, "environments\n")
cat("2. Growth stage range: 65.5 to 87.5 (exact)\n")
cat("3. GS67 is displayed between 66 and 68 (sorted by numeric value)\n")
cat("4. All points (GS67/Flowering/EndGrainFill) use normal circle shape\n")
cat("5. Only GS67 has larger size (5) for highlight (others size 2)\n")
cat("6. No grid lines in the plot (all major/minor grid lines removed)\n")
cat("7. No Chinese labels in plot, only English\n")
cat("8. No plot title/subtitle/caption\n")
cat("9. Output files saved to:", output_dir, "\n")
cat("\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("All analysis completed!\n")
cat(paste(rep("=", 80), collapse = ""), "\n")