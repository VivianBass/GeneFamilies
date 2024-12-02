
# Load environment variables to define directories for output data and results
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution.R")

# Libraries for efficient data handling
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)

# Libraries for data visualization
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# ---------------------------------------------------------------------------

process_regular_angles <- function(data_pattern, loaded_objects) {
    data_names <- loaded_objects[grepl(data_pattern, loaded_objects)]
    
    data_names %>%
        map_df(function(data_name) {
            current_data <- get(data_name)
            type_name <- sub(data_pattern, "", data_name)
            
            tibble(
                Type = type_name,
                Angle = unlist(current_data)
            )
        })
}

# Usage:
# For regular angles
load(file.path(output_data_dir, "exp.prof.angles.RData"))
loaded_objects <- ls()
df_complete_angles <- process_regular_angles(".lst_cos_angles_dists$", loaded_objects)

# For log2 angles
load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))
loaded_objects <- ls()
df_complete_angles_log2 <- process_regular_angles(".lst_cos_angles_dists_log2$", loaded_objects)





create_angle_boxplot <- function(data, type_combinations, test_type, y_breaks, title_prefix = "Complete Expression Angles") {
    ggplot(data, aes(x = Type, y = Angle, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = paste(title_prefix, paste0("(", test_type, ")")), y = "Angle") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
            test.args = list(alternative = "two.sided"),
            map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
            step_increase = 0.05,
            tip_length = 0.005,
            color = "black",
            size = 0.3,
            textsize = 2.5
        ) +
        theme(
            plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
            axis.text.x = element_text(size = 10),
            plot.margin = margin(r = 30)
        )
}

plot_and_save_combined_angles <- function(regular_data, log2_data, results_dir, filename_prefix) {
    # Process regular data
    df_regular_filtered <- regular_data %>% filter(is.finite(Angle))
    types_regular <- unique(df_regular_filtered$Type)
    type_combinations_regular <- combn(types_regular, 2, simplify = FALSE)
    y_breaks_regular <- seq(0, max(df_regular_filtered$Angle, na.rm = TRUE), by = 0.2)
    
    # Process log2 data
    df_log2_filtered <- log2_data %>% filter(is.finite(Angle))
    types_log2 <- unique(df_log2_filtered$Type)
    type_combinations_log2 <- combn(types_log2, 2, simplify = FALSE)
    y_breaks_log2 <- seq(0, max(df_log2_filtered$Angle, na.rm = TRUE), by = 0.2)
    
    # Create all plots
    plot_list <- list(
        create_angle_boxplot(df_regular_filtered, type_combinations_regular, "t.test", y_breaks_regular, "Regular Expression Angles"),
        create_angle_boxplot(df_log2_filtered, type_combinations_log2, "t.test", y_breaks_log2, "Log2 Expression Angles"),
        create_angle_boxplot(df_regular_filtered, type_combinations_regular, "wilcox.test", y_breaks_regular, "Regular Expression Angles"),
        create_angle_boxplot(df_log2_filtered, type_combinations_log2, "wilcox.test", y_breaks_log2, "Log2 Expression Angles")
    )
    
    # Save combined plots
    output_pdf <- file.path(results_dir, paste0(filename_prefix, "_combined.pdf"))
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=2, ncol=2, top=""),
           width = 16, height = 12,
           device = "pdf")
    
    return(plot_list)
}


combined_angle_plots <- plot_and_save_combined_angles(
    df_complete_angles,
    df_complete_angles_log2,
    results_dir,
    "complete_boxplots_expression_angles"
)
