require(GeneFamilies)
options(mc.cores = getMcCores())  
library(dotenv)  

# Load necessary libraries for data manipulation, plotting, and statistical analysis
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(rstatix)
library(ggpubr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R")

# Automatically select data objects that contain tissue-specific statistics
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names_tissue <- loaded_objects[grepl("_tissue_stats$", loaded_objects)]

# Initialize empty data frames to hold tissue-specific mean and median distance data
df_mean.dists_tissue <- data.frame()
df_median.dists_tissue <- data.frame()

# Process each tissue dataset to extract and organize mean and median distances by tissue and type
for (data_name in data_names_tissue) {
    current_data <- get(data_name)  # Retrieve the data object by name
    
    type_name <- sub(".lst_dists_tissue_stats", "", data_name)  # Remove suffix to get the type name
    
    # Create a temporary data frame for mean distances and add it to the main data frame
    temp_mean_df <- tibble(
        Type = type_name,
        Tissue = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists_tissue <- bind_rows(df_mean.dists_tissue, temp_mean_df)
    
    # Create a temporary data frame for median distances and add it to the main data frame
    temp_median_df <- tibble(
        Type = type_name,
        Tissue = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists_tissue <- bind_rows(df_median.dists_tissue, temp_median_df)
}

# ----------------------------------------------------------------------

# Perform t-tests for mean and median distances by tissue to compare between types (e.g., orthologs vs paralogs)
# Adjust p-values and apply the significance level function

# Function to assign significance level based on p-value
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# Check if sufficient data exists for t-tests by counting observations per tissue and type
valid_groups_tissue <- bind_rows(
    df_mean.dists_tissue %>% mutate(source = "mean"),
    df_median.dists_tissue %>% mutate(source = "median")
) %>%
    group_by(Tissue, Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~select(.x, Tissue, Type))

# Perform t-tests for mean and median distances with error handling to catch any issues
t_test_results_tissue <- tryCatch({
    # Median t-test per tissue
    if (nrow(valid_groups_tissue$median) >= 2) {
        t_test_median_tissue <- df_median.dists_tissue %>%
            semi_join(valid_groups_tissue$median, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median")
    } else {
        message("Not enough valid groups for median tissue t-test")
    }
    
    # Mean t-test per tissue
    if (nrow(valid_groups_tissue$mean) >= 2) {
        t_test_mean_tissue <- df_mean.dists_tissue %>%
            semi_join(valid_groups_tissue$mean, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean")
    } else {
        message("Not enough valid groups for mean tissue t-test")
    }
    
    # Combine results if both t-tests were successful
    if (exists("t_test_median_tissue") && exists("t_test_mean_tissue")) {
        t_test_summary_tissue <- bind_rows(t_test_median_tissue, t_test_mean_tissue)
        write.csv(t_test_summary_tissue, file.path(results_dir, "t_test_summary_tissue.csv"), row.names = FALSE)
        message("Tissue-specific t-test summary exported to CSV")
    }
    
    list(median = if(exists("t_test_median_tissue")) t_test_median_tissue else NULL,
         mean = if(exists("t_test_mean_tissue")) t_test_mean_tissue else NULL,
         summary = if(exists("t_test_summary_tissue")) t_test_summary_tissue else NULL)
}, error = function(e) {
    message("Error in tissue-specific t-tests: ", e$message)
    return(NULL)
})


# ---------------------------------------------------------------------------

# Generate boxplots for mean expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_mean_boxplot_combined.pdf")
tissue_types <- unique(df_mean.dists_tissue$Tissue)

pdf(output_pdf, width = 12, height = 8)

for (tissue in tissue_types) {
  df_tissue <- subset(df_mean.dists_tissue, Tissue == tissue)
  
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Mean Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    scale_y_continuous(breaks = seq(0, max(df_tissue$Distance), by = 0.1)) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10),  
      plot.margin = margin(r = 30)
    )
  
  # Add significance annotations only if there are multiple types
  types <- unique(df_tissue$Type)
  if (length(types) >= 2) {
    type_combinations <- combn(types, 2, simplify = FALSE)
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = type_combinations, map_signif_level = TRUE)
  }
  
  print(boxplot_tissue)  # Output the plot to the PDF
}

dev.off()


# ---------------------------------------------------------------------------

# Generate boxplots for median expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_median_boxplot_combined.pdf")
tissue_types <- unique(df_median.dists_tissue$Tissue)
pdf(output_pdf, width = 12, height = 8)

for (tissue in tissue_types) {
  df_tissue <- subset(df_median.dists_tissue, Tissue == tissue)
  
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Median Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    scale_y_continuous(breaks = seq(0, max(df_tissue$Distance), by = 0.1)) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10),  
      plot.margin = margin(r = 30)
    )
  
  # Add significance annotations only if there are multiple types
  types <- unique(df_tissue$Type)
  if (length(types) >= 2) {
    type_combinations <- combn(types, 2, simplify = FALSE)
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = type_combinations, map_signif_level = TRUE)
  }
  
  print(boxplot_tissue)  # Output the plot to the PDF
}
dev.off()
