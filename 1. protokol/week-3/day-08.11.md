**Date**: 08.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

### Tasks:

- Improved scripts, especially in computing Euclidean distances and statistics.
- Enhanced plotting scripts.
- Updated expression profiles to include `FBpp` protein sequences instead of `FBgn` names.
- Executed all scripts across test directories, generating plots in each results folder.

### Doubts and Issues:

- Should I compute Euclidean distances for gene families as well?
- Uncertainty about the accuracy of t-test results—they do not match the significance values shown in plots.
- Plots use `geom_signif()` from `ggsignif`, which directly incorporates significance levels, bypassing the t-test results; t-test results are saved separately in a CSV.
  
    ```R
    significance_level <- function(p) {
        if (p < 0.001) return("***")
        else if (p < 0.01) return("**")
        else if (p < 0.05) return("*")
        else return("ns")  
    }
    ```

- Current expression profiles lack data for con-orthologs etc.; only `specific in-paralogs` and `in-paralogs` are covered.
- therefor con-orthologs etc. are excluded due to the absence of matching `FBgn` names in the provided expression profiles.
- Tried using expression profiles from the diet paper, but they have even less overlap with the provided gene group proteins.
- Need expression profiles that better match our gene groups data.

### Next Steps:

- Adjust calculation of logarithmic distances for expression values.
- Complete unit tests for all functions and define test scenarios.
- Adjust t-test results to align with the displayed significance levels of the plot, .

---

**Code:**

- Added code in `test_compute_exp.prof.dists_statistics.R` to validate the availability of required files.

```R
valid_groups <- bind_rows(
    df_mean.dists %>% mutate(source = "mean"),
    df_median.dists %>% mutate(source = "median")
) %>%
    group_by(Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~pull(.x, Type))
```

- Current code in `plot_exp.prof.dists_distribution.R`.

```R
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

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution.R")

# Load expression profile distance statistics
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_dists_stats$", loaded_objects)]

# Initialize empty data frames to store the mean and median distances for each dataset
df_mean.dists <- data.frame()
df_median.dists <- data.frame()

# Process each dataset to extract and compile mean and median distance information
for (data_name in data_names) {
    current_data <- get(data_name)
    type_name <- sub(".lst_dists_stats", "", data_name)
    
    # Compile mean distance data into a tibble, adding columns for type, cluster, and distance
    temp_mean_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists <- bind_rows(df_mean.dists, temp_mean_df)
    
    # Compile median distance data into a tibble, adding columns for type, cluster, and distance
    temp_median_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists <- bind_rows(df_median.dists, temp_median_df)
}

# -----------------------------------------------------------------------------

# Function to label significance levels based on p-values
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# Check if there is enough data to perform t-tests on both mean and median distances
# Ensure each group has more than one observation for valid t-testing
valid_groups <- bind_rows(
    df_mean.dists %>% mutate(source = "mean"),
    df_median.dists %>% mutate(source = "median")
) %>%
    group_by(Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~pull(.x, Type))

# Perform t-tests with error handling to account for potential issues
t_test_results <- tryCatch({

    # Perform t-test on median distances if there are enough groups
    if (length(valid_groups$median) >= 2) {
        t_test_median <- df_median.dists %>%
            filter(Type %in% valid_groups$median) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median")
    } else {
        message("Not enough valid groups for median t-test")
    }
    
    # Perform t-test on mean distances if there are enough groups
    if (length(valid_groups$mean) >= 2) {
        t_test_mean <- df_mean.dists %>%
            filter(Type %in% valid_groups$mean) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean")
    } else {
        message("Not enough valid groups for mean t-test")
    }
    
    # Save results if both t-tests were completed successfully
    if (exists("t_test_median") && exists("t_test_mean")) {
        t_test_summary <- bind_rows(t_test_median, t_test_mean)
        write.csv(t_test_summary, file.path(results_dir, "t_test_summary.csv"), row.names = FALSE)
        message("T-test summary exported to CSV")
    }
    
    # Return results as a list for further reference
    list(median = if(exists("t_test_median")) t_test_median else NULL,
         mean = if(exists("t_test_mean")) t_test_mean else NULL,
         summary = if(exists("t_test_summary")) t_test_summary else NULL)
}, error = function(e) {
    message("Error in t-tests: ", e$message)
    return(NULL)
})

# -----------------------------------------------------------------------------

# Define combinations of types for pairwise significance annotations in boxplots
types <- unique(df_mean.dists$Type)
type_combinations <- combn(types, 2, simplify = FALSE)

# Plot the boxplot with significance annotations for mean distances
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_mean.dists$Distance), by = 0.2)) +
  geom_signif(comparisons = type_combinations, map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )

# Plot the boxplot with significance annotations for median distances
boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_median.dists$Distance), by = 0.2)) +
  geom_signif(comparisons = type_combinations, map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )

# Save both plots to a multi-page PDF
output_pdf <- file.path(results_dir, "boxplots_expression_distances.pdf")

# Combine plots into a list
plot_list <- list(boxplot_median, boxplot_mean)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)
```


- Current code in `plot_exp.prof.dists_distribution_tissue.R`


```R
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
```