require(GeneFamilies)
options(mc.cores = getMcCores())

library(dotenv)
library(dplyr)
library(purrr)
library(tidyr)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfExpressionProfileDistances.R path/2/GeneFamilies")

load("experiments/RPKM_flybase/data/ExpressionProfileDistances.RData")


# Function to calculate statistics (mean, median) for orthologs or paralogs expression profile distances
calculate_expression_profile_statistics <- function(data) {
 
  result <- map_dfr(names(data), function(name) {
    # Extract the distance matrix for current element
    dist_matrix <- as.matrix(data[[name]])
    
    tibble(
      Family = name,
      Mean = mean(dist_matrix, na.rm = TRUE),
      Median = median(dist_matrix, na.rm = TRUE),
      Max = max(dist_matrix, na.rm = TRUE),
      Min = min(dist_matrix[dist_matrix > 0], na.rm = TRUE)
    )
  })

  result <- result %>%
    filter_all(all_vars(!is.na(.) & !is.infinite(.)))
  
  return(result)
}

# orthologs
df_median_mean_orthologs <- calculate_expression_profile_statistics(input.args[[1]])

# paralogs
df_median_mean_paralogs <- calculate_expression_profile_statistics(input.args[[2]])


# Function to calculate statistics (mean, median) for orthologs or paralogs per tissue
calculate_expression_profile_statistics_per_tissue <- function(data) {
  
  result <- map_dfr(names(data), function(cluster) {
    map_dfr(names(data[[cluster]]), function(tissue) {
      tibble(
        Cluster = cluster,
        Tissue = tissue,
        Mean = mean(data[[cluster]][[tissue]], na.rm = TRUE),
        Median = median(data[[cluster]][[tissue]], na.rm = TRUE)
      )
    })
  })
  
  # Transform result output to horizontal (tissue names in Header)
  df_transformed <- result %>%
    pivot_wider(
      id_cols = Cluster,
      names_from = Tissue,
      values_from = c(Mean, Median),
      names_sep = "_"
    ) %>%
    select(Cluster, sort(colnames(.)))
  
  return(df_transformed)
}

# orthologs per tissue:
df_median_mean_orthologs_tissue <- calculate_expression_profile_statistics_per_tissue(input.args[[3]])

# paralogs per tissue:
df_median_mean_paralogs_tissue <- calculate_expression_profile_statistics_per_tissue(input.args[[4]])

# Save:
save(df_median_mean_paralogs, 
     df_median_mean_orthologs,
     df_median_mean_paralogs_tissue, 
     df_median_mean_orthologs_tissue,
     file = file.path(output_data_dir, "expression_profile_distances_statistics.RData"))

message("DONE")



