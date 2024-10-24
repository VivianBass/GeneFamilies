require(GeneFamilies)
options(mc.cores = getMcCores())

library(dotenv)
library(dplyr)
library(purrr)
library(tidyr)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfExpressionProfileDistances.R path/2/GeneFamilies")

load("experiments/RPKM_flybase/data/ExpressionProfileDistances.RData")

input.args <- list()
input.args[[1]] <- orthologs.exp.prof.dists
input.args[[2]] <- paralogs.exp.prof.dists
input.args[[3]] <- orthologs.exp.prof.dists.tissue
input.args[[4]] <- paralogs.exp.prof.dists.tissue

message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfExpressionProfileDistances.R path/2/GeneFamilies")

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

# ------------------------------------------------------------------------------------


# how exactly do the family statistics look like ??
# Families:
 families.exp.prof.dists.orth.dist <- setNames(mclapply(names(families.exp.prof.dists), 
     function(x) classSpecificExpressionProfileDists_test(families.exp.prof.dists[[x]], 
         gene.classes = gene.classes, fam.name = x)), names(families.exp.prof.dists))
 families.exp.prof.dists.orth.dist.df <- Reduce(rbind, mclapply(names(families.exp.prof.dists.orth.dist), 
     function(x) families.exp.prof.dists.orth.dist[[x]][["stats"]]))

# Orthologs:
 orthologs.exp.prof.dists.stats <- setNames(mclapply(names(orthologs.exp.prof.dists), 
     function(x) classSpecificExpressionProfileDists_test(orthologs.exp.prof.dists[[x]], 
         gene.classes = gene.classes["Orthologs"], fam.name = x)), names(orthologs.exp.prof.dists))
 orthologs.exp.prof.dists.stats.df <- Reduce(rbind, mclapply(names(orthologs.exp.prof.dists.stats), 
     function(x) orthologs.exp.prof.dists.stats[[x]][["stats"]]))

# Paralogs:
paralogs.exp.prof.dists.stats <- setNames(mclapply(names(paralogs.exp.prof.dists), 
     function(x) classSpecificExpressionProfileDists_test(paralogs.exp.prof.dists[[x]], 
         gene.classes = gene.classes["Paralogs"], fam.name = x)), names(paralogs.exp.prof.dists))
 paralogs.exp.prof.dists.stats.df <- Reduce(rbind, mclapply(names(paralogs.exp.prof.dists.stats), 
     function(x) paralogs.exp.prof.dists.stats[[x]][["stats"]]))


# Save results:
save(paralogs.exp.prof.dists.stats.df, families.exp.prof.dists.orth.dist, families.exp.prof.dists.orth.dist.df, 
    orthologs.exp.prof.dists.stats.df, file = file.path(output_data_dir, "ExpressionProfileDistanceDistributions.RData"))

