# Suponiendo que ya tienes listas las distancias de ortólogos y parálogos
load("data/ExpressionProfileDistances.RData")
load("data/orthologsTandems.RData")

# Orthologs, median, mean
 library(purrr)
 
library(dplyr)
library(tidyr)
library(tibble)

dist_matrix <- orthologs.exp.prof.dists
orthologs <- orthologs.lst
gene_classes <- list(orthologs = orthologs)
 
calculate_statistics <- function(dist_matrix, gene_classes, fam_name = NA, stats = c("median", "mean")) {
 
  dist_matrix <- as.matrix(dist_matrix)
  available_genes <- rownames(dist_matrix)
 
  results <- map_dfr(names(gene_classes$orthologs), function(cluster_name) {
   
    genes_in_class <- intersect(unique(gene_classes$orthologs[[cluster_name]]), available_genes)
   
    if (length(genes_in_class) > 1) {
     
      class_dist <- dist_matrix[genes_in_class, genes_in_class]
     
      stat_results <- map_dfr(stats, function(stat_name) {
        stat_value <- if (length(class_dist) > 0) {
          eval(call(stat_name, as.numeric(class_dist), na.rm = TRUE))
        } else NA      
        tibble(Family = fam_name, Statistic = stat_name, Gene_Class = cluster_name, Value = stat_value)
      })
 
      stat_results <- stat_results %>%
        filter(!is.na(Value) & is.finite(Value))
     
      return(stat_results)
    } else {
     
      return(tibble(Family = fam_name, Statistic = stats, Gene_Class = cluster_name, Value = NA))
    }
  })
 
  return(results)
}
 
 
ortholog_stats <- map_dfr(names(orthologs.exp.prof.dists), function(fam_name) {
  calculate_statistics(orthologs.exp.prof.dists[[fam_name]], gene_classes, fam_name)
})

 
df_cleaned <- ortholog_stats %>%
  filter(Value != "Invalid Number")
 
df_cleaned <- df_cleaned[ , c("Statistic", "Gene_Class", "Value")]
 
df_cleaned$Value <- as.numeric(df_cleaned$Value)
 
df_transformed <- df_cleaned %>%
  pivot_wider(names_from = Statistic, values_from = Value)
 
sorted_df <- df_transformed %>%
  arrange(as.numeric(gsub("cluster_", "", Family)))
 
orthologs.exp.prof.dists_stats_df <- sorted_df

save(orthologs.exp.prof.dists_stats_df, 
    file = "data/StatsExpressionProfileDistances.RData")


message("DONE")
