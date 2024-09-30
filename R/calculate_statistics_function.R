

# Function to calculate statistics for a given cluster, cluster Vector (mean, median)

calculate_statistics <- function(dist_matrix, gene_classes, fam_name = NA, stats = c("median", "mean")) {

  library(dplyr)
  library(purrr)
  library(tidyr)

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

# usage:
# ortholog_stats <- map_dfr(names(orthologs.exp.prof.dists), function(fam_name) {
#  calculate_statistics(orthologs.exp.prof.dists[[fam_name]], gene_classes, fam_name)
#})


