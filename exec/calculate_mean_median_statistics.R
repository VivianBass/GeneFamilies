

library(dplyr)
library(purrr)
library(tidyr)


# calculate orthologs and paralogs statistics (median, mean)

# input:
data <- paralogs.exp.prof.dists
data <- orthologs.exp.prof.dists

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

# transform result output to horizontal (tissue names in Header) 

df_transformed <- result %>%
  pivot_wider(
    id_cols = Cluster,
    names_from = Tissue,
    values_from = c(Mean, Median),
    names_sep = "_"
  ) %>%
  select(Cluster, sort(colnames(.)))

# Output, as used in this Package:
# df_median_mean_paralogs.exp.prof.dists
# df_median_mean_orthologs.exp.prof.dists

# ------------------------------------------------------------------------------

# calculate orthologs and paralogs statistics per tissue (median, mean)

# input:
data <- paralogs.exp.prof.dists.tissue
data <- orthologs.exp.prof.dists.tissue

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

# transform result output to horizontal (tissue names in Header) 

df_transformed <- result %>%
  pivot_wider(
    id_cols = Cluster,
    names_from = Tissue,
    values_from = c(Mean, Median),
    names_sep = "_"
  ) %>%
  select(Cluster, sort(colnames(.)))

# Output, as used in this Package:
# df_median_mean_paralogs.exp.prof.dists.tissue
# df_median_mean_orthologs.exp.prof.dists.tissue




