

library(dplyr)
library(purrr)
library(tidyr)



load("xxx-data.RData")
source("R/expression_funks.R")

# Orthologs, median, mean

dist_matrix <- orthologs.exp.prof.dists
orthologs <- vv_df3_ortho.lst
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

library(dplyr)
library(tidyr)

df_cleaned <- ortholog_stats %>% 
  filter(Value != "Invalid Number")

df_cleaned <- df_cleaned[ , c("Statistic", "Gene_Class", "Value")]

df_cleaned$Value <- as.numeric(df_cleaned$Value)

df_transformed <- df_cleaned %>%
  pivot_wider(names_from = Statistic, values_from = Value)

sorted_df <- df_transformed %>%
  arrange(as.numeric(gsub("cluster_", "", Family)))

orthologs.exp.prof.dists_stats_df <- sorted_df



library(dplyr)
library(purrr)

orthologs <- orthologs.exp.prof.dists

orthologs <- map(orthologs, function(cluster) {
  if (is.atomic(cluster) && length(cluster) == 1) {
    if (!is.na(cluster) && cluster != "NA") cluster else NULL
  } else {
    cluster[!is.na(cluster) & cluster != "NA"]
  }
})

# Remove any empty lists
orthologs <- compact(orthologs)

orthologs_df <- tibble(
  cluster = rep(names(orthologs), sapply(orthologs, length)),
  value = unlist(orthologs)
)

orthologs_df <- orthologs_df %>%
  arrange(as.numeric(gsub("orthologs_cluster_", "", cluster)))

 ---------------------- 
# Assuming orthologs is your list
orthologs_df <- tibble(cluster = names(orthologs),
                       values = orthologs) %>%
  mutate(values = map(values, ~if(is.atomic(.x)) list(.x) else .x)) %>%
  unnest(values) %>%
  group_by(cluster) %>%
  mutate(gene = row_number()) %>%
  ungroup() %>%
  rename(value = values)

df_ortho 

colnames(orthologs_df) <- c("cluster", "dist", "gene")

df <- orthologs_df %>%
  filter(dist != "Invalid Number")

df_cleaned <- orthologs_df %>%
  filter(dist != "Invalid Number")

df_cleaned <- orthologs_df %>%
  filter(!is.na(as.numeric(dist)))

-----------------------------------------------------

library(dplyr)
library(purrr)

paralogs <- paralogs.exp.prof.dists

paralogs <- map(paralogs, function(cluster) {
  if (is.atomic(cluster) && length(cluster) == 1) {
    if (!is.na(cluster) && cluster != "NA") cluster else NULL
  } else {
    cluster[!is.na(cluster) & cluster != "NA"]
  }
})

# Remove any empty lists
paralogs <- compact(paralogs)

paralogs_df <- tibble(
  cluster = rep(names(paralogs), sapply(paralogs, length)),
  value = unlist(paralogs)
)

# If you want to ensure the values are stored as a list
paralogs_df <- paralogs_df %>%
  mutate(values = map(values, ~if(is.atomic(.x)) list(.x) else .x))

paralogs_df <- paralogs_df %>%
  arrange(as.numeric(gsub("paralogs_cluster_", "", cluster)))



merged_df <- bind_cols(paralogs.exp.prof.dists_dataframe, orthologs.exp.prof.dists_dataframe)


paralogs.exp.prof.dists_dataframe
orthologs.exp.prof.dists_dataframe
plot_df_exp.prof.dists_ortho_para 

merged_df <- merged_df[ , c("Family...1", "paralogs.exp.prof.dists", "orthologs.exp.prof.dists")]
conames(merged_df) <- c("Family", "paralogs", "orthologs")
---------------------------------------------------

df_transformed <- paralogs_df %>%
  group_by(cluster) %>%
  mutate(id = row_number()) %>%
  pivot_wider(
    id_cols = cluster,
    names_from = id,
    values_from = value,
    names_prefix = "value_"
  ) %>%
  select(cluster, sort(colnames(.)))





paralogs_df <- orthologs_df %>%
  mutate(cluster = paste0("cluster_", cluster))

# Create dataframe from paralogs list
paralogs_df <- tibble(cluster = names(paralogs),
                      values = paralogs) %>%
  mutate(values = map(values, ~if(is.atomic(.x)) list(.x) else .x)) %>%
  unnest(values) %>%
  group_by(cluster) %>%
  mutate(gene = row_number()) %>%
  ungroup() %>%
  rename(value = values)





# Paralogs, median, mean

# inputs
dist_matrix <- paralogs.exp.prof.dists
paralogs <- vv_df3_para.lst
gene_classes <- list(paralogs = paralogs)

calculate_statistics <- function(dist_matrix, gene_classes, fam_name = NA, stats = c("median", "mean")) {

  dist_matrix <- as.matrix(dist_matrix)
  available_genes <- rownames(dist_matrix)

  results <- map_dfr(names(gene_classes$paralogs), function(cluster_name) {
    
    genes_in_class <- intersect(unique(gene_classes$paralogs[[cluster_name]]), available_genes)
    
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

paralog_stats <- map_dfr(names(paralogs.exp.prof.dists), function(fam_name) {
  calculate_statistics(paralogs.exp.prof.dists[[fam_name]], gene_classes, fam_name)
})

library(dplyr)
library(tidyr)

df_cleaned <- paralog_stats %>%
  filter(Value != "Invalid Number")

df_cleaned <- df_cleaned[ , c("Family", "Statistic", "Gene_Class", "Value")]

df_cleaned$Value <- as.numeric(df_cleaned$Value)

df_transformed <- df_cleaned %>%
  pivot_wider(names_from = Statistic, values_from = Value)

sorted_df <- df_transformed %>%
  arrange(as.numeric(gsub("cluster_", "", Family)))

paralogs.exp.prof.dists_stats_df <- sorted_df

paralogs.exp.prof.dists_stats_df <- paralogs.exp.prof.dists_stats_df[ , c("Gene_Class", "median", "mean")]

colnames(paralogs.exp.prof.dists_stats_df) <- c("Family", "median", "mean")

paralogs.exp.prof.dists_stats_df <- paralogs.exp.prof.dists_stats_df %>%
  arrange(as.numeric(gsub("cluster_", "", Family)))








----------------------------------------------------------------------
df <- merged_df %>% select(-median_paralogs.exp.prof.dists.y, -mean_paralogs.exp.prof.dists.y)

merged_df <- left_join(merged_df, filtered_df, by = "Family")


filtered_df <- orthologs.exp.prof.dists_stats_df %>% distinct(Family, .keep_all = TRUE)

colnames(orthologs.exp.prof.dists_stats_df) <- c("Family", "median_orthologs.exp.prof.dists", "mean_orthologs.exp.prof.dists")

orthologs.exp.prof.dists_stats_df

calculate_mean <- function(x) {
  return(mean(x, na.rm = TRUE))
}

# Apply the function to both columns
orthologs.exp.prof.dists_stats_df$median <- sapply(orthologs.exp.prof.dists_stats_df$median, calculate_mean)
orthologs.exp.prof.dists_stats_df$mean <- sapply(orthologs.exp.prof.dists_stats_df$mean, calculate_mean)



#---------------------------------------------------------------------

# Orthologs and Paralogs per Tissue

library(dplyr)
library(purrr)
library(tidyr)

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

# horizontal , tissue names in Header 

df_transformed <- result %>%
  pivot_wider(
    id_cols = Cluster,
    names_from = Tissue,
    values_from = c(Mean, Median),
    names_sep = "_"
  ) %>%
  select(Cluster, sort(colnames(.)))


df_median_mean_paralogs.exp.prof.dists.tissue_vertical
df_median_mean_paralogs.exp.prof.dists.tissue

df_median_mean_orthologs.exp.prof.dists.tissue_vertical
df_median_mean_orthologs.exp.prof.dists.tissue



# additional statistics

result <- map_dfr(names(data), function(cluster) {
  map_dfr(names(data[[cluster]]), function(tissue) {
    tibble(
      Cluster = cluster,
      Tissue = tissue,
      Mean = mean(data[[cluster]][[tissue]], na.rm = TRUE),
      Median = median(data[[cluster]][[tissue]], na.rm = TRUE),
      SD = sd(data[[cluster]][[tissue]], na.rm = TRUE),
      Variance = var(data[[cluster]][[tissue]], na.rm = TRUE),
      Min = min(data[[cluster]][[tissue]], na.rm = TRUE),
      Max = max(data[[cluster]][[tissue]], na.rm = TRUE),
      Range = max(data[[cluster]][[tissue]], na.rm = TRUE) - min(data[[cluster]][[tissue]], na.rm = TRUE),
      IQR = IQR(data[[cluster]][[tissue]], na.rm = TRUE),
      Q1 = quantile(data[[cluster]][[tissue]], 0.25, na.rm = TRUE),
      Q3 = quantile(data[[cluster]][[tissue]], 0.75, na.rm = TRUE),
      Skewness = skewness(data[[cluster]][[tissue]], na.rm = TRUE),
      Kurtosis = kurtosis(data[[cluster]][[tissue]], na.rm = TRUE),
      Count = sum(!is.na(data[[cluster]][[tissue]])),
      CV = sd(data[[cluster]][[tissue]], na.rm = TRUE) / mean(data[[cluster]][[tissue]], na.rm = TRUE)
    )
  })
})




