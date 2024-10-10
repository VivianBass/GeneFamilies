# Primero, asegúrate de que tienes los datos en el formato correcto.
# Suponiendo que tienes `gene_family_expression_dists_stats` como un dataframe
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")


load(file.path(output_data_dir, "/gene_family_expression_dists.RData"))

# Función para determinar el nivel de significancia
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  # Not significant
}
# Crear un dataframe con las distancias medianas
median_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$median_ortholog_exp_dists,
               gene_family_expression_dists_stats$median_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# t_test_medians <- t.test(distance ~ group, data = median_distances, alternative = 'greater')


# # Crear un dataframe con las distancias medias
mean_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$mean_ortholog_exp_dists,
               gene_family_expression_dists_stats$mean_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# Realizar el t-test para las distancias medias
# t_test_mean_results <- t.test(distance ~ group, data = mean_distances, alternative = 'greater')

# Guardar los resultados de los t-tests en un dataframe
# t_test_summary <- data.frame(
#   analysis = c("Median", "Mean"),
  
#   # Resultados sin corregir
#   statistic_original = c(
#                           t_test_medians$statistic, 
#                           t_test_mean_results$statistic),
#   p_value_original = c(
#                        t_test_medians$p.value, 
#                        t_test_mean_results$p.value),
#   conf_int_lower_original = c(
#                                t_test_medians$conf.int[1], 
#                                t_test_mean_results$conf.int[1]),
#   conf_int_upper_original = c(
#                                t_test_medians$conf.int[2], 
#                                t_test_mean_results$conf.int[2]),
#   # Usar las medias reportadas por los t-tests
#   mean_orthologs = c(
#                       t_test_medians$estimate[1],     
#                       t_test_mean_results$estimate[1]),  
#   mean_paralogs = c(
#                      t_test_medians$estimate[2],     
#                      t_test_mean_results$estimate[2]), 
  
#   # Significancia
#   significance_original = c(
#                     significance_level(t_test_medians$p.value),
#                     significance_level(t_test_mean_results$p.value)),
#   p_value_corrected = p.adjust(c(
#                                   t_test_medians$p.value, 
#                                   t_test_mean_results$p.value), 
#                                 method = "BH")
# )

# # # Añadir la columna de significancia basada en el p-valor corregido
# t_test_summary$significance_corrected <- sapply(t_test_summary$p_value_corrected, significance_level)

boxplot_median <- ggplot(median_distances, aes(x = group, y = distance, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  # Sin outliers en el boxplot para evitar superposición
  geom_jitter(width = 0.1, alpha = 0.3,size = 1) +  # Añadir puntos jitter
  labs(title = "Median Expression Distances",
       y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Orthologs", "Paralogs")),
              map_signif_level = TRUE)



boxplot_mean <- ggplot(mean_distances, aes(x = group, y = distance, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  # Sin outliers en el boxplot para evitar superposición
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +  # Añadir puntos jitter
  labs(title = "Mean Expression Distances",
       y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Orthologs", "Paralogs")),
              map_signif_level = TRUE)

# Guardar los tres boxplots en un solo PDF
pdf(file.path(results_dir,"boxplots_expression_distances_with_jitter.pdf"), height = 15, width = 7)
grid.arrange(boxplot_median, boxplot_mean, ncol = 1)
dev.off()


# # Guardar los resultados del t-test originales y corregidos en un archivo CSV
# write.csv(t_test_summary, file.path(results_dir, "t_test_results_original_and_corrected.csv"), row.names = FALSE)
