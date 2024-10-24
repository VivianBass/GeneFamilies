require(GeneFamilies)
options(mc.cores = getMcCores())

library(ggplot2)
library(ggsignif)
library(gridExtra)

library(dplyr)
library(purrr)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(file.path(output_data_dir, "expression_profile_distances_statistics.RData"))
load(file.path(output_data_dir, "ExpressionProfileDistances.RData"))

load(output_data_dir+"/gene_family_expression_dists.RData")
load(output_data_dir+"/boxplot_expression_distances.RData")

significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

selected_tissues <- c("mE_mRNA_1182-4H_cells", "mE_mRNA_A_1d_carcass", "mE_mRNA_A_1d_dig_sys","mE_mRNA_A_20d_carcass","mE_mRNA_A_20d_dig_sys")  

# Filtrar los datos para incluir solo los tissues seleccionados
filtered_data <- boxplot_data %>% filter(tissue %in% selected_tissues)

save(filtered_data, file = file.path(output_data_dir, "/filtered_data_expression_distances.RData"))


filtered_data <- na.omit(filtered_data)
# Realizar la prueba t-test y guardar los p-valores
t_test_results <- filtered_data %>%
  filter(!is.na(distance)) %>%
  group_by(tissue) %>%
  summarise(
    t_test = list(t.test(distance ~ type, data = ., alternative = 'greater')),
    .groups = "drop"  # Para evitar warnings
  ) %>%
  ungroup() %>%
  mutate(
    p_value = map_dbl(t_test, ~ .x$p.value),        # Extraer el p-value de la prueba t
    statistic = map_dbl(t_test, ~ .x$statistic)     # Extraer la estadística de la prueba t
  ) %>%
  select(-t_test)  # Eliminar la columna con los objetos de prueba t para limpiar los resultados

# Aplicar la corrección de p-valores
t_test_results <- t_test_results %>%
  mutate(
    p_value_adjusted = p.adjust(p_value, method = "BH"), 
    significance_original = sapply(p_value, significance_level), 
    significance_adjusted = sapply(p_value_adjusted, significance_level)  
  )

write.csv(t_test_results, results_dir+"/t_test_results_by_tissue.csv", row.names = FALSE)



# Comprobar que ambos grupos están presentes
p <- ggplot(filtered_data, aes(x = tissue, y = distance, fill = type)) +
  geom_boxplot() +
  labs(title = "Expression Distances by Tissue",
       x = "Tissue",
       y = "Distance") +
  scale_fill_manual(values = c("Orthologs" = "blue", "Paralogs" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

for (i in 1:nrow(t_test_results)) {
  p <- p + annotate("text", 
                    x = i,  
                    y = max(filtered_data$distance, na.rm = TRUE) + 0.02, 
                    label = t_test_results$significance_adjusted[i], 
                    size = 5, 
                    color = "black")  
}

pdf(results_dir+"/boxplot_expression_distances_by_tissue.pdf", height = 15, width = 7)
print(p)
dev.off()

