library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dplyr)
library(purrr)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(output_data_dir+"/gene_family_expression_dists.RData")
load(output_data_dir+"/boxplot_expression_distances.RData")

significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  # Not significant
}


selected_tissues <- c("mE_mRNA_1182-4H_cells", "mE_mRNA_A_1d_carcass", "mE_mRNA_A_1d_dig_sys","mE_mRNA_A_20d_carcass","mE_mRNA_A_20d_dig_sys")  

# Filtrar los datos para incluir solo los tissues seleccionados
filtered_data <- boxplot_data %>%
  filter(tissue %in% selected_tissues)

save(filtered_data, file = output_data_dir+"/filtered_data_expression_distances.RData")


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




#------------------------------------------------------------------------------


library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)



# Select relevant columns for plotting, 
col_names <- c("Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")

#if you want to dynamically filter columns based on a pattern in their names (e.g., selecting all columns #that contain "count"), you could use the dplyr::select() with helper functions like contains():

col_names <- names(df %>% select(contains("count")))
col_names <- names(df %>% select_if(is.numeric))



plot_title <- "df_median/mean_exp.prof.dists.tissue"
x_axis_label <- ""
y_axis_label <- "exp.prof.dists_per_cluster"

# Pivot the dataframe to long format
df <- your_df

df_long <- df %>% pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Count")

create_boxplot <- function(df_long, plot_title, y_axis_label, output_file = "median_mean_exp.prof.dists.tissue.pdf") {
  library(ggplot2)

  plot_vv <- ggplot(df_long, aes(x = Category, y = Count, fill = Category)) +
    geom_boxplot(color = "black", outlier.shape = 21, outlier.color = "red", outlier.fill = "red") +
    scale_fill_brewer(palette = "Pastel1") +
    labs(title = plot_title, x = " ", y = y_axis_label) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.40),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.32),
      axis.title.y = element_text(size = 14, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10, angle = 35, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "gray80"),
      legend.position = "none"
    )

  ggsave(output_file, plot = plot_vv, width = 8, height = 6)
  
  return(plot_vv)
}

-----------------------------------------------------------------------------------------








# ggplot 2 and ggrain, for Raincloud plots

library(ggplot2)
library(ggrain)
library(dplyr)
library(tidyr)


# Pivot the dataframe to long format

df <- plot_df_exp.prof.dists_ortho_para

df_long <- df %>% pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Value")

# Select relevant columns for plotting
col_names <- c("paralogs", "orthologs")

# Define titles and labels for the plot
plot_title <- "Raincloud_Plot_of_Expression_Profile_Distances_para/ortho"
x_axis_label <- ""
y_axis_label <- "Expression Profile Distances"

create_raincloud_plot <- function(df_long, plot_title, x_axis_label, y_axis_label, output_file = "tissue.pdf") {

  library(ggplot2)
  library(ggrain)

  plot_vv <- ggplot(df_long, aes(x = Category, y = Value, fill = Category)) +
    geom_rain(alpha = 0.6, width = 0.6, justification = 0.5) +
    geom_boxplot(color = "black", width = 0.1, position = position_nudge(x = 0), alpha = 0.4) +
    geom_jitter(aes(color = Category), width = 0.2, size = 0.6, alpha = 0.4) +  
    theme_minimal(base_size = 14) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = plot_title,
      x = x_axis_label,
      y = y_axis_label
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 10), hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text.x = element_text(size = 10, angle = 30, hjust = 0.5, vjust = 0.65),
      axis.text.y = element_text(size = 12),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  ggsave(output_file, plot = plot_vv, width = 8, height = 6)
  
  return(plot_vv)
}




