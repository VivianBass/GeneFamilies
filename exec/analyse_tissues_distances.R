# Cargar la librería dplyr para un manejo más sencillo de datos
library(dplyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Inicializar lista para almacenar los datos de boxplot
boxplot_data <- data.frame(tissue = character(),
                            distance = numeric(),
                            type = character(),
                            stringsAsFactors = FALSE)

# Iterar sobre las familias de genes en gene_family_expression_distances_by_tissue
for (family in names(gene_family_expression_distances_by_tissue)) {
    print(family)
  # Obtener las distancias por tejido
  tissue_data <- gene_family_expression_distances_by_tissue[[family]]
  
  # Iterar sobre cada tejido y recolectar distancias para ortólogos y parálogos
  for (tissue in names(tissue_data$ortholog_expression_dists_by_tissue)) {
    ortholog_dists <- tissue_data$ortholog_expression_dists_by_tissue[[tissue]]
    paralog_dists <- tissue_data$paralog_expression_dists_by_tissue[[tissue]]
    
    # Agregar datos de ortólogos si no está vacío
    if (length(ortholog_dists) > 0) {
      boxplot_data <- rbind(boxplot_data, data.frame(tissue = tissue,
                                                       distance = ortholog_dists,
                                                       type = "Orthologs",
                                                       stringsAsFactors = FALSE))
    }
    
    # Agregar datos de parálogos si no está vacío
    if (length(paralog_dists) > 0) {
      boxplot_data <- rbind(boxplot_data, data.frame(tissue = tissue,
                                                       distance = paralog_dists,
                                                       type = "Paralogs",
                                                       stringsAsFactors = FALSE))
    }
  }
}

# Cargar la librería dplyr para un manejo más sencillo de datos

# Filtrar datos de ortólogos donde la distancia no está vacía, y contar por tissue
tissue_ortholog_counts <- boxplot_data %>%
  filter(type == "Orthologs" & !is.na(distance)) %>%
  group_by(tissue) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Mostrar los tissues con más datos de ortólogos (donde la distancia no es vacía)
head(tissue_ortholog_counts)
tail(tissue_ortholog_counts)

save(boxplot_data, 
     file = output_data_dir+"/boxplot_expression_distances.RData")