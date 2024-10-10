library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")


# Cargar los archivos RData necesarios
load(file.path(output_data_dir, "families.RData"))      # contiene el objeto de familias con columnas Family y Gene
load(file.path(output_data_dir, "orthologsTandems.RData")) # contiene paralogs.lst y orthologs.lst
load(file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))      # contiene las expresiones genéticas por tissue, con los genes en la última columna

# Inicializar las estructuras que vamos a usar
gene_family_expression_distances <- list()  # Para distancias combinadas
gene_family_expression_distances_by_tissue <- list()  # Para distancias por tejido

# Estructura para estadísticas combinadas y por tejido
gene_family_expression_dists_stats <- data.frame(
  Gene_Fam_ID = character(),
  mean_ortholog_exp_dists = numeric(),
  median_ortholog_exp_dists = numeric(),
  mean_paralog_exp_dists = numeric(),
  median_paralog_exp_dists = numeric(),
  stringsAsFactors = FALSE
)

gene_family_expression_dists_stats_by_tissue <- list()  # Para estadísticas por tejido

# Función para medir las distancias de expresión por tejido
measure_expression_distances <- function(gene_ids, expression_data) {
  # Filtrar solo los genes que están en el conjunto de gene_ids
  gene_expression <- expression_data[expression_data$gene %in% gene_ids, ]
  
  # Inicializamos una lista para guardar las distancias por tejido
  tissue_distances <- list()
  
  # Iterar sobre todas las columnas menos la última (que es 'gene')
  for (tissue in colnames(gene_expression)[-ncol(gene_expression)]) {
    # Obtenemos las expresiones de los genes para el tejido actual
    tissue_data <- gene_expression[, tissue, drop = FALSE]
    
    # Calculamos la matriz de distancias entre genes en este tejido
    dist_matrix <- as.matrix(dist(tissue_data, method = "euclidean"))
    
    # Extraemos las distancias (vector) y las guardamos
    tissue_distances[[tissue]] <- dist_matrix[lower.tri(dist_matrix)]
  }
  
  # Unir todas las distancias de todos los tejidos en un solo vector
  all_distances <- unlist(tissue_distances)
  
  return(list(all_distances = all_distances, tissue_distances = tissue_distances))
}

# Iterar sobre cada familia de genes
# Limitar el bucle a las primeras 5 familias únicas
for (f_i in unique(families.genes.df$Family)) {
  
  # Obtener los genes que pertenecen a la familia actual f_i
  fam_gene_ids <- families.genes.df$Gene[families.genes.df$Family == f_i]
  
  # Asegúrate de que fam_gene_ids sea un vector de caracteres
  if (!is.character(fam_gene_ids)) {
    fam_gene_ids <- as.character(fam_gene_ids)
  }
  
  # Obtener los ortólogos de los genes de la familia
  fam_orthologs <- unique(unlist(
    lapply(orthologs.lst, function(cluster_genes) {
      # Verifica que cluster_genes no sea nulo o vacío
      if (!is.null(cluster_genes) && length(cluster_genes) > 0) {
        intersect(as.character(cluster_genes), fam_gene_ids)
      } else {
        return(NULL)
      }
    })
  ))

  # Obtener los parálogos de los genes de la familia
  fam_paralogs <- unique(unlist(
    lapply(paralogs.lst, function(cluster_genes) {
      # Verifica que cluster_genes no sea nulo o vacío
      if (!is.null(cluster_genes) && length(cluster_genes) > 0) {
        intersect(as.character(cluster_genes), fam_gene_ids)
      } else {
        return(NULL)
      }
    })
  ))

  fam_paralogs <- setdiff(fam_paralogs, fam_orthologs)
  
  # Medir distancias de expresión para los ortólogos y parálogos
  ortholog_distances <- measure_expression_distances(fam_orthologs, rna.seq.exp.profils)
  paralog_distances <- measure_expression_distances(fam_paralogs, rna.seq.exp.profils)
  
  # Guardar las distancias combinadas (de todos los tejidos)
  gene_family_expression_distances[[f_i]] <- list(
    ortholog_expression_dists = ortholog_distances$all_distances,
    paralog_expression_dists = paralog_distances$all_distances
  )
  
  # Guardar las distancias por tejido
  gene_family_expression_distances_by_tissue[[f_i]] <- list(
    ortholog_expression_dists_by_tissue = ortholog_distances$tissue_distances,
    paralog_expression_dists_by_tissue = paralog_distances$tissue_distances
  )
  
  # Calcular estadísticas combinadas (media y mediana)
  gene_family_expression_dists_stats <- rbind(gene_family_expression_dists_stats, data.frame(
    Gene_Fam_ID = f_i,
    mean_ortholog_exp_dists = mean(ortholog_distances$all_distances, na.rm = TRUE),
    median_ortholog_exp_dists = median(ortholog_distances$all_distances, na.rm = TRUE),
    mean_paralog_exp_dists = mean(paralog_distances$all_distances, na.rm = TRUE),
    median_paralog_exp_dists = median(paralog_distances$all_distances, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
  
  # Calcular estadísticas por tejido
  tissue_stats <- list()
  for (tissue in names(ortholog_distances$tissue_distances)) {
    tissue_stats[[tissue]] <- data.frame(
      Gene_Fam_ID = f_i,
      tissue = tissue,
      mean_ortholog_exp_dists = mean(ortholog_distances$tissue_distances[[tissue]], na.rm = TRUE),
      median_ortholog_exp_dists = median(ortholog_distances$tissue_distances[[tissue]], na.rm = TRUE),
      mean_paralog_exp_dists = mean(paralog_distances$tissue_distances[[tissue]], na.rm = TRUE),
      median_paralog_exp_dists = median(paralog_distances$tissue_distances[[tissue]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  # Guardar estadísticas por tejido en la lista
  gene_family_expression_dists_stats_by_tissue[[f_i]] <- do.call(rbind, tissue_stats)
}

# Combinar las estadísticas por tejido en un único data frame
all_tissue_stats <- do.call(rbind, gene_family_expression_dists_stats_by_tissue)

# Guardar los resultados en archivos RData
save(gene_family_expression_distances, gene_family_expression_distances_by_tissue,
     gene_family_expression_dists_stats, all_tissue_stats, 
     file = file.path(output_data_dir, "gene_family_expression_dists.RData"))
