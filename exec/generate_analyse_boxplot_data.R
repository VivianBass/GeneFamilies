# Load the dplyr library for easier data manipulation
library(dplyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Initialize a list to store boxplot data
boxplot_data <- data.frame(tissue = character(),
                           distance = numeric(),
                           type = character(),
                           stringsAsFactors = FALSE)

# Iterate over the gene families in gene_family_expression_distances_by_tissue
for (family in names(gene_family_expression_distances_by_tissue)) {
    print(family)
  # Get distances by tissue
  tissue_data <- gene_family_expression_distances_by_tissue[[family]]
  
  # Iterate over each tissue and collect distances for orthologs and paralogs
  for (tissue in names(tissue_data$ortholog_expression_dists_by_tissue)) {
    ortholog_dists <- tissue_data$ortholog_expression_dists_by_tissue[[tissue]]
    paralog_dists <- tissue_data$paralog_expression_dists_by_tissue[[tissue]]
    
    # Add ortholog data if not empty
    if (length(ortholog_dists) > 0) {
      boxplot_data <- rbind(boxplot_data, data.frame(tissue = tissue,
                                                       distance = ortholog_dists,
                                                       type = "Orthologs",
                                                       stringsAsFactors = FALSE))
    }
    
    # Add paralog data if not empty
    if (length(paralog_dists) > 0) {
      boxplot_data <- rbind(boxplot_data, data.frame(tissue = tissue,
                                                       distance = paralog_dists,
                                                       type = "Paralogs",
                                                       stringsAsFactors = FALSE))
    }
  }
}

# Load the dplyr library for easier data manipulation

# Filter ortholog data where distance is not empty, and count by tissue
tissue_ortholog_counts <- boxplot_data %>%
  filter(type == "Orthologs" & !is.na(distance)) %>%
  group_by(tissue) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Show the tissues with the most ortholog data (where distance is not empty)
head(tissue_ortholog_counts)
tail(tissue_ortholog_counts)

save(boxplot_data, file = file.path(output_data_dir,"/boxplot_expression_distances.RData"))
