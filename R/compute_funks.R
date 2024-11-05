
#' Compute Euclidean Distances Between Gene Expression Profiles
#'
#' @description Calculates pairwise Euclidean distances across gene expression profiles for a given set of genes, aggregated across all specified tissues.
#'
#' @param gene.accessions A list of gene accession identifiers for which distances should be computed.
#' @param expression.profiles A data frame or tibble containing gene expression data, with gene IDs in one column and expression values in other columns. Default is `rna.seq.exp.profils`.
#' @param expr.prof.gene.col The column name in `expression.profiles` containing gene IDs. Default is `"FBpp_ID"`.
#' @param tissues A character vector of column names in `expression.profiles` that represent tissue-specific expression values. By default, includes all columns except `expr.prof.gene.col`.
#' @param dist.method Distance calculation method. Default is `"euclidean"`.
#' @return A vector representing the distance matrix, or `NA` if there are fewer than two expression profiles.
#' @examples
#' # Compute overall distances for specified genes across all tissues
#' distances <- exp.prof.dists(gene.accessions = list("gene1", "gene2", "gene3"))
#'
exp.prof.dists <- function(gene.accessions, expression.profiles = rna.seq.exp.profils,
                  expr.prof.gene.col = "FBpp_ID", 
                  tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)), 
                  dist.method = "euclidean") {
    
    # Get all genes and filter expression profiles
    all_genes <- unlist(gene.accessions)
    exp.profs <- expression.profiles %>% 
                 filter(!!sym(expr.prof.gene.col) %in% all_genes) %>% 
                 column_to_rownames(expr.prof.gene.col)
    
    if (nrow(exp.profs) > 1) {

        dist_matrix <- exp.profs %>% 
                   select(all_of(tissues)) %>% 
                   as.matrix() %>%
                   dist(method = dist.method) %>%
                   as.vector()

        return(dist_matrix)

    } else {NA}
}

#' Compute Euclidean Distances Between Gene Expression Profiles by Tissue
#'
#' @description Computes pairwise Euclidean distances between expression profiles of specified genes, optionally by each tissue. Useful for analyzing gene similarity within individual conditions.
#'
#' @param gene.accessions A vector of gene accession identifiers for which distances should be computed.
#' @param expression.profiles A data frame or tibble containing gene expression data, with gene IDs in one column and expression values in other columns. Default is `rna.seq.exp.profils`.
#' @param expr.prof.gene.col The column name in `expression.profiles` containing gene IDs. Default is `"FBpp_ID"`.
#' @param tissues A character vector of column names in `expression.profiles` that represent tissue-specific expression values. By default, includes all columns except `expr.prof.gene.col`.
#' @param dist.method Distance calculation method. Default is `"euclidean"`.
#' @return A list of distance vectors, one per tissue, or `NA` if there are fewer than two expression profiles.
#' @examples
#' # Compute per-tissue distances for specified genes
#' tissue_distances <- exp.prof.dists_tissue(gene.accessions = c("gene1", "gene2", "gene3"))
#'
exp.prof.dists_tissue <- function(gene.accessions, expression.profiles = rna.seq.exp.profils,
                                 expr.prof.gene.col = "FBpp_ID",
                                 tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)), 
                                 dist.method = "euclidean") {
    
    all_genes <- unlist(gene.accessions)
    exp.profs <- expression.profiles %>% 
        filter(!!sym(expr.prof.gene.col) %in% all_genes) %>% 
        column_to_rownames(expr.prof.gene.col)
    
    if (nrow(exp.profs) > 1) {

        tissue_distances <- tissues %>% 
            set_names() %>%
            map(~{
                exp.profs %>% 
                    select(all_of(.x)) %>% 
                    dist(method = dist.method) %>%
                    as.vector()
            })
        
        tissue_distances
    } else {NA}
}


#' Add Object to List if Exists
#'
#' Checks if a variable exists in the environment and, if so, retrieves it,
#' adds it to the specified list with a given name, and logs a message about its creation.
#'
#' @param list_obj A list to which the variable should be added if it exists.
#' @param obj_name A character string indicating the name to assign to the object within the list.
#' @param var_name A character string indicating the name of the variable to check for in the environment.
#'
#' @return The updated list with the variable added if it exists; otherwise, returns the list unchanged.
#' @examples
#' my_list <- list()
#' my_list <- add_to_list_if_exists(my_list, "example_obj", "example_var")
add_to_list_if_exists <- function(list_obj, obj_name, var_name) {
    if (exists(var_name)) {
        list_obj[[obj_name]] <- get(var_name)
        cat(sprintf("Object '%s' created successfully.\n", obj_name))
    } else {
        cat(sprintf("Object '%s' does not exist and was not created.\n", obj_name))
    }
    return(list_obj)
}

#' Calculate Statistics for Expression Profile Distances
#'
#' Computes the mean and median of the expression profile distances for each gene family
#' within a given dataset.
#'
#' @param data A list of matrices where each element represents a gene family's expression profile distance matrix.
#'
#' @return A tibble with columns: `Family`, `Mean`, and `Median`, containing the mean and median
#' values of the distances for each family.
#' @examples
#' result <- calculate_exp.prof.dists.statistics(my_data)
calculate_exp.prof.dists.statistics <- function(data) {
  result <- map_dfr(names(data), function(name) {
    dist_matrix <- as.matrix(data[[name]])
    tibble(
      Family = name,
      Mean = mean(dist_matrix, na.rm = TRUE),
      Median = median(dist_matrix, na.rm = TRUE)
    )
  })
  result <- result %>%
    filter_all(all_vars(!is.na(.) & !is.infinite(.)))
  
  return(result)
}

#' Calculate Statistics for Expression Profile Distances by Tissue
#'
#' Computes the mean and median of the expression profile distances for each gene cluster
#' per tissue within a given dataset, transforming the output to have tissues in headers.
#'
#' @param data A nested list where each element represents a gene cluster, and each cluster
#' contains a list of tissues with associated expression profile distances.
#'
#' @return A tibble where each row represents a gene cluster, and columns include `Cluster` and 
#' `Tissue-specific` mean and median values.
#' @examples
#' result <- calculate_exp.prof.dists.tissue.statistics(my_data)
calculate_exp.prof.dists.tissue.statistics <- function(data) {
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
  
  df_transformed <- result %>%
    pivot_wider(
      id_cols = Cluster,
      names_from = Tissue,
      values_from = c(Mean, Median),
      names_sep = "_"
    ) %>%
    select(Cluster, sort(colnames(.)))
  
  return(df_transformed)
}
