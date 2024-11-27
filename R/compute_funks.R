
#' Compute Euclidean Distances Between Gene Expression Profiles
#'
#' @description Calculates pairwise Euclidean distances across gene expression profiles for a given set of genes, aggregated across all specified tissues. The distances are computed within each species, providing a measure of similarity in gene expression profiles.
#'
#' @param gene.accessions A list of gene accession identifiers for which distances should be computed.
#' @param expression.profiles A data frame or tibble containing gene expression data. This should include a column for gene IDs and additional columns for expression values. Default is `rna.seq.exp.profils`.
#' @param expr.prof.gene.col The column name in `expression.profiles` containing gene IDs. Default is `"FBpp_ID"`.
#' @param tissues A character vector of column names in `expression.profiles` that represent tissue-specific expression values. By default, includes all columns except `expr.prof.gene.col`, `"Parent_FBgn"`, and `"Species"`.
#' @param dist.method Distance calculation method. Default is `"euclidean"`. Other methods (e.g., `"manhattan"`, `"maximum"`, etc.) can be specified.
#'
#' @return A list of vectors, each representing the pairwise Euclidean distances between gene expression profiles for a species. Each vector contains the distances between the specified genes' expression profiles across the selected tissues. If there are fewer than two expression profiles for a species, `NA` is returned for that species.
#'
#' @examples
#' # Compute overall distances for specified genes across all tissues
#' distances <- exp.prof.dists(gene.accessions = list("gene1", "gene2", "gene3"))
#' "Species"
#' @export
exp.prof.dists <- function(gene.accessions,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, "Species")),
                          dist.method = "euclidean") {
    
    all_genes <- unlist(gene.accessions)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            dist_matrix <- as.vector(dist(species_data, method = dist.method))
            
            return(dist_matrix)
        }
        return(NA)
    })
    
    return(distances)
}


exp.prof.dists_log2 <- function(gene.accessions,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, "Species")),
                          dist.method = "euclidean") {
    

    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            dist_matrix <- abs(log2(as.vector(dist(species_data, method = dist.method))))
            
            return(dist_matrix)
        }
        return(NA)
    })
    
    return(distances)
}



exp.prof.dists_cosine <- function(gene.accessions,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), 
                                          c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            
            # Calculate pairwise cosine similarities
            n <- nrow(species_data)
            cosine_dists <- numeric()
            
            for(i in 1:(n-1)) {
                for(j in (i+1):n) {
                    # Cosine similarity calculation
                    dot_product <- sum(species_data[i,] * species_data[j,])
                    norm_i <- sqrt(sum(species_data[i,]^2))
                    norm_j <- sqrt(sum(species_data[j,]^2))
                    cosine_sim <- dot_product / (norm_i * norm_j)
                    cosine_dists <- c(cosine_dists, cosine_sim)
                }
            }
            
            return(cosine_dists)
        }
        return(NA)
    })
    
    return(distances)
}





#' Compute Euclidean Distances Between Gene Expression Profiles by Tissue
#'
#' @description Computes pairwise Euclidean distances between expression profiles of specified genes, optionally by each tissue. This function is useful for analyzing gene similarity within individual tissues or conditions.
#'
#' @param gene.accessions A vector of gene accession identifiers for which distances should be computed.
#' @param expression.profiles A data frame or tibble containing gene expression data. This should include a column for gene IDs and additional columns for expression values. Default is `rna.seq.exp.profils`.
#' @param expr.prof.gene.col The column name in `expression.profiles` containing gene IDs. Default is `"FBpp_ID"`.
#' @param tissues A character vector of column names in `expression.profiles` that represent tissue-specific expression values. By default, includes all columns except `expr.prof.gene.col`, `"Parent_FBgn"`, and `"Species"`.
#' @param dist.method Distance calculation method. Default is `"euclidean"`. Other methods can be specified if desired.
#'
#' @return A list of distance vectors, one per tissue. Each vector contains the pairwise Euclidean distances between the specified genes' expression profiles for that tissue. If there are fewer than two expression profiles for a tissue, `NA` is returned.
#'
#' @examples
#' # Compute per-tissue distances for specified genes
#' tissue_distances <- exp.prof.dists_tissue(gene.accessions = c("gene1", "gene2", "gene3"))
#'
#' @export
exp.prof.dists_tissue <- function(gene.accessions,
                                 expression.profiles = rna.seq.exp.profils,
                                 expr.prof.gene.col = "FBpp_ID",
                                 tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, "Species")),
                                 dist.method = "euclidean") {
    
    all_genes <- unlist(gene.accessions)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    rownames(exp.profs) <- exp.profs[[expr.prof.gene.col]]
    exp.profs <- exp.profs[, tissues]
    
    if (nrow(exp.profs) > 1) {
        tissue_distances <- tissues %>%
            set_names() %>%
            map(~{
                tissue_data <- exp.profs %>%
                    select(all_of(.x)) %>%
                    sapply(as.numeric)
                
                as.vector(dist(tissue_data, method = dist.method))
            })
        
        return(tissue_distances)
    }
    return(NA)
}

#' Calculate Statistics for Expression Profile Distances
#'
#' @description 
#' Computes the mean and median of expression profile distances for each gene family
#' within a given dataset of distance matrices. The function handles multiple gene families
#' and returns summary statistics for each.
#'
#' @param data A list where each element is a matrix representing the expression profile
#'   distance matrix for a gene family. Each matrix should contain numeric values.
#'
#' @return A tibble with three columns:
#'   \item{Family}{The name of the gene family (from the names of the list elements).}
#'   \item{Mean}{The mean of the distances for the gene family.}
#'   \item{Median}{The median of the distances for the gene family.}
#'   The function returns a tibble with these statistics, excluding any entries with `NA` or infinite values.
#'
#' @examples
#' # Calculate the mean and median of expression profile distances for each gene family
#' result <- calculate_exp.prof.dists.statistics(my_data)
#'
#' @importFrom dplyr filter_all
#' @importFrom purrr map_dfr
#' @export
calculate_exp.prof.dists.statistics <- function(data) {
    if (length(data) == 0) {
        return(tibble(Family = character(), Mean = numeric(), Median = numeric()))
    }
    
    result <- map_dfr(names(data), function(name) {
        dist_matrix <- data[[name]]
        # Flatten nested species lists and combine all numeric values
        values <- unlist(dist_matrix)
        values <- values[is.finite(values)]
        
        tibble(
            Family = name,
            Mean = if(length(values) > 0) mean(values, na.rm = TRUE) else NA_real_,
            Median = if(length(values) > 0) median(values, na.rm = TRUE) else NA_real_
        )
    })
    
    result <- result[!is.na(result$Mean) & !is.na(result$Median), ]
    return(result)
}

#' Calculate Statistics for Expression Profile Distances by Tissue
#'
#' @description
#' Computes the mean and median of the expression profile distances for each gene cluster, 
#' per tissue within a given dataset. The output is structured such that tissues are represented in the columns, 
#' with mean and median statistics for each tissue.
#'
#' @param data A nested list where each element represents a gene cluster. Each cluster contains a list of 
#'   tissues with associated expression profile distances. The structure should be such that each cluster (gene family) 
#'   contains named tissues as its sublist, and each tissue contains numeric distances.
#'
#' @return A tibble where each row corresponds to a gene cluster, and the columns include:
#'   - `Family`: The gene cluster name.
#'   - `Tissue`: The name of the tissue.
#'   - `Mean`: The mean of the expression profile distances for that tissue.
#'   - `Median`: The median of the expression profile distances for that tissue.
#'   The result contains one row for each tissue in each gene cluster, with statistics for each tissue.
#'
#' @examples
#' # Calculate the mean and median of expression profile distances for each gene cluster per tissue
#' result <- calculate_exp.prof.dists.tissue.statistics(my_data)
#'
#' @export
calculate_exp.prof.dists.tissue.statistics <- function(data) {
    if (length(data) == 0) {
        return(tibble(
            Family = character(),
            Tissue = character(),
            Mean = numeric(),
            Median = numeric()
        ))
    }

    result <- map_dfr(names(data), function(name) {
        tissue_data <- data[[name]]
        if (!is.null(tissue_data)) {
            map_dfr(names(tissue_data), function(tissue) {
                values <- unlist(tissue_data[[tissue]])
                values <- values[is.finite(values)]
                
                tibble(
                    Family = name,
                    Tissue = tissue,
                    Mean = if(length(values) > 0) mean(values, na.rm = TRUE) else NA_real_,
                    Median = if(length(values) > 0) median(values, na.rm = TRUE) else NA_real_
                )
            })
        }
    })

    result <- result %>%
        filter(if_all(c(Mean, Median), ~!is.na(.) & !is.infinite(.)))

    return(result)
}

#' Validate and Filter Loaded Data Objects
#'
#' This function filters a list of loaded data objects by a specified pattern and validates them to ensure they meet the required criteria. 
#' It checks that each object is a non-empty list or vector and contains valid data (i.e., it is not entirely `NA` values).
#'
#' @param loaded_objects A character vector containing the names of objects loaded in the R environment.
#' @param pattern A character string representing the regex pattern to identify the names of the target data objects.
#' @return A character vector of validated object names that match the specified pattern and contain valid data.
#'         Invalid data types (e.g., non-lists and non-vectors), empty objects, or objects with only `NA` values are excluded.
#' @examples
#' # Load data and then filter and validate the names of loaded gene expression profiles:
#' loaded_objects <- ls() # List all loaded objects
#' valid_data_names <- validate_data(loaded_objects, "(_dists$|\\.dists$)")
#' valid_data_names_tissue <- validate_data(loaded_objects, "(_dists_tissue$|\\.dists\\.tissue$)")
#' 
#' @export
validate_data <- function(loaded_objects, pattern) {
    data_names <- loaded_objects[grepl(pattern, loaded_objects)]
    valid_names <- character()
    
    for (name in data_names) {
        data_object <- get(name, envir = parent.frame())
        
        if (!is.list(data_object) && !is.vector(data_object)) {
            message(sprintf("Excluding %s - invalid data type", name))
            next
        }
        
        if (length(data_object) == 0 || all(is.na(unlist(data_object)))) {
            message(sprintf("Excluding %s - empty or NA-only data", name))
            next
        }
        
        valid_names <- c(valid_names, name)
    }
    
    return(valid_names)
}

#' Perform Tissue-Specific Statistical Tests
#'
#' Conducts t-tests and Wilcoxon tests for comparing gene distances between types within tissues. 
#' Results are adjusted for multiple comparisons using the Benjamini-Hochberg method.
#'
#' @param data A dataframe containing distance data with columns `Tissue`, `Type`, and `Distance`.
#' @param valid_groups A named list of dataframes with valid `Tissue` and `Type` combinations for analysis.
#' @param analysis_type A string indicating the analysis type ("mean" or "median").
#' @return A list containing t-test and Wilcoxon test results as dataframes, or `NULL` if no valid groups are available.
#' @examples
#' perform_tissue_tests(data, valid_groups, "mean")
perform_tests <- function(data, valid_groups, analysis_type) {
    if (length(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            t_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            wilcox_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}

#' Perform Tissue-Specific Statistical Tests
#'
#' Conducts t-tests and Wilcoxon tests for comparing gene distances between types within tissues. 
#' Results are adjusted for multiple comparisons using the Benjamini-Hochberg method.
#'
#' @param data A dataframe containing distance data with columns `Tissue`, `Type`, and `Distance`.
#' @param valid_groups A named list of dataframes with valid `Tissue` and `Type` combinations for analysis.
#' @param analysis_type A string indicating the analysis type ("mean" or "median").
#' @return A list containing t-test and Wilcoxon test results as dataframes, or `NULL` if no valid groups are available.
#' @examples
#' perform_tissue_tests(data, valid_groups, "mean")
perform_tissue_tests <- function(data, valid_groups, analysis_type) {
    if (nrow(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            wilcox_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}
