
#' Load Data Frame with Specific Header Type
#'
#' @description Loads a data frame from a file path with a predefined header format.
#' This function reads a tab-separated file, with each column as a character vector.
#' Typical usage involves files containing gene family data with columns for family,
#' gene, species, and ortholog or paralog information.
#'
#' @param file_path A string specifying the path to the file to be read.
#' @return A data frame with five columns, each read as character strings.
#' @examples
#' # Load a data frame from a file with a specific format
#' df <- load_data_frame("path/to/file.txt")
load_data_frame <- function(file_path) {
    read.table(file_path, header = TRUE, sep = "\t", 
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))
}

#' Create Nested List from Data Frame
#'
#' @description Creates a nested list from a data frame of gene family information.
#' The output structure depends on the `header_type` argument, which specifies
#' whether the data contains ortholog or paralog information.
#'
#' @param df A data frame with gene family data, including columns for family, gene, 
#'   species, and ortholog/paralog information.
#' @param header_type A string specifying the type of header format to use, either 
#'   `"Ortholog"` or `"Paralog"`.
#' @return A nested list organized by family, gene, and species, containing either 
#'   ortholog or paralog data based on the specified `header_type`.
#' @examples
#' # Create a nested list from a data frame with ortholog information
#' nested_list <- create_nested_list(df, header_type = "Ortholog")
create_nested_list <- function(df, header_type) {
    if (header_type == "Ortholog") {
        gene_col <- "Ortholog"
        species_col <- "Ortholog_species"
    } else if (header_type == "Paralog") {
        gene_col <- "Paralog"
        species_col <- "Paralog_species"
    } else {
        stop("Invalid header type provided.")
    }

    df %>%
        group_by(Family, Gene_species, Gene, !!sym(species_col)) %>%
        summarise(!!sym(gene_col) := list(!!sym(gene_col)), .groups = "drop") %>%
        group_by(Family, Gene_species, Gene) %>%
        summarise(nested_info = list(setNames(!!sym(gene_col), !!sym(species_col))), .groups = "drop") %>%
        group_by(Family) %>%
        summarise(gene_info = list(setNames(nested_info, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>%
        deframe()
}



#' Filter Genes Based on Expression Data
#'
#' @description Filters the input data frame by retaining only those rows where the gene in the `"Gene"` column has corresponding entries in the expression data. This function is intended to simplify the dataset by removing entries for genes with no available expression data.
#'
#' @param df A data frame containing gene data, with a column labeled `"Gene"`.
#' @param expression_data A data frame or tibble containing expression data, with a column `"FBpp_ID"` listing gene IDs with available expression data.
#' @return A filtered data frame containing only rows where the `"Gene"` column has matching entries in the `expression_data`.
#' @examples
#' # Filter genes based on available expression data
#' filtered_df <- filter_v1(df = gene_data, expression_data = rna_seq_data)
#'
filter_v1 <- function(df, expression_data) {
    # Filter genes that exist in the expression data
    genes_intersect <- intersect(df$Gene, expression_data$FBpp_ID)
    df <- df %>% filter(Gene %in% genes_intersect)
    
    return(df)
}


#' Filter Gene Pairs Based on Expression Data and Pair Type
#'
#' @description Filters the input data frame by retaining only those rows where both the gene in the `"Gene"` column and its corresponding ortholog or paralog (as specified in the `type` parameter) have entries in the expression data. This function is useful for reducing datasets to only include pairs where both genes have expression data available.
#'
#' @param df A data frame containing gene pairs, with one column labeled `"Gene"` and another column containing either `"Ortholog"` or `"Paralog"` gene identifiers, as specified by the `type` argument.
#' @param type A character string specifying the type of gene pair to filter by: `"Ortholog"` or `"Paralog"`. Default is `"Ortholog"`.
#' @param expression_data A data frame or tibble containing expression data, with a column `"FBpp_ID"` listing gene IDs for which expression data is available.
#' @return A filtered data frame containing only the gene pairs where both the `"Gene"` and the specified `type` column have corresponding entries in `expression_data`.
#' @examples
#' # Filter ortholog pairs based on available expression data
#' filtered_df <- filter_v2(df = gene_pairs, type = "Ortholog", expression_data = rna_seq_data)
#'
#' # Filter paralog pairs based on available expression data
#' filtered_df <- filter_v2(df = gene_pairs, type = "Paralog", expression_data = rna_seq_data)
#'
filter_v2 <- function(df, type = c("Ortholog", "Paralog"), expression_data) {
    type <- match.arg(type)
    
    # Filter genes that exist in expression data
    genes_intersect <- intersect(df$Gene, expression_data$FBpp_ID)
    df <- df %>% filter(Gene %in% genes_intersect)
    
    # Filter orthologs/paralogs that exist in expression data
    pairs_intersect <- intersect(df[[type]], expression_data$FBpp_ID)
    df <- df %>% filter(!!sym(type) %in% pairs_intersect)
    return(df)
}


