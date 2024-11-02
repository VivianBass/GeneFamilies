
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
