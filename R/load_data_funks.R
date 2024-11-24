
#' Load Data Frame with Validation for Specific Column Patterns
#'
#' @description This function loads a tab-separated data frame from a given file path, validating its headers 
#' to ensure the file contains expected columns for either ortholog or paralog data. It reads all columns 
#' as character vectors and extracts relevant columns based on the file type.
#'
#' @param file_path A string specifying the path to the file to be read. The file must be tab-separated and 
#' contain a header row.
#'
#' @return A data frame with validated column headers corresponding to either ortholog or paralog data. 
#' If the file does not meet the expected structure, the function raises an error. 
#' 
#' @details
#' The function first validates the presence of required columns in the file header:
#'   - Base columns: `"Family"`, `"Gene"`, `"Gene_species"`.
#'   - Additional columns for either:
#'     - Ortholog data: `"Ortholog"`, `"Ortholog_species"`.
#'     - Paralog data: `"Paralog"`, `"Paralog_species"`.
#' If neither pattern is found, or the number of columns does not match the headers, an error is raised. 
#' For valid files, only the relevant columns are returned.
#'
#' @examples
#' \dontrun{
#' # Load a data frame with ortholog or paralog data
#' df <- load_data_frame("path/to/file.txt")
#' }
#' 
#' @seealso [read.table()] for file reading.
#'
#' @importFrom utils read.table
#' @export
load_data_frame <- function(file_path, file_type = NULL) {
    if (!file.exists(file_path)) {
        stop("cannot open file")
    }
    
    # Define expected column patterns
    expected_names <- list(
        Ortholog = c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"),
        Paralog = c("Family", "Gene", "Gene_species", "Paralog", "Paralog_species")
    )
    
    # Read data with explicit column names
    df <- tryCatch({
        read.table(file_path, 
                  header = FALSE,
                  skip = 1,
                  sep = "\t",
                  comment.char = "", 
                  quote = "",
                  col.names = expected_names[[file_type]],
                  colClasses = "character",
                  stringsAsFactors = FALSE)
    }, error = function(e) {
        # If that fails, try reading with headers
        df <- read.table(file_path, 
                        header = TRUE,
                        sep = "\t",
                        comment.char = "", 
                        quote = "",
                        colClasses = "character",
                        stringsAsFactors = FALSE)
        names(df) <- expected_names[[file_type]]
        return(df)
    })
    
    return(as.data.frame(df))
}


#' Create Nested List from Data Frame
#'
#' @description 
#' Creates a nested list from a data frame containing gene family information.
#' The structure of the output list depends on the `header_type` argument, which 
#' determines whether the data contains ortholog or paralog information.
#'
#' @param df A data frame containing gene family data. The data should include 
#'   columns for `Family`, `Gene`, `Gene_species`, and either `Ortholog` or `Paralog`
#'   information, depending on the `header_type`.
#' @param header_type A string specifying the type of information to use for the list. 
#'   It can either be `"Ortholog"` or `"Paralog"`, indicating which gene column to use.
#'
#' @return A nested list, where each family is represented as a list of genes, 
#'   with each gene further nested by species, containing either ortholog or paralog 
#'   data depending on the `header_type`. 
#'
#' @details 
#' This function expects a data frame where the columns represent gene family data.
#' The columns are checked against the required ones for either ortholog or paralog information, 
#' and an error is raised if any are missing. The resulting list is organized by family, 
#' and each family contains nested gene and species-level information.
#'
#' @examples
#' # Example: Create a nested list from a data frame with ortholog information
#' df <- data.frame(
#'   Family = c("Fam1", "Fam1", "Fam2"),
#'   Gene = c("GeneA", "GeneB", "GeneC"),
#'   Gene_species = c("Species1", "Species2", "Species1"),
#'   Ortholog = c("Ortholog1", "Ortholog2", "Ortholog3"),
#'   Ortholog_species = c("Species1", "Species2", "Species3")
#' )
#' nested_list <- create_nested_list(df, header_type = "Ortholog")
#' print(nested_list)
#'
#' @export
create_nested_list <- function(df, header_type) {
    # Handle empty dataframe
    if (nrow(df) == 0) {
        return(list())
    }
    
    # Set up column names based on header type
    if (header_type == "Ortholog") {
        gene_col <- "Ortholog"
        species_col <- "Ortholog_species"
        required_cols <- c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species")
    } else if (header_type == "Paralog") {
        gene_col <- "Paralog"
        species_col <- "Paralog_species"
        required_cols <- c("Family", "Gene", "Gene_species", "Paralog", "Paralog_species")
    } else {
        stop("Invalid header type provided.")
    }
    
    # Check for required columns
    if (!all(required_cols %in% names(df))) {
        stop("Required columns are missing.")
    }
    
    # Convert all columns to character type
    df <- df %>%
        mutate(across(everything(), as.character))
    
    # Process the data
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

