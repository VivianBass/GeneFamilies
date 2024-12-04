
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

#' Create Nested List from Gene Family Data
#'
#' @param df Data frame with gene family information (Family, Gene, Gene_species, and Ortholog/Paralog data)
#' @param header_type String specifying data type ("Ortholog" or "Paralog")
#'
#' @return Nested list organized by family, containing gene and species-level information
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Family = c("Fam1", "Fam1"),
#'   Gene = c("GeneA", "GeneB"),
#'   Gene_species = c("Sp1", "Sp2"),
#'   Ortholog = c("Orth1", "Orth2"),
#'   Ortholog_species = c("Sp1", "Sp2")
#' )
#' nested_list <- create_nested_list(df, "Ortholog")
#' }
#'
#' @importFrom dplyr group_by summarise across mutate
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

