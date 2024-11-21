
#' Calculate Angles for a Gene Group
#'
#' This function calculates diagnostic angles for a given set of genes based on RNA-seq expression profiles.
#'
#' @param genes A list of gene groups where each element is a vector of gene identifiers.
#' @param rna.seq.exp.profils A data frame containing RNA-seq expression profiles. 
#'   Must include a column `FBpp_ID` for gene identifiers and additional columns for tissue-specific expression levels.
#' @param tissues A vector of column names representing tissues in `rna.seq.exp.profils` for angle calculation.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{FBpp_ID}{Gene identifiers matching the input `genes`.}
#'   \item{angle.diag}{Calculated diagnostic angle values for each gene.}
#' }
#' If no matching genes are found, a warning is issued and an empty data frame is returned.
#'
#' @details
#' The function computes angles using the `cosDiag` function applied to tissue expression levels for each gene, normalized by \eqn{\sqrt{2}}.
#' Only genes present in both the `genes` list and the `FBpp_ID` column of `rna.seq.exp.profils` are processed.
#' NA or invalid values are filtered out from the resulting data frame.
#'
#' @note Ensure that the `cosDiag` function is defined in your environment and accepts the appropriate input format.
#'
#' @importFrom dplyr filter
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' # Example data
#' genes <- list(group1 = c("gene1", "gene2"), group2 = c("gene3", "gene4"))
#' rna.seq.exp.profils <- data.frame(
#'   FBpp_ID = c("gene1", "gene2", "gene3"),
#'   tissue1 = c(1.2, 0.8, 1.1),
#'   tissue2 = c(0.9, 1.0, 0.7)
#' )
#' tissues <- c("tissue1", "tissue2")
#'
#' # Compute angles
#' calculate_angles(genes, rna.seq.exp.profils, tissues)
#' }
#'
#' @export
calculate_angles <- function(genes, rna.seq.exp.profils, tissues) {
  genes.expr <- intersect(unlist(genes), rna.seq.exp.profils$FBpp_ID)
  
  if(length(genes.expr) == 0) {
    warning(paste("No matching genes found for group:", group))
    return(data.frame())
  }

  expr.angle.diag.df <- data.frame(
    FBpp_ID = genes.expr,
    angle.diag = as.numeric(mclapply(genes.expr, function(x) {
      cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == x), tissues])/sqrt(2)
    })),
    stringsAsFactors = FALSE
  )

  expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")
}


#' Validate Data Frames and Retain Original Names
#'
#' This function validates a list of data frames, keeping only those that are non-null and have at least one row. 
#' The names of the valid data frames in the input list are preserved in the returned list.
#'
#' @param df_list A named list of data frames to be validated.
#'
#' @return A named list of validated data frames. 
#'   Only data frames that are non-null and contain at least one row are included in the output list.
#'   If a data frame is null or empty, a warning message is displayed.
#'
#' @details
#' This function iterates over the input list of data frames, checks each data frame for validity, 
#' and retains the original name of each valid data frame in the output list.
#'
#' @examples
#' # Example usage
#' df_list <- list(
#'   df1 = data.frame(a = 1:3, b = 4:6),
#'   df2 = NULL,
#'   df3 = data.frame()
#' )
#'
#' validate_angle_dataframes(df_list)
#' # Output: List containing only 'df1', with a warning for 'df2' and 'df3'.
#'
#' @export
validate_angle_dataframes <- function(df_list) {
    valid_dfs <- list()
    for (name in names(df_list)) {
        df <- df_list[[name]]
        if (!is.null(df) && nrow(df) > 0) {
            valid_dfs[[name]] <- df
        } else {
            message(sprintf("Warning: %s dataframe is empty or null", name))
        }
    }
    return(valid_dfs)
}