
#' Compute Euclidean Distances Between Pairs of Gene Expression Values
#'
#' @description Calculates pairwise Euclidean distances between gene expression profiles 
#'   for a given set of genes. Optionally, distances can be computed separately for each 
#'   tissue.
#'
#' @param gene.accessions A vector of gene accession identifiers for which distances are 
#'   to be computed.
#' @param expression.profiles A data frame or tibble containing gene expression profiles 
#'   with gene IDs in one column and expression values in other columns. Default is 
#'   `rna.seq.exp.profils`.
#' @param expr.prof.gene.col A string specifying the column name in `expression.profiles` 
#'   that contains the gene IDs. Default is `"id"`.
#' @param tissues A character vector of column names representing tissues or conditions 
#'   to be used for calculating distances. Default is all columns in `expression.profiles` 
#'   except `expr.prof.gene.col`.
#' @param dist.method A string specifying the distance calculation method. Default is 
#'   `"euclidean"`.
#' @param per.tissue A logical value. If `TRUE`, computes distances separately for each 
#'   tissue; if `FALSE`, computes distances across all tissues together. Default is `FALSE`.
#' @return A distance matrix if `per.tissue = FALSE`, or a list of distance objects 
#'   (one per tissue) if `per.tissue = TRUE`. Returns `NA` if there are fewer than 
#'   two expression profiles.
#' @examples
#' # Compute Euclidean distances for a set of genes across all tissues
#' distances <- exp.prof.dists(gene.accessions = c("gene1", "gene2", "gene3"))
#'
#' # Compute distances per tissue
#' tissue_distances <- exp.prof.dists(gene.accessions = c("gene1", "gene2", "gene3"), per.tissue = TRUE)
exp.prof.dists <- function(gene.accessions, 
                  expression.profiles = rna.seq.exp.profils,
                  expr.prof.gene.col = "FBpp_ID", 
                  tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)),
                  dist.method = "euclidean", per.tissue = FALSE) {
    
    # Convert to tibble if not already
    expression.profiles <- as_tibble(expression.profiles)
    
    # Get all genes and filter expression profiles
    all_genes <- unlist(gene.accessions)
    exp.profs <- expression.profiles %>% filter(!!sym(expr.prof.gene.col) %in% all_genes) %>% column_to_rownames(expr.prof.gene.col)
    
    # Check that there are sufficient profiles to compute distances
    if (nrow(exp.profs) > 1) { 
        if (per.tissue) {
            # Calculate distances separately for each tissue
            tissue_distances <- tissues %>% set_names() %>%
                map(~{
                    dist_matrix <- exp.profs %>% select(all_of(.x)) %>% dist(method = dist.method)
                    as.vector(dist_matrix)  # Convert dist object to vector of distances
                })
                
            return(tissue_distances)
        } else {
            # Calculate distances across all tissues combined
            dist_matrix <- exp.profs %>% select(all_of(tissues)) %>% dist(method = dist.method)
            return(as.vector(dist_matrix))  # Convert dist object to vector of distances
        }
    } else {
        # Not enough profiles to compute distances
        return(NA)
    }
}


