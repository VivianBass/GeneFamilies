
# function to compute euclidean distances pairs of gene expression values
# This function is used in compute_exp.prof.dists.R
exp.prof.dists <- function(gene.accessions, 
                  expression.profiles = rna.seq.exp.profils,
                  expr.prof.gene.col = "id", 
                  tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)),
                  dist.method = "euclidean", per.tissue = FALSE) {
    
    # Convert to tibble if not already
    expression.profiles <- as_tibble(expression.profiles)
    
    # Get all genes and filter expression profiles
    all_genes <- unlist(gene.accessions)
    exp.profs <- expression.profiles %>% filter(!!sym(expr.prof.gene.col) %in% all_genes) %>%
                 column_to_rownames(expr.prof.gene.col)
    
    # Calculate distances per tissue
    if (nrow(exp.profs) > 1) { if (per.tissue) {
                tissue_distances <- tissues %>% set_names() %>%
                map(~{ exp.profs %>% select(all_of(.x)) %>% dist(method = dist.method)})
                
                return(tissue_distances)

        } else {
                return(exp.profs %>% select(all_of(tissues)) %>% dist(method = dist.method))
        }
    } else {
        return(NA)
    }
}
