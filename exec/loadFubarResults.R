require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/loadFubarResults.R path/2/families/working_directory path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[1]] <- "data/phylogeny/"
input.args[[2]] <- "./"

readFubarTable2 <- function(path.2.fubar.csv, family.name, sign.post.prob = 0.9) {
    # Leer el archivo CSV
    df <- read.csv(path.2.fubar.csv, header = FALSE)
    
    # Filtrar filas donde la columna 6 es mayor o igual al umbral
    filtered_df <- df[df$V6 >= sign.post.prob, ]
    print(filtered_df)
    
    # Seleccionar las columnas relevantes
    if (nrow(filtered_df) > 0) {
        result <- data.frame(
            Codon = filtered_df$V1,
            Post.Prob = filtered_df$V6,
            BayesFactor = filtered_df$V7
        )
        result$family <- family.name
        return(result)
    } else {
        return(NULL)
    }
}



#' Find all FUBAR result files and load them as tables:
fubar.fls <- system(paste("find", input.args[[1]], "-name 'cluster_*_ml_tree_no_node_labels.newick.fubar.csv' -type f"), 
    intern = TRUE)
names(fubar.fls) <- regmatches(fubar.fls, regexpr("cluster_\\d+", fubar.fls))

# fubar.fls <- head(fubar.fls, 2)  # Cambia el número a 2 para limitar a los primeros dos

fubar.tbl <- Reduce(rbind, mclapply(names(fubar.fls), function(fam) {
    readFubarTable2(fubar.fls[[fam]], fam)
}))


#' Identify families with decisive evidence for positive selection at at least
#' a single codon [https://en.wikipedia.org/wiki/Bayes_factor]:
fubar.fams.decisive.evidence <- sort(unique(fubar.tbl[which(fubar.tbl$BayesFactor > 
    100), "family"]))


#' Save results:
save(fubar.tbl, fubar.fams.decisive.evidence, file = file.path(input.args[[2]], "data", 
    "fubar_results.RData"))

message("DONE")
