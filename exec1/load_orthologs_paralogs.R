

# load orthologs and paralogs

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[2]] <- 
input.args[[3]] <- 

orthologs <- read.table(input.args[[2]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 12))
orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) unlist(orthologs[x, 
    ])), orths.nms)
orthologs.genes <- unlist(orthologs.lst)


paralogs <- read.table(input.args[[3]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 2))
paralogs.nms <- unique(paralogs$Family)
paralogs.lst <- setNames(mclapply(paralogs.nms, function(x) paralogs[which(paralogs$Family == 
    x), "Gene"]), paralogs.nms)
paralogs.genes <- unlist(paralogs.lst)

save(paralogs, paralogs.lst, orhologs, orthologs.lst, file = file.path("ortho_para.RData"))

# --------------------------------------------------------------------------------------

# alternativ, get the orthologs and paralogs list from a Dataframe

orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) unlist(orthologs[x, ])), orths.nms)
orthologs.genes <- unlist(orthologs.lst)


paras.nms <- paste("paralog_cluster_", 1:nrow(paralogs), sep = "")
paralogs.lst <- setNames(mclapply(1:nrow(paralogs), function(x) unlist(paralogs[x, ])), paras.nms)
paralogs.genes <- unlist(paralogs.lst)


# ----------------------------------------------------------------------------------------

# load Genfamilies data

readMclOutputTest <- function(path.2.mcl.out, family.name.prefix = "cluster_") {
    sys.cmd <- paste("awk 'NR==1 {cluster=\"", family.name.prefix, "\" NR; for (i = 1; i <= NF; i++) { print cluster \"\\t\" $i }} NR>1 {cluster=\"", family.name.prefix, "\" NR; for (i = 2; i <= NF; i++) { print cluster \"\\t\" $i }}' ", path.2.mcl.out, sep = "")

    fread(sys.cmd, stringsAsFactors = FALSE, sep = "\t",
        colClasses = rep("character", 2),
        col.names = c("Family", "Gene"), na.strings = "",
        data.table = FALSE)
}

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[2]] <- "mcl_output.txt"
input.args[[3]] <- 

families.genes.df <- readMclOutputTest(input.args[[2]])
families.lst <- mclDataFrameAsList(families.genes.df)

families.df <- read.table(input.args[[3]], sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    comment.char = "", quote = "", na.strings = "", colClasses = c("character", rep("numeric",12)))

families.df$size <- apply(families.df[, 2:12], 1, sum)


save(families.genes.df, families.lst, families.df, file = file.path("families.RData"))





load("x/data/orthologs/orthologs_modified.tsv")

orthologs <- read.table("x/data/orthologs/orthologs_mod", header = TRUE, sep = "\t", comment.char = "", 
quote = "", na.strings = "", colClasses = rep("character", 12))

paralogs <- read.table("x/data/paralogs/paralogs_modify.txt", header = TRUE, sep = "\t", comment.char = "", quote = "", 
na.strings = "", colClasses = rep("character", 2))


save(list = c("orthologs", "paralogs"), file = "orthologs_paralogs.RData")