require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/generateFubarBatchFiles.R path/2/family_phylogenies_directory")

input.args <- commandArgs(trailingOnly = TRUE)
hyphy.fubar.bf2 <- "inputRedirect = {};\n\ninputRedirect[\"alignment\"]=\"<%= fam.cds.msa.path %>\";\ninputRedirect[\"tree\"]=\"<%= fam.tree.no.node.lables.path %>\";\nExecuteAFile (\"<%= hyphy.batch.files.dir %>/FUBAR.bf\", inputRedirect);"

input.args[[1]] <- "data/phylogeny/"
hyphy.batch.files.dir <- "/home/vivianbass/Documents/hyphy/res/TemplateBatchFiles/SelectionAnalyses"
fam.dirs <- system(paste("find ", normalizePath(input.args[[1]]), " -maxdepth 1 -type d -regex '.*cluster_[0-9]*'", 
    sep = ""), intern = TRUE)
fam.nms <- sub("^.*/", "", as.character(fam.dirs), perl = TRUE)
names(fam.dirs) <- fam.nms

for (fam in names(fam.dirs)) {
    fam.d <- fam.dirs[[fam]]
    fam.tree.no.node.lables.path <- file.path(fam.d, paste(fam, "_ml_tree_no_node_labels.newick", 
        sep = ""))
    fam.cds.msa.path <- file.path(fam.d, paste(fam, "_CDS_MSA.fasta", sep = ""))
    fam.hyphy.fubar.output.path <- file.path(fam.d, paste(fam, "_hyphy_fubar_output.txt", 
        sep = ""))
    fam.hyphy.fubar.batch.file.path <- file.path(fam.d, paste(fam, "_hyphy_fubar_input.bf", 
        sep = ""))
    brew(text = hyphy.fubar.bf2, output = fam.hyphy.fubar.batch.file.path)
}

message("DONE")
