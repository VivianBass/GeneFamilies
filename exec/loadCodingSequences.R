require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadCodingSequences.R path/2/GeneFamilies/data aet.fa aly.fa ath.fa bra.fa chi.fa cru.fa esa.fa tpa.fa")

input.args <- commandArgs(trailingOnly = TRUE)
input.args[[1]] <- "data/"
#' Load coding sequences in DNA alphabetically:
dana.cds <- loadDnaFasta("inputs/fastas/dana.fasta")
dere.cds <- loadDnaFasta("inputs/fastas/dere.fasta")
dgri.cds <- loadDnaFasta("inputs/fastas/dgri.fasta")
dmel.cds <- loadDnaFasta("inputs/fastas/dmel.fasta")
dmoj.cds <- loadDnaFasta("inputs/fastas/dper.fasta")
dper.cds <- loadDnaFasta("inputs/fastas/dper.fasta")
dpse.cds <- loadDnaFasta("inputs/fastas/dpse.fasta")
dsec.cds <- loadDnaFasta("inputs/fastas/dsec.fasta")
dsim.cds <- loadDnaFasta("inputs/fastas/dsim.fasta")
dvir.cds <- loadDnaFasta("inputs/fastas/dvir.fasta")
dwil.cds <- loadDnaFasta("inputs/fastas/dwil.fasta")
dyak.cds <- loadDnaFasta("inputs/fastas/dyak.fasta")
 
 
#' Compile a single list holding ALL CDS of all eight Brassicaceaen genomes:
all.cds <- Reduce(append, list(dana.cds, dere.cds, dgri.cds, dmel.cds, dmoj.cds, dper.cds, dpse.cds, dsec.cds, dsim.cds, dvir.cds, dwil.cds, dyak.cds ))
# all.cds <- Reduce(append, list(dmel.cds))
 
#' Save results:
save(dana.cds, dere.cds, dgri.cds, dmel.cds, dmoj.cds, dper.cds, dpse.cds, dsec.cds, dsim.cds, dvir.cds, dwil.cds, dyak.cds, all.cds, file = file.path(input.args[[1]], "codingSequences.RData"))
# save(dmel.cds, all.cds, file = file.path(input.args[[1]], "codingSequences.RData"))

message("DONE")
