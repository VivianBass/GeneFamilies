



#' Load coding sequences in DNA alphabetically:
dana.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dana-all-CDS-r1.06.fasta")
dere.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dere-all-CDS-r1.05.fasta")
dgri.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dgri-all-CDS-r1.05.fasta")
dmel.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dmel-all-CDS-r6.59.fasta")
dmoj.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dper-all-CDS-r1.3.fasta")
dper.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dper-all-CDS-r1.3.fasta")
dpse.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dpse-all-CDS-r3.04.fasta")
dsec.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dsec-all-CDS-r1.3.fasta")
dsim.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dsim-all-CDS-r2.02.fasta")
dvir.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dvir-all-CDS-r1.2.fasta")
dwil.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dwil-all-CDS-r1.3.fasta")
dyak.cds <- loadDnaFasta("data-drosophila/Genome-Fasta/dyak-all-CDS-r1.3.fasta")


#' Compile a single list holding ALL CDS of all eight Brassicaceaen genomes:
all.cds <- Reduce(append, list(dana.cds, dere.cds, dgri.cds, dmel.cds, dmoj.cds, dper.cds, dpse.cds, dsec.cds, dsim.cds, dvir.cds, dwil.cds, dyak.cds ))

#' Save results:
save(dana.cds, dere.cds, dgri.cds, dmel.cds, dmoj.cds, dper.cds, dpse.cds, dsec.cds, dsim.cds, dvir.cds, dwil.cds, dyak.cds, all.cds, file = file.path("data-drosophila/Genome-Fasta", "codingSequences-drosophila.RData"))
message("DONE")
