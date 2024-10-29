
# Overview of Project Files and Rscripts 

input.args[[1]] <- 
input.args[[2]] <- 
input.args[[3]] <- 


## Section-1 - Loading Data

- 1.  `load_expression_data.R`	          
- 2.  `load_paralogs_orthologs_tandems_data.R`    
- 3.  `load_genefamilies_data.R`	                                        

- 1.  `load_expression_data.R` :

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")
message("Note:   <RPKM_counts_table.tsv> is expected to be TAB-Delimited and have a header line:\n", "id | tissue | rank")

- input.args[[1]] <- <RPKM_counts_table.tsv>

- RNA_Seq_RPKM_and_profiles.RData -> rna.seq.exp.profils, rpkm.rna.seq.counts
- RNA_Seq_RPKM_and_profiles.tsv





## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`

 

## Section-3 - Plotting Distributions


















