

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`	              

- 1. example: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

    Rscript exec/load_gene_expression_data.R experiments/test/RPKM.tsv

    input.args[[1]] <- "experiments/test/RPKM.tsv"

- 2. example: Rscript exec/load_gene_groups_data.R 
    <conserved_orthologs.tsv> <in_paralogs.tsv> <special_in_paralogs.tsv> <out_paralogs.tsv> <special_out_paralogs.tsv>

    Rscript exec/load_gene_groups_data.R ...

    input.args[[1]] <- "experiments/test/conserved_orthologs_test.txt"
    input.args[[2]] <- "experiments/test/in_paralogs_test.tsv"
    input.args[[3]] <- "experiments/test/special_in_paralogs_test.tsv"
    input.args[[4]] <- "experiments/test/out_paralogs_test.tsv"
    input.args[[5]] <- "experiments/test/special_out_paralogs_test.tsv"

- 3. example: Rscript exec/load_gene_families_data.R <families_file> 

    Rscript exec/load_gene_families_data.R experiments/test/Orthogroups.tsv

    input.args[[1]] <- "experiments/test/Orthogroups.tsv"


## Section-2 - Computing Distances, Statistics, T-tests 
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`

## Section-3 - Plotting Distributions




