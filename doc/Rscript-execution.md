

## Section-1 - Loading Data

- 1.  `load_expression_data.R`	          
- 2.  `load_paralogs_orthologs_tandems_data.R`    
- 3.  `load_genefamilies_data.R`	              

- 1. example: Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>

    Rscript exec/load_expression_data.R experiments/RPKM_flybase/RPKM.tsv

    input.args[[1]] <- "experiments/RPKM_flybase/RPKM.tsv"

- 2. example: Rscript exec/load_paralogs_orthologs_tandems_data.R <all_vs_all_file> <orthologs_file> <paralogs_file>

    Rscript exec/load_paralogs_orthologs_tandems_data.R experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv

    input.args[[1]] <- "all vs all"
    input.args[[2]] <- "experiments/RPKM_flybase/orthologs.csv"
    input.args[[3]] <- "experiments/RPKM_flybase/paralogs.csv"

- 3. example: Rscript exec/load_genefamilies_data.R <families_file> <counts_file>

    Rscript exec/load_genefamilies_data.R experiments/RPKM_flybase/families.tsv experiments/RPKM_flybase/counts.txt

    input.args[[1]] <- "experiments/RPKM_flybase/families.tsv"
    input.args[[2]] <- "experiments/RPKM_flybase/counts.txt"


## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`
- 4.  `compute_exp.prof.dists_angles.R`  


- 1.  Rscript exec/compute_exp.prof.dists.R 

- 2.  Rscript exec/compute_exp.prof.dists_statistics.R

    Rscript exec/compute_exp.prof.dists_statistics.R <orthologs.exp.prof.dists> <paralogs.exp.prof.dists> <orthologs.exp.prof.dists.tissue> <paralogs.exp.prof.dists.tissue>

- 3.  Rscript exec/compute_t-tests_wilcox-tests.R

- 4.  Rscript exec/compute_exp.prof.dists_angles.R 



## Section-3 - Plotting Distributions

- 1. `plot_expression_angles.R`    

- 1.  Rscript exec/plot_expression_angles.R



