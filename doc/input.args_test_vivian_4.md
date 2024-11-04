

**Note:**  
The input arguments and file paths provided in this document are set to the `experiments/test` 
directory for testing purposes. Additionally, the output directory specified in the `.env` file 
is configured to save generated files in the designated output path within the test environment.

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`	              

1. example: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

    - input.args[[1]] <- "experiments/test_vivian_4/RPKM_expression.tsv"

2. example: Rscript exec/load_gene_groups_data.R 
    <conserved_orthologs.tsv> <in_paralogs.tsv> <out_paralogs.tsv> <special_in_paralogs.tsv>  <special_out_paralogs.tsv>

    - input.args[[1]] <- "experiments/test_vivian_4/conserved_orthologs_90.tsv"
    - input.args[[2]] <- "experiments/test_vivian_4/in_paralogs_90.tsv"
    - input.args[[3]] <- "experiments/test_vivian_4/out_paralogs_90.tsv"
    - input.args[[4]] <- "experiments/test_vivian_4/special_in_paralogs_90.tsv"
    - input.args[[5]] <- "experiments/test_vivian_4/special_out_paralogs_90.tsv"

3. example: Rscript exec/load_gene_families_data.R <families_file> 

    - input.args[[1]] <- "experiments/test_vivian_4/Orthogroups.tsv"


## Section-2 - Computing Distances & Statistics
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`

### 1. `compute_exp.prof.dists.R`




### 2. `compute_exp.prof.dists_statistics.R`

## Section-3 - Plotting Distributions

- 1. `plot_exp.prof.dists_distribution.R`
- 2. `plot_exp.prof.dists_distribution_tissue.R`

### 1. `plot_exp.prof.dists_distribution.R`
### 2. `plot_exp.prof.dists_distribution_tissue.R`


