

**Note:**  
The input arguments and file paths provided in this document are set to the `experiments/test` 
directory for testing purposes. Additionally, the output directory specified in the `.env` file 
is configured to save generated files in the designated output path within the test environment.

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`	              

1. example: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

    input.args[[1]] <- "experiments/test_diet_M_new_1/tpm_values_largest.tsv"

2. example: Rscript exec/load_gene_groups_data.R 
    <conserved_orthologs.tsv> <in_paralogs.tsv> <out_paralogs.tsv> <special_in_paralogs.tsv>  <special_out_paralogs.tsv>

    input.args[[1]] <- "experiments/test_diet_M_new_1/dmel_dsec_conserved_orthologs.tsv"
    input.args[[2]] <- "experiments/test_diet_M_new_1/dmel_dsec_in_paralogs.tsv"
    input.args[[3]] <- "experiments/test_diet_M_new_1/dmel_dsec_out_paralogs.tsv"
    input.args[[4]] <- "experiments/test_diet_M_new_1/dmel_dsec_special_in_paralogs.tsv"
    input.args[[5]] <- "experiments/test_diet_M_new_1/dmel_dsec_special_out_paralogs.tsv"

3. example: Rscript exec/load_gene_families_data.R <families_file> 

    input.args[[1]] <- "experiments/test_diet_M_new_1/orthogroups_dmel_dsec.csv"


## Section-2 - Computing Distances & Statistics

- for this section you will not need to provide any input arguments. just run the rscripts.
			                                      
### 1. `compute_exp.prof.dists.R`

Rscript exec/compute_exp.prof.dists.R

### 2. `compute_exp.prof.dists_statistics.R`

USAGE: Rscript exec/compute_exp.prof.dists_statistics.R"

## Section-3 - Plotting Distributions

- for this section you will not need to provide any input arguments. just run the rscripts.

### 1. `plot_exp.prof.dists_distribution.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution.R

### 2. `plot_exp.prof.dists_distribution_tissue.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R

