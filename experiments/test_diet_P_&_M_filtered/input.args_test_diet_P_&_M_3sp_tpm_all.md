

**Note:**  
The input arguments and file paths provided in this document are set to the `experiments/test` 
directory for testing purposes. Additionally, the output directory specified in the `.env` file 
is configured to save generated files in the designated output path within the test environment.

## Section-1 - Loading Data

<br>

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`

<br>

1. example: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>
```R
    input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/3sp_tpm_all.tsv"
```
2. example: Rscript exec/load_gene_groups_data.R 
    <conserved_orthologs.tsv> <in_paralogs.tsv> <out_paralogs.tsv> <special_in_paralogs.tsv>  <special_out_paralogs.tsv>
```R
    input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_con_orthologs.tsv"
    input.args[[2]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_in_paralogs.tsv"
    input.args[[3]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_out_paralogs.tsv"
    input.args[[4]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_special_in_paralogs.tsv"
    input.args[[5]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_special_out_paralogs.tsv"
```

3. example: Rscript exec/load_gene_families_data.R <families_file> 
```R
    input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/orthogroups_filtered_more_than_5.tsv"
```

<br>

## Section-2 - Computing Distances & Statistics

**no input.args required**

### 1. `compute_exp.prof.dists.R`

USAGE: Rscript exec/compute_exp.prof.dists.R

### 2. `compute_exp.prof.dists_statistics.R`

USAGE: Rscript exec/compute_exp.prof.dists_statistics.R"

### 3. `compute_exp.prof.dists_angles.R`

USAGE: Rscript exec/compute_exp.prof.dists_angles.R

<br>

## Section-3 - Plotting Distributions

**no input.args required**

### 1. `plot_exp.prof.dists_distribution.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution.R

### 2. `plot_exp.prof.dists_distribution_tissue.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R

### 3. `plot_exp.prof.dists_angles.R`

USAGE: Rscript exec/plot_exp.prof.dists_angles.R

### 4. `generate_t-test_wilcox_test.R`

USAGE: Rscript exec/generate_t-test_wilcox_test.R

### 5. `generate_t-test_wilcox_test_tissue.R`

USAGE: Rscript exec/generate_t-test_wilcox_test_tissue.R
