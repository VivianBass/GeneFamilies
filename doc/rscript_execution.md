

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`	              

### 1. `load_gene_expression_data.R`
  
USAGE: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

Format Description of the Input files:

- <RPKM_counts_table.tsv> :     Header: id | tissue	| expression

### 2. `load_gene_groups_data.R`  

USAGE: Rscript exec/load_gene_groups_data.R <conserved_orthologs.tsv> 
<in_paralogs.tsv> <out_paralogs.tsv> <special_in_paralogs.tsv> <special_out_paralogs.tsv>

Format Description of the Input files:

- <conserved_orthologs.tsv> :   Header: Family | Gene | Gene_species | Ortholog | Ortholog_species
- <in_paralogs.tsv> :           Header: Family | Gene | Gene_species | Paralog  | Paralog_species
- <out_paralogs.tsv> :          Header: Family | Gene | Gene_species | Paralog  | Paralog_species
- <special_in_paralogs.tsv> :   Header: Family | Gene | Gene_species | Paralog  | Paralog_species 
- <special_out_paralogs.tsv> :  Header: Family | Gene | Gene_species | Paralog  | Paralog_species

### 3. `load_gene_families_data.R`

USAGE: Rscript exec/load_gene_families_data.R <families_file> 

Format Description of the Input files:

- <families_file> :             Header: Family | ...



## Section-2 - Computing Distances, Statistics, T-tests 
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`

### 1. `compute_exp.prof.dists.R`
### 2. `compute_exp.prof.dists_statistics.R`

## Section-3 - Plotting Distributions

- 1. `plot_exp.prof.dists_distribution.R`
- 2. `plot_exp.prof.dists_distribution_tissue.R`

### 1. `plot_exp.prof.dists_distribution.R`
### 2. `plot_exp.prof.dists_distribution_tissue.R`