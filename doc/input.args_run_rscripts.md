

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`	              

### 1. `load_gene_expression_data.R`
  
USAGE: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

Format Description of the Input files:

- <RPKM_counts_table.tsv> :     Header: FBpp_ID | tissue | expression

```R
	id	    tissue		expression
	gene1	tissue1		####
	gene2	tissue1		####
	gene1	tissue2		####
	gene2	tissue2		####
```

### 2. `load_gene_groups_data.R`  

USAGE: Rscript exec/load_gene_groups_data.R <conserved_orthologs.tsv> 
<in_paralogs.tsv> <out_paralogs.tsv> <special_in_paralogs.tsv> <special_out_paralogs.tsv>

Format Description of the Input files:

- <conserved_orthologs.tsv> :   Header: Family | Gene | Gene_species | Ortholog | Ortholog_species
- <in_paralogs.tsv> :           Header: Family | Gene | Gene_species | Paralog  | Paralog_species
- <out_paralogs.tsv> :          Header: Family | Gene | Gene_species | Paralog  | Paralog_species
- <special_in_paralogs.tsv> :   Header: Family | Gene | Gene_species | Paralog  | Paralog_species 
- <special_out_paralogs.tsv> :  Header: Family | Gene | Gene_species | Paralog  | Paralog_species

```R
	con_orthologs:
	Family		Gene		Gene_species	Ortholog	Ortholog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj

	in_paralogs:
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	special_in_paralogs:
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	out_paralog_
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	special_out_paralog:
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
```


### 3. `load_gene_families_data.R`

USAGE: Rscript exec/load_gene_families_data.R <families_file> 

Format Description of the Input files:

- <families_file> :             Header: Family | species1 | species2 ...

```R
	Family		species1		species2		species3
	family_1	gene1,gene2,gene3	gene4,gene5,gene6	gene7,gene8,gene9
```


## Section-2 - Computing Distances & Statistics

- no input.args required

### 1. `compute_exp.prof.dists.R`

Rscript exec/compute_exp.prof.dists.R

### 2. `compute_exp.prof.dists_statistics.R`

USAGE: Rscript exec/compute_exp.prof.dists_statistics.R"

## Section-3 - Plotting Distributions

- no input.args required

### 1. `plot_exp.prof.dists_distribution.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution.R

### 2. `plot_exp.prof.dists_distribution_tissue.R`

USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R