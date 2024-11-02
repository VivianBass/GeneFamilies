
# Overview Input files required in the Loading Section

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`                                        

### 1.  `load_gene_expression_data.R`	 

**expression data:**

- <RPKM_counts_table.tsv> : 

- is expected is expected to be TAB-Delimited and have the Header: id | tissue | expression

```R
	id	tissue		expression
	gene1	tissue1		####
	gene2	tissue1		####
	gene1	tissue2		####
	gene2	tissue2		####
```

	
### 2. `load_gene_groups_data.R` 

**gene-groups data:**

- <conserved_orthologs.tsv> :
- <in_paralogs.tsv> :
- <special_in_paralogs.tsv> :
- <out_paralogs.tsv> :
- <special_out_paralogs.tsv> :	

```R
	conserved_orthologs:
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

**families data:** 

- <families_file> :

```R
	Family		species1		species2		species3
	family_1	gene1,gene2,gene3	gene4,gene5,gene6	gene7,gene8,gene9
```






