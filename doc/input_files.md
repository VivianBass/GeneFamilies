
# Overview of Project Files and Rscripts 

## Section-1 - Loading Data

- 1.  `load_gene_expression_data.R`	          
- 2.  `load_gene_groups_data.R`    
- 3.  `load_gene_families_data.R`                                        

- 1.  `load_gene_expression_data.R`	 :

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")
message("Note:   <RPKM_counts_table.tsv> is expected to be TAB-Delimited and have a header line:\n", "id | tissue | rank")

- input.args[[1]] <- <RPKM_counts_table.tsv>

- RNA_Seq_RPKM_and_profiles.RData -> rna.seq.exp.profils, rpkm.rna.seq.counts
- RNA_Seq_RPKM_and_profiles.tsv


- Expression values  (the same as we are using now only delete rank and variance) Expected input:
	id	tissue		expression
	gene1	tissue1		####
	gene2	tissue1		####
	gene1	tissue2		####
	gene2	tissue2		####


	
- 2.  `load_gene_groups_data.R` 






- 3.  `load_gene_families_data.R` 

- Families with respective genes. Expected input:
	Family		species1		species2		species3
	family_1	gene1,gene2,gene3	gene4,gene5,gene6	gene7,gene8,gene9





## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`

## Section-3 - Plotting Distributions


















