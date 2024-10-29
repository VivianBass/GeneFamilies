In the new EasyVectorOmics we will need from the user:

- Expression values  (the same as we are using now only delete rank and variance) Expected input:
	id	tissue		expression
	gene1	tissue1		####
	gene2	tissue1		####
	gene1	tissue2		####
	gene2	tissue2		####
	

- Families with respective genes. Expected input:
	Family		species1		species2		species3
	family_1	gene1,gene2,gene3	gene4,gene5,gene6	gene7,gene8,gene9
	

- Load conserved_orthologs, in_paralogs, out_paralogs, special_in_paralogs, special_out_paralogs:
	Expected command for this script:
	Rscript loadOrthologsParalogs.R <in_paralogs.tsv> <special_in_paralogs.tsv> <out_paralogs.tsv> <special_out_paralogs.tsv> 		<conserved_orthologs.tsv>

	Expected input for conserved_orthologs (the best pair of orthologs):
	Family		Gene		Gene_species	Ortholog	Ortholog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj

	Expected input for in_paralogs (paralogs in the same species not related to a conserved ortholog):
	
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	special_in_paralogs (paralogs in the same species related to a conserved ortholog)
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	out_paralog (paralogs in other species not related to a conserved ortholog)
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	special_out_paralog (paralogs in other species related to a conserved ortholog)
	Family		Gene		Gene_species	Paralog		Paralog_species
	OG0000000	FBpp0117097	dana		FBpp0172663	dmoj
	
	
	
	
