
# `Section-1` - Loading

- 1. `loadRpkmRnaSeqCounts.R`	-> Load gene expressions.	   X                                  
- 2. `loadcodingSequences.R`       dont need ??
- 3. `LoadOrthologsAndTandems.R`   sim, all vs all fehlt noch  X
- 4. `loadGeneFamilies.R`           X

- basically already creates the expression profiles , do we need the script generateRpmkExpressionProfils ???
1.  Rscript exec/loadRpkmRnaSeqCounts.R <RPKM_counts_table.tsv>
-   Rscript exec/loadRpkmRnaSeqCounts.R experiments/RPKM_flybase/RPKM.tsv

- do we need all.cds at all ??
2.  loadcodingSequences.R`


- would have to add Tandems to, calculate  the sim ??, need all vs all blast/Diamond
3.  Rscript exec/loadOrthologsAndTandems.R <all_vs_all_file> <orthologs_file> <paralogs_file>
-   Rscript exec/loadOrthologsAndTandems.R inputs/all_vs_all_results.txt experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv
-   Rscript exec/loadOrthologsAndTandems.R experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv

4.  Rscript exec/loadGeneFamilies.R <families_file> <counts_file>,
-   Rscript exec/loadGeneFamilies.R path/2/GeneFamilies/data mcl_output.txt mcl_table.tsv"
-   Rscript exec/loadGeneFamilies.R experiments/RPKM_flybase/families.tsv experiments/RPKM_flybase/counts.txt







# `Section-2` - computing

- 1. `generateRpkmExpressionProfiles.R`		                                       
- 2. `computeExpressionProfileDistances.R` 	 
- 3. `investigateDistributionsOfExpressionProfileDistances.R`

1. Rscript exec/generateRpkmExpressionProfiles.R
- Rscript path/2/GeneFamilies/exec/generateRpkmExpressionProfiles.R path/2/GeneFamilies")



- load("data/data.RData")
2. Rscript exec/computeExpressionProfileDistances.R

- Rscript exec/investigateDistributionsOfExpressionProfileDistances.R <input_file>



# `Section-3` - Plotting

	                                  
- 2.  `plotDistributionsOfExpressionProfileDistances.R`	
- 2.  `expressionAngles.R`	 



# t tests
