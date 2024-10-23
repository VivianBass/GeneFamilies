

## Section-1 - Loading Data

- 1.  `loadRpkmRnaSeqCounts.R`	          -> `load_expression_data`
- 1.  `LoadOrthologsParalogsTandems.R`    -> `load_paralogs_orthologs_tandems_data`
- 3.  `loadGeneFamilies.R`	              -> `load_genefamilies_data`


- basically already creates the expression profiles , do we need the script generateRpmkExpressionProfils ???
1.  Rscript exec/loadRpkmRnaSeqCounts.R <RPKM_counts_table.tsv>
-   Rscript exec/loadRpkmRnaSeqCounts.R experiments/RPKM_flybase/RPKM.tsv

- would have to add Tandems to, calculate  the sim ??, need all vs all blast/Diamond
3.  Rscript exec/loadOrthologsAndTandems.R <all_vs_all_file> <orthologs_file> <paralogs_file>
-   Rscript exec/loadOrthologsAndTandems.R inputs/all_vs_all_results.txt experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv
-   Rscript exec/loadOrthologsAndTandems.R experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv

4.  Rscript exec/loadGeneFamilies.R <families_file> <counts_file>,
-   Rscript exec/loadGeneFamilies.R path/2/GeneFamilies/data mcl_output.txt mcl_table.tsv"
-   Rscript exec/loadGeneFamilies.R experiments/RPKM_flybase/families.tsv experiments/RPKM_flybase/counts.txt



## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `computeExpressionProfileDistances.R` 	      
- 2.  `investigateDistributionsOfExpressionProfileDistances.R`
- 3.  `t-tests`
- 4.  `angles`  expressionAngles.R


- Rscript exec/investigateDistributionsOfExpressionProfileDistances.R <input_file>



## Section-3 - Plotting Distributions

                                       
- 1.  `plotDistributionsOfExpressionProfileDistances.R`	
- plot Distances
- plot Distances Tissue

- 2.  `plot_expressionAngles.R`	