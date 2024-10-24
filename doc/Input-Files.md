
# Overview of Project Files and Rscripts
-----------------------------------------------------------------------------------------


## Section-1 - Loading Data

- 1.  `load_expression_data.R`	          
- 2.  `load_paralogs_orthologs_tandems_data.R`    
- 3.  `load_genefamilies_data.R`	                                        



- 1.  `load_expression_data.R` :


    - `RPKM_counts_table.tsv`             
    - `RNA_Seq_RPKM_and_profiles.tsv`  
    --> `RNA_Seq_RPKM_and_profiles.RData`

    - `rpkm.rna.seq.counts` (as `RPKM_counts_table.tsv`)         
    - `rpkm.expr.profiles.df` 

    --> `rpkmExpressionProfiles.RData` 


    - `eight_brassicaceae_tandems.txt`, `eight_brassicaceae_orthologs.txt`                      X
    - `orthologs` & `tandems` 										                
    - `orthologs.lst` & `tandems.lst` 								             

    - `mcl_output.txt` , `mcl_table.tsv`    
    - `families.lst` & `families.genes.df, families.df` 	

    --> `orthologsTandems.RData` 	 
    --> `families.RData`     




## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`
- 4.  `compute_exp.prof.dists_angles.R` 





    - `families.exp.prof.dists` & `families.exp.prof.dists.tissue`			  
    - `orthologs.exp.prof.dists` & `orthologs.exp.prof.dists.tissue`			
    - `tandems.exp.prof.dists` & `tandems.exp.prof.dists.tissue`

    - `tandems.exp.prof.dists.orth.dist`   
    - `tandems.exp.prof.dists.orth.dist.df` 	    
    - `families.exp.prof.dists.orth.dist`
    - `families.exp.prof.dists.orth.dist.df` 	 
    - `orthologs.exp.prof.dists.stats`
    - `orthologs.exp.prof.dists.stats.df`		

    --> `ExpressionProfileDistances.RData` 	
    --> `ExpressionProfileDistanceDistributions.RData`
    --> `families.exp.RData` 





## Section-3 - Plotting Distributions

- 1. `plot_expression_angles.R`                                                 



    - 1. `expressionAngleToDiagonalBoxplot.pdf`					        
    - 2. `relativeExpressionVersatilityBoxplot.pdf`	

    - 3. `medianExpressionProfileDistancesPerGeneGroupClassesBoxplot.pdf`
    - 4. `medianExpressionProfileDistancesPerGeneGroupClassesBoxplot_NoPosSel.pdf`
    - 5. `medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot.pdf`
    - 6. `medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot_NoPosSel.pdf`
    - 7. `medianExpressionDistsTandemNonOrthsAndExpandedNonOrthsHist.pdf`



















