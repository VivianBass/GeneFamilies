


-----------------------------------------------------------------------------------------
# Overview of Project Files and Rscripts
-----------------------------------------------------------------------------------------

- Section-1: `1. Load tissue specific expression data`  
- Section-2: `2. Compute pairwise gene expression profile distances within gene groups`
- Section-3: `3. Expression Profile based Function Diversity`
- Section-4: `4. Expression vector space analysis`

-----------------------------------------------------------------------------------------
## Rscripts: 
-----------------------------------------------------------------------------------------

- 1.  `loadRpkmRnaSeqCounts.R`			                                      (Section-1)
- 2.  `loadcodingSequences.R`			                                        (Section-1)
- 1.  `LoadOrthologsAndTandems.R`			                                    (Section-2)
- 3.  `loadGeneFamilies.R`				                                        (Section-2)

- 4.  `computeExpressionProfileDistances.R` 	                            (Section-2)

- 1.  `investigateDistributionsOfExpressionProfileDistances.R`		        (Section-3)                          
- 2.  `plotDistributionsOfExpressionProfileDistances.R`				            (Section-3)


- 1.  `generateRpkmExpressionProfiles.R` 		                              (Section-4) 
- 2.  `expressionAngles.R`	                                              (Section-4) 



-----------------------------------------------------------------------------------------
## Files:
-----------------------------------------------------------------------------------------

  Section-1:

- `RPKM_counts_table.tsv`             
- `RNA_Seq_RPKM_and_profiles.tsv`     
- `aet.fa aly.fa ath.fa bra.fa chi.fa cru.fa esa.fa tpa.fa` 
- `all.cds`

--> `RNA_Seq_RPKM_and_profiles.RData`
--> `codingSequences.RData`

-----------------------------------------------------------------------------------------

  Section-2:

- `eight_brassicaceae_tandems.txt`, `eight_brassicaceae_orthologs.txt`                      X
- `orthologs` & `tandems` 										                
- `orthologs.lst` & `tandems.lst` 								             

- `mcl_output.txt` , `mcl_table.tsv`    
- `families.lst` & `families.genes.df, families.df` 					

- `families.exp.prof.dists` & `families.exp.prof.dists.tissue`			  
- `orthologs.exp.prof.dists` & `orthologs.exp.prof.dists.tissue`			
- `tandems.exp.prof.dists` & `tandems.exp.prof.dists.tissue`

--> `orthologsTandems.RData` 
--> `GeneGroups.RData` 				 
--> `families.RData`       
--> `ExpressionProfileDistances.RData`

-----------------------------------------------------------------------------------------

  Section-3:

- `families.exp.prof.dists` 			    
- `families.exp.prof.dists.tissue`		
- `orthologs.exp.prof.dists` 			    
- `orthologs.exp.prof.dists.tissue`		
- `tandems.exp.prof.dists` 				    
- `tandems.exp.prof.dists.tissue`		 

- `tandems.exp.prof.dists.orth.dist`   
- `tandems.exp.prof.dists.orth.dist.df` 	    
- `families.exp.prof.dists.orth.dist`
- `families.exp.prof.dists.orth.dist.df` 	 
- `orthologs.exp.prof.dists.stats`
- `orthologs.exp.prof.dists.stats.df`		

--> `ExpressionProfileDistances.RData` 	
--> `ExpressionProfileDistanceDistributions.RData`
--> `families.exp.RData`    
--> `BUSTED_Results.RData`  
--> `GeneGroups.RData`

-----------------------------------------------------------------------------------------

  Section-4:

- `rpkm.rna.seq.counts` (as `RPKM_counts_table.tsv`) 
- `all.cds` 	         
- `rpkm.expr.profiles.df` 					             

- `tands.expr, orths.expr, dupl.expr, psel.expr` 	
- `tands.w.orths` 							                 
- `dupl.w.orths`							                    

- `tands.w.orths.angles.df`					             
- `tands.psel.w.orths.angles.df`				         
- `dupl.w.orths.angles.df` 					              
- `dupl.psel.w.orths.angles.df`	

--> `rpkmExpressionProfiles.RData` 
--> `GeneGroups.RData`


-----------------------------------------------------------------------------------------
## Plots:
-----------------------------------------------------------------------------------------

      Section-3:
- 1. `expressionAngleToDiagonalBoxplot.pdf`					        
- 2. `relativeExpressionVersatilityBoxplot.pdf`	
- 3. `medianExpressionProfileDistancesPerGeneGroupClassesBoxplot.pdf`
- 4. `medianExpressionProfileDistancesPerGeneGroupClassesBoxplot_NoPosSel.pdf`
- 5. `medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot.pdf`
- 6. `medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot_NoPosSel.pdf`
- 7. `medianExpressionDistsTandemNonOrthsAndExpandedNonOrthsHist.pdf`

      Section-4:
- 1. `meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf`	 
- 2. `meanTissueVersatilityDiffsAfterDuplicationHistograms.pdf`		
- 3. `afterDuplicationAngleBetweenOrth2DiagVecsBoxplot.pdf`			  
- 4. `afterDuplicationAngleBetweenOrth2DiagVecsHistograms.pdf`		
- 5. `meanChangeInAngleBetweenOrthOnDiagsBoxplot.pdf`


## GeneSets
-----------------------------------------------------------------------------------------

  - Fubar-, Busted-, Cafe Results are used to create the following GeneSets:

  - exec/`defineGeneSets.R`

  `GeneSet 1` :   sets of genes and gene families showing signs of positive selection:
  `GeneSet 2` : 	genes and gene families showing species specific sign of
                  significant expansion or contraction:
  `GeneSet 3` : 	Conserved families are those that have an identical number of genes 
                  within each species
  `GeneSet 4` : 	Generate data containers required for the computation of annotation based
                  function diversity (Shannon-Entropies).
  `GeneSet 5` : 	Define Duplicated Gene Sets with Orthologs:
  `GeneSet 6` : 	Define lists that can be used to classify above groups into subgroups:
  `GeneSet 7` :   Define sets of expressed genes:

## Plots created by exec/`defineGeneSets.R`
-----------------------------------------------------------------------------------------
                      
  - 1. `GeneGroupsVenn.pdf`
  - 2. `ExpressedGeneGroupsVenn.pdf`

- `Venn Diagram` 

- This visualization illustrates gene families present across various species, 
  highlighting both shared and unique families. The central area shows gene families 
  common to all species, indicating conserved genes, while the outer segments represent 
  species-specific families. In the Brassicaceae family, 12,654 gene families are shared 
  among eight species, demonstrating conservation within the group. Unique gene families 
  are also noted, with C. hirsuta having 694, A. thaliana 1,020, and C. rubella 541. 
  This analysis helps explore evolutionary relationships and functional differences among these species.