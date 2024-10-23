
# EasyVectorOmics Package

In this profect the goal is to reproduce the vector space analyses as done in
the _Cardamine hirsuta_ genome project [1]. The R-Code required to do so shall
be isolated and made usabel with _any_ data. Next, this code shall be used to
analyse public open access data of model species to evaluate whether the method
works for other data-sets as well.

# References

[1] Gan, X., Hay, A., Kwantes, M., Haberer, G., Hallab, A., Ioio, R. D.,
  Hofhuis, H., Pieper, B., Cartolano, M., Neumann, U., Nikolov, L. A., Song,
  B., Hajheidari, M., Briskine, R., Kougioumoutzi, E., Vlad, D., Broholm, S.,
  Hein, J., Meksem, K., … Tsiantis, M. (2016). The Cardamine hirsuta genome
  offers insight into the evolution of morphological diversity. Nature Plants,
  2, 16167. https://doi.org/10.1038/nplants.2016.167


# Expression-Vector-Space-Analysis

- This process measures and compares `gene expression levels` under different conditions, 
  such as across tissues, developmental stages, or environmental contexts. 
  Differential expression analysis provides insights into the genetic mechanisms 
  underlying leaf shape and structure, identifying key genes involved in leaf development.

- main question:
- Is there a difference in expression diversity between Orthologs and Paralogs?



The steps and data required to carry out the original analysis are devided into 3 Sections 

- `1. Loading Data`                                                             (Section-1) 
- `2. Computing Distances, Statistics, T-tests and Angles`                      (Section-2) 
- `3. Plotting Distributions`                                                   (Section-3) 

-------------------------------------------------------------------------------------------

## Section-1 - Loading Data

- 1.  `loadRpkmRnaSeqCounts.R`	          -> `load_expression_data`
- 1.  `LoadOrthologsParalogsTandems.R`    -> `load_paralogs_orthologs_tandems_data`
- 3.  `loadGeneFamilies.R`	              -> `load_genefamilies_data`

- Tools used, to get paralogs, orthologs, tandems, genefamilies: 

- `BLAST` (Basic Local Alignment Search Tool) or `Diamond`: Identifies regions of similarity between 
  biological sequences (DNA, RNA, proteins) using local alignment and scores alignments based on sequence similarity.
  Helps find homologous sequences (orthologs, paralogs) by comparing query sequences against a database.

- `Markov Clustering (MCL)`: MCL Algorithm Groups genes into clusters (for example genefamily Clusters)
- by using pairwise sequence similarities from BLAST results in a similarity matrix.


- 1. `load_expression_data`:
- load expression data for different sets of genes (orthologs, paralogs, tandems)
- expression data represents count numbers (expressed RNA molecules per unit time for each   
  gene. showing  whether genes are active (ON), inactive (OFF))
- expresion data is used later on to generate expression profiles 
->  High expression indicates active genes producing a lot of RNA, 
->  while low expression indicates inactive or repressed genes producing little to no RNA.
- generate Expression Profils
- Gene Expression as Vectors: Gene expression profiles (e.g., RPKM normalized values) across tissues are treated as vectors in a multi-dimensional space where each axis represents a tissue.

- `Normalization of the Expression` (Transcriptome) Data 
- loads the species DNA-Fasta-Files, store them into Objects/Variables
- These files represent the genomic sequences of each species in Fasta format.
- Normalization Purpose: Adjusts for technical variations (e.g., sequencing depth, RNA quality) 
to ensure fair comparison of gene expression levels across different samples.
- Normalization Methods: TPM (Transcripts Per Million), FPKM and RPKM are commonly used methods


- 2. `load_paralogs_orthologs_tandems_data`:
- Gene-Groups
- load Information about the Orthologs, Paralogs and Tandems identified within the species Genomes
- the Files contain Information about Pairwise Sequence Similarities (Orthologs & Paralogs)

- 3. `load_genefamilies_data`:
- load Information about Gene Families . MCL output files provide the required information
- MCL generates Gene-cluster (Gene-Families) a cluster refers to a group of genes that are related 
  by sequence similarity


-------------------------------------------------------------------------------------------

## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `computeExpressionProfileDistances.R` 	      
- 2.  `investigateDistributionsOfExpressionProfileDistances.R`
- 3.  `t-tests`
- 4.  `angles`

- Ausgangspunkt quasi Gene-Familys Dataframe with genefamily clusters

- 1. Euclidean Distances for Gene Expression Analysis:
- purpose: to Measures and identify similarities or differences in gene activity. representing changes in gene activity within the cluster.
- Euclidean distances are calculated within each gene-family cluster 
- therefor an All vs. All comparison is performed within each gene-family-cluster (basically the expression count number of each gene)
- the gene family clusters have an array of genes associated with theire specific expression profiles (count number of expression)
- resulting array of distances within each gene-family cluster indicates how similar or different the gene-activity / gene-expression is 
- basically resulting in an comparison matrix with distances for all vs. all gene comparison in the specific cluster. 
- Small distances indicate similar gene expression levels, suggesting similar gene activity, while large distances suggest divergent activity levels.

- 2. statistical Analysis, investigate Distributions of ExpressionProfile Distances :
- statististics within or to the genefamily clusters and their euclidean distances
-  The R script summarizes distances per gene- group with statistics 
  including mean, median, max, and max-minus-min for both Orthologs and Paralogs.
- - for each family and each subset (Orthologs, Paralogs) measure the mean expression distance
- statistics: hauptsächlich Mean / Median - Euclidean Distances 


- 3. t test and wilcox test: show significance level of the calculated statistics
- Tests Across Distributions: Analyze distributions of distances (mean/median), expression angles, and tissue specificity/versatility, and apply appropriate statistical tests (t-test or Wilcoxon) to determine significant differences between orthologs and paralogs.

- Purpose of T-Tests: Statistical tests used to determine if there are significant differences between the means of two groups, 
(such as means in gene expression levels in different tissues or gene classes.)
-  Use a parametric t-test to compare the means of orthologs and paralogs
- Null hypothesis: the means are equal; alternative hypothesis: the means differ.
- Perform pairwise t-tests to compare orthologs vs paralogs for different distributions (mean, median distances, tissue-specific distances, expression angles, and versatility),
- assessing whether the means and overall distributions are significantly different.

- Wilcoxon Rank-Sum Test: Use this non-parametric test when assumptions for the t-test are not met. 
- It compares ranks to determine if one group tends to have larger values than the other.
- P-Value Adjustment: After obtaining p-values from multiple tests, apply Benjamini-Hochberg (BH) correction using p.adjust() to control the False Discovery Rate (FDR) and reduce the likelihood of false positives.

- 4. Angles , Angle Analysis:
- Angle Calculation: The script computes angles between expression vectors (of gene expression profiles in each gene family) and diagonals (function cosDiag),
- Diagonal Reference: The diagonal in this space refers to the line connecting the origin to the point (1,1,…), representing equal expression across tissues. 
Angles to the diagonal quantify the deviation from balanced expression. So Deviation from this Diagonal indicates genes which expression is more tissue specific
- Diagonal as equal distibution of a genes Activity across tissue. meaning expression level is the same in each tissue
- the expression vector basically represents a gene-family cluster with the respective genes and their gene expression profile across the tissues
- the mean expression vectors are used to calculte angles (statistics), 
- so the angles measuring rotation around the diagonal to assess functional diversification. 

- angles indicate: Tissue specifity and Tissue Versatility
the analysis investigates shifts in tissue specificity post-duplication.
Tissue Specificity: Angles between expression vectors and the diagonal help assess how specific or versatile a gene's expression is across tissues. Smaller angles suggest balanced expression; larger angles indicate tissue-specific expression.
Tissue Versatility: This concept measures how broadly a gene is expressed across multiple tissues. A low angle to the diagonal indicates high versatility, while a large angle indicates expression in specific tissues
- creating data frames to evaluate tissue specificity and the distinctiveness of gene expression profiles in each gene family

- comparing ortholog and paralog:
- Mean expression vectors for orthologs and paralogs are compared within each gene family.
- Vector Space for Orthologs and Paralogs: By calculating the angles between the diagonal and mean expression vectors of orthologs and paralogs, (within gene-family clusters, )


-------------------------------------------------------------------------------------------

## Section-3 - Plotting Distributions

                                       
- 1.  `plotDistributionsOfExpressionProfileDistances.R`	
- plot Distances
- plot Distances Tissue

- 2.  `plot_expressionAngles.R`				          

- Create two scientific plots, both boxplots, or raincloud plots (package grain)
  - Each plot should visualize the distribution of the expression distances in orthologs and paralogs
  - the first plot shows all measured distances the second only the mean distances

Plotting Angles: For gene families with a positive distance between ortholog and paralog expression vectors, the calculated angles are visualized using boxplots or raincloud plots to show the distribution of expression diversity.

Plot distributions of the above calculated angles, please, using boxplots and/or raincloud plots.



