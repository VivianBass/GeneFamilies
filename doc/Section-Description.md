
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

- what it does 

- Is there a difference in expression diversity between Orthologs and Paralogs?
- Measure for each gene family the pairwise expression distances of the subset of Paralogs and the other subset of Orthologs
- for each family and each subset (Orthologs, Paralogs) measure the mean expression distance
- Create two scientific plots, both boxplots, or raincloud plots (package grain)
- Each plot should visualize the distribution of the expression distances in orthologs and paralogs
- the first plot shows all measured distances
- the second only the mean distances



The steps and data required to carry out the original analysis are devided into 4 Sections 

- `1. Load tissue specific expression data`                                     (Section-1) 
- `2. Compute pairwise gene expression profile distances within gene groups`    (Section-2) 
- `3. Expression Profile based Function Diversity`                              (Section-3) 
- `4. Expression vector space analysis`                                         (Section-4)

-------------------------------------------------------------------------------------------

## Section-1 - Loading Data

- 1.  `loadRpkmRnaSeqCounts.R`
- 1.  `generateRpkmExpressionProfiles.R`	
- 2.  `loadcodingSequences.R`	
- 1.  `LoadOrthologsAndTandems.R`
- 3.  `loadGeneFamilies.R`	

- load expression data for different sets of genes (orthologs, paralogs, tandems)
- expression data represents count numbers (expressed RNA molecules per unit time for each   
  gene. showing  whether genes are active (ON), inactive (OFF))
- expresion data is used later on to generate expression profiles 

->  High expression indicates active genes producing a lot of RNA, 
->  while low expression indicates inactive or repressed genes producing little to no RNA.

- loads the species DNA-Fasta-Files, store them into Objects/Variables
- These files represent the genomic sequences of each species in Fasta format.

- `Normalization of the Expression` (Transcriptome) Data 

- Normalization Purpose: Adjusts for technical variations (e.g., sequencing depth, RNA quality) 
to ensure fair comparison of gene expression levels across different samples.

- Normalization Methods: TPM (Transcripts Per Million), FPKM and RPKM are commonly used methods

Identifies regions of similarity between biological sequences (DNA, RNA, proteins) using local alignment and scores alignments based on sequence similarity.
Helps find homologous sequences (orthologs, paralogs) by comparing query sequences against a database.

- Tools used in Section 2: Blast and Markov Clustering (MCL)
- BLAST (Basic Local Alignment Search Tool):

- Markov Clustering (MCL):

Groups genes into clusters (e.g., orthologs, paralogs) by using pairwise sequence similarities from BLAST results in a similarity matrix.
The MCL algorithm processes this matrix to produce gene families, outputting clusters, which can be further formatted for analysis.

- load Information about the Orthologs Paralogs and Tandems identified within the species Genomes
- the Files contain Information about Pairwise Sequence Similarities (Orthologs & Paralogs)

- load Information about Gene Families . MCL output files provide the required information
- MCL generates Gene-cluster (Gene-Families) a cluster refers to a group of genes that are related 
  by sequence similarity

-------------------------------------------------------------------------------------------

## Section-2 - computing distances, statistics and angles
			                                      
- 4.  `computeExpressionProfileDistances.R` 	      
- 1.  `investigateDistributionsOfExpressionProfileDistances.R`

Euclidean Distance for Gene Expression Analysis:

Measures pairwise similarities or differences in gene activity by calculating Euclidean distances between expression profiles within gene-family clusters.
Small distances indicate similar gene expression levels, suggesting similar gene activity, while large distances suggest divergent activity levels.
All-vs-all comparisons are performed within each gene cluster, treating the expression profiles as vectors to assess functional diversification.


-------------------------------------------------------------------------------------------

## Section-3 - plotting

                                       
- 2.  `plotDistributionsOfExpressionProfileDistances.R`				          


Statistical Analysis: The R script summarizes distances per gene group with statistics 
including mean, median, max, and max-minus-min for both Orthologs and Paralogs.

- for each family and each subset (Orthologs, Paralogs) measure the mean expression distance
- statistics: Mean / Median - Euclidean Distances 

- for each family and each subset (Orthologs, Paralogs) measure the mean expression distance
- statistics: Mean / Median - Euclidean Distances 

Class-Specific Distances: Distances are categorized into Orthologs and Non-Orthologs, 
with statistics computed separately for each class within gene clusters.

Gene Group Lists: The analysis involves three main gene group lists—orthologs, tandems, 
and families—each containing clusters of genes within specific gene families.

Distance Calculations: For each cluster, Euclidean distances between gene expression 
profiles are calculated pairwise (all vs all), representing changes in gene activity within the cluster.

In order to investigate which classes of gene groups mostly contribute to developmental diversity, 
measure and subsequently compare the distributions of expression profile distances 
within the different tissues.

Purpose of T-Tests: Statistical tests used to determine if there are significant differences 
between the means of two groups, such as gene expression levels in different tissues or gene classes.


- Create two scientific plots, both boxplots, or raincloud plots (package grain)
  - Each plot should visualize the distribution of the expression distances in orthologs and paralogs
  - the first plot shows all measured distances the second only the mean distances


-------------------------------------------------------------------------------------------

## Section-4

 		                              
- 2.  `expressionAngles.R`	                                                  


Gene Expression as Vectors: Gene expression profiles (e.g., RPKM normalized values) across tissues are treated as vectors in a multi-dimensional space where each axis represents a tissue.

Angle Analysis: The study explores angles between mean expression vectors of two sub-groups (orthologs and paralogs) within gene-family clusters, measuring rotation around the diagonal to assess functional diversification.

Diagonal Reference: The diagonal in this space refers to the line connecting the origin to the point (1,1,…), representing equal expression across tissues. Angles to the diagonal quantify the deviation from balanced expression.

Tissue Specificity: Angles between expression vectors and the diagonal help assess how specific or versatile a gene's expression is across tissues. Smaller angles suggest balanced expression; larger angles indicate tissue-specific expression.

Ortholog vs. Paralog Comparison: Mean expression vectors for orthologs and paralogs are compared within each gene family. Non-overlapping standard deviation regions between these vectors suggest functional changes after gene duplication.

Vector Space for Orthologs and Paralogs: By calculating the angles between the diagonal and mean expression vectors of orthologs and paralogs, the analysis investigates shifts in tissue specificity post-duplication.

Tissue Versatility: This concept measures how broadly a gene is expressed across multiple tissues. A low angle to the diagonal indicates high versatility, while a large angle indicates expression in specific tissues.

Angle Calculation: The script computes angles between expression vectors and diagonals (cosDiag), creating data frames to evaluate tissue specificity and the distinctiveness of gene expression profiles in each gene family.

Plotting Angles: For gene families with a positive distance between ortholog and paralog expression vectors, the calculated angles are visualized using boxplots or raincloud plots to show the distribution of expression diversity.

Functional Diversification: The analysis of angles between orthologs and paralogs highlights post-duplication functional diversification, reflected by changes in tissue-specific gene expression.

Plot distributions of the above calculated angles, please, using boxplots and/or raincloud plots.
-------------------------------------------------------------------------------------------

- t test and wilcox test

Pairwise Comparison of Distributions: Perform pairwise tests to compare orthologs vs paralogs for different distributions (mean, median distances, tissue-specific distances, expression angles, and versatility), assessing whether the means and overall distributions are significantly different.

T-Test: Use a parametric t-test to compare the means of orthologs and paralogs if data assumptions are met (normality and equal variance). Null hypothesis: the means are equal; alternative hypothesis: the means differ.

Wilcoxon Rank-Sum Test: Use this non-parametric test when assumptions for the t-test are not met. It compares ranks to determine if one group tends to have larger values than the other.

P-Value Adjustment: After obtaining p-values from multiple tests, apply Benjamini-Hochberg (BH) correction using p.adjust() to control the False Discovery Rate (FDR) and reduce the likelihood of false positives.

Tests Across Distributions: Analyze distributions of distances (mean/median), expression angles, and tissue specificity/versatility, and apply appropriate statistical tests (t-test or Wilcoxon) to determine significant differences between orthologs and paralogs.