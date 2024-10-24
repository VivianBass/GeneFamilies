
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

- 1.  `load_expression_data.R`	          
- 2.  `load_paralogs_orthologs_tandems_data.R`    
- 3.  `load_genefamilies_data.R`

Tools for identifying paralogs, orthologs, tandems, and gene families:

- **BLAST** or **Diamond**: 
These tools identify regions of similarity between 
biological sequences (DNA, RNA, or proteins) using local alignment. 
They score alignments based on sequence similarity to find 
homologous sequences (orthologs, paralogs) by comparing query sequences to a database.

- **Markov Clustering (MCL)**: 
MCL groups genes into clusters (e.g., gene families) 
by creating a similarity matrix from pairwise sequence similarities derived from BLAST results.

- 1. `load_expression_data`:

Load expression data for different gene sets (e.g., orthologs, paralogs, tandem genes). 
Expression data, represented as RNA count numbers per unit time, indicates gene activity: 
high expression shows genes are active (ON), producing RNA, 
while low expression shows inactive or repressed genes (OFF), producing little to no RNA.

Generate gene expression profiles. Expression values (e.g., RPKM) across tissues 
are treated as vectors in multi-dimensional space, where each axis represents a tissue, 
capturing gene activity across different biological contexts.

- Normalization of expression (transcriptome) data:

Load species DNA Fasta files, representing genomic sequences, into objects/variables.
  
**Normalization Purpose**: Adjusts for technical variations like sequencing depth 
and RNA quality to allow fair comparison of gene expression levels across different samples.

**Normalization Methods**: TPM (Transcripts Per Million), FPKM, and RPKM are commonly 
used to normalize gene expression data for accurate analysis.

- 2. `load_paralogs_orthologs_tandems_data`:

Gene Groups: Load information on orthologs, paralogs, and tandem genes 
identified within species genomes.

The files contain data on pairwise sequence similarities for orthologs and paralogs, 
used to identify these relationships.

- 3. `load_genefamilies_data`:

Load information on gene families from MCL output files.
MCL generates gene clusters (gene families), 
where each cluster groups genes that are related by sequence similarity.


## Section-2 - Computing Distances, Statistics, T-tests and Angles
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`
- 3.  `compute_t-tests_wilcox-tests.R`
- 4.  `compute_exp.prof.dists_angles.R`  

- Starting point: a gene family dataframe containing gene family clusters.

- 1.  `compute_exp.prof.dists.R` 

**Euclidean Distances for Gene Expression Analysis**:

**Purpose**: Measure similarities or differences in gene activity, 
highlighting changes within gene-family clusters.
  
Euclidean distances are calculated through an all-vs-all comparison 
of expression counts within each gene-family cluster.
  
Each gene-family cluster consists of genes with associated expression profiles (count numbers), 
forming an array of distances.

The resulting comparison matrix shows distances for all gene pairs within the cluster. 
Small distances indicate similar gene expression (similar activity), 
while larger distances suggest divergent activity levels.


- 2.  `compute_exp.prof.dists_statistics.R`

**Statistical Analysis of Expression Profile Distances**:

**Purpose**: Investigate distributions of Euclidean distances within gene-family clusters.

The R script summarizes distances per gene group (Orthologs, Paralogs) with key statistics, 
including mean, median, max, and max-minus-min.

For each gene family and subset (Orthologs, Paralogs), the mean and median Euclidean distances 
are calculated to assess similarities in gene expression levels.


- 3.  `compute_t-tests_wilcox-tests.R`

**T-Test and Wilcoxon Test: Significance of Statistical Analysis**

**Tests Across Distributions**: Analyze distributions of distances (mean/median), 
expression angles, and tissue specificity/versatility. Apply t-tests or Wilcoxon tests 
to determine significant differences between orthologs and paralogs.

**Purpose of T-Tests**: Parametric tests used to check if the means of two groups 
(e.g., gene expression levels across tissues or gene classes) are significantly different. 
Compare the means of orthologs vs. paralogs using pairwise t-tests across metrics 
like mean/median distances, tissue-specific distances, and expression angles.  
  - **Null Hypothesis**:        The means are equal.
  - **Alternative Hypothesis**: The means differ.

**Wilcoxon Rank-Sum Test**: A non-parametric alternative when t-test assumptions aren't met, 
comparing ranks to assess if one group tends to have larger values.

**P-Value Adjustment**: After running multiple tests, apply Benjamini-Hochberg (BH) 
correction with p.adjust() to control the False Discovery Rate (FDR) and minimize false positives.


- 4.  `compute_exp.prof.dists_angles.R`


**Angle Analysis in Gene Expression**:

**Angle Calculation**: The script computes angles between expression vectors (gene expression profiles) 
within gene families and the diagonal using the `cosDiag` function.

**Diagonal Reference**: The diagonal represents equal expression across all tissues 
(a line from the origin to the point (1,1,...)). Angles to this diagonal measure 
deviations from balanced expression, indicating tissue specificity.

**Interpretation**: 
  - **Small Angles**: Suggest balanced gene expression across tissues (high tissue versatility).
  - **Large Angles**: Indicate tissue-specific expression.

**Expression Vectors**: Represent gene-family clusters with gene expression profiles across tissues. 
Mean expression vectors are used to calculate angles, 
assessing functional diversification through rotation around the diagonal.

**Tissue Specificity and Versatility**:
Angles measure how specific or versatile a gene's expression is. 
Small angles indicate high versatility (broad expression), 
while larger angles show more tissue-specific expression.

**Ortholog vs. Paralog Comparison**:
Mean expression vectors for orthologs and paralogs are compared within 
each gene family by calculating angles between the diagonal and their mean expression vectors. 
This helps evaluate tissue specificity and functional shifts post-duplication.


## Section-3 - Plotting Distributions

- 1. `plot_expression_angles.R`    


- **Scientific Plots Overview**:

**Boxplots of Expression Distances**:
**First Plot**: ...
**Second Plot**: ...

**Boxplot of Expression Angles**:
Visualizes the distribution of calculated angles between expression vectors and 
the diagonal for gene families with a positive distance between ortholog and paralog expression profiles.