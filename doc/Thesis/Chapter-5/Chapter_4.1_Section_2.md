



#### **Section 2: Distance and Statistical Computations**
# ------------------------------------------------------------------------------------------------
`Section - 2: Compute pairwise gene expression profile distances within gene groups`

- **Scripts**:
  - `compute_exp.prof.dists.R`: Calculates Euclidean distances.
  - `compute_exp.prof.dists_statistics.R`: Summarizes distance statistics.
  - `compute_exp.prof.dists_angles.R`: Computes cosine angles for expression vectors.
- **Output**: Results stored in R objects for further analysis (e.g., `ExpressionProfileDistances.RData`).

Section-2 - Computing Distances, Statistics 

### Compute pairwise gene expression profile distances within gene groups

> Differences in expression are correlated with function diversification. In
> order to detect such diversification of function we measure euclidean
> distances between expression profiles of genes within gene groups. Distances
> are computed between whole expression profiles as well as on each euclidean
> dimension (tissue).

### Expression Profile based Function Diversity

In this step the euclidean distances in the expression vector space are
invetigated in order to find out whether belonging to a certain gene group,
i.e. gene family, orthologs, paralogs, or tandem duplicates has an influence on
how gene expression changes. These changes are measured in distances between
mean expression vectors, e.g. between the mean expression of orthologs and the
mean expression of paralogs, all of which belong to the same gene family.

So, the above, among others, tests whether within each gene family the gene
expression of the orthologs is different from the gene expression of the
paralogs by comparing the vector clouds of orthologs and paralogs,
respectively. See slide fifteen (15) of `./inst/ExpressedGeneGroupsVenn.pdf`.

Results of this step are stored as plots and tables of t-Test outcomes. See the
same section in Vignette `./vignettes/GeneFamilies.Rmd` and the respective
Rscripts, respectively.

I think the plots are those seen on slide fourteen (14).

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

