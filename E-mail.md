

# `Is there a difference in expression diversity between Orthologs and Paralogs? Drosophila dataset:`

## 1. Gene-Family Datatable ---> 2 Subsets (Paralogs / Orthologs) 
---------------------------------------------------------------------------

- Drosophila Datasets:

- Two datasets were used: one for orthologs and one for paralogs, both belonging to the same gene family.
Data were obtained from FlyBase and supplemented with Orthofinder analysis to retrieve ortholog data. 

- A gene-family data table was created, containing genes clustered into gene family clusters.
Paralogs and orthologs were mapped to their respective gene family clusters.
Rows with fewer than two genes per gene family cluster were removed for further analysis.
The resulting dataframe, vv_df2_merged_449_rows, was used to create distributions of expression profiles and other analyses.
vv_df2_merged_449_rows Details:

The dataframe contains 449 rows, 1,244 gene family clusters in total, 
and two subsets (paralogs and orthologs).

From this, a dataframe was created to count the number of genes in each gene family cluster.
The distribution of counts per cluster was then plotted using the plot_distribution_functions.R script.

- File locations:

- R script:   R/plot_distribution_functions.R
- Plot:       plots/distribution_of_count.pdf


## 2. Expression-Profiles 
---------------------------------------------------------------------------

`Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive
 Systems among Drosophila Species`

from: https://www.cell.com/cell-reports

The paper investigates how different Drosophila species exhibit distinct responses 
to carbohydrate-rich diets, highlighting variations in gene regulation and metabolic pathways.

Generalist Species Diets: Generalist species like D. melanogaster and D. simulans 
thrive on a broad range of rotting fruits, vegetables, and plant matter, 
which contain varying carbohydrate levels.
Specialist Species Diets: Specialist species such as D. sechellia feed exclusively 
on specific plants like Morinda citrifolia, which may have a more specialized 
nutritional composition, requiring unique metabolic adaptations

Three diet types were analyzed (M, P, C).
Distinct Diets (M and P): The M diet is carbohydrate-rich and protein-poor, while the P diet is protein-rich with low carbohydrates. These diets help explore how Drosophila species respond to different nutrient balances. C Diet (Control): A balanced diet with moderate levels of carbohydrates and proteins, used as a reference in studies

- Expression Profiles Data:

The dataset contained expression values for four Drosophila species 
(D. melanogaster, D. simulans, etc.) --> directory: data/mmc7.xlsx
From this, the expression.profiles dataframe was derived.
For the analysis, we focused only on one diet and D. melanogaster, 
selecting expression data from four tissues: Fat Body, Muscle, Gut, and Whole Body.
For simplicity, only D. melanogaster data were used in this analysis.

Expression Dataframes: Based on the above, expression dataframes were created:
- df14_P_dmel
- df14_M_dmel
- df14_P_dmel represents RNA-seq expression profiles for one of the diets.
- for the 3 samples taken for each tissue , the mean value was calculated 



## 3. Measure Euclidean Distances 
---------------------------------------------------------------------------

- Measure for each gene family the pairwise expression distances of the subset of Paralogs and Orthologs
were measured 
- basically calculating the euclidean distance between the expression profiles of each gene in the gene-family cluster
- so each geneFamily cluster was treated as a vector of Gene names belonging to that cluster, in the subset of orthologs and paralogs genes
- within each cluster the expression profiles of each gene were compared to each other, basically all vs all gene expression values gathered from the expression profiles were compared and the distances were calculated pairwise for each gene matchup in the Gene Family cluster

--> GeneFamilies dataframe : df14_P_dmel   (P diet and dmel melanogaster)
--> rna.seq.exp.profiles.df : df14_P_dmel as (rna.seq.exp.profils)

- Input files used in this script for caluclating the euclidean distances:
- RScript: exec/calculate_exp.prof.dists.R togerther with R/calculate_expression_profile_dists_function.R
-> vv_df3_ortho.lst as (orthologs.lst)
-> vv_df3_para.lst as (paralogs.lst)
-> df14_P_dmel as (rna.seq.exp.profils) 

- within 2 Subgroups (Paralogs / Orthologs)
- also per Tissue: (Fat_body, Muscle, Gut, Whole_body)

Output Files created:
orthologs.exp.prof.dists 
orthologs.exp.prof.dists.tissue 

paralogs.exp.prof.dists 
paralogs.exp.prof.dists.tissue 


## 4. statistics: Mean / Median - Euclidean Distances 
---------------------------------------------------------------------------


- for each family and each subset (Orthologs, Paralogs) measure the mean expression distance

- Pseudocode provided:

for each gene-family f_i in drosophila families, do:
  fam_gene_ids = get_family_gene_member_ids( f_i )
    fam_orthologs = get_orthologs_of_genes( fam_gene_ids )
    fam_paralogs = get_paralogs_of_genes( fam_gene_ids )
    fam_ortholog_expression_distances = measure_expression_distances( fam_orthologs )
    fam_paralogs_expression_distances = measure_expression_distances( fam_paralogs ) 

  gene_family_expression_distances.append(// key (family-ID) = value (list with two keys)
    f_i = list(
       ortholog_expression_dists = fam_ortholog_expression_distances,
       paralog_expression_dists = fam_paralogs_expression_distances
    )
  gene_family_expression_dists_stats.add_column(
    mean-ortholog-exp-dist = mean( fam_ortholog_expression_distances ),
    median-ortholog-exp-dist = median( fam_ortholog_expression_distances )
    // the same for mean and median of paralog expression dists
  )
  
 
- exec\calculate_mean_median_statistics.R was used to calculate the mean and median of the expression distances 
for the 2 Subsets each

# df_median_mean_paralogs.exp.prof.dists
# df_median_mean_orthologs.exp.prof.dists

# df_median_mean_paralogs.exp.prof.dists.tissue
# df_median_mean_orthologs.exp.prof.dists.tissue
  									

## 5. Plot the results (Scientific Plots)
------------------------------------------------------------------------------------

- Create two scientific plots, both boxplots, or raincloud plots (package grain)
  - Each plot should visualize the distribution of the expression distances in orthologs and paralogs
  - the first plot shows all measured distances the second only the mean distances

- exec\calculate_mean_median_statistics.R
- R/plot_distribution_functions.R

-> siehe : plots directory

## 6. Angles, Change in Tissue Versatility
------------------------------------------------------------------------------------

Tissue Versatility: In the context of gene expression, tissue versatility refers to 
how similar or distinct the gene expression profiles are across different tissues. 
Greater versatility would mean that a gene is expressed in many tissues, 
while less versatility suggests that the gene expression is more restricted to certain tissues.

Euclidean Distances Between Expression Profiles: The Euclidean distance 
between two points (representing expression levels in tissue X and Y) reflects 
how different the gene expression profiles are between the two tissues. In the plot:

 ) appear to represent the Euclidean distances between gene expression profiles 
 of two groups (e.g., gene sets) in tissues X and Y.
Angles in the Plot: The angles in the plot seem to represent the degree 
of separation in gene expression profiles between the tissues:

A smaller angle between two points in the gene expression space (such as for 

 ) suggests greater similarity in expression profiles.
A larger angle indicates more divergent gene expression patterns between the two tissues.

Then, we want the angles -> for all families, that have a positive distance between 
their ortholog and paralog expression vector clouds (slide 15).

- Change in tissue versatility (slide 16) 
- Change in tissue specificity (slide 18)

Pseudo-Code:

fams_with_separated_orth_para_expr_clouds =
  gene_family_expression_dists_stats[ which(
    dist_between_ortholog_and_paralog_expression_clouds > 0
  ) ].get_column( Gene-Family-ID )

for each gene-family f_i in fams_with_separated_orth_paral_expr_clouds
  fam_gene_ids = get_family_gene_member_ids( f_i )
    fam_orthologs = get_orthologs_of_genes( fam_gene_ids )
    fam_paralogs = get_paralogs_of_genes( fam_gene_ids )
      // delta angle from slide 16
 
    delta_angle_orthologs = cosDiag( fam_orthologs )
    delta_angle_paralogs = cosDiag( fam_paralogs )



    exprVecSpaceEvolutionAfterDupl <- function(genes, family.name, classifier.lst, 
  



  !!!
  cosAngleVec <- function(a, b)
    phi_angle_orthologs_paralogs =
    // store in a table with columns Family-ID, delta_angle_orths, delta_angle_paral, phi_orth_paral
    ...


## 7. Plot distributions of the above calculated angles, please, using boxplots and/or raincloud plots.
------------------------------------------------------------------------------------








## 8. Pairwise t-tests and wilcox-tests
------------------------------------------------------------------------------------

For all distributions that we generate do pairwise tests, 
to see whether the respective empirical distributions differ significantly:


- Are the mean values significantly different? 
  - Compare orthologs vs paralogs
- Are the overall distributions different, i.e. the values of the first above those of the other?
  use wilcox ranked sum test:
- After obtaining all p-values correct them for multiple hypothesis testing with
  p.adjust( vector-of-p_values, method="BH")


- Distributions of distances of paralogs and orthologs
Distribution of mean distances
Distribution of median distances
- Distributions of distances per tissue
Distribution of mean distances per tissue
Distribution of median distances per tissue
- Expression angles
- Expression versatility



T-test Overview: 

-  A parametric test used to compare the means of two groups to determine 
if they are significantly different. Assumes normal distribution and equal variances.
Assumes Normality of data distribution and Homogeneity of variance (in both groups).

Null hypothesis (H₀): The means of the two groups are equal.
Alternative hypothesis (H₁): The means of the two groups are different.


Wilcoxon Test (Wilcoxon Rank-Sum Test) Overview: 

- A non-parametric alternative to the t-test, used when the assumptions 
of normality or equal variances are not met. It compares ranks instead of means.
The data can be ordinal or continuous but doesn’t require normality.
It assesses whether one group tends to have larger values than the other.

Null hypothesis (H₀): The distributions of the two groups are the same.
Alternative hypothesis (H₁): The distributions of the two groups are different.

Test statistic (W or U): Measures the sum of ranks; large or small values indicate differences between groups.
p-value: If p < 0.05, reject the null hypothesis, indicating significant differences.
Effect size (r): Indicates the magnitude of the difference; small (≈0.1), medium (≈0.3), and large (≈0.5).


P-Value Adjustment (p.adjust): 

After obtaining p-values from multiple comparisons (t-tests and Wilcoxon tests), 
they are adjusted to control for multiple hypothesis testing using p.adjust(). 
This helps to reduce false positives.


Benjamini-Hochberg (BH) Correction:

This method is used to adjust p-values for multiple hypothesis testing by 
controlling the False Discovery Rate (FDR). It ensures a balance between detecting 
true effects and minimizing false discoveries.




## Open questions and issues: `Why do we divide the cosine by sqrt(2)?`
------------------------------------------------------------------------------------

- Try out the functions with given angles and see what comes out. 
- Maybe then we can see why we do this division. (angles to try could be 0,45,90,..,180,..,360)


- Functions: 

- Angle for each row Vector group (a and b) -> Vector a and b (contain numeric dist values)
- a as x Vector and b as diagonal Vector  d.v
- Vektor b ein Diagonalvektor ist, das heißt, er ist ein Vektor, in dem alle Einträge gleich 1 sind
- b=[1,1,1,…] , n ist die Länge des Vektors 𝑏 also die Anzahl der Elemente in 𝑎
- Ein höheres 𝑛 (mayne as the number of tissues present in the rpkm.expr.profiles for the vector)
  führt zu einer kleineren Größe des Terms sqrt(2) / sqrt(n), was bedeutet, dass der Wert von cosDiag
  cosDiag verringert wird, wenn mehr Datenpunkte betrachtet werden.


cosDiag(x = rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == x), expr.cols], )/`sqrt(2)`

a = x = = rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == x), expr.cols]
b = d.v = rep(1, length(x))

cosDiag --> {`sum(a * b)/(sqrt(sum(a^2)) * sqrt(sum(b^2)))`} /`sqrt(2)`

--> `(sum(a * b) * sqrt(2)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))`

--> `(sum(a) * sqrt(2)) / (sqrt(n) * sqrt(sum(a^2)))`

-->  `sqrt(2) / sqrt(n)`





- 
Dividing by 2 scales distances and angles uniformly, preserving the shape of vector clouds but does not fully normalize the data.
In high-dimensional spaces like gene expression data, dividing by 2 helps maintain consistent scaling for comparisons of distances and angles, especially between orthologs and paralogs.
This rescaling aids in gene expression analysis by standardizing the interpretation of tissue-specific and versatile expression patterns across datasets.

The expression `sqrt(2) / sqrt(n)`  is a form of standardization that accounts for the influence of the number of elements on the measure of 
cosDiag
cosDiag. This can help in comparing results when dealing with different data set sizes or varying contexts.

Correlation and Variation:
The term also impacts the interpretation of 
cosDiag
cosDiag in the context of correlation or variation of the data, as it adjusts the result relative to the number of data points.


At 0° (or 360°), the cosine value is 1, indicating maximum alignment, where vectors point in the same direction, resulting in perfect similarity. At 45°, the cosine is approximately 0.707, suggesting partial alignment and moderate similarity between the vectors. At 90°, the cosine value is 0, meaning the vectors are perpendicular, showing no similarity or projection between them.

At 180°, the cosine value is -1, indicating that the vectors are pointing in opposite directions, representing maximum dissimilarity. Similarly, at 270°, the cosine value is 0 again, indicating perpendicular vectors but with a reversed orientation. Finally, at 360° (or back to 0°), the cosine returns to 1, indicating perfect alignment once more after completing a full rotation.

without sqrt(2) 
As n increases (more data points or higher dimensions), the overall result becomes smaller because the denominator increases. For large n, this would reduce the impact of the cosine similarity value.

For smaller n, the scaling would have less of an impact, and the cosine similarity would remain closer to its unscaled value.

This reflects how larger datasets or higher-dimensional data may require additional normalization to account for complexity.