

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

- Expression Profiles Data:

Expression profile data were gathered from a shared drive.
The dataset contained expression values for four Drosophila species 
(D. melanogaster, D. simulans, D. yakuba, and D. erecta).
Three diet types were analyzed (M, P, etc.).
From this, the expression.profiles dataframe was derived.
For the analysis, we focused only on one diet and D. melanogaster, 
selecting expression data from four tissues: Fat Body, Muscle, Gut, and Whole Body.

Expression Dataframes: Based on the above, expression dataframes were created:
- df14_P_dmel
- df14_M_dmel

Details:

"P" and "M" represent the different diets used in the experiment.
For simplicity, only D. melanogaster data were used in this analysis.
df14_P_dmel represents RNA-seq expression profiles for one of the diets.


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


Then, we want the angles -> for all families, that have a positive distance between their ortholog and paralog expression vector clouds (slide 15).

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
  
  cosAngleVec <- function(a, b)
    phi_angle_orthologs_paralogs =
    // store in a table with columns Family-ID, delta_angle_orths, delta_angle_paral, phi_orth_paral
    ...


## 7. Plot distributions of the above calculated angles, please, using boxplots and/or raincloud plots.
------------------------------------------------------------------------------------
















# Open questions and issues:

- Why do we divide the cosine by sqrt(2)? Try out the functions with given angles and see what comes out. Maybe then we can see why we do this division.
  - angles to try could be 0,45,90,..,180,..,360
  - maybe some kind of normalization?
  - check the function bodies. How does Asis calculate these angles?

find out why is he dividing the cosine by sqrt(2)?
-----------------------------------

cosAngleVec <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'cosAngleVec(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    sum(a * b)/(euclNorm(a) * euclNorm(b))
}

cosDiag <- function(x, d.v = rep(1, length(x))) {
    cosAngleVec(x, d.v)
}

euclNorm <- function(x) {
    if (hasInvalidVecSpaceComponents(x)) 
        return(NaN)
    sqrt(sum(x^2))
}

paralog.expr.angle.diag.df <- data.frame(gene = para.expr, angle.diag = as.numeric(mclapply(para.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)


## 8. Pairwise t-tests and wilcox-tests
------------------------------------------------------------------------------------

For all distributions that we generate do pairwise tests, to see whether the respective empirical distributions
differ significantly:
- Are the mean values significantly different? Use a t-test; in R
  See ?t.test for details and the attached snippet
  - Compare orthologs vs paralogs
  - you can consider putting the comparisons into the plots - try it out and / or write all comparisons into a table
- Are the overall distributions different, i.e. the values of the first above those of the other?
  use wilcox ranked sum test:
  ?wilcox.test
- After obtaining all p-values correct them for multiple hypothesis testing with
  ?p.adjust( vector-of-p_values, method="BH")


Distributions of distances of paralogs and orthologs
Distribution of mean distances
Distribution of median distances
Distributions of distances per tissue
Distribution of mean distances per tissue
Distribution of mean distances per tissue
Expression angles
Expression versatility


find out why is he dividing the cosine by sqrt(2)?





 [8] "df14_M_dmel"
[10] "df14_P_dmel"

[13] "orthologs.exp.prof.dists"
[17] "orthologs.exp.prof.dists_stats_df"
[16] "orthologs.exp.prof.dists_dataframe"
[14] "orthologs.exp.prof.dists.tissue"
[15] "orthologs.exp.prof.dists.tissue_dataframe"
[1]  "df_median_mean_orthologs.exp.prof.dists.tissue"

[22] "paralogs.exp.prof.dists"
[25] "paralogs.exp.prof.dists_dataframe"
[26] "paralogs.exp.prof.dists_stats_df"
[23] "paralogs.exp.prof.dists.tissue"
[24] "paralogs.exp.prof.dists.tissue_dataframe"
 [3] "df_median_mean_paralogs.exp.prof.dists.tissue"


[28] "vv_df2_merged_449_rows"
[30] "vv_df3_ortho"
[31] "vv_df3_ortho.lst"
[32] "vv_df3_para"
[33] "vv_df3_para.lst"

[18] "orths.expr.angle.diag.df"
[21] "paralog.expr.angle.diag.df"