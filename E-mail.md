

# `Is there a difference in expression diversity between Orthologs and Paralogs? Drosophila dataset:`


load("data/x-data.RData")

## 1. Gene-Family Datatable ---> 2 Subsets (Paralogs / Orthologs) 
---------------------------------------------------------------------------


- Gene-Family Datatable   ---> 2 Subsets (Paralogs / Orthologs) 

- cluster_1 ...
- cluster_2 ...
- cluster_3 ...
- 1244 Clusters (in total)

- 1. combined Families.df -> Clusters / Ortholog Count / Paralog Count /dmel Count


--> vv_df2_merged_449_rows
--> vv_df1_Family_id_more_then_1_both_449_rows

save(vv_df1,
     vv_df_1_top,
     vv_df1_Family_id_more_then_1_both_449_rows,
     vv_df1_Family_id_more_then_1_ortholog_count_657_rows,
     vv_df1_Family_id_more_then_1_paralog_count_781_rows,
     vv_df1_gene_family_ortholog_paralog_count_1244_rows,
     vv_df2_merged_449_rows,
     file = "data/vv_data.RData")

load("data/vv_data.RData")


- 2. combined Families.lst -> clustered ortholog genes and paralog genes 


## 2. Expression-Profiles 
---------------------------------------------------------------------------

load("data/expression_data.RData")

--> df14_P / df14_M 
--> df14_P_append / df14_P_append

df14_M
vv_df2_merged_449_rows
vv_df1_Family_id_more_then_1_both_449_rows
save(df14_M, vv_df2_merged_449_rows, vv_df1_Family_id_more_then_1_both_449_rows , file = "data.RData")


## 3. Measure Euclidean Distances 
---------------------------------------------------------------------------

- Measure for each gene family the pairwise expression distances of the subset of Paralogs and Orthologs

--> FamiliesExpressionProfileDistances --> orthologs / Paralogs

- gene_family_expression_distances = list()





## 4. statistics: Mean / Median - Euclidean Distances 
---------------------------------------------------------------------------


- for each family and each subset (Orthologs, Paralogs) measure the mean expression distance


for each gene-family f_i in drosophila families, do:
  fam_gene_ids = get_family_gene_member_ids( f_i )
    fam_orthologs = get_orthologs_of_genes( fam_gene_ids )
    fam_paralogs = get_paralogs_of_genes( fam_gene_ids )
    fam_ortholog_expression_distances = measure_expression_distances( fam_orthologs )
    fam_paralogs_expression_distances = measure_expression_distances( fam_paralogs )   
  gene_family_expression_distances.append(
    // key (family-ID) = value (list with two keys)
    f_i = list(
       ortholog_expression_dists = fam_ortholog_expression_distances,
       paralog_expression_dists = fam_paralogs_expression_distances
    )
  gene_family_expression_dists_stats.add_column(
    mean-ortholog-exp-dist = mean( fam_ortholog_expression_distances ),
    median-ortholog-exp-dist = median( fam_ortholog_expression_distances )
    // the same for mean and median of paralog expression dists
  )
  
 
  
- gene_family_expression_dists_stats  --> table( columns=c(Gene-Fam-ID, mean-ortholog-exp-dists, median-ortholog-exp-dists,
										    mean-paralog-exp-dists, median-paralog-exp-dists))
										



## 5. Plot the results (Scientific Plots)
------------------------------------------------------------------------------------

- Create two scientific plots, both boxplots, or raincloud plots (package grain)
  - Each plot should visualize the distribution of the expression distances in orthologs and paralogs
  - the first plot shows all measured distances
  - the second only the mean distances

Plots:

// median expression distances
boxplot( list(
  median-ortholog-expression-dists = gene_family_expression_dists_stats.get_column( median-ortholog-exp-dists ),
  median-paralogs-expression-dists = gene_family_expression_dists_stats.get_column( median-paralogs-exp-dists ),
))

// the same for the means

// example for all measured distances:
 
boxplot( list(
  all-ortholog-expression-dists = unlist(lapply( gene_family_expression_distances, function(x) {
     get_list_entry(x, "ortholog_expression_dists")
  })) ,
  all-paralogs-expression-dists = unlist(lapply( gene_family_expression_distances, function(x) {
     get_list_entry(x, "paralogs_expression_dists")
  }))
------------------------------------------------------------------------------------

## 6. Angles, Change in Tissue Versatility
------------------------------------------------------------------------------------


Then, we want the angles -> for all families, that have a positive distance between their ortholog and paralog expression vector clouds (slide 15).

- Change in tissue versatility (slide 16) 
- Change in tissue specificity (slide 18)

- Angles
- expression vector clouds

df14_P_dmel
paralog.lst
ortholog.lst


df_median_mean_paralogs.exp.prof.dists.tissue

df_median_mean_orthologs.exp.prof.dists.tissue


tandems.exp.prof.dists.orth.dist.df
 `orthologs.exp.prof.dists.stats.df`




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

Open questions and issues:
- Why do we divide the cosine by sqrt(2)? Try out the functions with given angles and see what comes out. Maybe then we can see why we do this division.
  - angles to try could be 0,45,90,..,180,..,360
  - maybe some kind of normalization?
  - check the function bodies. How does Asis calculate these angles?

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


---------------------------------------------------

- Can we get access to information about tandems in Drosophila
  - investigate the online resources
  - feel free to ask Jan-Lukas for help (put Asis in Cc)
  - feel free to ask Gemini or ChatGPT



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


a <- rnorm(10000, 100, 5)
b <- rnorm(10000, 1000, 50)
summary(a)
summary(b)
boxplot(list(a=a, b=b))
t.test(a, b, alternative='greater')
t.test(b, a, alternative='greater')
?t.test
t.test(b, a, alternative='greater')
savehistory('~/Desktop/t_test_examples.R')
