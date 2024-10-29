

# - two factor data (ST-Exp)

- tissue factor >= & species factor >= 4 (at least for family members)
- paralogs, orthologs, tandems, (maybe trans duplicates)
- everything in one distribution plot: boxplot with gitter points or raincloud plot
- measure the expression distances separately for each species, and plot the distances for each species:

  So, run another experiment with our EasyVectorOmics code:

    - Measure expression distances among orthologs using the above species_tissue two factor data (ST-Exp)

    - Measure the same for paralogs in the ST-Exp

    - Measure the same for tandem and trans duplicates in the ST-Exp

    - put everything in one distribution plot: boxplot with gitter points or raincloud plot


    Identify groups of orthologs that have at least four members across the species:

    - measure the expression distances separately for each species, and plot the distances for each species:

    one distribution (e.g. boxplot) per species; see
    https://docs.google.com/drawings/d/1OY8pnMThwUQJgPdPNQCKAOSiJVwH6bo3TWeL5zYMWo0/edit?usp=sharing
    - measure the same paralogs, tandems, trans duplicates 
    - the latter only if, we have enough groups that actually have more (>=) than three, better four members.


    For all distributions that we generate do pairwise tests, to see whether the respective empirical distributions differ significantly?

    - Are the mean values significantly different? Use a t-test; in R- See ?t.test for details and the attached snippet

    - Compare orthologs vs paralogs, orthologs vs tandem, ortholog vs trans duplicated

    - and all pairs of species in the other above experiment

    - you can consider putting the comparisons into the plots - try it out and / or write all comparisons into a table

    - Are the overall distributions different, i.e. the values of the first above those of the other?
    use wilcox ranked sum test:

    ?wilcox.test

    - After obtaining all p-values correct them for multiple hypothesis testing with
    ?p.adjust( vector-of-p_values, method="BH")


To Do:

  - species and tissue two factor analysis, species~tissue setting:
  - Create expression vectors in an expression vector space JUST defined by tissues, i.e. the axes are tissues

     - Using gene identity we now can measure distances between species, i.e. dist( Gene-A~species-1, Gene-A~species-2)
     - Now we can categorize these distances by whether a gene is an ortholog or an paralog and investigate whether gene expression appears to be more conserved in orthologs or in paralogs. For this, we do an ANOVA, comparing distances between orthologs with distances between paralogs. 
       https://docs.google.com/drawings/d/112ljytgADb1bCjDSfnmwJJgfjFurd6XyaJCocj2KMIU/edit?usp=sharing

After we have all our analyses clear and done, we will redo them using logarithmic transformation of RPKM values.

### --------------------------------------------------

two-factor data (ST-Exp) should contain and its purpose:

1. **Species and Tissue as Primary Factors:** ST-Exp should include gene expression data across multiple species and tissues, treating both species and tissue types as two independent factors.

2. **Minimum Family Members:** Each family group must have at least four members across different species to allow for meaningful comparisons.

3. **Gene Duplication Types:** Measure expression distances for orthologs, paralogs, 5 groups

4. **Expression Distance Metrics:** Calculate expression distances per gene type (ortholog, paralog, etc.) across species, with one distribution plot per species (boxplot/raincloud) and pairwise statistical testing (t-tests, Wilcoxon tests).

5. **ANOVA for Conservation Analysis:** Compare expression conservation across species by measuring distances in a tissue-defined vector space, assessing if orthologs retain expression more than paralogs.