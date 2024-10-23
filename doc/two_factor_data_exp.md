

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

- open issues from the last meeting:

  - do the expression distances and angle measures for tandems and trans-duplicates separately

  - species and tissue two factor analysis

- If the biological question is whether after duplication there was a change in gene expression we might want to measure the distances between pairs of one ortholog and one paralog. The assumption being that the ortholog was the ancestor (conserved expression) of the paralog (duplicated and evolved expression). For this let's do an ANOVA like ChatGPT suggested:
  https://chatgpt.com/share/66fc1cb1-6c4c-800b-a190-5d0d12972b41

  - This might even work better for our multi-species setup. Do one additional ANOVA in this species~tissue setting:

     - Create expression vectors in an expression vector space JUST defined by tissues, i.e. the axes are tissues

     - Using gene identity we now can measure distances between species, i.e. dist( Gene-A~species-1, Gene-A~species-2)
     - Now we can categorize these distances by whether a gene is an ortholog or an paralog and investigate whether gene expression appears to be more conserved in orthologs or in paralogs. For this, we do an ANOVA, comparing distances between orthologs with distances between paralogs. 
       https://docs.google.com/drawings/d/112ljytgADb1bCjDSfnmwJJgfjFurd6XyaJCocj2KMIU/edit?usp=sharing

After we have all our analyses clear and done, we will redo them using logarithmic transformation of RPKM values.

Soon, we will collect all results and interprete them before more experiments.
