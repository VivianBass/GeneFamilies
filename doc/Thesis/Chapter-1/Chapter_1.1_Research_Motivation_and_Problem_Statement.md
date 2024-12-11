

## 1.1. **Research Motivation and Problem Statement**



  A. Research Objectives
      - Clear statement of research goals
      - Specific research questions and hypotheses
      - Scope and limitations of the study

  B. Problem Definition
      - Description of the research problem
      - Significance and relevance
      - Current challenges in the field

  C. Conceptual Framework
      - Overview of ortholog vs paralog comparison
      - Key terminology and concepts
      - Theoretical background


- **Objective of the EasyVectorOmics Analysis**

In this project the goal is to reproduce the vector space analyses as done in
the _Cardamine hirsuta_ genome project [1]. The R-Code required to do so shall
be isolated and made usabel with _any_ data. Next, this code shall be used to
analyse public open access data of model species to evaluate whether the method
works for other data-sets as well.

- **Vector Space Analysis**
Vector space analysis is a computational approach in bioinformatics that evaluates gene expression diversity across conditions, such as tissues, species, or developmental stages. By representing gene expression levels as vectors in a multidimensional space, this method helps uncover genetic mechanisms that drive biological structures and functions.


- **Goals**

- Understand and reproduce the original R code and its results from the 
Cardamine hirsuta genome project

- Create EasyVectorOmics, a new R package that generalizes the vector space analysis methods 
to work with any dataset, structured as a pipeline of executable scripts

- [DOI: https://doi.org/10.1038/nplants.2016.167](https://doi.org/10.1038/nplants.2016.167)


- **Research Question / Hypothesis**

> How does the mean distance between orthologous genes compare to that of paralogous genes, 
> and what implications does this have for understanding evolutionary relationships and 
> functional divergence in gene families?
  

# ----------------------------------------------------------------------------------

 # References

  [1] Gan, X., Hay, A., Kwantes, M., Haberer, G., Hallab, A., Ioio, R. D.,
    Hofhuis, H., Pieper, B., Cartolano, M., Neumann, U., Nikolov, L. A., Song,
    B., Hajheidari, M., Briskine, R., Kougioumoutzi, E., Vlad, D., Broholm, S.,
    Hein, J., Meksem, K., … Tsiantis, M. (2016). The Cardamine hirsuta genome
    offers insight into the evolution of morphological diversity. Nature Plants,
    2, 16167. https://doi.org/10.1038/nplants.2016.167


- See the namesake section in Vignette `./vignettes/GeneFamilies.Rmd` for
details. Basically, here the angles between mean expression vectors of two
sub-groups within the same super-group of genes are explored. For example, for
each gene family the angle of rotation around the diagonal between the mean
orthologs and mean paralogs expression. See slides sixteen (16) to nineteen
(19).
