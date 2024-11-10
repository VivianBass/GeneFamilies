
# EasyVectorOmics Package

In this project the goal is to reproduce the vector space analyses as done in
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


**Homology and Gene Relationships**

Homology refers to shared ancestry or functional similarity among genes or structures, categorized as follows:

1. **Historical Homology**: Traits inherited from a common ancestor, like a family heirloom.
2. **Proximal-Cause Homology**: Traits shaped by similar developmental mechanisms, akin to different recipes using shared ingredients.
3. **Factorial Homology**: Homology with layered relationships across biological levels, like a patchwork quilt.
4. **Context-Dependent Homology**: Definition varies by context (morphology, genetics), similar to words with multiple meanings.
5. **Continuity of Evolution**: Homology as a spectrum, where traits gradually evolve.

**Gene Types**:

- **Orthologs**: Genes in different species from a common ancestor, maintaining similar functions (e.g., human and mouse hemoglobin).
- **Paralogs**: Genes duplicated within a lineage, allowing for new functions. Includes:
  - **In-Paralogs**: Duplication within one species post-speciation (like specialized siblings).
  - **Out-Paralogs**: Duplication pre-speciation, resulting in related genes across species (like distant cousins).
  - **Special In/Out-Paralogs**: Unique to one species (in-paralogs) or cross-species (out-paralogs) with conserved functions, resembling orthologs due to evolutionary pressure.

**Key Differences**:
- **Orthologs**: Directly related genes across species, maintaining similar roles.
- **Paralogs**: Duplicated genes within or across species, often diverging in function. 
- **Special Out-Paralogs**: Cross-species paralogs that retain conserved functions, closely resembling orthologs but from ancient duplications.


# Algorythm, distinguishing between gene groups 
  
- key parameters the algorythm evaluates or considers:
- similratity score
- overlapp

- phylogenic trees, newick trees (for example from Orthofinder)

- dublication or speciation events

- Flybase syntenic relationchips

- Cactus Bioinformatics tool

- harmonic means on similarity scores


- could we basically also apply all vs all ?? 
- the protein genes can belong to multiple groups (from 5 groups + tandems)
- startting from a specific leaf/node

- distribution of gene pairs, and overall gene pairs 

# Harmonic mean

The harmonic mean is a statistical measure used to average rates or ratios, giving more weight to smaller values, making it suitable for data where each point's contribution is part of a whole.
In bioinformatics, it is applied in sequence alignment scoring, estimating effective population size, assessing microbial diversity, calculating the F1 score in machine learning, and averaging evolutionary or mutation rates.
It is preferred over the arithmetic mean in contexts where low values significantly impact the overall result, providing a more balanced and realistic estimate.

**Expression-Vector-Space Analysis Summary**

This analysis evaluates gene expression diversity across conditions (e.g., tissues, species, developmental stages) to uncover genetic mechanisms that shape biological structure. The two-factor dataset (ST-Exp) contains gene expression data with *species* and *tissue* as primary factors, requiring at least four family members per species for robust comparisons. Expression distances are calculated for gene types (orthologs, paralogs, etc.), with species-specific distribution plots and statistical tests to assess differences. An ANOVA-based conservation analysis further examines whether orthologs retain expression levels more consistently than paralogs, addressing the key question: *Is there a difference in expression diversity between orthologs and paralogs?*

In the new EasyVectorOmics we will need from the user:

The steps and data required to carry out the original analysis are devided into 3 Sections 

- `1. Loading Data`                                                             (Section-1) 
- `2. Computing Distances, Statistics, T-tests`                                 (Section-2) 
- `3. Plotting Distributions`                                                   (Section-3) 

## Section 1 - Data Loading

1. **`load_gene_expression_data.R`**  
2. **`load_gene_groups_data.R`**  
3. **`load_gene_families_data.R`**  

### Tools for Identifying Gene Families and Homologs

- **BLAST** / **Diamond**: 
  Identify similar regions in DNA, RNA, or protein sequences to find homologous genes (orthologs, paralogs).
- **Markov Clustering (MCL)**: 
  Groups genes into clusters (gene families) using similarity matrices from BLAST results.

### 1. `load_gene_expression_data.R`

Loads expression data for gene groups, typically in DNA FASTA format. Expression levels (e.g., RNA counts, RPKM) reflect gene activity: higher counts indicate active genes (ON), while lower counts indicate inactive genes (OFF). Profiles are multi-dimensional vectors representing gene activity across tissues. Normalization methods (e.g., TPM, FPKM) adjust for technical variations, ensuring accurate comparisons between samples.

### 2. `load_gene_groups_data.R`

Loads ortholog and paralog relationships across species, with files categorized by group type. Each file follows a specific header convention:

- **Conserved Orthologs**: 
  One ortholog per species-gene pair; represents highly conserved, reliable relationships.
- **In Paralogs**: 
  Species-specific paralogs from recent duplications with no conserved ortholog.
- **Special In Paralogs**: 
  Species-specific paralogs with a conserved ortholog relationship, indicating potential conservation.
- **Out Paralogs**: 
  Paralogs from older duplications, not tied to a conserved ortholog.
- **Special Out Paralogs**: 
  Cross-species paralogs with a conserved ortholog relationship, suggesting retained functional ties.

**File Naming Convention**:
- Orthologs Header: `Family | Gene | Gene_species | Ortholog | Ortholog_species`
- Paralogs Header:  `Family | Gene | Gene_species | Paralog  | Paralog_species`

### 3. `load_gene_families_data.R`

Loads gene family clusters created by MCL. Each family clusters related genes based on sequence similarity, organizing them by shared evolutionary lineage.


## Section-2 - Computing Distances, Statistics, T-tests 
			                                      
- 1.  `compute_exp.prof.dists.R` 	                        
- 2.  `compute_exp.prof.dists_statistics.R`

### 1. `compute_exp.prof.dists.R` 

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


### 2. `compute_exp.prof.dists_statistics.R`

**Statistical Analysis of Expression Profile Distances**:

**Purpose**: Investigate distributions of Euclidean distances within gene-family clusters.

The R script summarizes distances per gene group (Orthologs, Paralogs) with key statistics, 
including mean, median, max, and max-minus-min.

For each gene family and subset (Orthologs, Paralogs), the mean and median Euclidean distances 
are calculated to assess similarities in gene expression levels.


## Section-3 - Plotting Distributions

- 1. `plot_exp.prof.dists_distribution.R`
- 2. `plot_exp.prof.dists_distribution_tissue.R`




- **Scientific Plots Overview**:

**Boxplots of Expression Distances**:
**First Plot**: ...
**Second Plot**: ...

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