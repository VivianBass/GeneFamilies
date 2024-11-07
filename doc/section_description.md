
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


# Homologues

**Types of Homology:**

- **Historical Homology**:  
  - Defined by shared ancestry; structures are homologous if inherited from a common ancestor.
  - *Metaphor*: A family watch passed down, used differently by each sibling but with shared origins.

- **Proximal-Cause Homology**:  
  - Focuses on developmental processes and similar genetic mechanisms that shape homologous traits.
  - *Metaphor*: Different recipes (pizza vs. flatbread) with common ingredients result in similar dishes.

- **Factorial/Combinatorial Homology**:  
  - Homology is seen as a layered, relative concept with partial relationships across biological levels.
  - *Metaphor*: A patchwork quilt with some squares identical by origin, while others are added later.

- **Context-Dependent Homology**:  
  - Homology definitions vary based on perspective (morphological, genetic, developmental).
  - *Metaphor*: The word “bank,” which changes meaning based on context (e.g., financial vs. riverbank).

- **Continuity of Evolution**:  
  - Homology as a spectrum, where traits evolve gradually rather than fitting binary categories.
  - *Metaphor*: Evolution is like a road trip with no single dividing line between start and destination.



---

**Orthologs**: (similar genes across different species, from common ancestor)
- Orthologs are genes in different species that evolved from a common ancestral gene, 
  typically retaining similar functions across species.
- Example: Hemoglobin genes in humans and mice are orthologs, originating from a shared ancestor.
- *Metaphor*: Cousins inheriting a distinctive family trait (like eye color) from a shared grandparent, even though they’re in different families.

**Paralogs Overview:**
- Paralogs are genes that originate from duplication within a species, 
  allowing them to evolve new functions over time.
- Types include in-paralogs (pre-speciation, unique to one species) 
- and out-paralogs (past-speciation, found in multiple species).
- Broad category of genes that arise from duplication events within a species or across species.
- Can evolve new or modified functions over time due to gene duplication within the same lineage.

---




**Conserved Orthologs**: (like Monet Picture)
- A single version of a gene maintained across species with minimal change due to its essential role.
- *Metaphor*: A family heirloom, like a unique painting passed down unchanged to one family member in each generation, symbolizing continuity across lineages.

**In-Paralogs**: (like gene siblings in a family)
- Result from gene duplication within a single species after it diverged from others, often adapting to distinct functions within that species.
- *Metaphor*: Siblings in a family business with specialized roles, sharing a close origin but taking on different responsibilities.

**Out-Paralogs**: (like gene cousins in a family)
- Arise from gene duplication events in a common ancestor before species split, creating related genes across different species with unique roles.
- *Metaphor*: Distant cousins whose family branches split generations ago, sharing an ancestral link but evolving separately.

**Special In-Paralogs Overview:**
- Species-specific duplicates with a conserved ortholog relationship, indicating that while they are unique to one species, they retain a function linked to an original ortholog across species.
- *Metaphor*: Identical cousins—distinct to one family branch but closely resembling relatives from other branches, reflecting a shared ancestral trait.

**Special Out-Paralogs Overview:**
- because out-Paralogs could be mistaken for Orthologs
- Cross-species paralogs that retain a conserved ortholog function, suggesting evolutionary pressure to preserve a similar role across species despite divergence.
- *Metaphor*: A foundational recipe adapted across different cultures, maintaining a core purpose and ingredients, showing its original, shared importance.

- difference
- **Orthologs**: Genes in different species that evolved from a single ancestral gene after speciation, retaining similar functions and representing a direct lineage.
- **Special Out-Paralogs**: Genes duplicated before speciation, existing across species and retaining conserved functions similar to orthologs despite originating from an ancient duplication.

**Key Differences**:
   - **Paralogs** include any duplicated genes (regardless of species or time).
   - **In-paralogs** are paralogs specific to one species, emerging post-speciation.
   - **Special out-paralogs** are cross-species paralogs with conserved functions, resembling orthologs but originating from a pre-speciation duplication event.



# Expression-Vector-Space-Analysis 

- This process measures and compares `gene expression levels` under different conditions, 
  such as across tissues, developmental stages, or environmental contexts. 
  Differential expression analysis provides insights into the genetic mechanisms 
  underlying leaf shape and structure, identifying key genes involved in leaf development.


two-factor data (ST-Exp) should contain and its purpose:

1. **Species and Tissue as Primary Factors:** ST-Exp should include gene expression data across multiple species and tissues, treating both species and tissue types as two independent factors.

2. **Minimum Family Members:** Each family group must have at least four members across different species to allow for meaningful comparisons.

3. **Gene Duplication Types:** Measure expression distances for orthologs, paralogs, 5 groups

4. **Expression Distance Metrics:** Calculate expression distances per gene type (ortholog, paralog, etc.) across species, with one distribution plot per species (boxplot/raincloud) and pairwise statistical testing (t-tests, Wilcoxon tests).

5. **ANOVA for Conservation Analysis:** Compare expression conservation across species by measuring distances in a tissue-defined vector space, assessing if orthologs retain expression more than paralogs.


- main question: Is there a difference in expression diversity between Orthologs and Paralogs?

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

### 3. `compute_t-tests_wilcox-tests.R`

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