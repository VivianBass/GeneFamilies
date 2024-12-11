



#### **Section 3: Visualization and Statistical Significance**
# ------------------------------------------------------------------------------------------------  

- **Scripts**:
  - `plot_exp.prof.dists_distribution.R`: Generates distribution plots.
  - `plot_exp.prof.dists_distribution_tissue.R`: Tissue-specific distance visualizations.
  - `plot_exp.prof.dists_angles.R`: Visualizes cosine angle distributions.

- **Statistical Tests**:
  - Parametric t-tests and non-parametric Wilcoxon tests compare orthologous and paralogous gene expression distributions.
  - P-values adjusted using Benjamini-Hochberg correction.


**Statistical Analysis of Expression Profile Distances**:

**Purpose**: Investigate distributions of Euclidean distances within gene-family clusters.

The R script summarizes distances per gene group (Orthologs, Paralogs) with key statistics, 
including mean, median, max, and max-minus-min.

For each gene family and subset (Orthologs, Paralogs), the mean and median Euclidean distances 
are calculated to assess similarities in gene expression levels.


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


#### **4.5 Outputs and Applications**

- **Tables and Statistics**:
  - Summaries of mean/median distances and statistical test results.
- **Plots**:
  - High-resolution boxplots, distribution histograms, and annotated results for publication-ready figures.
- **Functional Insights**:
  - Identification of functional divergence in gene families.
  - Comparison of evolutionary trends in gene expression among orthologs and paralogs.

**Scientific Plots Overview**:

**Boxplots of Expression Distances**:
**First Plot**: ...
**Second Plot**: ...