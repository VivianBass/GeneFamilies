


# ------------------------------------------------------------------------------------------------
### **4.1 Purpose and Overview**
# ------------------------------------------------------------------------------------------------

#### **Purpose**

The **EasyVectorOmics** R-package facilitates bioinformatics analysis by enabling streamlined data handling, statistical computation, and visualization for comparative gene expression studies. It plays a crucial role in transitioning from data preprocessing (as discussed in Chapter 3) to results interpretation and visualization, setting up the foundation for the results (Chapter 5) and discussion (Chapter 6).

The package is designed to:
- Load and normalize gene expression data.
- Compute pairwise distances in gene expression profiles for comparative analysis.
- Perform statistical tests to identify significant differences in gene expression patterns.
- Generate high-quality visualizations for data interpretation and reporting.

#### **Overview** 

Key features of **EasyVectorOmics** include:
- Comprehensive data loading capabilities for expression profiles, gene groups, and gene families.
- Advanced distance computation, including Euclidean distances and cosine angles.
- Statistical analysis, including t-tests, Wilcoxon tests, and descriptive statistics.
- Versatile plotting functions for distribution visualizations and annotated statistical results.


#### **Core Functions** 
# ------------------------------------------------------------------------------------------------

1. **Data Loading**  
   Scripts like `load_gene_expression_data.R`, `load_gene_groups_data.R`, and `load_gene_families_data.R` handle preprocessing and loading of raw data into R-friendly structures.

2. **Distance Computation**  
   - Pairwise **Euclidean Distances** within gene groups and tissues.
   - **Cosine Angles** for similarity measures.
   - Angle-to-diagonal metrics for assessing data alignment.

3. **Statistical Analysis**  
   - Descriptive statistics (mean, median, range) of expression distances.
   - Comparative tests (t-tests and Wilcoxon tests) for significant differences between orthologs and paralogs.

4. **Visualization**  
   - Distribution plots of gene expression distances.
   - Annotated boxplots summarizing statistical results.
   - Tissue-specific plots for versatile data exploration.



# ------------------------------------------------------------------------------------------------
### **4.4 Detailed Breakdown of Analysis Workflow**
# ------------------------------------------------------------------------------------------------

The **EasyVectorOmics** package serves as the computational core of the broader bioinformatics pipeline, facilitating:
1. **Automated Analysis**: Standardized scripts for reproducible computations.
2. **Exploratory Data Analysis**: Interactive plots and summaries to uncover patterns in gene expression.
3. **Cross-Species Comparisons**: Insights into functional diversification through orthologous and paralogous relationships.

#### **Section 1: Data Loading**

- **Scripts**:
  - `load_gene_expression_data.R`: Handles RPKM data loading.
  - `load_gene_groups_data.R`: Imports ortholog and paralog relationships.
  - `load_gene_families_data.R`: Processes MCL output for gene family clustering.
- **Normalization**: Adjusts raw expression data using RPKM methodology to ensure comparability.

#### **Section 2: Distance Computations**

- **Scripts**:
  - `compute_exp.prof.dists.R`: Calculates Euclidean distances.
  - `compute_exp.prof.dists_statistics.R`: Summarizes distance statistics.
  - `compute_exp.prof.dists_angles.R`: Computes cosine angles for expression vectors.
- **Output**: Results stored in R objects for further analysis (e.g., `ExpressionProfileDistances.RData`).

#### **Section 3: Visualization and Statistical Computations**

- **Scripts**:
  - `plot_exp.prof.dists_distribution.R`: Generates distribution plots.
  - `plot_exp.prof.dists_distribution_tissue.R`: Tissue-specific distance visualizations.
  - `plot_exp.prof.dists_angles.R`: Visualizes cosine angle distributions.

- **Statistical Tests**:
  - Parametric t-tests and non-parametric Wilcoxon tests compare orthologous and paralogous gene expression distributions.
  - P-values adjusted using Benjamini-Hochberg correction.