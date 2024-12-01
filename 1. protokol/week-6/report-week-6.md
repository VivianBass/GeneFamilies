### **Weekly Protocol (25.11.2024 - 29.11.2024)**

**Git Branch**: `Gene-Families-tests-Andre`

---

#### **Overview**
The week focused on validating statistical tests, refining distance methods (Euclidean and Angular), and analyzing gene expression data across multiple species. New experiments were conducted with filtered orthogroups, and various distance calculations were applied to assess gene expression patterns.

#### **Key Tasks & Progress**
- **Validation of Tests**: Verified t-tests, added Wilcoxon test calculations, and switched to "two.sided" method for better statistical accuracy.
- **Distance Methods**: Explored Euclidean and angular distances (including log transformations) for gene expression data.
- **Expression Profile Analysis**: Normalized RNA-seq data and performed distance calculations using raw and log-transformed values.
- **Plotting & Statistical Testing**: Created plots for t-tests and Wilcoxon tests, visualizing significance levels for different distance methods.
- **Transcript Experiments**: Conducted experiments with three species (dmel, dsec, dsim) and filtered orthogroups, focusing on gene families with at least 6 gene pairs.

#### **Doubts & Issues**
   - How should distance methods be applied to gene families?  
   - Clarification needed for interpreting cosine angular distance results, as the range seems limited to [0, 1.570796].

#### **Next Steps**
   - Refine and document the workflow's four steps.
   - Add background information on the experiments, data availability, and context.
   - Streamline the code by removing redundancies.
   - Provide a detailed report of the data, including tissue types, species, experimental conditions, and gene family counts.
   - Further investigate the use of vector space analysis and expression ratio calculations.
   - Calculate fold changes for experimental conditions (e.g., diet vs. control, tissue vs. whole body).
   - Continue testing and plotting results for both Euclidean and angular distances.