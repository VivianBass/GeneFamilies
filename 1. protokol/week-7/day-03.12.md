**Date**: 03.12.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Completed Tasks**:

1. **Statistical Test Consolidation**:

   - Removed redundant statistical test code.
   - Unified test functions (t-test/Wilcoxon) into a single implementation.
   - Integrated statistical tests directly into plotting scripts.

2. **Visualization Enhancements**:

   - Removed redundant code and organized functions into dedicated modules.
   - Created 4-in-1 boxplot layouts (regular Euclidean/log2 with t-test/Wilcox).
   - Added gene group count annotations for all plots.
   - Added mean lines and values for each box in the boxplots.
   - Added horizontal reference line showing overall mean. of the 5 groups in the plots
   - Optimized the annotations in the plots and the code to handle them.


3. **Affected Files**:

   - `exec/plot_exp.prof.dists_angles_distributions.R`
   - `exec/plot_exp.prof.dists_euclidean_distributions.R`
   - `exec/plot_exp.prof.dists_euclidean_tissue_distributions.R`

---

**Open Questions**:

1. Validation needed for mean/median annotations in plots to ensure accuracy.
2. Clarification needed: Should angle-to-diagonal calculations be applied to:

   - Euclidean distances?
   - Log2 distances?
   - Angle distances?

---

**Next Steps**:

   - Implement angle-to-diagonal statistics.
   - Add angle-to-diagonal calculations for log2 data.
   - Add angle-to-diagonal calculations for angle distances.

---

**Code**:

- The functions for boxplots and statistical tests are also organized into dedicated files.

     ```
     R/
     ├── plot_distribution_funks.R         # General plotting utilities
     ├── plot_distribution_tissue_funks.R  # Tissue-specific plotting
     ├── statistical_tests_funks.R         # Core statistical functions
     └── statistical_tests_tissue_funks.R  # Tissue-specific statistics
     ```

