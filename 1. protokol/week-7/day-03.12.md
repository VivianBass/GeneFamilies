
**Date**: 03.12.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Completed Tasks**:

1. Statistical Test Consolidation
   - Removed redundant statistical test code
   - Unified test functions (t-test/Wilcoxon) into single implementation
   - Integrated statistical tests directly into plotting scripts

2. Visualization Enhancements
   - Removed redundant code and organized functions into dedicated modules 
   - Created 4-in-1 boxplot layouts (regular euclidean/log2 with t-test/Wilcox)
   - Added gene group count annotations for all plots
   - Implemented mean/median value annotations
   - Added horizontal reference line showing overall mean

3. Affected Files
   - exec/plot_exp.prof.dists_angles_distributions.R
   - exec/plot_exp.prof.dists_euclidean_distributions.R
   - exec/plot_exp.prof.dists_euclidean_tissue_distributions.R

**Open Questions**:

1. Validation needed for mean/median annotations in plots. To reassure they are right
2. Clarification needed: Should angle-to-diagonal calculations be applied to:
   - Euclidean distances
   - Log2 distances
   - Angle distances

**Next Steps**:

1. Documentation
   - Enhance roxygen2 documentation for all functions
   - Improve lab book documentation
   - Document SLURM scripts and tool usage

2. Analysis Extensions
   - Implement angle-to-diagonal statistics
   - Add angle-to-diagonal calculations for log2 data
   - Add angle-to-diagonal calculations for angle distances

3. Presentation Preparation
   - Outline key sections: intro, methods, results, discussion
   - Select and prepare figures
   - Create bullet-point overview for thesis structure


---

**Code:**

- The functions for boxplots and statistical tests are also organized into dedicated files.

R/
├── plot_distribution_funks.R         # Core plotting functions for distributions
├── plot_distribution_tissue_funks.R  # Tissue-specific plotting functions
├── statistical_tests_funks.R         # General statistical test implementations
└── statistical_tests_tissue_funks.R  # Tissue-specific statistical analyses