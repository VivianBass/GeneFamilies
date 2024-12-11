**Date**: 04.12.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Completed Tasks**:

1. **Documentation**:

   - Completed roxygen2 documentation for each function.

2. **Visualization Enhancements**:

   - Improved annotations and adjusted code to handle them more effectively.
   - Added horizontal lines representing mean values for each boxplot:
      - Lines are drawn within the box, with the mean value displayed to the left of the line.
      - Included an overall mean of the means of the 5 groups, represented by a dotted horizontal line.

3. **Statistical Analyses**:

   - Implemented angle-to-diagonal statistics (Wilcoxon test, t-test).
   - Generated angles-to-diagonal data for log2-transformed values.

4. **Cosine Angle Distances**:

   - Calculated both radians and degrees for cosine angle distances.
   - Plots now include both angles and degrees.

     ```R
     degrees = radians * (180 / math.pi)
     ```

     $$\text{degrees} = \text{radians} \times \frac{180^{\circ}}{\pi}$$

---

**Doubts and Issues**:


---

**Next Steps**:



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

