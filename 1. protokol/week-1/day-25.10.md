

**Date**: 25.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

1. **Plot Script Improvements**:
- imroved plot scripts, including `plot_exp.prof.dists_distributions.R` and `plot_expression_angles.R`.
- Integrated t-tests directly into these plotting scripts to streamline statistical analysis alongside visualization.

2. **Script for Expression Vector Space**:
- created `defining_exp.vector.space_tissue.R` to define the expression vector space and select specific tissues for analysis.

3. **Plot and T-Test Validation**:
- Ran trial plots and t-tests with mixed data to test plot functionality and appearance, - located in the experiments/RPKM_flybase/results folder.

**Next Steps**:
- would need a naming convention for files across the package to maintain consistent and standardized names.
- Define input/output file formats, including data structure and required formats, to ensure compatibility throughout the pipeline.

---

**Code**

- Code for the t-tests

```R
t_test_median <- df_median.dists %>%
  t_test(Distance ~ Type, alternative = "greater") %>% 
  adjust_pvalue(method = "BH") %>%                    
  mutate(significance = sapply(p, significance_level),
         analysis = "Median")  

t_test_mean <- df_median.dists %>%
  t_test(Distance ~ Type, alternative = "greater") %>% 
  adjust_pvalue(method = "BH") %>%                    
  mutate(significance = sapply(p, significance_level),
         analysis = "Mean")  

# Combine the two results into one summary dataframe
t_test_summary <- bind_rows(t_test_median, t_test_mean)

write.csv(t_test_summary, file.path(results_dir, "t_test_summary.csv"), row.names = FALSE)
```

