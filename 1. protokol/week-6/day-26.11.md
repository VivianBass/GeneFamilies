
**Date**: 26.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

### **Tasks**:

1. **Validation of Tests**:
   - Verified the correctness of t-test calculations using reference DataFrames.
   - Added Wilcoxon test calculations and switched from the "greater" method to "two.sided":
     - **"two.sided"**: Tests if the means differ in either direction (currently used).
     - **"less"**: Tests if the mean of \(x\) is less than \(y\).
     - **"greater"**: Tests if the mean of \(x\) is greater than \(y\).

2. **Angle Calculations**:
   - Incorporated t-test and Wilcoxon tests for angular distance calculations.
   - Created separate plots for t-test and Wilcoxon annotations.

3. **Significance Levels**:
   - Added significance levels directly in the `perform_t_test` function:
     ```R
     p.adj.signif = case_when(
         p >= 0.05 ~ "ns",
         p < 0.001 ~ "***",
         p < 0.01 ~ "**",
         p < 0.05 ~ "*"
     )
     ```

4. **Input Order**:
   - Verified that changing the input order does not affect the results.

5. **Transcript Experiments**:
   - Conducted experiments with complete transcript values and the largest transcript per gene.
   - Tested two new directories with TPM values, focusing on three species (dmel, dsec, dsim):
     - `experiments/test_diet_P_&_M_new_3_species_tpm_all`
     - `experiments/test_diet_P_&_M_new_3_species_tpm_largest_prot`

---

**Execution Examples**:

1. **Gene Expression Data**:
```R
   Rscript exec/load_gene_expression_data.R experiments/test_diet_P_&_M_new_3_species_tpm_all/3sp_tpm_all.tsv
```

2. **Gene Groups Data**:
```R
   Rscript exec/load_gene_groups_data.R \
       experiments/test_diet_P_&_M_new_3_species_tpm_all/dmel_dsec_dsim_conserved_orthologs.tsv \
       experiments/test_diet_P_&_M_new_3_species_tpm_all/dmel_dsec_dsim_in_paralogs.tsv \
       experiments/test_diet_P_&_M_new_3_species_tpm_all/dmel_dsec_dsim_out_paralogs.tsv \
       experiments/test_diet_P_&_M_new_3_species_tpm_all/dmel_dsec_dsim_special_in_paralogs.tsv \
       experiments/test_diet_P_&_M_new_3_species_tpm_all/dmel_dsec_dsim_special_out_paralogs.tsv
```

3. **Gene Families Data**:
```R
   Rscript exec/load_gene_families_data.R experiments/test_diet_P_&_M_new_3_species_tpm_all/orthogroups.tsv
```

---

### **Doubts/Issues**:

- 

---

### **Next Steps**:

- 

---

### **Code Overview**:

- **General Tests**:  
  Located in `R/compute_funks.R`. Performs t-tests and Wilcoxon tests for analysis types, adjusts p-values, and incorporates significance levels.  
```R
perform_tests <- function(data, valid_groups, analysis_type) {
    if (length(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            t_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            wilcox_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}
```

- **Tissue-Specific Tests**:  
  Similar structure, specific to tissue-based data.  
```R
  perform_tissue_tests <- function(data, valid_groups, analysis_type) {
    if (nrow(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            wilcox_test(Distance ~ Type, alternative = "two.sided") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}
```
