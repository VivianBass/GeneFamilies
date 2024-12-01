
**Date**: 29.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

### **Summary**:

- Conducted a new experiment with **filtered orthogroups** (at least 6 gene pairs per group) from 3 species, focusing only on distinct families/orthogroups.
- Created a new test directory: `experiments/test_diet_P_&_M_filtered`.
  
### **Examples**:

1. **Load Gene Expression Data**:
   ```R
   input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/3sp_tpm_all.tsv"
   ```

2. **Load Gene Groups Data**:
   ```R
   input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_con_orthologs.tsv"
   input.args[[2]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_in_paralogs.tsv"
   input.args[[3]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_out_paralogs.tsv"
   input.args[[4]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_special_in_paralogs.tsv"
   input.args[[5]] <- "experiments/test_diet_P_&_M_filtered/dmel_dsec_dsim_filtered_special_out_paralogs.tsv"
   ```

3. **Load Gene Families Data**:
   ```R
   input.args[[1]] <- "experiments/test_diet_P_&_M_filtered/orthogroups_filtered_more_than_5.tsv"
   ```



**Doubts and Issues**:


- how to accuratly interprete all those results ??

- cos_angles_dists are only in range [1] 0.000000 1.570796 is that right ? 

> all_values <- unlist(con_orthologs_v.lst_cos_angles_dists)
> range(all_values, na.rm = TRUE)
[1] 0.000000 1.570796


**Next Steps**:

- Refine and document the workflow’s four steps in detail.
- Add background information on the experiments, data availability, and the research context.
- Remove redundant code to streamline the process.
- Report data descriptions: include tissue types, species, experimental conditions, and the number of gene families used.
- Explain the vector space analysis for comparing gene expression.
- Use spreadsheets to calculate expression ratios, e.g., log (te/tc) = log(te) - log(tc), to assess gene fold changes.
- Specify the log base (2, 10, natural log) used for calculations.
- Analyze gene expression differences between paralogs and orthologs in vector space.
- Define a control group (referencing the relevant paper).
- Calculate fold changes for all experimental conditions: diet vs. control, tissue vs. whole body.


