
**Date**: 28.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

### **Tasks**:

- Determine how to plot complete distances without calculating mean and median; create plot scripts for "all distances."
- Plan how to conduct t-tests and Wilcoxon tests for the complete distances.

**Data Loaded**:
```R
load(file.path(output_data_dir, "exp.prof.dists.RData"))
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))

load(file.path(output_data_dir, "exp.prof.angles.RData"))
load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))
```

Loaded datasets:

- **Distance Data**:  
  - `con_orthologs_v.lst_dists`, `con_orthologs_v.lst_dists_tissue`
  - `in_paralogs_v.lst_dists`, `in_paralogs_v.lst_dists_tissue`
  - `out_paralogs_v.lst_dists`, `out_paralogs_v.lst_dists_tissue`
  - `special_in_paralogs_v.lst_dists`, `special_in_paralogs_v.lst_dists_tissue`
  - `special_out_paralogs_v.lst_dists`, `special_out_paralogs_v.lst_dists_tissue`

- **Log2-Transformed Distance Data**:  
  - `con_orthologs_v.lst_dists_log2`, `con_orthologs_v.lst_dists_tissue_log2`
  - `in_paralogs_v.lst_dists_log2`, `in_paralogs_v.lst_dists_tissue_log2`
  - `out_paralogs_v.lst_dists_log2`, `out_paralogs_v.lst_dists_tissue_log2`
  - `special_in_paralogs_v.lst_dists_log2`, `special_in_paralogs_v.lst_dists_tissue_log2`
  - `special_out_paralogs_v.lst_dists_log2`, `special_out_paralogs_v.lst_dists_tissue_log2`

- **Cosine Angles Data**:  
  - `con_orthologs_v.lst_cos_angles_dists`
  - `in_paralogs_v.lst_cos_angles_dists`
  - `out_paralogs_v.lst_cos_angles_dists`
  - `special_in_paralogs_v.lst_cos_angles_dists`
  - `special_out_paralogs_v.lst_cos_angles_dists`

- **Log2-Transformed Cosine Angles Data**:  
  - `con_orthologs_v.lst_cos_angles_dists_log2`
  - `in_paralogs_v.lst_cos_angles_dists_log2`
  - `out_paralogs_v.lst_cos_angles_dists_log2`
  - `special_in_paralogs_v.lst_cos_angles_dists_log2`
  - `special_out_paralogs_v.lst_cos_angles_dists_log2`

---

### **Key Tasks**:

- **Plot Complete Distances**:
  - We need to create plots for **all distances** (without calculating mean or median).
  - The script should handle plotting of complete distance matrices for both **Euclidean** and **Angular** distance measures, with and without log transformation.

- **t-Test and Wilcoxon Test**:
  - Conduct t-tests and Wilcoxon tests on the **complete distances**.
  - Separate plots should be created for **t-test** and **Wilcoxon test** results, as combining both annotations in a single plot is not feasible.

- **Distance Methods**:
  - The function should be flexible to handle both **Euclidean distance** and **Angular (Cosine) distance** measures.
  - Both **raw** and **log2-transformed** distance measures should be included for the tests.
  
  **Required tests**:
  - **All distances (Euclidean)**
  - **Mean (Euclidean)**
  - **Median (Euclidean)**
  - **All distances (Euclidean) with Log2**
  - **Mean (Euclidean) with Log2**
  - **Median (Euclidean) with Log2**
  - **All distances (Angular)**
  - **Mean (Angular)**
  - **Median (Angular)**
  - **All distances (Angular) with Log2**
  - **Mean (Angular) with Log2**
  - **Median (Angular) with Log2**

- **Separate Plots for t-test and Wilcoxon test**:
  - Generate individual plots for t-tests and Wilcoxon tests for each experiment (i.e., for mean, median, and all distances).

---

### **Next Steps**:



---

**Doubts and Issues**:


