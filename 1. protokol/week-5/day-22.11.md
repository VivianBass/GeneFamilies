
**Date**:  22.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`  

---

**Tasks**:  

- **New Data Integration**:  
  - Incorporated new data from raw reads (`gffread`, `Trimmomatic`, `kallisto` workflow, etc.).  
  - Recreated results using the new data and provisional TSVs (expressions for *D. melanogaster* and *D. sechellia*).  
  - Set up a new experiment directory and regenerated results and plots.  

- **Statistical Validation**:  
  - Verified t-test results and their alignment with significance levels in the plots.  
  - Resolved issues with t-tests and Wilcoxon tests.  

- **Gene Group Analysis**:  
  - Confirmed distance calculations for each gene group are calculated separately.  

---

**Doubts and Issues**:  

- Validation of calculations: Are angles, distances, etc., computed correctly?  
- Clarification on plots:  
  - What specific plots are required?  
  - Should gene families be included in the distribution plots, beyond gene group distributions?  
- so t test for means and wilcox for medians ?? (when exactly to use t-test and when wilcox ?)

---

**Next Steps**:  

- Refine Lab Book descriptions for the workflow’s four steps and provide detailed documentation for each step.  

