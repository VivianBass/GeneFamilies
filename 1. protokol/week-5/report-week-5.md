
### **Weekly Protocol (18.11.2024 – 22.11.2024)**  

**Git Branch**: `Gene-Families-tests-Andre`  

---

#### **Overview**  
This week focused on integrating new data, debugging R scripts, conducting statistical analyses, and validating methods for EasyVectorOmics analysis. Work included the processing of *D. melanogaster* and *D. sechellia* datasets, documentation updates, and plotting results.  

---

#### **Key Tasks & Progress**  

- **Data Integration & Preparation**:  
  - Incorporated new data from raw reads, including workflows with `gffread`, `Trimmomatic`, and `Kallisto`.  
  - Recreated results using provisional TSVs for expression data of *Dmel* and *Dsec*.  
  - Established a new experiment directory and regenerated results and plots.  

- **Analysis & Validation**:  
  - Debugged R scripts, ensuring accurate calculation of angles and distances for gene groups.  
  - Verified t-test and Wilcoxon test results, aligning them with significance levels in plots.  
  - Ensured distance calculations for each gene group are performed separately.  

- **Statistical Scripts**:  
  - Enhanced functionality in `angles_funks.R` with two key functions:  
    - `calculate_angles`: Computes angles for gene groups.  
    - `validate_angle_dataframes`: Ensures dataframes are non-empty for tissue versatility plots.  
  - Created dedicated scripts for statistical tests:  
    - `generate_t-test_wilcox_test_tissue.R`  
    - `generate_t-test_wilcox_test.R`  

- **Documentation**:  
  - Continued refining Lab Book entries to document workflows and methods comprehensively.  
  - Documented unit test scenarios for angle functions in `unit_tests_angles_functions.md`.  

- **R Script Refinements**:  
  - Removed redundant code from `compute_exp.prof.dists_angles.R` by replacing repetitive sections with a function and loop.  
  - Added `roxygen2` documentation for improved function clarity.  

---

#### **Doubts & Issues**  

- **Validation**:  
  - Are calculations (e.g., angles, distances) implemented correctly?  
- **Plots**:  
  - What specific plots are required?  
  - Should gene families be included in angle distribution plots alongside gene groups?  

---

#### **Next Steps**  

- Continue refining Lab Book documentation for workflows and steps.  
- Address remaining doubts about calculations and required plots.  
 

