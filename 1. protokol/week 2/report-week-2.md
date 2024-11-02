

**Weekly Summary (28.10.2024 - 30.10.2024):**

**Git Branch**: `Gene-Families-tests-Andre`

---

### Overview
This week focused on updating loading data scripts and processing scripts to handle gene-pair data structures and accommodate five distinct gene groups in a pairwise analysis format.

### Key Tasks & Progress

1. **Script Updates for Pairwise Data**:
- Modified scripts to load 5 gene groups (one conserved orthologs group and 4 distingt paralogs gourps) and enable **pairwise data  processing**.
   

2. **Gene Group Classification**:
- Established and implemented five distinct gene groups:
    - **Conserved Orthologs**
    - **In-Paralogs with Orthologs**
    - **In-Paralogs without Orthologs**
    - **Out-Paralogs with Orthologs**
    - **Out-Paralogs without Orthologs**
- Defined headers and structure for each group file.

3. **Data Format and Loading Functions**:
- Refined data structure for nested gene-pair information, aligning with   
  OrthoFinder output format. basically how to display the data in nested. 
- Created reusable loading functions in `load_data_funks.R` and adjusted the   
  main scripts to source from this module. Containing 2 functions load_data_frame() & create_nested_list()

4. **Documentation & Testing**:
- Updated documentation to specify the format for gene group input files.
- Created dummy datasets for testing and organized these in `experiments/dummy-datasets`.
- Documented unit testing procedures using the `testthat` package.

5. **Distance Calculations**:
- Initiated adjustments for log-transformed distance calculations across the five gene groups.

### Doubts & Issues
- Awaiting detailed distinctions among gene groups for more accurate classification.

### Next Steps
- Finalize distance and statistical computations for the five gene groups.
- Continue refining input formats and adjusting for log-transformed 
  calculations in `compute_exp.prof.dists.R`.
- Continue implementing Unit Tests for the loading functions.

