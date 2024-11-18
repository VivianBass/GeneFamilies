
**Date**: 12.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

### **Tasks**  

- Worked on unit tests and documented test scenarios for each function in `compute_funks.R` and `load_data_funks.R`.  
  - Prepared multiple test cases per function.  
  - Created test scripts for each function using the `testthat` package.  
  - Set up a `tests` directory as required by `testthat`, following naming conventions for the test scripts.  
  - Documented test scenarios in the documentation folder and saved initial test results in `tests/test_results`.  
  - Example code used to run tests:  
    ```R
    sink("tests/test_results/test_results_exp.prof.dists.md")
    test_file("tests/testthat/test-exp.prof.dists.R")
    sink()
    ```  
  - Initial results can be found in `tests/test_results`.  

- Paused work on unit tests as per discussion with Vivian, to focus on gathering data for expression profiles this week.  
  - Thus current results are limited to *Drosophila melanogaster* (dmel) due to the lack of mapping data for other species.  

- Attended a meeting at 4 PM in TH-Bingen:  
  - **Expression Profiles**:  
    - FlyBase gene identifiers available only for dmel, limiting analysis.  
    - Need data from at least three (preferably four) species and four tissues to create a viable expression vector space.  

  - **Assignments**:  
    - Review papers on yeast and *Drosophila* in the shared Drive.  
    - Explore sources for obtaining expression data for other species.  

  - **Current Results**:  
    - Results are restricted to special cases (e.g., *in_paralogs* and *special_in_paralogs*) using dmel gene identifiers.  
    - Found in `experiments/test_vivian_1/results`.  

  - **Naming Conventions**:  
    - Finalized and standardized.  

---

### **Next Steps**  

- Review the first 25 papers on yeast for relevant data in Drive.  

