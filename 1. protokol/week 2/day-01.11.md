
**Date**:  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- check statistics rscript with new input data , 5 gene groups
- check distances rscript with new input data , 5 gene groups
- see how to change dists rscript for logarythmic data


- weekly report

- unit tests for the functions using testthat package
- All the functions inside your load_data_funks.R should have unit tests.
- for the functions load_data_frame() & create_nested_list() sourced from: source("R/load_data_funks.R")

- created doc\testing_loading_functions.md for testing scenarios for the loading functions 

- need to do roxygen2 on the functions


- VPN, sophos, work on server



**Doubts and Issues**:

- problem in computing distances with compute_exp.prof.dists.R ,  because gene-id in expression profils is FBgn0000003 and we are using FBpp proteinsequences for the distances
-> so we would need eiter to get the expressionprofils by proteinsequences ids or map the proteinsequences ids first back to the corresponding gene sequences
-> maybe a dummy dataset for trying
-> do we have the mapped data already `?

- for now i will consider the FBgn names for calculating distances

- need data with protein sequences rather then gene sequences for further testing,
the scriprts are already adjusted for gene sequences like FBgn names
- need data for FBpp names 

- are we working with protein sequences or gene sequences, for expression profils ect . or interchangeable?


**Next Steps**:


---

**Code:**