
**Date**:  06.11.2024 
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:


- Refactor `compute_exp.prof.dists.R` and `compute_exp.prof.dists_statistics.R` to reduce redundancy, 
  simplify code, and improve readability. Potentially create a reusable function.

- Update descriptions in `roxygen2` documentation for clarity. And uptdate the Documentation overall

- as i had to low data from the current expression profils, try to use other expression profils
  but with the same mappping as done with the fasta file from the server FBgn to FBpp mapping
- try with different expression profils and see if the inteersection and data amount differs
- also tzry to filter the genegroups dataframes just by 1 column first and see if difference

- continue to do unit tests for each function , create different scenarious to test for 

- send Asis the weekly report per email

- would have to add the families also for distance calculation !?

- did a v1 and a v2 filter to preprocess data , sothat it matches the provided expression profiles;
- check if v1 filetered has the same results/lenght as v2 filtered


**Doubts and Issues**:


**Next Steps**:


---

**Code:**