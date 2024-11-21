
**Date**:  21.11.2024
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Debug R scripts for angle calculations and plots. 
- added 2 more functions in R\angles_funks.R
- calculate_angles , function to calculate angles for gene groups
- and a validate_angle_dataframes  function to validate if the created angles, dataframes for each genegroup
are available or empty, this was important for the plotting process of tissue versatility etc. 
- did a basic roxygen2 documentation for the functions 
- also removed redundant code from exec\compute_exp.prof.dists_angles.R , to basically use a function and for loop instead of 
doing the same code for each gene group for angles calculation 


- Implemented Wilcoxon tests in the analysis. 
- also included wilcox tests and created 2 seperate rscripts just for t-tests and wilcox tests
- exec\generate_t-test_wilcox_test_tissue.R
- exec\generate_t-test_wilcox_test.R


- Questions i need to clarify 
- calculate the angles for the gene groups or what ??
- include gene families in plotting ?? distances etc. or what ?? and how and why ??
- # would need to merge all 5 gene groups angles dataframes
# 5 different types ?? or just ortholog and paralogs ??
# merge paralog.expr.angle.diag.df and orths.expr.angle.diag.df and plot results:
- check if the distances are calculated for each gene group separately
- 


- Perform unit tests for each function.   
- Finalize and validate unit tests for all scripts and Functions.  
- created doc\unit_tests_angles_functions.md to document test sceneraios for angles_funks functions




- improve Lab-Book descriptions for the 4 steps, and precisly document how we did the 4 steps
- check if the distances are calculated for each gene group separately


**Doubts and Issues**:


**Next Steps**:


---

**Code:**