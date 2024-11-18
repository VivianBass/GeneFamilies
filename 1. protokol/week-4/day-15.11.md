
**Date**: 15.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`  

---

### **Tasks**  

- Recreated dataframes to include both FBgn and FBpp identifiers using the mapping table.  
- Updated the two test directories (`test_diet_M` and `test_diet_P`) to include results for both diets.  
- Adjusted R scripts to fit the updated data structure and resolved issues during result generation.  
- Generated results and boxplots for both directories.  
- Expanded results to include con_orthologs and other categories, in addition to the previously analyzed *special_in_paralogs*.  
- Confirmed partial success with mapping but noted insufficient data for complete coverage.  

---

### **Doubts and Issues**  

- Should the boxplots differentiate between species (dmel/dsec) or focus solely on gene groups?  
- Discrepancy between significance levels in boxplots and CSV files—are the t-tests being applied correctly?  

---

### **Next Steps**  

- Investigate and correct t-test implementation to resolve significance level discrepancies.  
