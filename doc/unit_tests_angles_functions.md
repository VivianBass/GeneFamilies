
# Testing scenarios for loading functions 

**Function Under Test**: `calculate_angles()`
**Function Location**: `R/angles_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_calculate_angles.md")
test_file("tests/testthat/test-calculate_angles.R")
sink()
```

---

### **Test Scenarios**

1. **"Computes diagnostic angles for valid genes and tissues"**  
   - **Objective**: Verify that the function calculates diagnostic angles correctly when given a valid list of genes, RNA-seq expression profiles, and tissue columns.  
   - **Expected Outcome**: The result is a data frame containing accurate `angle.diag` values for each gene in the list, and it includes only the matching gene identifiers.  

2. **"Handles empty gene groups gracefully"**  
   - **Objective**: Ensure the function behaves appropriately when the input `genes` list is empty.  
   - **Expected Outcome**: The function returns an empty data frame and raises a warning indicating no matching genes were found.  

3. **"Ignores genes with no match in RNA-seq profiles"**  
   - **Objective**: Verify the function's behavior when no genes in the `genes` list match the `FBpp_ID` column in `rna.seq.exp.profils`.  
   - **Expected Outcome**: The function returns an empty data frame and raises a warning indicating no matching genes were found.  

4. **"Handles invalid tissue columns"**  
   - **Objective**: Test the function's response when the `tissues` vector includes column names that do not exist in `rna.seq.exp.profils`.  
   - **Expected Outcome**: The function raises an error indicating the tissue columns are invalid or missing.  

5. **"Handles empty RNA-seq expression profile data"**  
   - **Objective**: Verify the function's behavior when `rna.seq.exp.profils` is an empty data frame.  
   - **Expected Outcome**: The function returns an empty data frame and raises a warning about the lack of input data.  

6. **"Processes duplicate gene identifiers correctly"**  
   - **Objective**: Ensure the function handles duplicate gene identifiers in the input `genes` list appropriately.  
   - **Expected Outcome**: The function computes angles for each unique gene identifier and includes them in the result without duplication.  

7. **"Handles non-numeric tissue data"**  
   - **Objective**: Verify the function's response when `rna.seq.exp.profils` contains non-numeric data in the tissue columns.  
   - **Expected Outcome**: The function raises an error indicating the presence of invalid data types.  

8. **"Processes large datasets efficiently"**  
   - **Objective**: Ensure the function handles a large number of genes and tissues in a reasonable amount of time without errors or performance degradation.  
   - **Expected Outcome**: The function successfully computes angles for all valid genes and tissues in the input data, with a result of appropriate size.  

9. **"Returns correctly structured output"**  
   - **Objective**: Verify that the function's output has the expected structure (a data frame with `FBpp_ID` and `angle.diag` columns).  
   - **Expected Outcome**: The output is a data frame with the specified column names and no missing or invalid data in the `angle.diag` column.  

10. **"Handles genes without expression values"**  
    - **Objective**: Test the function's behavior when some genes in the `genes` list have no expression values in the tissue columns of `rna.seq.exp.profils`.  
    - **Expected Outcome**: The function excludes genes with missing or NA values in tissue columns and computes angles for the rest.  

---

**Function Under Test**: `validate_angle_dataframes()`
**Function Location**: `R/angles_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_validate_angle_dataframes.md")
test_file("tests/testthat/test-validate_angle_dataframes.R")
sink()
```

---

### **Test Scenarios**

1. **"Validates and retains non-empty data frames"**  
   - **Objective**: Ensure the function retains data frames with at least one row.  
   - **Expected Outcome**: The returned list contains only the non-empty data frames with their names preserved.

2. **"Excludes null data frames"**  
   - **Objective**: Ensure the function excludes `NULL` entries from the input list.  
   - **Expected Outcome**: The returned list does not include any `NULL` entries, and a warning message is displayed for each excluded entry.

3. **"Excludes empty data frames"**  
   - **Objective**: Ensure the function excludes data frames with zero rows.  
   - **Expected Outcome**: The returned list does not include any empty data frames, and a warning message is displayed for each excluded entry.

4. **"Handles a list with only null or empty data frames"**  
   - **Objective**: Test the function's behavior when all entries in the input list are either `NULL` or empty.  
   - **Expected Outcome**: The function returns an empty list and displays warnings for all entries.

5. **"Processes a list with mixed valid and invalid entries"**  
   - **Objective**: Verify the function can handle a mix of valid, empty, and `NULL` entries in the input list.  
   - **Expected Outcome**: The returned list contains only the valid entries, and warnings are displayed for invalid ones.

6. **"Handles an empty input list"**  
   - **Objective**: Ensure the function handles an empty input list without errors.  
   - **Expected Outcome**: The function returns an empty list without any warnings.

