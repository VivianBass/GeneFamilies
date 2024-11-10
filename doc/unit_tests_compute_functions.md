

# Testing scenarios for loading functions 

**Function Under Test**: `exp.prof.dists()`
**Function Location**: `R/compute_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_exp.prof.dists.md")
test_file("tests/testthat/test-exp.prof.dists.R")
sink()
```

---

### Test Scenarios

1. **"Calculates Euclidean distances for multiple genes"**  
   - **Objective**: Verify that Euclidean distances are correctly calculated for multiple genes across tissues.
   - **Expected Outcome**: Returns a vector of distances, with the expected length (for three genes, there are three pairwise comparisons).

2. **"Returns NA if only one gene is provided"**  
   - **Objective**: Ensure that if only one gene is available, the function returns `NA`.
   - **Expected Outcome**: `NA` is returned since pairwise distances require at least two genes.

3. **"Handles missing columns gracefully"**  
   - **Objective**: Confirm that an error is raised if specified tissue columns are missing in `expression.profiles`.
   - **Expected Outcome**: An error indicating that one or more columns are undefined.

4. **"Computes distances using specified distance method"**  
   - **Objective**: Validate that the function can use a non-default distance calculation method (`manhattan`).
   - **Expected Outcome**: Returns a vector of distances calculated with the specified method.

5. **"Returns empty vector if no genes match"**  
   - **Objective**: Test the function’s behavior when there are no matching genes in `expression.profiles`.
   - **Expected Outcome**: The function should return `NA` as no matches were found.

---

**Function Under Test**: `exp.prof.dists_tissue()`
**Function Location**: `R/compute_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_exp.prof.dists_tissue.md")
test_file("tests/testthat/test-exp.prof.dists_tissue.R")
sink()
```

---

### Test Scenarios

1. **"Calculates distances per tissue for multiple genes"**  
   - **Objective**: Verify that the function computes distances for each tissue individually when multiple genes are provided.
   - **Expected Outcome**: A list with an entry per tissue, each containing a vector of distances with the correct length (for three genes, three pairwise comparisons).

2. **"Returns NA if only one gene is provided"**  
   - **Objective**: Ensure that `NA` is returned when there is only one gene, as pairwise distances require at least two genes.
   - **Expected Outcome**: `NA` since there are not enough profiles for distance calculation.

3. **"Handles non-matching genes gracefully"**  
   - **Objective**: Check that the function returns `NA` if there are no matching genes between `gene.accessions` and `expression.profiles`.
   - **Expected Outcome**: `NA` as there are no profiles to calculate distances.

4. **"Handles missing tissue columns gracefully"**  
   - **Objective**: Confirm that an error is raised if specified tissue columns are not present in `expression.profiles`.
   - **Expected Outcome**: An error indicating missing columns.

5. **"Computes distances using specified distance method"**  
   - **Objective**: Validate that the function can use a non-default distance calculation method, such as `"manhattan"`.
   - **Expected Outcome**: A list of vectors with distances calculated using the specified method.



---

**Function Under Test**: `calculate_exp.prof.dists.statistics()`
**Function Location**: `R/compute_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_calculate_exp.prof.dists.statistics.md")
test_file("tests/testthat/test-calculate_exp.prof.dists.statistics.R")
sink()
```

---

### Test Scenarios

1. **"Computes mean and median for each family"**  
   - **Objective**: Verify the function correctly calculates the mean and median for each gene family.
   - **Expected Outcome**: The result contains accurate mean and median values for each family based on provided data.

2. **"Handles matrices with NA values"**  
   - **Objective**: Ensure the function computes statistics while ignoring `NA` values.
   - **Expected Outcome**: Mean and median values should be calculated by ignoring `NA`s.

3. **"Returns empty tibble for empty data"**  
   - **Objective**: Check that an empty list input results in an empty tibble with the correct column names.
   - **Expected Outcome**: An empty tibble with columns `Family`, `Mean`, and `Median`.

4. **"Filters out infinite values"**  
   - **Objective**: Confirm that `Inf` values are filtered out and not included in calculations.
   - **Expected Outcome**: Mean and median are computed after filtering out `Inf` values.

5. **"Handles single-value matrices"**  
   - **Objective**: Validate that the function handles matrices with only a single value, returning that value as both the mean and median.
   - **Expected Outcome**: Mean and median are equal to the single matrix value for each family.


---

**Function Under Test**: `calculate_exp.prof.dists.tissue.statistics()`
**Function Location**: `R/compute_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_calculate_exp.prof.dists.tissue.statistics.md")
test_file("tests/testthat/test-calculate_exp.prof.dists.tissue.statistics.R")
sink()
```

---

### Test Scenarios

1. **"Computes mean and median per tissue for each family"**  
   - **Objective**: Verify the function accurately computes mean and median values for each gene cluster and tissue when valid data is provided.
   - **Expected Outcome**: `Mean` and `Median` values should match expected calculations, and `Family` and `Tissue` columns should label each entry correctly.

2. **"Handles tissues with NA values"**  
   - **Objective**: Confirm that `NA` values are ignored in mean and median calculations.
   - **Expected Outcome**: Results should only include non-`NA` values in calculations.

3. **"Returns empty tibble for empty data"**  
   - **Objective**: Ensure that an empty list input results in an empty tibble with the correct columns.
   - **Expected Outcome**: An empty tibble with columns `Family`, `Tissue`, `Mean`, and `Median`.

4. **"Filters out infinite values"**  
   - **Objective**: Confirm that the function ignores `Inf` values during calculations.
   - **Expected Outcome**: Results should include only finite values in mean and median calculations.

5. **"Handles clusters with single-value tissues"**  
   - **Objective**: Validate the function’s ability to handle tissues with a single distance value, where mean and median are the same.
   - **Expected Outcome**: `Mean` and `Median` values should equal the single value for each tissue.


---

**Function Under Test**: `validate_data()`
**Function Location**: `R/compute_funks.R`

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_validate_data.md")
test_file("tests/testthat/test-validate_data.R")
sink()
```

---

### Test Scenarios

1. **"Correctly identifies valid lists and vectors"**  
   - **Objective**: Ensure that the function includes only valid lists and vectors.
   - **Expected Outcome**: `test_list_valid` and `test_vector_valid` should be included in `valid_data_names`.

2. **"Excludes objects with only NA values"**  
   - **Objective**: Confirm that the function excludes lists or vectors containing only `NA` values.
   - **Expected Outcome**: `test_list_na` and `test_vector_na` should be excluded from `valid_data_names`.

3. **"Excludes empty lists"**  
   - **Objective**: Verify that empty lists are excluded from the results.
   - **Expected Outcome**: `test_empty_list` should not appear in `valid_data_names`.

4. **"Excludes invalid data types"**  
   - **Objective**: Ensure that objects that aren’t lists or vectors (like data frames) are excluded.
   - **Expected Outcome**: `test_invalid_type` should not appear in `valid_data_names`.

5. **"Returns empty vector if no objects match pattern"**  
   - **Objective**: Test that the function correctly returns an empty vector when no object names match the specified pattern.
   - **Expected Outcome**: `valid_data_names` should be an empty character vector.

6. **"Returns empty vector if no valid objects"**  
   - **Objective**: Test that the function returns an empty vector if all matched objects are invalid, empty, or contain only `NA`.
   - **Expected Outcome**: `valid_data_names` should be an empty character vector when matching only invalid or empty objects.

