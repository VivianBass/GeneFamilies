
# Testing scenarios for compute functions 

#### **Function Under Test**: `exp.prof.dists()`
#### **Function Location**: `R/compute_funks.R`
#### **used in rscript**: `exec/compute_exp.prof.dists.R`

The `exp.prof.dists` function calculates pairwise Euclidean distances between gene expression profiles for a specified list of genes. It filters the `expression.profiles` data to include only the specified genes, then computes the distances across the selected tissue columns. The function splits the data by species and calculates the distance matrix for each species group. It returns a vector of distance values or `NA` if there are fewer than two expression profiles for a species.

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

#### **Function Under Test**: `exp.prof.dists_tissue()`
#### **Function Location**: `R/compute_funks.R`
#### **used in rscript**: `exec/compute_exp.prof.dists.R`

The `exp.prof.dists_tissue` function computes pairwise Euclidean distances between gene expression profiles for a list of specified genes, calculated separately for each tissue. It filters the expression data to include only the specified genes and the relevant tissue columns, then calculates the distance for each tissue individually. The function returns a list of distance vectors, one for each tissue, or `NA` if there are fewer than two expression profiles.

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

#### **Function Under Test**: `calculate_exp.prof.dists.statistics()`
#### **Function Location**: `R/compute_funks.R`
#### **used in rscript**: `exec/compute_exp.prof.dists_statistics.R`

The `calculate_exp.prof.dists.statistics` function computes the mean and median of expression profile distances for each gene family. It takes a list of matrices, where each matrix represents the distance matrix for a gene family's expression profiles. The function calculates the mean and median distances for each family, handling missing or infinite values by excluding them from the calculations. It returns a tibble with the gene family's name, along with the mean and median values of the distances.

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

#### **Function Under Test**: `calculate_exp.prof.dists.tissue.statistics()`
#### **Function Location**: `R/compute_funks.R`
#### **used in rscript**: `exec/compute_exp.prof.dists_statistics.R`

The `calculate_exp.prof.dists.tissue.statistics` function calculates the mean and median of expression profile distances for each gene cluster, separated by tissue type. 
It takes a nested list where each gene cluster contains expression profile distances for multiple tissues. For each gene cluster, the function computes the mean and median of the distances for each tissue, and returns the results in a tibble. Each row in the tibble corresponds to a specific gene cluster and tissue, with columns for the gene cluster name (`Family`), the tissue name (`Tissue`), the mean of the distances (`Mean`), and the median of the distances (`Median`).


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

#### **Function Under Test**: `validate_data()`
#### **Function Location**: `R/compute_funks.R`
#### **used in rscript**: `exec/compute_exp.prof.dists_statistics.R`

The `validate_data` function filters and validates object names from a list (`loaded_objects`) based on a regular expression or predefined naming `pattern`. It retrieves matching objects using `get(name)` and excludes those that:  

1. Are not lists or vectors.  
2. Are empty or contain only `NA` values.  

The function returns a vector of valid object names and logs messages for excluded objects, specifying the reason for exclusion. basically catching objects that have the matching name according to pattern and also verifing the objects are not empty  

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_validate_data.md")
test_file("tests/testthat/test-validate_data.R")
sink()
```
---

### Test Scenarios  

1. **Identifies Valid Lists and Vectors**  
   - **Objective**: Verify that the function includes only valid lists and vectors that match the pattern.  
   - **Expected Outcome**: `test_list_valid` and `test_vector_valid` should be included in `valid_data_names`.  

2. **Excludes Objects Containing Only NA Values**  
   - **Objective**: Ensure that lists or vectors consisting entirely of `NA` values are excluded.  
   - **Expected Outcome**: `test_list_na` and `test_vector_na` should not appear in `valid_data_names`.  

3. **Excludes Empty Lists**  
   - **Objective**: Confirm that the function excludes lists that are empty.  
   - **Expected Outcome**: `test_empty_list` should not be present in `valid_data_names`.  

4. **Excludes Invalid Data Types**  
   - **Objective**: Ensure that objects that are neither lists nor vectors (e.g., data frames) are excluded.  
   - **Expected Outcome**: `test_invalid_type` should not be included in `valid_data_names`.  

5. **Handles Patterns with No Matches**  
   - **Objective**: Verify that the function returns an empty character vector when no object names match the specified pattern.  
   - **Expected Outcome**: `valid_data_names` should be an empty character vector.  

6. **Handles Scenarios with No Valid Objects**  
   - **Objective**: Ensure that the function returns an empty vector if all objects matching the pattern are invalid, empty, or consist only of `NA`.  
   - **Expected Outcome**: `valid_data_names` should be an empty character vector.  