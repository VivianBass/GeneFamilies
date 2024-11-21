
# Testing scenarios for loading functions 

#### **Function Under Test**: `load_data_frame()`
#### **Function Location**: `R/load_data_funks.R`
#### **used in rscript**: `exec/load_gene_groups_data.R`

The `load_data_frame` function reads a tab-separated file from a specified file path and loads its content into a data frame. The file is expected to have five columns, all of which are read as character strings. The function ensures the file has a header row and avoids interpreting special characters, quotes, or empty fields as data.

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_load_data_frame.md")
test_file("tests/testthat/test-load_data_frame.R")
sink()
```
---

### **Testing Scenarios**

1. **Szenario 1: Basic Case: Valid File with Correct Header and Data**
   - **Objective**: Ensure the function correctly loads a valid file with the expected header and data structure.
   - **Input**: A sample file with correct columns (`Family`, `Gene`, `Gene_species`, `Ortholog`, `Ortholog_species`).
   - **Expected Outcome**: The data frame should load without errors, have five columns with the expected names, and contain the test data accurately.
   
2. **Szenario 2: Missing Columns**
   - **Objective**: Check the function’s response to files missing required columns.
   - **Input**: A file missing one or more columns (e.g., only `Family`, `Gene`, `Gene_species`, and `Ortholog`).
   - **Expected Outcome**: An error or clear message indicating missing columns.

3. **Szenario 3: Incorrect Column Names**
   - **Objective**: Verify that the function correctly identifies files with unexpected column names.
   - **Input**: A file with incorrect column names (e.g., `Family`, `Gene`, `Gene_species`, `WrongName`, `Ortholog_species`).
   - **Expected Outcome**: An error or message indicating column names do not match the expected schema.

4. **Szenario 4: Empty File**
   - **Objective**: Confirm that the function can handle files with headers but no data rows.
   - **Input**: A file with the correct headers but zero data rows.
   - **Expected Outcome**: A data frame with the correct columns but zero rows.

5. **Szenario 5: Non-Tab-Delimited File**
   - **Objective**: Test the function's handling of files with different delimiters, such as commas or spaces.
   - **Input**: A file with columns separated by commas instead of tabs.
   - **Expected Outcome**: A warning or improperly loaded data frame, as the file does not meet the expected tab-delimited structure.

6. **Szenario 6: File with Extra Columns**
   - **Objective**: Check if the function ignores or handles files with additional, unexpected columns.
   - **Input**: A file with extra columns beyond the expected five.
   - **Expected Outcome**: The function should load only the specified columns or raise a warning about extra columns.

7. **Szenario 7: File with Incorrect Data Types**
   - **Objective**: Confirm that all columns are read as character data, even if the file contains different data types (e.g., numeric values).
   - **Input**: A file with numeric values in one or more columns.
   - **Expected Outcome**: The function should load all columns as strings without errors.

8. **Szenario 8: File Does Not Exist**
   - **Objective**: Ensure that the function handles missing file paths gracefully.
   - **Input**: A non-existent file path.
   - **Expected Outcome**: An error or clear message indicating the file could not be found.

---

#### **Function Under Test**: `create_nested_list()`
#### **Function Location**: `R/load_data_funks.R`
#### **used in rscript**: `exec/load_gene_groups_data.R`

The `create_nested_list` function transforms a data frame containing gene family information into a hierarchical nested list based on the specified `header_type` ("Ortholog" or "Paralog"). It groups the data by family, gene, and species, then organizes the corresponding ortholog or paralog information into a structured list. Each gene is associated with its species-specific data, and the final output is a list of families, each containing nested gene-specific information.

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_create_nested_list.md")
test_file("tests/testthat/test-create_nested_list.R")
sink()
```
---

### **Testing Scenarios**

1. **Basic Functionality: Nested List for Orthologs**
   - **Description**: Verifies that `create_nested_list` can correctly create a nested list from a data frame containing ortholog gene data. Checks that the output matches the expected nested structure when the header type is `"Ortholog"`.

2. **Basic Functionality: Nested List for Paralogs**
   - **Description**: Tests that `create_nested_list` correctly generates a nested list for paralog gene data. Confirms the expected structure when using the `"Paralog"` header type.

3. **Empty Data Frame Handling**
   - **Description**: Ensures that `create_nested_list` can handle an empty data frame input without error, and returns an empty list as expected.

4. **Multiple Families Handling**
   - **Description**: Verifies that `create_nested_list` correctly separates and nests genes within distinct families in the output structure. Tests the function’s ability to handle multiple families within the same data frame.

5. **Error Handling: Missing Required Columns**
   - **Description**: Confirms that `create_nested_list` throws an informative error when required columns (e.g., `Gene`, `Ortholog`) are missing from the data frame. This validates the function’s input validation.

6. **Header Type Variability**
   - **Description**: Tests that `create_nested_list` can handle different header types (e.g., `"Ortholog"` and `"Paralog"`) and produce the correct structure based on the header type specified.

7. **Duplicate Genes within the Same Family**
   - **Description**: Ensures that `create_nested_list` handles duplicate genes within the same family by merging their associated data under the same gene key in the output list.

8. **Non-Character Data Type Handling**
   - **Description**: Checks that `create_nested_list` can process non-character data types (e.g., factors and integers) and converts them to characters as needed, ensuring consistent data types in the output list.

9. **Error Handling: Invalid Header Type**
   - **Description**: Tests that `create_nested_list` throws an error if an invalid header type is provided, confirming that only valid types (like `"Ortholog"` or `"Paralog"`) are accepted.

---

#### **Function Under Test**: `filter_v1`
#### **Function Location**: `R/load_data_funks.R`
#### **used in rscript**: `exec/load_gene_groups_data.R`

The `filter_v1` function filters a data frame of gene data by retaining only those rows where the gene listed in the `"Gene"` column has a corresponding entry in the `expression_data` data frame (which contains gene identifiers in the `"FBpp_ID"` column). The function returns a new data frame containing only the genes with available expression data, removing any rows for genes not present in the expression dataset.

Test the file in the R terminal using the following command:
```R
sink("tests/test_results/test_results_filter_v1.md")
test_file("tests/testthat/test-filter_v1.R")
sink()
```
---

### Testing Scenarios for `filter_v1`

The `filter_v1` function filters rows from `df` based on matching gene identifiers in `expression_data`. The following test scenarios validate its functionality under typical and edge-case conditions:

1. **Basic Filtering with Matching Genes**  
   - **Objective**: Verify that only rows in `df` with genes present in `expression_data` are retained.
   - **Expected Outcome**: The filtered data frame should contain only rows with genes found in both data frames.

2. **No Matching Genes**  
   - **Objective**: Confirm that if no genes in `df` are found in `expression_data`, an empty data frame is returned.
   - **Expected Outcome**: A data frame with zero rows, indicating no matches.

3. **All Genes Match**  
   - **Objective**: Ensure that if all genes in `df` are present in `expression_data`, all rows in `df` are retained.
   - **Expected Outcome**: The output should be identical to `df`, containing all original rows.

4. **Empty `df` Input**  
   - **Objective**: Test function handling when `df` is empty (no rows).
   - **Expected Outcome**: Return an empty data frame, indicating correct handling of empty input.

5. **Empty `expression_data` Input**  
   - **Objective**: Validate that an empty `expression_data` results in no matches and thus an empty data frame.
   - **Expected Outcome**: Return an empty data frame, as there are no genes to match against.



