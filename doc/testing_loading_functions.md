

# Testing scenarios for loading functions 

**Function Under Test**: `load_data_frame()`
**Function Location**: `R/load_data_funks.R`

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

**Function Under Test**: `load_data_frame()`
**Function Location**: `R/load_data_funks.R`

---

### **Testing Scenarios**