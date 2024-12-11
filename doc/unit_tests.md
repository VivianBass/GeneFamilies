
## Definition of a Unit-Test

Unit tests verify individual code components (like functions) in isolation to ensure correct functionality. They're automated, fast, and give quick feedback on code changes, catching bugs early and preventing regressions. Each test is independent and repeatable, often using mock data to simulate inputs. Unit tests also serve as documentation by defining expected behavior, which supports future maintenance. They follow an "Arrange, Act, Assert" structure: set up inputs, perform actions, and check outcomes.

## `testthat` Package Workflow for Unit Testing in R

### **1. Install and Load `testthat`:**
- Install the package with `install.packages("testthat")`.
- Load it into your script or session using `library(testthat)`.

### **2. Set Up Testing Infrastructure and Directory:**
- Run `usethis::use_testthat()` from the R terminal. This command creates a `tests/testthat` directory in your package and sets up a `tests/testthat.R` file, which is the main test file. The `testthat.R` file should be placed directly in the `tests` directory, while all individual test files will go into `tests/testthat/`.
- Ensure that a `tests/testthat` directory exists in your package directory, along with the `testthat.R` file in the `tests` folder. This setup is essential for starting unit tests with the `testthat` package.

- In `tests/testthat`, create test files for each function, 
  naming them with a `test-` prefix (e.g., `test-myfunction.R`).

**Example Directory Structure:**

```R
GeneFamilies/
├── R/
│   └── load_data_funks.R
├── tests/
│   ├── testthat/
│   │   └── test-load_data_frame.R
│   └── testthat.R
```

This basic setup provides the necessary structure to organize and execute unit tests using the `testthat` package.


### **3. create dummy Datasets for testing purposes**

#### Approaches for Adding Testing Data in R Packages

1. **`inst/extdata` Folder**: 
   For large or external data files, place `.csv`, `.rds`, or `.rda` files in `inst/extdata`. Load them in tests with:
   ```R
   file_path <- system.file("extdata", "your_data_file.csv", package = "GeneFamilies")
   your_data <- read.csv(file_path)
   ```

2. **`data/` Folder for Package Data**:
   For smaller datasets you want accessible via `data()`, save `.rda` files in `data/`. Add `LazyData: true` in `DESCRIPTION` to enable:
   ```R
   data("families", package = "GeneFamilies")
   ```

3. **`tests/testthat/` Folder**:
   For test-specific data, place files in `tests/testthat/` (e.g., `tests/testthat/data/your_data_file.csv`). Load with:
   ```R
   your_data <- read.csv("tests/testthat/data/your_data_file.csv")
   ``` 

Each option is suited to different needs.
Choose based on file size, package accessibility, or test-specific use.


### **4. Write Test Cases / Test Scenarios Using `test_that`:**

- Each test file contains `test_that()` blocks to define individual test scenarios. 
- Each `test_that()` block specifies a test scenario tailored to verify a particular aspect of the function's behavior under various conditions, helping ensure the function implementation works as intended.
- A good testing scenario should be clear, specific, and focused on a single condition or behavior.
- You can have multiple test scenarios within a single test file.

- Use `expect_*` functions such as `expect_equal()`, `expect_true()`, and `expect_error()` to confirm expected behavior.

    ```R
    test_that("myfunction returns the square of input", {
      expect_equal(myfunction(4), 16)
      expect_error(myfunction("text"))
    })
    ```

- This setup follows the standard unit testing pattern:

  - **Arrange**: Set up test data.
  - **Act**: Call the function being tested.
  - **Assert**: Verify that the results match expectations.

**Example Scenario 1:**

```R
test_that("load_data_frame handles extra columns gracefully", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\tExtraColumn\nF1\tG1\tSp1\tO1\tSpO1\textra", temp_file)
    df <- load_data_frame(temp_file)
    expect_equal(ncol(df), 5) # Only 5 columns should be loaded
    unlink(temp_file)
})
```

**Explanation of Example Scenario 1**:
- Creates a temporary test file with sample data.
- Calls `load_data_frame()` with this test data.
- Checks if the output matches expectations, specifically that only five columns are loaded, using `expect_equal()`.
- Cleans up the temporary file afterward.

This setup provides a reliable structure for testing various cases, ensuring that functions perform as expected across different scenarios.


### **5. Run Tests:**

- Use `devtools::test()` command from R-Terminal to run all tests in the `tests/testthat` directory.
- **All Tests**: Run all tests with `devtools::test()` from the package root.
- **Single Test File**: Run a specific test file with `test_file("tests/testthat/test-myfunction.R")`.
- Note: you dont need to link the functions you are testing. explicetly in the test file. 
they are called by the package. just have to match the names !!

- Capture test results to a file:

```R
sink("test_results.md")
devtools::test()
sink()
```

```R
cat("# Test Results\n\n```\n", file="test_results.md")
capture.output(devtools::test(), file="test_results.md", append=TRUE)
cat("\n```", file="test_results.md", append=TRUE)
```

This will create a test_results.md file in your current working directory with the test output formatted as a code block in markdown.


══ Results-File ════════════════════════════════════════════════════════════════════

This `testthat` output provides a structured summary of the test results for the `load_data_frame` function, highlighting successes and failures across test cases. Here’s a breakdown of its structure and what each part means:

1. **Summary Line** in the Header of the Results File:

```R
    ✔ | F W  S  OK | Context
    ⠏ |          0 |                                                                                                  
    ⠹ | 5 1      7 |                                                                                                   
    ✖ | 5 1      7 | load_data_frame
```

- The symbols indicate the following:
- **✔ (OK)**: Number of successful tests.
- **F (FAIL)**: Number of failed tests.
- **W (WARN)**: Number of warnings encountered during testing.
- **S (SKIP)**: Number of skipped tests.
- **OK**: Total count of successful tests at the end.
- The context (`load_data_frame`) shows which function is under test.

2. **Detailed Error Messages by Test**

The following lines detail errors and failures from specific tests. 
Each error or failure block includes:
- The **test name** or description.
- **Error message** and **backtrace** to help identify the source of the error.

3. **Summary Section (Results)**

The final summary aggregates the test outcomes:

```
[ FAIL 5 | WARN 1 | SKIP 0 | PASS 7 ]
```

- **FAIL 5**: Five tests failed due to issues like missing columns, incorrect column names, or other discrepancies.
- **WARN 1**: One test generated a warning (e.g., a `file not found` scenario).
- **SKIP 0**: No tests were skipped.
- **PASS 7**: Seven tests passed successfully.

This breakdown highlights areas where `load_data_frame` may need improvement, 
particularly in error handling for incorrect or unexpected file formats.