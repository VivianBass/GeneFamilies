
## Definition of a Unit-Test

Unit tests verify individual code components (like functions) in isolation to ensure correct functionality. They're automated, fast, and give quick feedback on code changes, catching bugs early and preventing regressions. Each test is independent and repeatable, often using mock data to simulate inputs. Unit tests also serve as documentation by defining expected behavior, which supports future maintenance. They follow an "Arrange, Act, Assert" structure: set up inputs, perform actions, and check outcomes.

## `testthat` Package Workflow for Unit Testing in R

**1. Install and Load `testthat`:**
- Install with `install.packages("testthat")`.
- Load using `library(testthat)`.

**2. Set Up Testing Infrastructure:**
- Create a `tests/testthat` directory within your package by running `usethis::use_testthat()`. This also sets up a `tests/testthat.R` file to recognize tests.

**3. Create Test Files:**
- In `tests/testthat`, create test files for each function, naming them with a `test-` prefix (e.g., `test-myfunction.R`).

**4. Write Test Cases Using `test_that`:**
- Each test file uses `test_that()` blocks to define individual tests:
    ```R
    test_that("myfunction returns the square of input", {
    expect_equal(myfunction(4), 16)
    expect_error(myfunction("text"))
    })
    ```
- Use `expect_*` functions like `expect_equal()`, `expect_true()`, and `expect_error()` to assert correct behavior.

**5. Run Tests:**
- **All Tests**: Run all tests with `devtools::test()` from the package root.
- **Single Test File**: Run a specific test file with `test_file("tests/testthat/test-myfunction.R")`.

**6. Review Test Output:**
- `testthat` will display any failures or errors, making it easy to spot and resolve issues quickly.