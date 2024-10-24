

# R-Package-Structure

- **Package Creation**: Develop the package named `EasyVectorOmics`. 

- **Package Skeleton**: 

  - Reference: [package.skeleton](https://search.r-project.org/R/refmans/utils/html/package.skeleton.html) 
  
  - Structure of `EasyVectorOmics`:
    ```
    EasyVectorOmics (Package):
    ├── DESCRIPTION      # Contains title, author, description, license, dependencies, etc.
    ├── NAMESPACE        # For exports and imports declarations
    ├── R/               # Directory for R code (use roxygen2 comments for documentation)
    ├── man/             # Directory for documentation files (.Rd)
    ├── data/            # For included datasets
    ├── inst/            # For additional files to be installed with the package
    ├── tests/           # For unit tests (using testthat)
    ├── exec/            # For executable scripts
    ├── vignettes/       # For long-form documentation
    ├── LICENSE          # Appropriate license
    └── README.md        # Information about the package (optional)
    ```

- **Essential Directories and Files**:
  - `R/`: Contains R code files (.R).
  - `man/`: For documentation files (.Rd).
  - `DESCRIPTION`: Metadata about the package; list all required packages.
  - `NAMESPACE`: Defines exports and imports.

- **Optional but Recommended Directories**:
  - `data/`: For included datasets.
  - `inst/`: For additional files to be installed.
  - `tests/`: For unit tests.
  - `vignettes/`: For detailed documentation.

- **Creating the Package Skeleton**:
Use `package.skeleton(name = "EasyVectorOmics")` or `devtools::create_package("path/to/mypackage")`.
Example command: `devtools::create_package("C:/Users/andre/Desktop/EasyVectorOmics")`.

- **Building and Checking**:
After creating the package structure, populate the `DESCRIPTION`, `NAMESPACE`, 
and `R/` directories with metadata, namespace information, and functions.
Use roxygen2 comments to automatically generate documentation files in the `man/` directory.
Add a proper `LICENSE` file and `README.md`.

- **Use devtools for Package Management**:
Build the package with `devtools::build()`.
Run checks with `devtools::check()`.

- **Documentation**:
Use roxygen2 for function-level documentation.
For package-level docs, create a `package.R` file with roxygen comments.
Create vignettes using `usethis::use_vignette()` for detailed tutorials and examples.

- **Testing**:
Write unit tests for individual functions and integration tests 
for how functions work together using the `testthat` package.
Set up testing infrastructure with `usethis::use_testthat()`.
  
- **Error Handling**: 
Implement proper error checking in functions.

- **CRAN Policies**: 
Review policies if submitting to CRAN.

- **Continuous Integration**: 
Consider setting up CI/CD (e.g., GitHub Actions).

- **Namespace Documentation**: 
Refer to `vignette("namespace")` for generating 
a `NAMESPACE` file and understanding namespacing in R with roxygen2.
Here’s a summary of the most important functions for managing an R package:


## usefull Functions for Package Development

### Package Setup Functions (Typically Called Once)

- **`create_package()`**: Initializes a new R package structure.
- **`use_git()`**: Sets up Git for version control.
- **`use_mit_license()`**: Adds an MIT license file to the package.
- **`use_testthat()`**: Integrates the `testthat` package for unit testing.
- **`use_github()`**: Connects the package to a GitHub repository.
- **`use_readme_rmd()`**: Creates a README file in RMarkdown format.

### Development Functions (Called Regularly)

- **`use_r()`**: Adds a new R script to the package.
- **`use_test()`**: Creates a new test script for functions.
- **`use_package()`**: Adds a package dependency to the `DESCRIPTION` file.

### Frequent Development Functions (Called Multiple Times Per Day)

- **`load_all()`**: Loads all functions and scripts in the package for immediate use.
- **`document()`**: Generates documentation files from roxygen2 comments.
- **`test()`**: Runs all unit tests to verify code functionality.
- **`check()`**: Validates the package, ensuring it meets R package standards.
- **`install()`**: Installs or reinstalls the package, similar to an update.

These functions facilitate the development, testing, and maintenance of R packages, 
ensuring efficient workflows and adherence to best practices.
