
## a quick 5-step workflow for using **Roxygen2** to document functions in an R package:

1. **Set Up Roxygen2 in Your Package**:

- Ensure your package has a `DESCRIPTION` file with `Roxygen2` listed as a dependency 
  under `Imports` or `Suggests`.
- Set `Roxygen` options in the `DESCRIPTION` file, typically like this:
    ```plaintext
    Roxygen: list(markdown = TRUE)
    ```

2. **Add Documentation Tags to Functions**:

- Above each function definition, add a `#'` line with Roxygen2 tags 
  to document the function. Common tags include:

    - `@title` (function name or title)
    - `@description` (what the function does)
    - `@param` (parameter details)
    - `@return` (what the function returns)
    - `@examples` (usage examples)

- Example:
    ```r
    #' Title: Add Two Numbers
    #' @description Adds two numeric values.
    #' @param x Numeric value.
    #' @param y Numeric value.
    #' @return The sum of x and y.
    #' @examples
    #' add_numbers(3, 5)
    add_numbers <- function(x, y) {
        x + y
    }
    ```

3. **Generate Documentation Files**:

- Run `devtools::document()` or `roxygen2::roxygenise()` to automatically generate `.Rd` files in the `man` directory. This command parses your Roxygen2 comments and creates corresponding documentation files.

4. **Check and Build the Package**:

- Use `devtools::check()` to ensure your package passes all checks, including documentation accuracy. This step helps verify that there are no missing parameters or mismatches in the documentation.

5. **Update and Maintain Documentation**:

- Every time you modify a function or add a new one, update the Roxygen2 comments accordingly. Re-run `devtools::document()` to keep the `.Rd` files in sync with your code changes. 

