

**Date**:  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:


1. **Define Gene Groups:**  
   - Research and define the five gene groups precisely, sourcing papers and PDFs for reference.
   - Update the documentation in the `doc` directory to include detailed descriptions of each gene group and specify the format for input files.

2. **Documentation Updates:**  
   - Modify the *Section description* and *input_files* documentation to cover the five gene groups, ensuring they are well-described.

3. **Create Dummy Test Files:**  
   - Generate dummy files for quick testing of R scripts, including input files for the loading section.
   - Organize these in a `dummy-datasets` folder under `experiments/dummy-datasets` with mixed data types for diverse testing.

4. **Refactor Loading Functions:**  
   - Develop a function for loading gene groups, organizing all relevant loading functions in a file called `load_data_funks.R` within the `R` directory.
   - Remove these functions from the main scripts and instead source them from `load_data_funks.R`.

- added description in doc directory on how to per form unit-tests in R using the testthat package


**Doubts and Issues**:

**Next Steps**:


---

**Code:**

```R
# Function to load data frames with a specific header type.
# Header Type for Orthologs: Family | Gene | Gene_species | Ortholog | Ortholog_species
# Header Type for Paralogs:  Family | Gene | Gene_species | Paralog  | Paralog_species
# This function is used in load_gene_groups_data.R
load_data_frame <- function(file_path) {
    read.table(file_path, header = TRUE, sep = "\t", 
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))
}


# Function to create nested lists from data frames
# This function is used in load_gene_groups_data.R
create_nested_list <- function(df, header_type) {
    # Determine the appropriate column names based on header_type
    if (header_type == "Ortholog") {
        gene_col <- "Ortholog"
        species_col <- "Ortholog_species"
    } else if (header_type == "Paralog") {
        gene_col <- "Paralog"
        species_col <- "Paralog_species"
    } else {
        stop("Invalid header type provided.")
    }

    df %>%
        group_by(Family, Gene_species, Gene, !!sym(species_col)) %>%
        summarise(!!sym(gene_col) := list(!!sym(gene_col)), .groups = "drop") %>%
        group_by(Family, Gene_species, Gene) %>%
        summarise(nested_info = list(setNames(!!sym(gene_col), !!sym(species_col))), .groups = "drop") %>%
        group_by(Family) %>%
        summarise(gene_info = list(setNames(nested_info, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>%
        deframe()
}
```