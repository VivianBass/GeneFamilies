

**Date**:  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:


- exact definition of the 5 gene groups, get some papers , pdfs etc on that 
- 

- have to update the Section descrition file and the input_files in documentation to 
account for 5 gene groups and also include them in the description 

- VPN, sophos, work on server 



- Function for genegroups 
- unit tests for the functions and do a functions file for the loading script functions
- Add all functions that we need to a special function file called "load_data_funks.R" inside R directory and delete them from your scripts. Call the file in your scripts so we don't have functions inside the exec scripts anymore.
- All the functions inside your load_data_funks.R should have unit tests. Please read:
	https://smartbear.com/learn/automated-testing/what-is-unit-testing/
	https://www.geeksforgeeks.org/unit-testing-in-r-programming/



- improved the the documentation files inside the doc directory, to also include the 
descriptions for the 5 genegroups and the format of the input files

- see how to change dists rscript for logarythmic data

- construct some dumy files for quickly testing , 
- created a folder called dummy-datasets in experiments/dummy-datasets used to located
mixed data of different kind for testing purposes

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
# used in load_gene_groups_data.R
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