**Date**: 05.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Improved `compute_exp.prof.dists.R` and `compute_exp.prof.dists_statistics.R`.
- Created functions to handle nested lists for computing distances and statistics.
- Filtered dataframes from `genegroups` to retain only genes intersecting with provided expression profiles, reducing computation time and improving results by minimizing NA/null values.
- Used checks (via `exists()`) to ensure required objects/lists (e.g., `genegroups`) are available, allowing scripts to run smoothly even with partial data.
- Created `add_to_list_if_exists()` in `compute_exp.prof.dists_statistics.R` to add only available objects to lists.

**Doubts and Issues**:

- Insufficient intersecting data between expression profiles and current gene data (FBpp), especially for `con_orthologs`.
- Filtering often results in too little data to compute meaningful distances/statistics or generate plots.
- Only `in_paralogs` and `special_in_paralogs` have adequate data for computations.
- Should we work with filtered data or retain all data, even if redundant?

**Next Steps**:

- Refactor `compute_exp.prof.dists.R` and `compute_exp.prof.dists_statistics.R` to reduce redundancy, simplify code, and improve readability. Potentially create a reusable function.
- Update descriptions in `roxygen2` documentation for clarity.
- try with other exp

---

**Code Samples**:

1. **Object Existence Check Function**:
    ```R
    add_to_list_if_exists <- function(list_obj, obj_name, var_name) {
        if (exists(var_name)) {
            list_obj[[obj_name]] <- get(var_name)
            cat(sprintf("Object '%s' created successfully.\n", obj_name))
        } else {
            cat(sprintf("Object '%s' does not exist and was not created.\n", obj_name))
        }
        return(list_obj)
    }
    ```

2. **Current `compute_exp.prof.dists.R` Implementation**:
    ```R
    require(GeneFamilies)
    options(mc.cores = getMcCores())
    library(dotenv)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(tibble)

    output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
    message("USAGE: Rscript exec/compute_exp.prof.dists.R")

    input.args <- commandArgs(trailingOnly = TRUE)
    load(file.path(output_data_dir, "gene_families.RData"))                       
    load(file.path(output_data_dir, "gene_groups_filtered.RData"))           
    load(file.path(output_data_dir, "gene_expression.RData")) 

    # Sourced function for distance calculations
    source("R/compute_funks.R")

    # Track created objects
    created_objects <- character()

    # Calculate distances and add created objects if they exist
    if (exists("con_orthologs_filtered.lst")) {
        con_orthologs.filtered.dists <- mclapply(con_orthologs_filtered.lst, exp.prof.dists)
        con_orthologs.filtered.dists.tissue <- mclapply(con_orthologs_filtered.lst, exp.prof.dists_tissue)
        created_objects <- c(created_objects, "con_orthologs.filtered.dists", "con_orthologs.filtered.dists.tissue")
    }

    if (exists("in_paralogs_filtered.lst")) {
        in_paralogs.filtered.dists <- mclapply(in_paralogs_filtered.lst, exp.prof.dists)
        in_paralogs.filtered.dists.tissue <- mclapply(in_paralogs_filtered.lst, exp.prof.dists_tissue)
        created_objects <- c(created_objects, "in_paralogs.filtered.dists", "in_paralogs.filtered.dists.tissue")
    }

    if (exists("out_paralogs_filtered.lst")) {
        out_paralogs.filtered.dists <- mclapply(out_paralogs_filtered.lst, exp.prof.dists)
        out_paralogs.filtered.dists.tissue <- mclapply(out_paralogs_filtered.lst, exp.prof.dists_tissue)
        created_objects <- c(created_objects, "out_paralogs.filtered.dists", "out_paralogs.filtered.dists.tissue")
    }

    if (exists("special_in_paralogs_filtered.lst")) {
        special_in_paralogs.filtered.dists <- mclapply(special_in_paralogs_filtered.lst, exp.prof.dists)
        special_in_paralogs.filtered.dists.tissue <- mclapply(special_in_paralogs_filtered.lst, exp.prof.dists_tissue)
        created_objects <- c(created_objects, "special_in_paralogs.filtered.dists", "special_in_paralogs.filtered.dists.tissue")
    }

    if (exists("special_out_paralogs_filtered.lst")) {
        special_out_paralogs.filtered.dists <- mclapply(special_out_paralogs_filtered.lst, exp.prof.dists)
        special_out_paralogs.filtered.dists.tissue <- mclapply(special_out_paralogs_filtered.lst, exp.prof.dists_tissue)
        created_objects <- c(created_objects, "special_out_paralogs.filtered.dists", "special_out_paralogs.filtered.dists.tissue")
    }

    # Save only created objects
    if (length(created_objects) > 0) {
        save(list = created_objects, file = file.path(output_data_dir, "exp.prof.dists_filtered.RData"))
        message("Objects saved to ", file.path(output_data_dir, "exp.prof.dists_filtered.RData"))
    } else {
        message("No objects were created, so nothing was saved.")
    }
    ```