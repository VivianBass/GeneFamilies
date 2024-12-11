

**Date**: 29.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- **Goal Definition**: Established primary objectives for the coming weeks.
  
- **Data Loading and Distance Calculation**: 
Set up to load all five files related to gene families and expression values. 
Distances will be computed between gene pairs.

- **Gene Group Classification**: Created five distinct gene groups:
  1. **Conserved Orthologs**
  2. **In-Paralogs with Orthologs**
  3. **In-Paralogs without Orthologs**
  4. **Out-Paralogs with Orthologs**
  5. **Out-Paralogs without Orthologs**

- **File Structure and Headers**:
  - Defined headers for paralog files: 
    `Family`/ `Gene`/ `Gene_species`/ `Paralog`/ `Paralog_species`

  - Example header structure for **Conserved Orthologs**:  
    ```
    Family      Gene        Gene_species    Ortholog      Ortholog_species
    OG0000000   FBpp0117097 dana            FBpp0172663   dmoj
    ```

- **Script Adjustments**:
  - **Data Loading**:
    - Updated `exec/3.load_gene_groups_data.R` to load the five gene groups 
    (conserved orthologs, in-paralogs, out-paralogs, special_in-paralogs, special_out-paralogs).
    - Command for this script:  
        ```
        Rscript load_gene_groups_data.R <in_paralogs.tsv> <special_in_paralogs.tsv> <out_paralogs.tsv> <special_out_paralogs.tsv> <conserved_orthologs.tsv>
        ```
  - **Gene Families**:
    - Adapted `exec/2.load_gene_families_data.R` to the new OrthoFinder input format:
      ```
      Family      species1        species2        species3
      family_1    gene1,gene2,gene3  gene4,gene5,gene6  gene7,gene8,gene9
      ```
  - **Gene Expression**:
    - Revised `exec/1.load_gene_expression_data.R`.

- **Improved messages**: improved messages across all three data loading scripts.


**Doubts and Issues**:

- Clarify the specific differences among the five gene groups and how to categorize them accordingly.
- Need to clarify the method for performing log-transformed distance calculations. What is the   
  exact process for applying logarithmic transformation to distance measurements?


**Next Steps**:

- Adapt distance and statistical computations to include all five gene groups.
- Finalize the structure and required inputs for each of the five gene group files.
- perform log-transformed distance calculations (exec/5.compute_exp.prof.dists.R)

---

**Code**


- code for loading conserved orthologs data
```R
library(dplyr)
library(tidyr)
library(tibble)

con_orthologs <- read.table(input.args[[1]], header = TRUE, sep = "\t", 
                comment.char = "", quote = "", na.strings = "", 
                colClasses = rep("character", 5))

con_orthologs.lst <- con_orthologs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
                mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
                select(cluster_name, Gene) %>% deframe()
```



