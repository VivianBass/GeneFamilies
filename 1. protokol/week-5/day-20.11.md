
**Date**:  20.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Continued documenting in the Lab Book.  
- Documented the current dataset, paper, and data used for analysis.  
- **Key Paper**:  
  - *"Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species"*  
  - DOI: [https://doi.org/10.1016/j.celrep.2019.08.011](https://doi.org/10.1016/j.celrep.2019.08.011)  
  - Investigates dietary adaptability in *Drosophila*, comparing generalists (*D. melanogaster*, *D. simulans*) to specialists (*D. sechellia*).  

- Gathered data for *D. melanogaster* (*dmel*) and *D. sechellia* (*dsec*) from the referenced study.  
  - Obtained GFF and FASTA files from FlyBase.  
  - Downloaded FASTQ files from [NCBI BioProject PRJDB4481](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB4481).  
  - Plan to create reference transcriptomes using tools like GFFread, Trimmomatic, and Kallisto.  

- Created R scripts:  
  - **`exec/compute_exp.prof.dists_angles.R`**: Calculates angles for gene expression profiles.  
  - **`exec/plot_exp.prof.dists_angles.R`**: Plots these angles.  

---

**Doubts and Issues**:  

1. Need to verify and fix minor issues in the R scripts for angle calculation and plotting.   

---

**Next Steps**:

- Debug R scripts for angle calculations and plots.  
- Implement and test t-tests and Wilcoxon tests in the analysis. 
- Perform unit tests for each function.   
- Finalize and validate unit tests for all scripts.  

---

**Code:**

- **Function and Code to Calculate Angles for Each Gene Group**  

```R
calculate_angles <- function(genes, rna.seq.exp.profils, tissues) {
  genes.expr <- intersect(unlist(genes), rna.seq.exp.profils$FBpp_ID)
  
  if (length(genes.expr) == 0) {
    warning(paste("No matching genes found for group:", group))
    return(data.frame())
  }

  expr.angle.diag.df <- data.frame(
    FBpp_ID = genes.expr,
    angle.diag = as.numeric(mclapply(genes.expr, function(x) {
      cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == x), tissues]) / sqrt(2)
    })),
    stringsAsFactors = FALSE
  )

  expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")
}

angle_results <- list()

# Loop through each gene group
for (group in gene_groups) {
  gene_list <- get(group)
  df_name <- paste0(gsub("_v\\.lst$", "", group), ".expr.angle.diag.df")
  angle_results[[df_name]] <- calculate_angles(gene_list, rna.seq.exp.profils, tissues)
  assign(df_name, angle_results[[df_name]])
}

# Save the results
if (length(angle_results) == 0) {
  stop("No results to save")
}

save(list = names(angle_results), file = file.path(output_data_dir, "exp.prof.dists_angles.RData"))
message("DONE")
```