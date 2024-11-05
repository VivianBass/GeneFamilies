**Date**: 04.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Started a new assignment with four experiments based on datasets for different gene groups:
  - Gene groups were defined using a harmonic mean cut-off for conserved orthologs.
  - Datasets categorized by cut-off levels: 00% (no cut-off), 70%, 80%, and 90%.
  
- Generated a new expression profiles dataframe to match protein sequence names (FGpp) instead of gene IDs (FGgn):
  - Used a server-downloaded mapping table (FGgn to FGpp).
  - Selected only the FGpp with the longest matching sequence to avoid multiple proteins mapping to a single gene.

- Created four directories in the `experiments` folder for each dataset and began generating data with scripts from the `exec` folder, especially using the loading rscripts.

**Doubts and Issues**:

- Encountered issues with the nested structure in gene group lists, requiring adjustments to loading functions and scripts.
- Need to add logarithmic calculations for distances after resolving formatting issues.

**Next Steps**:

1. Adjust function to account for different gene group list formats for accurate distance calculations.
2. Integrate logarithmic calculations into the scripts.
3. Start statistical analysis of distances using the R statistics script.