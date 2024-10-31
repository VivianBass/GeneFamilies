
**Date**: 22.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

1. **Script Execution and Testing**:
   - Ran the R scripts in the `exec` directory up to the distance calculation step to confirm functionality.  
   - **Result**: Scripts function as expected up to this point.

2. **Script Review**:
   - Next step: Evaluate if all scripts are necessary or if they can be consolidated.

3. **Project Structure**:
   - Consider implementing a folder structure to categorize code sections into "Loading," "Computing," "Plotting," etc.

**Doubts & Issues**:
- **all.cds File**: Is the `all.cds` file necessary? It was previously used for generating expression profile data.
- **All-vs-All Similarity Metrics**: Is an all-vs-all BLAST/Diamond run required for similarity metrics? If so, should we use BLAST or DIAMOND?
- One requires `pairwiseSimilarities` data, which depends on the `All vs All` file, potentially generated from `all.cds`.

**Next Steps**:
- Implement the proposed structural adjustments.