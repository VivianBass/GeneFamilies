
**Date**: 28.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Received updated instructions on how to proceed.
- Update scripts to load **gene pairs** (both orthologs and paralogs) 
  instead of individual gene lists.
- Adjust scripts to support **pairwise data processing**.
- Measure distances as before, adding **log-transformed distance calculations**.

- Maintain five **Gene Groups**:

    - Conserved orthologs  
    - In-paralogs with orthologs  
    - In-paralogs without orthologs  
    - Out-paralogs with orthologs  
    - Out-paralogs without orthologs 

- Ensure all five groups are loaded: conserved orthologs, in-paralogs with orthologs,   
  in-paralogs without orthologs, out-paralogs with orthologs, and out-paralogs without orthologs.

**Doubts and Issues**:

- Need to clarify the method for performing log-transformed distance calculations. What is the   
  exact process for applying logarithmic transformation to distance measurements?

**Next Steps**:

- Adjust loading scripts to accommodate the five gene groups, updated data formats, and the OrthoFinder input format.