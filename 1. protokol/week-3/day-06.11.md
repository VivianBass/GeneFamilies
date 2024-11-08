
**Date**:  06.11.2024 
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

- Streamline `compute_exp.prof.dists.R` and `compute_exp.prof.dists_statistics.R` to reduce redundancy and improve readability.
- Consider creating a reusable function for common tasks.

- Improve clarity of `roxygen2` descriptions and update overall documentation for consistency.

- Since the current expression profiles have limited data, try using alternative profiles.
- Ensure they use the same FBgn-to-FBpp mapping from the server’s FASTA file.
- Test different expression profiles to see if they yield a higher intersection and data volume.

- Try filtering `genegroups` dataframes by only one column to observe differences.

- Applyed two filtering methods (`v1` and `v2`) to align data with expression profiles.
- Verify if `v1` and `v2` yield the same results and length.


**Doubts and Issues**:

**Next Steps**:


---

**Code:**