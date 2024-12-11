
## **Statistical Analysis of Expression Distances**

### Purpose

- Analyze distribution of Euclidean distances within gene families
- Compare expression patterns between orthologs and paralogs
- Assess tissue-specific expression differences

### Analysis Levels

1. **Overall Expression Distances**
   - Distribution of pairwise distances
   - Mean distances per gene group
   - Median distances per gene group

2. **Tissue-Specific Analysis**
   - Per-tissue distance distributions
   - Tissue-specific mean distances
   - Tissue-specific median distances

3. **Additional Metrics**
   - Expression angles
   - Expression versatility
   

### Statistical Testing

1. **Mean Comparison**
   ```R
   # Compare means between orthologs and paralogs
   t.test(ortholog_distances, paralog_distances)
   ```

2. **Distribution Comparison**
   ```R
   # Compare overall distributions
   wilcox.test(ortholog_distances, paralog_distances)
   ```

3. **Multiple Testing Correction**
   ```R
   # Adjust p-values using BH method
   p.adjust(p_values, method="BH")
   ```


