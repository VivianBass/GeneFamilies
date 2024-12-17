
## **Chapter 7: Discussion**

### 5.1 **Key Insights**

- Insights from the analysis of gene families and homologs.
- Comparison of results between statistical tests and algorithms.

### 5.2 **Limitations**

- Challenges encountered (e.g., data size, assumptions of statistical tests).

### 5.3 **Applications**

- Potential applications of EasyVectorOmics in other biological studies.



    - what did we find ?
    - what does it mean ?
    - compare with state of Knowledge
    - Confidence of the results
    - whats next and what questions arise, what could be done next or differently


◦ Do we have similar results in other public data as in the paper?









      - **Distance Interpretation**:

        - **Small Distances**: Indicate similar gene activity patterns, suggesting functional conservation
        - **Large Distances**: Reveal significant variations, suggesting functional divergence





- log2 slightly higher singnif levels
- different count numbers of genes in the gene groups , 
- which test better to use for which purpose ? mean median ?



To determine the most significant result for the t-test, you typically focus on the p-value (p) and the adjusted p-value (p.adj). The adjusted p-value accounts for multiple comparisons and is a more conservative measure of significance. Here's how you can use these metrics to identify the most significant result:

    p-value (p): This is the raw p-value from the t-test, indicating the probability of observing the test results under the null hypothesis.
    Adjusted p-value (p.adj): This is the p-value adjusted for multiple comparisons, often using methods like Bonferroni correction or False Discovery Rate (FDR). It provides a more stringent measure of significance.
    Significance Level (p.adj.signif): This indicates the level of significance (e.g., *, **, ***) based on the adjusted p-value.

Steps to Identify the Most Significant Result:

    Sort by Adjusted p-value (p.adj): The most significant result will have the smallest adjusted p-value.
    Check Significance Level (p.adj.signif): Ensure that the result with the smallest adjusted p-value also has the highest significance level (e.g., ***).
