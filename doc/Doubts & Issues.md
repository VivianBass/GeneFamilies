
- **Distinguishing Orthologs and Paralogs**: 
How can we effectively differentiate between orthologs and paralogs? 
Can a gene be classified as both, depending on context?

- **Finding Tandem Duplicates**: 
What is the best approach for identifying tandem duplicates? 
Is it feasible to write a short script to parse the GFF or GTF file for *Drosophila melanogaster*, 
considering the proximity of genes within five positions?

- **Naming Convention for Drosophila Genes**: 
How do we address the naming conventions of genes across different *Drosophila* species 
that share the same identifiers? 

- **Plotting Considerations**: 
Which statistics should we plot—means, medians, or both? 
Should we use ortholog/paralog groups or gene family clusters for generating distribution plots?


- negative angles möglich ?? invalid numbers ?? bei compute angles R







#### Open questions and issues: `Why do we divide the cosine by sqrt(2)?`

**Function**: **cosDiag** 

- what cosDiag is doing?
Calculates the angle for each vector group (a and b), where vector **a** 
contains numeric distance values and vector **b** is the diagonal vector 
defined as \( b = [1, 1, 1, \ldots] \), with \( n \) representing 
the length of vector **b** (the number of elements in **a**).

- **Cosine Values**: 
At specific angles, the cosine values indicate similarity: 
0° (1) and 360° (1) show maximum alignment; 45° (0.707) suggests moderate similarity; 
90° (0) indicates no similarity; 180° (-1) shows maximum dissimilarity; 
270° (0) reflects reversed perpendicular orientation.

- **Normalization Impact**: 
As \( n \) increases, the impact of cosine similarity 
diminishes due to the increasing denominator. For smaller \( n \), scaling has less effect, 
highlighting the need for normalization in larger or higher-dimensional datasets.