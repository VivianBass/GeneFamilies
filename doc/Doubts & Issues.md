



- results
Looking at our current results, it seems that the orthologs are quite tissue specific, i.e. the are expressed in a few, maybe even just a single tissue while the paralogs
obtained more tissues in which they are expressed, i.e. they are less tissue specific.

Evidence for this:

- intra-ortholog distances compared to intra-paralog distances
([mean|median]_exp.prof.dists_boxplot-with-pvalues.pdf)

- the angles between mean ortholog expression to the diag
are signif. greater than those between mean paralog expression to the diag (Angles_boxplot-with-pvalues.pdf)

- the intra-ortholog rotation angles (change in tissue specificity) are greater than the intra-paralog rotation angles
(relativeExpressionVersatilityBoxplot.pdf)





-------------------------------------------------------------------------------------------







- #### Doubt: regarding Genesets, Subgroups

- See whether we can find other meaningful sub-groups in the families to investigate what happened to gene expression during the evolution of a particular gene family? 
- E.g. tandem clusters or conserved domain architecture or something else.



## Open questions and issues: `Why do we divide the cosine by sqrt(2)?`

- Try out the functions with given angles and see what comes out. 
- Maybe then we can see why we do this division. (angles to try could be 0,45,90,..,180,..,360)


- Functions: 

- Angle for each row Vector group (a and b) -> Vector a and b (contain numeric dist values)
- a as x Vector and b as diagonal Vector  d.v
- Vektor b ein Diagonalvektor ist, das heißt, er ist ein Vektor, in dem alle Einträge gleich 1 sind
- b=[1,1,1,…] , n ist die Länge des Vektors 𝑏 also die Anzahl der Elemente in 𝑎
- Ein höheres 𝑛 (mayne as the number of tissues present in the rpkm.expr.profiles for the vector)
  führt zu einer kleineren Größe des Terms sqrt(2) / sqrt(n), was bedeutet, dass der Wert von cosDiag
  cosDiag verringert wird, wenn mehr Datenpunkte betrachtet werden.


cosDiag(x = rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == x), expr.cols], )/`sqrt(2)`

a = x = = rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == x), expr.cols]
b = d.v = rep(1, length(x))

cosDiag --> {`sum(a * b)/(sqrt(sum(a^2)) * sqrt(sum(b^2)))`} /`sqrt(2)`

--> `(sum(a * b) * sqrt(2)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))`

--> `(sum(a) * sqrt(2)) / (sqrt(n) * sqrt(sum(a^2)))`

-->  `sqrt(2) / sqrt(n)`





- 
Dividing by 2 scales distances and angles uniformly, preserving the shape of vector clouds but does not fully normalize the data.
In high-dimensional spaces like gene expression data, dividing by 2 helps maintain consistent scaling for comparisons of distances and angles, especially between orthologs and paralogs.
This rescaling aids in gene expression analysis by standardizing the interpretation of tissue-specific and versatile expression patterns across datasets.

The expression `sqrt(2) / sqrt(n)`  is a form of standardization that accounts for the influence of the number of elements on the measure of 
cosDiag
cosDiag. This can help in comparing results when dealing with different data set sizes or varying contexts.

Correlation and Variation:
The term also impacts the interpretation of 
cosDiag
cosDiag in the context of correlation or variation of the data, as it adjusts the result relative to the number of data points.


At 0° (or 360°), the cosine value is 1, indicating maximum alignment, where vectors point in the same direction, resulting in perfect similarity. At 45°, the cosine is approximately 0.707, suggesting partial alignment and moderate similarity between the vectors. At 90°, the cosine value is 0, meaning the vectors are perpendicular, showing no similarity or projection between them.

At 180°, the cosine value is -1, indicating that the vectors are pointing in opposite directions, representing maximum dissimilarity. Similarly, at 270°, the cosine value is 0 again, indicating perpendicular vectors but with a reversed orientation. Finally, at 360° (or back to 0°), the cosine returns to 1, indicating perfect alignment once more after completing a full rotation.

without sqrt(2) 
As n increases (more data points or higher dimensions), the overall result becomes smaller because the denominator increases. For large n, this would reduce the impact of the cosine similarity value.

For smaller n, the scaling would have less of an impact, and the cosine similarity would remain closer to its unscaled value.

This reflects how larger datasets or higher-dimensional data may require additional normalization to account for complexity.