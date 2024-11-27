
# ----------------------------------------------------------------------
# euclidean Distances


```R
exp.prof.dists <- function(special_in_paralogs_v.lst,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), 
                                          c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            dist_matrix <- as.vector(dist(species_data, method = "euclidean"))
            
            return(dist_matrix)
        }
        return(NA)
    })
    
    return(distances)
}
```

Yes! Here's the mathematical formula for euclidean distances used in the code:

For two expression profiles A and B with n tissues, the euclidean distance is:

√(Σ(Ai - Bi)²) where i goes from 1 to n

In our specific case with 4 tissues (fat, gut, muscle, whole_body), it expands to:

√((A_fat - B_fat)² + (A_gut - B_gut)² + (A_muscle - B_muscle)² + (A_wholebody - B_wholebody)²)



# ----------------------------------------------------------------------
# euclidean log2 Distances


```R
exp.prof.dists_log2 <- function(gene.accessions,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, "Parent_FBgn", "Species")),
                          dist.method = "euclidean") {
    

    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            # abs (absolute) because negative values for log between 0 - 1 (log2(0.5) = -1)
            # logarythmus with base 2
            species_data <- sapply(species_data, as.numeric)
            dist_matrix <- abs(log2(as.vector(dist(species_data, method = dist.method))))
            
            return(dist_matrix)
        }
        return(NA)
    })
    
    return(distances)
}
```

# ----------------------------------------------------------------------
# cosine angles Distances


```R
exp.prof.angles <- function(special_in_paralogs_v.lst,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles), 
                                          c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    angles <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            
            # Calculate pairwise angles
            angle_matrix <- matrix(0, nrow=nrow(species_data), ncol=nrow(species_data))
            for(i in 1:(nrow(species_data)-1)) {
                for(j in (i+1):nrow(species_data)) {
                    # Cosine similarity
                    cos_sim <- sum(species_data[i,] * species_data[j,]) / 
                              (sqrt(sum(species_data[i,]^2)) * sqrt(sum(species_data[j,]^2)))
                    # Convert to angle in radians
                    angle_matrix[i,j] <- acos(cos_sim)
                }
            }
            return(as.vector(angle_matrix[upper.tri(angle_matrix)]))
        }
        return(NA)
    })
    
    return(angles)
}

distances_result <- exp.prof.angles(
    special_in_paralogs_v.lst = special_in_paralogs_v.lst,
    expression.profiles = rna.seq.exp.profils
)




exp.prof.dists.angles <- function(special_in_paralogs_v.lst,
                                  expression.profiles = rna.seq.exp.profils,
                                  expr.prof.gene.col = "FBpp_ID",
                                  tissues = setdiff(colnames(expression.profiles), 
                                                    c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    # Retrieve all genes
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            
            # Initialize an empty matrix for pairwise angles
            num_genes <- nrow(species_data)
            angle_matrix <- matrix(0, nrow = num_genes, ncol = num_genes, 
                                   dimnames = list(rownames(species_data), rownames(species_data)))
            
            # Compute pairwise angles
            for (i in 1:(num_genes - 1)) {
                for (j in (i + 1):num_genes) {
                    u <- as.numeric(species_data[i, ])
                    v <- as.numeric(species_data[j, ])
                    
                    # Calculate cosine similarity
                    cosine_similarity <- sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
                    
                    # Compute angle in radians and convert to degrees
                    angle <- acos(cosine_similarity) * (180 / pi)
                    
                    # Fill the matrix symmetrically
                    angle_matrix[i, j] <- angle
                    angle_matrix[j, i] <- angle
                }
            }
            
            return(as.vector(angle_matrix[upper.tri(angle_matrix)]))
        }
        return(NA)
    })
    
    return(distances)
}


distances_result <- exp.prof.dists.angles(
    special_in_paralogs_v.lst = special_in_paralogs_v.lst,
    expression.profiles = rna.seq.exp.profils
)








To calculate distances based on angles between vectors instead of Euclidean distances, you can use the **cosine similarity** or directly compute the **angles** using the arccosine function. Heres an approach for adapting your function to calculate distances using angles:

### Mathematical Concept:
The **cosine similarity** between two vectors \( u \) and \( v \) is:
\[
\text{cosine\_similarity}(u, v) = \frac{u \cdot v}{\|u\| \|v\|}
\]
Where:
- \( u \cdot v \) is the dot product of the two vectors.
- \( \|u\| \) and \( \|v\| \) are the magnitudes (norms) of the vectors.

The angle \( \theta \) between two vectors is derived as:
\[
\theta = \arccos\left(\frac{u \cdot v}{\|u\| \|v\|}\right)
\]
- Convert \( \theta \) from radians to degrees if needed.




```R
exp.prof.dists.angles <- function(special_in_paralogs_v.lst,
                                  expression.profiles = rna.seq.exp.profils,
                                  expr.prof.gene.col = "FBpp_ID",
                                  tissues = setdiff(colnames(expression.profiles), 
                                                    c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    # Retrieve all genes
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            
            # Initialize an empty matrix for pairwise angles
            num_genes <- nrow(species_data)
            angle_matrix <- matrix(0, nrow = num_genes, ncol = num_genes, 
                                   dimnames = list(rownames(species_data), rownames(species_data)))
            
            # Compute pairwise angles
            for (i in 1:(num_genes - 1)) {
                for (j in (i + 1):num_genes) {
                    u <- as.numeric(species_data[i, ])
                    v <- as.numeric(species_data[j, ])
                    
                    # Calculate cosine similarity
                    cosine_similarity <- sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
                    
                    # Compute angle in radians and convert to degrees
                    angle <- acos(cosine_similarity) * (180 / pi)
                    
                    # Fill the matrix symmetrically
                    angle_matrix[i, j] <- angle
                    angle_matrix[j, i] <- angle
                }
            }
            
            return(as.vector(angle_matrix[upper.tri(angle_matrix)]))
        }
        return(NA)
    })
    
    return(distances)
}
```

---

### Explanation of Changes:
1. **Numerical Conversion**:
   - Ensure all values in the tissue columns are numeric for vector computations.

2. **Angle Computation**:
   - For each pair of rows (genes), calculate the cosine similarity and then derive the angle using \( \arccos \).

3. **Pairwise Angles**:
   - Store the angles in a symmetric matrix (`angle_matrix`) and return the upper triangular values as the result.

---

### Output
The output will be a list of pairwise angles (in degrees) for each species group. If only one gene is present in a group, `NA` is returned.

---

### Advantages
- Angles provide a measure of the directional similarity between gene expression profiles, useful for studying correlations.
- Unlike Euclidean distances, they are scale-invariant and focus on vector orientation. 

You can run this function similarly to your original code but interpret the output in terms of angular distances rather than direct Euclidean measures.





If you only have single scalar values (for example, just one expression level per gene, instead of vectors), you can't directly compute the **cosine similarity** as it requires vectors with multiple dimensions (i.e., multiple expression values across tissues or conditions). However, you can think of it in a simplified way by comparing the **relative difference** between the two values. Cosine similarity is typically used for vectors (with multiple dimensions), but here's an approach to handle scalar values:

### Scalar Cosine Similarity

In the case of two scalar values \(x_1\) and \(x_2\), the **cosine similarity** would theoretically be 1 if both values are the same direction (positive or negative), and -1 if they are in exactly opposite directions, and 0 if they are orthogonal (independent). However, for scalar values, this doesn't directly apply unless you're considering the values as part of a broader feature.

To handle scalar comparisons, you might use **simple subtraction** or **relative distance** between the two values.

### Cosine-like measure for scalars

One approach to compare scalar values could be:

\[
\text{Cosine\_Similarity}(x_1, x_2) = \frac{x_1 \cdot x_2}{|x_1| \cdot |x_2|}
\]
This reduces to:
\[
\text{Cosine\_Similarity}(x_1, x_2) = \frac{x_1 \cdot x_2}{\sqrt{x_1^2} \cdot \sqrt{x_2^2}} = \frac{x_1 \cdot x_2}{|x_1| \cdot |x_2|}
\]
For scalar values, this is equivalent to calculating the **sign of the product**:

- If both \(x_1\) and \(x_2\) are positive or both negative, the result is **positive**.
- If one is positive and the other negative, the result is **negative**.

### Implementation of Scalar Cosine-like Measure

You can adjust your function to compute this kind of similarity:

```r
exp.prof.dists.scalar_angles <- function(special_in_paralogs_v.lst,
                                         expression.profiles = rna.seq.exp.profils,
                                         expr.prof.gene.col = "FBpp_ID",
                                         tissues = setdiff(colnames(expression.profiles), 
                                                           c(expr.prof.gene.col, "Parent_FBgn", "Species"))) {
    # Retrieve all genes
    all_genes <- unlist(special_in_paralogs_v.lst)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    distances <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Ensure numeric conversion
            species_data <- sapply(species_data, as.numeric)
            
            # Initialize an empty matrix for pairwise cosine similarity
            num_genes <- nrow(species_data)
            similarity_matrix <- matrix(0, nrow = num_genes, ncol = num_genes, 
                                       dimnames = list(rownames(species_data), rownames(species_data)))
            
            # Compute pairwise cosine-like similarities for scalars
            for (i in 1:(num_genes - 1)) {
                for (j in (i + 1):num_genes) {
                    x1 <- species_data[i, ]
                    x2 <- species_data[j, ]
                    
                    # Compute the cosine-like similarity
                    similarity <- sum(x1 * x2) / (sqrt(sum(x1^2)) * sqrt(sum(x2^2)))
                    
                    # Store the similarity in the matrix (angle similarity)
                    similarity_matrix[i, j] <- similarity
                    similarity_matrix[j, i] <- similarity
                }
            }
            
            return(as.vector(similarity_matrix[upper.tri(similarity_matrix)]))
        }
        return(NA)
    })
    
    return(distances)
}
```

distances_result <- exp.prof.dists.scalar_angles(
    special_in_paralogs_v.lst = special_in_paralogs_v.lst,
    expression.profiles = rna.seq.exp.profils
)



### Notes:
1. **Scalar values comparison**:
    - For scalar comparisons, we are computing a "cosine-like" measure that treats the values as vectors of length 1.
    - For each pair of scalar values (expression values), the **cosine-like similarity** compares how similarly aligned they are (based on their signs and magnitudes).
   
2. **Angle measure**:
    - You can then use this cosine-like measure to interpret how aligned (or opposite) two expression values are in your analysis.
   
3. **What if you want actual angles?**:
    - If you want to compute actual angles (in degrees or radians), you'd apply the `acos()` function to the cosine similarity:
      ```r
      angle <- acos(similarity) * (180 / pi)  # For degrees
      angle <- acos(similarity)  # For radians
      ```
    - This would give you the angle between the two values as a measure of similarity.

4. **Practical Use**:
    - This approach is most useful when you are comparing the "directionality" of gene expression in a simplified context (e.g., comparing how two genes in the same tissue are correlated).

### Summary:

If you're comparing scalar values (single gene expression values), you can use a **cosine-like** approach that treats each value as a vector of length 1. This measure compares how aligned the values are, and you can compute the **angle** between them by applying the arccosine function.