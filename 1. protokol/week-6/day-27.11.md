
**Date**: 27.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

### **Tasks**:

- We tested different distance methods: Euclidean, Angular, with and without log transformation.

### **1. Distance Calculations:**

#### **Euclidean Distance**
- A standard method for comparing gene vectors:
  \[
  d = \sqrt{\sum (gene1[i] - gene2[i])^2}
  \]

#### **Log₂ Transformation**
- Applied after distance calculations to normalize values:
  - Avoided negative outputs by adding a constant:  
    \[
    \log_2(\text{value} + 1)
    \]

**Transformation Effects:**

| **Original Range** | **Transformed Range**          | **Key Features**                     |
| ------------------ | ------------------------------ | ------------------------------------ |
| 0 < x < 1          | Negative values (e.g., -1, -2)  | Smaller inputs result in larger negatives. |
| 1 ≤ x < 2          | 0 to 1                         | Inputs closer to 2 approach 1.       |
| 2 ≤ x ≤ 10         | 1 to ~3.32                     | Larger values yield progressively higher outputs. |

**Summary:**
- Values between 0 and 1 shrink and become negative.
- Values between 1 and 2 map to a small positive range (0–1).
- Values greater than 2 scale to progressively larger positive values.

We discussed applying the log₂ transformation **before** calculating the distances. Specifically, RNA-seq expression profiles should be log₂-transformed directly, and the transformed values should then be used for distance calculations.

---

### **2. Expression Profiles**
- Normalized RPKM counts using log₂ transformation:
  ```R
  rpkm.rna.seq.counts_log2 <- rpkm.rna.seq.counts %>%
      mutate(expression = log2(expression + 1))
  ```
- Ensured the log₂-transformed values avoided negative outputs.

---

### **3. Cosine Similarity and Angular Distance**

Gene expression profiles are organized by tissues, e.g.:

| **gene_id** | **tissue1** | **tissue2** | **tissue3** | **tissue4** |  
|-------------|-------------|-------------|-------------|-------------|  
| fbPP        | 0.9         | 3.4         | 5.4         | 0.2         |  

Example vector for a gene: `[0.9, 3.4, 5.4, 0.2]`.  

#### **Cosine Similarity**  
- Measures the directional alignment of vectors.  
- Formula:  
  \[
  \text{Cosine Distance}(u, v) = 1 - \frac{u \cdot v}{\|u\| \|v\|}
  \]
- **Cosine similarity** captures the directional alignment of vectors, independent of magnitude.
- **Angular distance** (arccosine of cosine similarity) measures divergence between vectors.
- Useful for analyzing relative directionality in gene expression profiles.

#### **Angular Distance (Cosine Angle)**  
1. **Steps**:
   - Compute the **dot product**:
     \[
     \text{dot\_product} = \sum(gene1[i] \times gene2[i])
     \]
   - Calculate **magnitudes**:
     \[
     \text{magnitude\_gene1} = \sqrt{\sum(gene1[i]^2)} \quad \text{and} \quad \text{magnitude\_gene2} = \sqrt{\sum(gene2[i]^2)}
     \]
   - Find **cosine similarity**:
     \[
     \text{cosine} = \frac{\text{dot\_product}}{\text{magnitude\_gene1} \times \text{magnitude\_gene2}}
     \]
   - Convert cosine similarity to **angular distance**:
     \[
     \theta = \arccos(\text{cosine})
     \]

2. **Implementation**:
   ```R
   # Example
   gene1 <- c(23.19, 12.34, 8.45)
   gene2 <- c(2.50, 9.87, 11.23)
   
   dot_product <- sum(gene1 * gene2)
   magnitude_gene1 <- sqrt(sum(gene1^2))
   magnitude_gene2 <- sqrt(sum(gene2^2))
   
   cosine_angle <- dot_product / (magnitude_gene1 * magnitude_gene2)
   ```

---

### **4. Distance Computation Blueprint**

Distances were calculated for both raw and log-transformed data using the following methods:

#### **Euclidean Distances:**
- **All Distances**
- **Mean**
- **Median**

#### **Angular Distances:**
- **All Distances**
- **Mean**
- **Median**

For each method, both raw values and log-transformed values were used for the calculations.

---

### **5. Statistical Analysis**  
- Conducted t-tests and Wilcoxon tests for all distance methods.  
- Created separate plots to visualize the results for each statistical method.  

---

---

### **Doubts/Issues**:

**Gene Families**:
   - How should distance methods be applied to gene families?

---

### **Next Steps**:

1. Confirm the following:
   - Order of log transformation (before or after distance calculation).
   - Required statistical analyses for cosine/angular distances.
2. Refine plots to clearly compare methods (Euclidean vs. cosine/angular).
3. Test cosine and angular methods on normalized and log-transformed profiles.

---

### **Code Overview**:

1. **Cosine Angle Calculation**:
   - Computes the angle between two vectors:
     ```R
     calculate_cosine_angle <- function(vec1, vec2) {
         dot_product <- sum(vec1 * vec2)
         magnitude1 <- sqrt(sum(vec1^2))
         magnitude2 <- sqrt(sum(vec2^2))
         
         # Handle zero vectors
         if (magnitude1 == 0 || magnitude2 == 0) return(NA)
         
         cosine <- dot_product / (magnitude1 * magnitude2)
         cosine <- min(max(cosine, -1), 1)  # Correct numerical precision
         return(acos(cosine))
     }
     ```

2. **Angular Distance Profiles**:
   - Calculates pairwise angles for gene expression across tissues:
     ```R
     exp.prof.angles <- function(gene.accessions,
                          expression.profiles = rna.seq.exp.profils,
                          expr.prof.gene.col = "FBpp_ID",
                          tissues = setdiff(colnames(expression.profiles),
                                          c(expr.prof.gene.col, "Species"))) {
    
    all_genes <- unlist(gene.accessions)
    exp.profs <- as.data.frame(expression.profiles[expression.profiles[[expr.prof.gene.col]] %in% all_genes, ])
    species_groups <- split(exp.profs, exp.profs$Species)
    
    angles <- lapply(species_groups, function(species_data) {
        if (nrow(species_data) > 1) {
            rownames(species_data) <- species_data[[expr.prof.gene.col]]
            species_data <- species_data[, tissues]
            
            # Convert to numeric matrix
            expr_matrix <- sapply(species_data, as.numeric)
            
            # Calculate pairwise angles
            n <- nrow(expr_matrix)
            angle_vector <- numeric()
            
            for (i in 1:(n-1)) {
                for (j in (i+1):n) {
                    angle <- calculate_cosine_angle(expr_matrix[i,], expr_matrix[j,])
                    angle_vector <- c(angle_vector, angle)
                }
            }
            
            return(angle_vector)
        }
        return(NA)
    })
    
    return(angles)	
	}
     ```


3. **Log2 Transformation**:
   - Normalizes RPKM values, avoiding negative outputs:
     ```R
     rpkm.rna.seq.counts_log2 <- rpkm.rna.seq.counts %>%
         mutate(expression = log2(expression + 1))
     ```








