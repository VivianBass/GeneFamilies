
## **Expression Vector Space**:

- **Objective**

  - To quantitatively measure and compare gene expression patterns across various conditions, revealing functional diversification patterns through vector space analysis.


### **Gene Expression**  

- **Gene Activity Levels**

  - **High Expression**: Represents actively transcribed genes with substantial RNA production
  - **Low Expression**: Indicates repressed or inactive genes with minimal RNA production


### **Distance Measurements**

  - Gene expression profiles are compared using Euclidean distances across:
    - Complete expression profiles
    - Specific dimensions (e.g., individual tissue types)
  
  - **Distance Interpretation**:

    - **Small Distances**: Indicate similar gene activity patterns, suggesting functional conservation
    - **Large Distances**: Reveal significant variations, suggesting functional divergence
    

### **Dimensionality**

  - Each tissue type/condition represents one dimension in the vector space
  - Total dimensions equal the number of tissues/conditions analyzed

- **Vector Construction**

  - Individual gene expression patterns form vectors in multi-dimensional space
  - Expression values are normalized (RPKM/TPM) for accurate comparisons


### **Expression Change Analysis**

  - Changes are quantified by calculating distances between mean expression vectors
  - Provides robust measurement of gene activity variation across conditions

- **Biological Interpretation**

  - Increased distances suggest functional divergence
  - Vector similarity indicates potential conservation of regulatory mechanisms
  - Vector orientation changes reveal tissue-specific adaptations
