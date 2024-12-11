
## **Chapter 5: Experiments and Observations**

### 5.1 **Overview of Experiments**

- **Objective**:

  - Outline the goals of the experiments conducted using the R package (e.g., validating the package’s functionality, analyzing specific gene families).

- **Experiments Performed**:

  - List the key experiments:

    - Analysis of **Conserved Orthologs** vs. **In Paralogs**.
    - Comparison of **Special Out Paralogs** and other groups.
    - Statistical analysis using T-tests and Wilcoxon tests.
    - Generating distance and angle distributions.

### 5.2 **Experiment Setup**

- **Input Data**:
  - Specify the datasets used (e.g., gene groups, expression profiles, processed data from the lab-book).

- **Implementation Steps**:
  - Describe the workflow followed for each experiment using the EasyVectorOmics package.

- **Tools Used**:
  - Reiterate the tools (e.g., GFFread, Trimmomatic, Kallisto) and their role in preparing data for the experiments.



## **Data Used in Analysis**

- **Species Under Investigation**  
  - Generalists: *D. melanogaster* (Drosophila) (dmel)
  - Specialists: *D. sechellia* (Drosophila) (dsec)

- **Conditions**  
  - Diets: M (medium), P (protein-rich), C (carbohydrate-rich)

- **Tissues:**  
  - Muscle  
  - Fat Body  
  - Whole Body  
  - Gut  




- **Experiment: test_diet_P_&_M_filtered**

  - RNA-seq data from two Drosophila species (melanogaster, sechellia, simulans)

  - Four tissue types with three biological replicates each
  - 2 diet conditions (medium, protein-rich, carbohydrate-rich)

  - `Quantifikation of Gene-Expression`
  - Quality filtering of raw RNA-seq data
  - Reference transcriptome generation