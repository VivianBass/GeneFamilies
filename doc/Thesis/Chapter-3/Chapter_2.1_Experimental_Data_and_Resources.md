
## **Experimental Data and Resources**

_"Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species."_

### **Drosophila as Model Organism**

  - Family: *Drosophilidae* (fruit flies)
  - Genome: 13,600 protein-coding genes, 180M base pairs
  - Chromosomes: 2n = 8 (3 autosome pairs, 1 sex chromosome pair)

- **Research Value**

  - Well-annotated genome
  - Conserved genetic pathways shared with humans
  - Important model for genetics, development, and evolution studies



### **Experimental Conditions**

- **Species Comparison**

  - **Generalist**: *D. melanogaster* (dmel) - dietary flexibility
  - **Specialist**: *D. sechellia* (dsec) - Morinda fruit specific

- **Dietary Treatments**

  - M Diet: Medium composition (balanced control)
  - P Diet: Protein-enriched formulation
  - C Diet: Carbohydrate-enriched formulation

- **Tissue Samples Analyzed**

  - **Muscle tissue**: Metabolically active tissue
  - **Fat body**: Primary metabolic and storage organ
  - **Gut**: Main site of nutrient absorption
  - **Whole body**: Complete organismal response


# ------------------------------------------------------------------------------------------------


## **Data Resources**

### **1. RNA Sequencing Data**

We analyzed data from two **Drosophila** species. **Drosophila sechellia (dsec)** and **Drosophila melanogaster (dmel)** across multiple tissue types. Each tissue was represented by **three biological replicates**, with the respective accession numbers listed below:

All RNA-sequencing data have been deposited in the **DDBJ Sequence Read Archive** under BioProject accession number **PRJDB4481**. 

- **Repository**: DDBJ Sequence Read Archive
- **BioProject**: [PRJDB4481](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB4481)
- **Dataset Accessions**: DRA004295, DRA006831, DRA007810
- **Data Format**: FASTQ


| **Species**               | **Tissue**     | **Replicate 1** | **Replicate 2** | **Replicate 3** |
|---------------------------|----------------|-----------------|-----------------|-----------------|
| **D. sechellia**          | Whole Body     | DRR051498       | DRR051500       | DRR051502       |
|                           | Gut            | DRR129481       | DRR129483       | DRR129485       |
|                           | Fat Body       | DRR129493       | DRR129495       | DRR129497       |
|                           | Muscle         | DRR129535       | DRR129537       | DRR129539       |
| **D. melanogaster**       | Whole Body     | DRR051486       | DRR051488       | DRR051490       |
|                           | Gut            | DRR129475       | DRR129477       | DRR129479       |
|                           | Fat Body       | DRR129487       | DRR129489       | DRR129491       |
|                           | Muscle         | DRR129511       | DRR129513       | DRR129515       |



### **2. Reference Genome Resources**

To create _"reference Transcriptomes"_ for each species, we use **GFFread**, a tool that requires both GFF and FASTA files.

#### Drosophila melanogaster (FlyBase FB2024_05)
- **GFF Annotation File**: dmel-all-r6.60.gff.gz (Flybase)
- **Genome FASTA File**: dmel-all-chromosome-r6.60.fasta.gz (Flybase)
- **Local Path (Server)**: `/media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dmel`

#### Drosophila sechellia (FlyBase FB2018_05)
- **GFF Annotation File**: dsec-all-r1.3.gff.gz (Flybase)
- **Genome FASTA File**: dsec-all-chromosome-r1.3.fasta.gz (Flybase)
- **Local Path (Server)**: `/media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dsec`

These files are stored on the server at:  
`/media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dsec`




### **Supplementary Resources**

- Additional data, including spreadsheets and expression profiles, can be downloaded from the publication's supporting information page:  

- Comprehensive supplementary data available through Cell Reports:
- [Supporting Materials](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31064-2)

- Expression profiles
- Data spreadsheets
- Analytical workflows







# ------------------------------------------------------------------------------------------------

## **References**

1. Watanabe, Kaori, Yasutetsu Kanaoka, Shoko Mizutani, Hironobu Uchiyama, Shunsuke Yajima, Masayoshi  
   Watada, Tadashi Uemura, and Yukako Hattori.  
   _"Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species."_*Cell Reports* 28, no. 10 (2019): 2594–2607.e7.  
- [DOI: https://doi.org/10.1016/j.celrep.2019.08.011](https://doi.org/10.1016/j.celrep.2019.08.011)
- [Supporting Materials](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31064-2)








