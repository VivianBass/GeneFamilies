
## 3.2 **Genome Assembly and Annotation**

- Genome assembly and annotation form the critical foundation for downstream analyses
- High-quality assembly and annotation are essential for accurate bioinformatics studies

### **Assembly Pipeline**

- **DNA Sequencing Process**


- `RNA Sequencing for Gene Expressions`

- RNA-Seq for Gene Expression: RNA sequencing (RNA-Seq) measures gene expression 
by sequencing the complete set of RNA transcripts (the transcriptome) in a sample, 
allowing comparison of expression levels across different conditions.

. Process: RNA is converted to complementary DNA (cDNA) through reverse transcription 
and then sequenced using high-throughput technologies. 
The resulting sequences are mapped to a reference genome.

. Gene Expression Counts: The number of reads mapped to each gene indicates 
its expression level, reflecting how often the gene was transcribed in the sample.


- Dataset Structure: Each row represents a gene's expression vector across different plant tissues, 
  with axes for tissue types (cotyledon, developing leaf, seedling, flower stage 16, flower stage 9)            X
  and columns for gene ID and expression levels.

- Expression Measurement: Expression levels reflect gene activity, 
  indicating how much of the gene’s DNA is transcribed into RNA, varying by tissue.

- Gene Activity: Expression counts show whether genes are 
  active (ON), inactive (OFF), or repressed in each tissue.


  

- **Shotgun Sequencing & Modern Technologies**
  - Random DNA fragmentation and PCR amplification
  - Utilizes platforms like Illumina and PacBio
  - Generates millions of short DNA sequence reads
  - Enables comprehensive genome coverage
  - Foundation for complete genome assembly

- **Assembly Process**

  - **Contig Assembly**
    - Uses SOAPdenovo for read alignment
    - Merges overlapping reads into contigs
    - Creates continuous DNA sequences

  - **Superscaffold Construction**
    - BAMLINK algorithm links contigs
    - Uses Bayesian framework
    - Produces larger genome segments

  - **Scaffold Anchoring**
    - Utilizes 8,249 bacterial artificial chromosomes (BAC)
    - Ensures accurate orientation
    - Validates assembly quality

### **Genome Annotation**

  - Structural annotation
  - Functional annotation
  - Gene prediction







