

## **Data Mining in R (DATR) exercises**


### **Overview**

- **Purpose**: Processing raw RNA-seq data into quantified gene expression values

- Ensures accurate gene mapping and expression quantification across species
- Enables reliable cross-species expression comparisons

### **Analysis Tools**

1. **GFFread**

   - **Purpose**: Creates reference transcriptomes
   - **Input**: GFF annotations + FASTA genome sequences
   - **Output**: Transcriptome FASTA file
   - **Function**: Combines genomic annotations with sequences for RNA-seq alignment

2. **Trimmomatic**

   - **Purpose**: Quality control of raw reads
   - **Input**: Raw FASTQ files
   - **Output**: Filtered FASTQ files
   - **Function**: 

        - Removes adapter sequences
        - Trims low-quality bases
        - Filters short reads
        - Ensures high-quality input for quantification

3. **Kallisto**

   - **Purpose**: Expression quantification
   - **Input**: Filtered reads + Reference transcriptome
   - **Output**: Gene expression values (TPM)
   - **Function**: 

        - Performs pseudoalignment of reads
        - Quantifies transcript abundance
        - Generates normalized expression values


### **Workflow Implementation**

- _Lab-book:_ 

- Sequential processing through Slurm scripts
- Automated pipeline execution on bioserver
- Standardized output formats for downstream analysis


