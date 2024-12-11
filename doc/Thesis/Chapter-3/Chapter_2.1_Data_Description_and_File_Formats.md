
## **Data Description and File Formats**


- **FASTA Format** (.fa, .fasta)

  ```R
  >sequence_name
  ATGCTAGCTAGCTAGCTGATCGATCG
  GCTAGCTAGCTAGCTGATCGATCGAT
  ```
  - Used for storing nucleotide or protein sequences in plain text format

- Sourced from: _Flybase_




- **FASTQ Format** (.fastq, .fq)

  ```R
  @SRR001666.1
  GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  +
  !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
  ```

- Stores raw sequencing reads with quality scores for each base; essential input for RNA-seq analysis




- **GFF Format** (General Feature Format)

  ```R
  chr1  source  gene  1000  2000  .  +  .  ID=gene001;Name=BRCA1
  ```
  - Stores genomic features and their locations, including genes, exons, and regulatory elements

- Sourced from: _Flybase_, _Ensembl_




- **Kallisto Output**

  ```R
  target_id  length  eff_length  est_counts  tpm
  gene001    1200    1000        450         23.5
  ```
  - Contains RNA-seq quantification results, including transcript abundances and TPM values




- **Newick Format**

  ```R
  (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
  ```
  - Represents phylogenetic trees in a text-based format, showing evolutionary relationships









