
## **Step 3: Gene Expression Quantification with Kallisto**

### **Overview**

- Fast, alignment-free RNA-seq quantification tool
- Uses pseudoalignment for transcript abundance estimation
- Provides accurate expression measurements in TPM (Transcripts Per Million)


### **1. Index Creation**

```bash
# Build Kallisto index from reference transcriptome
kallisto index \
  -i results/transcriptome.idx \
  results/gffread_output.fa
```
- **Input**: Reference transcriptome (FASTA)
- **Output**: Indexed transcriptome file (.idx)
- **Purpose**: Enables rapid pseudoalignment



### **2. Expression Quantification**

```bash
# Quantify expression for paired-end reads
kallisto quant \
  -i results/transcriptome.idx \
  -o results/condition_name \
  -b 100 \
  reads_1T.fq.gz \
  reads_2T.fq.gz
```

#### Parameters

- `-i`: Path to index file
- `-o`: Output directory
- `-b`: Bootstrap samples (default: 100)
- Input files: Paired-end trimmed FASTQ files


#### Output

- Transcript-level abundances (TPM)
- Bootstrap estimates for uncertainty
- Organized by experimental condition


### **Important Notes**

- Run for all conditions/replicates
- Maintain consistent naming conventions
- Monitor computational resources
- Verify input file quality



# -------------------------------------------------------------------------------------


## Step 4: Extract TPM Values

### **Description**  

Kallisto provides TPM (Transcripts Per Million) values directly in the output file `abundance.tsv`. 
These values represent transcript-level quantification. For gene-level analysis, TPM values must be aggregated across all transcripts belonging to the same gene.
-> basically we take the outpu file from kallisto and extract the TPM values and aggregate them to the gene level

This step focuses on extracting and processing TPM values from Kallisto output to calculate overall gene expression levels.

---

### **Code**

```R
# Load required libraries
library(dplyr)

# Function to aggregate TPM values at the gene level
aggregate_TPM <- function(kallisto_dir, mapping_file, output_file = "gene_level_tpm.txt") {
  
  # Step 1: Load transcript-to-gene mapping file
  # Mapping file should have two columns: Transcript_ID and Gene_ID
  mapping <- read.table(mapping_file, header = TRUE, sep = "\t")
  
  # Step 2: Load Kallisto output
  abundance <- read.table(file.path(kallisto_dir, "abundance.tsv"), header = TRUE, sep = "\t")
  
  # Step 3: Merge TPM values with transcript-to-gene mapping
  abundance <- merge(abundance, mapping, by.x = "target_id", by.y = "Transcript_ID")
  
  # Step 4: Aggregate TPM values at the gene level
  gene_tpm <- abundance %>%
    group_by(Gene_ID) %>%
    summarize(TPM = sum(TPM, na.rm = TRUE))
  
  # Step 5: Save the aggregated results
  write.table(gene_tpm, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(gene_tpm) # Return the gene-level TPM values
}

# Example: Generate gene-level TPM values
aggregate_TPM("path/to/kallisto/output", "path/to/transcript_to_gene_map.txt", "results/gene_level_tpm.txt")
```

---

### **What It Does**  

1. **Loads Mapping File:**  
   - A transcript-to-gene mapping file links each transcript to its corresponding gene.  

2. **Reads Kallisto Output:**  
   - Parses the `abundance.tsv` file, which contains transcript-level TPM values.  

3. **Aggregates TPM Values:**  
   - Sums TPM values for all transcripts belonging to the same gene to calculate gene-level TPM.  

4. **Saves Results:**  
   - Outputs a file (`gene_level_tpm.txt`) with gene-level TPM values for downstream analysis.


