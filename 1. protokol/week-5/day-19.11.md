
**Date**: 19.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`  

---

### **Tasks**  

- **Assignment:** Perform the exercise from `doc/DATR_6_Uebung.pdf` (Data Mining in R).  
- **Purpose:** Process RNA-seq data to quantify gene expression and compare gene activity across conditions.  
- **Tools Used:** Trimmomatic, gffread, Kallisto.  
- **Steps Involved:**  

  1. Create Reference Transcriptome using GFFread.  
  2. Perform Quality Filtering with Trimmomatic.  
  3. Quantify Gene Expression using Kallisto & Bowtie2.  
  4. Extract TPM values.  

- **Data:**  
  - RNA-seq data: Forward and reverse reads (available on the server).  
  - Reference Transcriptome: Includes genomic sequence (FASTA) and annotations (GFF), obtained from SGD or FlyBase.  
  - Data for the exercise is available on the server (`/media/BioNAS/ag_hallab/DATR/material`).  

---

### **Directories Created**  

- Tests: `/media/BioNAS/ag_hallab/EasyVectorOmics/Tests_AK`  
  - Methods: `methods/`  
  - Results: `results/`  

---

### **Execution**  

- Verified the scripts provided in `doc/DATR_6_Uebung.pdf` are functional.  
- Successfully recreated the results and understood the functionality of the tools used.  

---

### **Current Work**  

- Recreating results using data from the Diet paper (*Drosophila*).  
- Need to gather complete data (previously only partial data was available) to proceed with the EasyVectorOmics package and related results.  
- Results will be stored in:  
  - `/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietM`  
  - `/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietP`  

- **Reference Transcriptome for *dmel* (D. melanogaster):**  
  - Using GFF and FASTA files from FlyBase (release FB2024_05):  
    - GFF: [FlyBase GFF](https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/gff/dmel-all-r6.60.gff.gz)  
    - FASTA: [FlyBase FASTA](https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz)  

---

### **Lab Book**  

- Documenting the process in `doc/Lab-book-EasyVectorOmics.md`.  

---

### **Next Steps**  

- Continue documenting in the Lab Book.  
- Gather data for *dmel* and *dsec* species from the Diet paper.  
- Fix t-tests and add Wilcox tests for analysis.  
- Create R script to calculate angles.  

---

### **Code**  

**Slurm Script for Trimming FASTQ Files with Trimmomatic:**  

```bash
#!/bin/bash

# Base directory for dietM
base_dir="/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietM"

# Adapter file
adapter_file="/usr/local/bin/trimmomatic/adapters/TruSeq3-SE.fa"

# Loop through species
for species in dsec; do
    species_dir="$base_dir/$species"

    if [ -d "$species_dir" ]; then
        # Loop through tissues
        for tissue in fat muscle whole_body gut; do
            tissue_dir="$species_dir/$tissue"

            if [ -d "$tissue_dir" ]; then
                # Create output directory for trimmed files
                output_dir="$tissue_dir/trimmed"
                mkdir -p "$output_dir"

                # Process each FASTQ file
                for fastq in "$tissue_dir"/*.fastq.gz; do
                    base_name=$(basename "$fastq" .fastq.gz)
                    java -jar /usr/local/bin/trimmomatic/trimmomatic-0.39.jar SE -threads 4 \
                        "$fastq" \
                        "$output_dir/${base_name}_trimmed.fastq" \
                        ILLUMINACLIP:$adapter_file:2:30:10:2 \
                        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70
                done
            else
                echo "Error: $tissue_dir does not exist."
            fi
        done
    else
        echo "Error: $species_dir does not exist."
    fi
done
```  