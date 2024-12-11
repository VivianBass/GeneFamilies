

## **Step 1: Reference Transcriptome Creation with GFFread**

### **Purpose**

- Creates a comprehensive collection of transcript sequences (mRNA, RNA)
- Essential for accurate RNA-seq quantification and differential expression analysis
- Combines genomic sequence with annotation information

### **Required Input Files**

1. **Genomic Sequence** (FASTA format)
   - Contains complete genome nucleotide sequences
   - Example: `/path/to/reference.fsa`

2. **Genomic Annotations** (GFF format)
   - Defines gene structures and features
   - Includes exons, introns, CDS locations
   - Example: `/path/to/reference.gff`


### **GFFread Command**

```bash
    gffread /path/to/reference.gff \
    -g /path/to/reference.fsa \
    -w /path/to/output_transcriptome.fa
```


### **Parameters**

- `-g`: Input genome FASTA file
- `-w`: Output transcriptome FASTA file
- Additional options available for specific feature extraction


### **Output**

- FASTA file containing transcript sequences
- Used for:

  - Alignment-free quantification
  - Gene expression analysis
  - Differential expression studies


# -------------------------------------------------------------------------------------


### ⚠️ **File Preparation Guidelines**


### 1. **GFF File Preparation**

- **Decompress GFF Files**

```bash
# Unzip GFF file (required for gffread)
gunzip /path/to/references/<species>/<gfffile.gff.gz>
```

- Ensure the GFF file is unzipped before running the `gffread` command, as compressed files are not supported. 


### **FASTA File Processing**

- **Step 1: Remove Empty Lines**

```bash
# Extract and clean FASTA content
zcat /path/to/references/<species>/<fastafile.fasta.gz> | \
awk 'NF' > cleaned.fasta
```

- Clean and adjust FASTA files for both Drosophila melanogaster (dmel) and Drosophila sechellia (dsec) to ensure compatibility 


- **Step 2: Simplify Headers**

```bash
# Keep only chromosome identifiers
awk '/^>/ {print $1; next} {print}' cleaned.fasta > cleaned_final.fasta
```


### Required for Both Species

- Process files for:

  - *D. melanogaster* (dmel)
  - *D. sechellia* (dsec)

- Store cleaned files in respective species directories
- Verify file integrity after processing


