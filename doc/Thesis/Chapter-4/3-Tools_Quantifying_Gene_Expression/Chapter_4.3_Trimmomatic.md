
## **Step 2: Quality Control with Trimmomatic**

### **Overview**

- Tool for cleaning and filtering raw RNA-seq data
- Processes paired-end FASTQ files
- Improves data quality for downstream analysis

### Input

- **Raw FASTQ Files**

  - Forward reads (`*_1.fastq`)
  - Reverse reads (`*_2.fastq`)
  - Contains sequences and quality scores

### Usage

```bash
java -jar trimmomatic-0.39.jar PE -threads 4 \
  input_1.fastq input_2.fastq \
  output_1T.fq.gz output_1U.fq.gz \
  output_2T.fq.gz output_2U.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
  MINLEN:70
```

### Key Parameters

- **PE**: Paired-end mode
- **-threads**: Number of processing threads
- **ILLUMINACLIP**: Adapter removal settings
- **MINLEN**: Minimum read length threshold

### Processing Steps

1. **Adapter Removal**
   - Uses specified adapter sequences
   - Maintains read pairing when possible

2. **Quality Filtering**
   - Removes low-quality bases
   - Filters short reads (<70 bp)
   - Preserves paired-end information

3. **Output Generation**
   - Compressed FASTQ format
   - Separate files for paired reads
   - Discards unpaired reads

### Data Source

- Available through NCBI
- Project ID: PRJDB4481
- [Access Link](https://www.ncbi.nlm.nih.gov/Traces/study/)
