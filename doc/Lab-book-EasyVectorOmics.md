
# EasyVectorOmics Analysis - Lab Book

## **Project Overview**

The goal of this project is to reproduce the vector space analyses conducted in the _Cardamine hirsuta_ genome project. The R code used for the analysis will be refined and made applicable to any dataset. This code will then be used to analyze publicly available data from model species to evaluate the generalizability of the method to other datasets. [DOI: https://doi.org/10.1038/nplants.2016.167](https://doi.org/10.1038/nplants.2016.167)


## **Current Work**

1. **Defining Expression Vector Spaces**
   - To construct a viable expression vector space, we require data from at least three, preferably four, species. 
   - Additionally, we need expression measurements from at least four different tissues, as these tissues serve as the axes of the expression vector space.

2. **Recreating Results from the Paper**
   - We are currently replicating the results with data from the study:  
     **"Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species"**  
     [DOI: https://doi.org/10.1016/j.celrep.2019.08.011](https://doi.org/10.1016/j.celrep.2019.08.011)


## **Study Overview**

This paper investigates dietary adaptability in *Drosophila* species, comparing generalists such as *D. melanogaster* and *D. simulans* to specialists like *D. sechellia*. Key findings include:

- **Generalists** (*D. melanogaster* and *D. simulans*) possess carbohydrate-responsive gene regulation through the TGF-β/Activin signaling pathway, enabling adaptation to diets with varying nutrient balances and maintaining metabolic balance.

- **Specialists** (*D. sechellia*) lack such regulation, resulting in poor adaptation to carbohydrate-rich diets, metabolite accumulation, and reduced dietary flexibility. 

These findings underscore the role of gene-environment interactions in shaping nutritional adaptability and evolutionary differences between generalist and specialist species.


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

## **References**

1. Watanabe, Kaori, Yasutetsu Kanaoka, Shoko Mizutani, Hironobu Uchiyama, Shunsuke Yajima, Masayoshi Watada, Tadashi Uemura, and Yukako Hattori.  
   _"Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species."_  
   *Cell Reports* 28, no. 10 (2019): 2594–2607.e7.  
   [DOI: https://doi.org/10.1016/j.celrep.2019.08.011](https://doi.org/10.1016/j.celrep.2019.08.011)

## Data Availability

All RNA-sequencing data have been deposited in the **DDBJ Sequence Read Archive** under BioProject accession number **PRJDB4481**. 
The specific accession numbers for the datasets are:

- **DRA004295**
- **DRA006831**
- **DRA007810**

### Supplementary Materials

Additional data, including spreadsheets and expression profiles, can be downloaded from the publication's supporting information page:  
[Cell Reports - Supporting Materials](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31064-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124719310642%3Fshowall%3Dtrue)

### Reference Transcriptomes from FlyBase

To create reference transcriptomes, we use **GFFread**, a tool that requires both GFF and FASTA files.

For **Drosophila melanogaster**, the files required for reference transcriptome generation are sourced from **FlyBase release FB2024_05** (current as of November 19, 2024). These include:

- **GFF File**:  
  [dmel-all-r6.60.gff.gz](https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/gff/dmel-all-r6.60.gff.gz)

- **FASTA File**:  
  [dmel-all-chromosome-r6.60.fasta.gz](https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz)

These files are stored on the server at:  
`/media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dmel`

For **Drosophila sechellia**, the corresponding files for reference transcriptome generation are sourced from **FlyBase release FB2018_05**. These include:

- **GFF File**:  
  [dsec-all-r1.3.gff.gz](https://ftp.flybase.net/releases/FB2018_05/dsec_r1.3/gff/dsec-all-r1.3.gff.gz)

- **FASTA File**:  
  [dsec-all-chromosome-r1.3.fasta.gz](https://ftp.flybase.net/releases/FB2018_05/dsec_r1.3/fasta/dsec-all-chromosome-r1.3.fasta.gz)

These files are stored on the server at:  
`/media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dsec`


### Quality Filtering

For quality filtering of FASTQ files, we used **Trimmomatic**, a tool designed to trim low-quality bases and adapter sequences from RNA-Seq data. The datasets used in this analysis are publicly available and can be accessed through the following link:  
[NCBI FASTQ Quality Filtering - Trimmomatic](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&query_key=3&WebEnv=MCID_673684b340e26547c60a0f01&o=organism_s%3Aa%253Bacc_s%3Bacc_s%3Aa)

**BioProject**: [PRJDB4481](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB4481)

We analyzed data from two **Drosophila** species. **Drosophila sechellia (dsec)** and **Drosophila melanogaster (dmel)** across multiple tissue types. Each tissue was represented by **three biological replicates**, with the respective accession numbers listed below:

#### Accession Numbers for RNA-Seq Data

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

<br>

## **Data Mining in R (DATR) exercises**

In this data mining exercise, as outlined in **doc/DATR_6_Uebung.pdf**, we will perform bioinformatic analyses of gene expression data using powerful tools in the field of transcriptomics. The goal is to gain hands-on experience with processing, quantification, and analysis of RNA-Seq data to understand gene expression levels across different conditions.

We will work with tools like **Trimmomatic**, **gffread**, **Kallisto** which are commonly used in modern bioinformatics pipelines. By the end of this exercise, we will have practical skills in filtering raw sequencing data, creating reference transcriptomes, quantifying gene expression, and performing downstream statistical analysis.

Previously, we relied on gene expression data extracted from spreadsheets linked in the paper. However, we faced issues accurately mapping gene expression data to the respective species, as the expression profiles only referenced **dmel** gene names. This led to boxplots that represented only **dmel** genes and excluded other species present in the profiles. Although we attempted to map genes using symbols provided in the spreadsheets, this approach yielded incomplete results, as not all genes in the expression profiles could be mapped.

To resolve this, we decided to use raw data from the paper to ensure precise gene mapping and represent all genes from all species present in the expression profiles. Therefore, we are utilizing the RNA-seq workflow to efficiently transform raw sequencing reads into accurate gene-level expression data.

This RNA-seq workflow efficiently transforms raw sequencing reads into precise gene-level expression data, ensuring reliable results. It enables condition-specific comparisons, identifies differentially expressed genes, and provides insights into gene regulation.

<br>

## Step 1: Create the Reference Transcriptome with GFFread

A **reference transcriptome** is a collection of all transcript sequences (mRNA, RNA, etc.) for a given organism, derived from the genomic sequence and annotated gene structures. It provides a comprehensive map of known transcripts, essential for RNA-seq analysis to quantify gene expression and identify differentially expressed genes.

To create the **reference transcriptome**, you need:

1. **Genomic Sequence (FASTA format)**: Contains nucleotide sequences of the entire genome.  
2. **Genomic Annotations (GFF format)**: Metadata about gene structures, specifying positions of exons, introns, coding sequences (CDS), and genes.

```bash
gffread /media/BioNAS/ag_hallab/DATR/material/reference.gff \
-g /media/BioNAS/ag_hallab/DATR/material/reference.fsa \
-w /media/BioNAS/ag_hallab/EasyVectorOmics/Tests_AK/results/output_transcriptome.fa

# `-g`: Specifies the genome FASTA file.
# `-w`: Defines the output FASTA file for the transcriptome.
```

### **What It Does:**
- Combines the genomic annotations (GFF) and the genomic sequence (FASTA) to extract the transcript sequences for each annotated gene.  
- Outputs a FASTA file containing all transcript sequences for downstream analyses, such as alignment-free quantification, gene quantification, and differential expression analysis.  

<br>

### ⚠️ **Important Note**

- Ensure the GFF file is unzipped before running the `gffread` command, as compressed files are not supported. 

```bash
gunzip /media/BioNAS/ag_hallab/EasyVectorOmics/material/references/<species>/<gfffile.gff.gz>
```

- Clean and adjust FASTA files for both Drosophila melanogaster (dmel) and Drosophila sechellia (dsec) to ensure compatibility: 

1. **Remove Line Breaks**:

```bash
zcat /media/BioNAS/ag_hallab/EasyVectorOmics/material/references/<species>/<fastafile.fasta.gz> | awk 'NF' > cleaned.fasta
```

2. **Adjust Headers**:  
- Simplify headers to retain only the chromosome names. 

```bash
awk '/^>/ {print $1; next} {print}' cleaned.fasta > cleaned_final.fasta
```
- Perform these steps for both species and ensure the cleaned files are saved in their respective directories.

<br>

## Step 2: Quality Filtering with Trimmomatic

- **Trimmomatic** is used to filter the quality of raw paired-end reads (`FASTQ` files), improving data quality before downstream analysis.

- **FASTQ Files** are text-based files that store raw sequencing data. They contain nucleotide sequences (reads) along with quality scores for each base, indicating the confidence of the sequencing results.

```bash
    # This runs Trimmomatic (`trimmomatic-0.39.jar`) using Java.
    # The tool is located in the specified path `/usr/local/bin/trimmomatic/`
    # PE indicates **Paired-End mode**, meaning both forward (`_1.fastq`) and reverse (`_2.fastq`) reads are processed.
    # Allocates 4 threads for parallel processing to speed up the trimming process.
    java -jar /usr/local/bin/trimmomatic/trimmomatic-0.39.jar PE -threads 4 \
    # Raw paired-end input FASTQ files: `_1.fastq` is the forward read, `_2.fastq` is the reverse read.
    /media/BioNAS/ag_hallab/DATR/material/control/SRR9929273_1.fastq \
    /media/BioNAS/ag_hallab/DATR/material/control/SRR9929273_2.fastq \
    # Trimmomatic produces separate output files for forward and reverse reads after trimming.
    # `/dev/null` is used for discarded or unpaired reads (not saved).
    # Output files: Trimmed reads stored in compressed FASTQ format (`*.fq.gz`).
    <your_dir>/material/control/SRR9929273_1T.fq.gz /dev/null \
    <your_dir>/material/control/SRR9929273_2T.fq.gz /dev/null \
    # **ILLUMINACLIP**: Removes adapter sequences using the specified adapter file.
    # `keepBothReads`: Retains both paired-end reads, even if one is trimmed due to low quality.
    ILLUMINACLIP:/usr/local/bin/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    # Discards reads shorter than 70 bases after trimming.
    MINLEN:70
```
### **What It Does:**
- **Removes adapter sequences** (using the specified adapter file, TruSeq3-PE.fa).
- **Filters low-quality reads** and bases below the defined quality threshold.
- **Retains paired-end reads** if possible, even if one is trimmed due to low quality.
- **Outputs high-quality trimmed reads**, discarding any shorter than 70 bases, in compressed FASTQ format, ready for further analysis.

- Trimmomatic trims raw paired-end FASTQ reads by removing adapter sequences, low-quality bases, and contaminants, ensuring high-quality reads for downstream quantification.

<br>

## Step 3: Quantify Gene Expression with Kallisto

- Kallisto is an efficient, alignment-free tool for RNA-seq analysis. It quantifies gene expression by pseudoaligning reads to a reference transcriptome, enabling fast and accurate transcript abundance estimation.

### **1. Index the Reference Transcriptome**

```bash
kallisto index -i <your_dir>/results/transcriptome.idx <your_dir>/results/<gffread_output.fa>

# Input: FASTA file generated using GFFread (Step 1).
# Output: An indexed transcriptome file used for pseudoalignment.
```

### **What It Does:**  
- Creates an index file (`transcriptome.idx`) from the reference transcriptome.

### **2. Quantify Gene Expression**

```bash
kallisto quant -i <your_dir>/results/transcriptome.idx \
-o results/<condition> -b 100 \
<your_dir>/material/control/SRR9929273_1T.fq.gz \
<your_dir>/material/control/SRR9929273_2T.fq.gz

# `-i`: Specifies the indexed transcriptome file created in the previous step.
# `-o`: Defines the output folder for results, named based on the condition (e.g., control, treated).
# `-b 100`: Runs 100 bootstraps to estimate variability in abundance estimates.
# Input: Paired-end trimmed FASTQ files (`_1T.fq.gz` and `_2T.fq.gz`) from Step 2.
# Output: TPM values and variability metrics for each transcript.
```

### **What It Does:**  
- Pseudoaligns trimmed reads to the indexed transcriptome.  
- Estimates transcript abundances (TPM) and calculates variability using bootstraps.  
- Outputs results to a folder specific to the condition (e.g., control, treated).  

**Note**: Run `kallisto quant` for all conditions and replicates to analyze the complete dataset.


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






# ---------------------------------------------------------------------








 

 
## Step 1: Create the Reference Transcriptome with GFFread
## Step 2: Quality Filtering with Trimmomatic
## Step 3: Quantify Gene Expression with Kallisto
## Step 4: Extract TPM values

We ran kallisto with 
kallisto index -i dmel_transcriptome.idx dmel_transcriptome.fa
 
But we detected some transcriptome duplicates, so we deleted everything after the first blank with
awk '/^>/{print $1; next} {print}' dmel_transcriptome.fa > dmel_transcriptome_clean.fa

and then we checked for duplicates:
grep "^>" dmel_transcriptome_clean.fa | sort | uniq -d

we got 2 duplicated transcriptomes, so we deleted the duplicates using:
awk '/^>/{if(seen[$0]++) next} {print}' dmel_transcriptome_clean.fa > dmel_transcriptome_deduplicated.fa

 
Then we could execute
kallisto index -i dmel_transcriptome.idx dmel_transcriptome_deduplicated.fa

we did that only for dmel because we didnt have problems with dsec for kallisto index

also please add in the documentation that we had to compress every trimmomatic result because the result files are so big
 