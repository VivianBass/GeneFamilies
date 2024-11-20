
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

For **Drosophila sechellia**, the corresponding files for reference transcriptome generation are also available on FlyBase. 

...

### Quality Filtering

For quality filtering of FASTQ files, we utilize **Trimmomatic**. The datasets can be accessed through the following link:  
[NCBI FASTQ Quality Filtering - Trimmomatic](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&query_key=3&WebEnv=MCID_673684b340e26547c60a0f01&o=organism_s%3Aa%253Bacc_s%3Bacc_s%3Aa)


## **Data Mining in R (DATR) exercises**

In this data mining exercise, as outlined in **doc/DATR_6_Uebung.pdf**, we will perform bioinformatic analyses of gene expression data using powerful tools in the field of transcriptomics. The goal is to gain hands-on experience with processing, quantification, and analysis of RNA-Seq data to understand gene expression levels across different conditions.

We will work with tools like **Trimmomatic**, **gffread**, **Kallisto** which are commonly used in modern bioinformatics pipelines. By the end of this exercise, we will have practical skills in filtering raw sequencing data, creating reference transcriptomes, quantifying gene expression, and performing downstream statistical analysis.

Previously, we relied on gene expression data extracted from spreadsheets linked in the paper. However, we faced issues accurately mapping gene expression data to the respective species, as the expression profiles only referenced **dmel** gene names. This led to boxplots that represented only **dmel** genes and excluded other species present in the profiles. Although we attempted to map genes using symbols provided in the spreadsheets, this approach yielded incomplete results, as not all genes in the expression profiles could be mapped.

To resolve this, we decided to use raw data from the paper to ensure precise gene mapping and represent all genes from all species present in the expression profiles. Therefore, we are utilizing the RNA-seq workflow to efficiently transform raw sequencing reads into accurate gene-level expression data.

This RNA-seq workflow efficiently transforms raw sequencing reads into precise gene-level expression data, ensuring reliable results. It enables condition-specific comparisons, identifies differentially expressed genes, and provides insights into gene regulation.

---

### Step 1: Create the Reference Transcriptome with GFFread

- Extract annotated transcript sequences for pseudoalignment.

To create the `reference Transcriptome`, you need:

1. **Genomic Sequence (FASTA format)**: Contains nucleotide sequences of the entire genome.  
2. **Genomic Annotations (GFF format)**: Metadata about gene structures, specifying positions of exons, introns, CDS, and genes. 

- `gffread` combines the annotations (GFF) with the sequence (FASTA) to extract transcript sequences, 
   producing a FASTA file of all annotated transcripts.

```bash
    gffread /media/BioNAS/ag_hallab/DATR/material/reference.gff \
    -g /media/BioNAS/ag_hallab/DATR/material/reference.fsa \
    -w /media/BioNAS/ag_hallab/EasyVectorOmics/Tests_AK/results/output_transcriptome.fa

    # `-g`: Specifies the genome FASTA file.  
    # `-w`: Defines the output FASTA file for the transcriptome.
```

#### **What It Does:**



- which information do the gff and fasta gz exactly contain ?

---

### Step 2: Quality Filtering with Trimmomatic

- Use Trimmomatic to Filter the quality of the raw reads (`FASTQ files`)
- Trims raw paired-end FASTQ reads using quality thresholds and adapter files.
- removes adapter sequences, low-quality bases, and contaminants from raw paired-end reads (FASTQ files).  
- Adapter sequences are short, synthetic DNA sequences added to the ends of DNA fragments during library preparation for sequencing.
- Output: High-quality trimmed reads ready for quantification.

```bash
    # This runs Trimmomatic (`trimmomatic-0.39.jar`) using Java.
    # The tool is located in the specified path `/usr/local/bin/trimmomatic/`
    # PE Indicates **Paired-End mode**, meaning both forward (`_1.fastq`) and reverse (`_2.fastq`) reads are processed.
    # Allocates 4 threads for parallel processing to speed up the trimming process.
    java -jar /usr/local/bin/trimmomatic/trimmomatic-0.39.jar PE -threads 4 \
    # These are the raw paired-end input FASTQ files: `_1.fastq`: Forward reads. and `_2.fastq`: Reverse reads.
    /media/BioNAS/ag_hallab/DATR/material/control/SRR9929273_1.fastq \
    /media/BioNAS/ag_hallab/DATR/material/control/SRR9929273_2.fastq \
    # Trimmomatic produces separate output files: one for Forward reads after trimming and one for Reverse reads after trimming.
    # `/dev/null`: Specifies where discarded or unpaired reads are sent (not stored).
    # Output files: Trimmed reads are stored in compressed FASTQ format (`*.fq.gz`).
    <your_dir>/material/control/SRR9929273_1T.fq.gz /dev/null \
    <your_dir>/material/control/SRR9929273_2T.fq.gz /dev/null \
    # **ILLUMINACLIP**: Removes adapter sequences using the `TruSeq3-PE.fa` adapter file.
    # **`keepBothReads`**: Ensures both paired-end reads are retained, even if one is trimmed due to low quality.
    ILLUMINACLIP:/usr/local/bin/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    # Discards reads shorter than 70 bases after trimming.
    MINLEN:70
```

#### **What It Does:**
1. Removes adapter sequences using the specified adapter file.  
2. Filters low-quality reads and bases.  
3. Ensures both reads in a pair are retained if possible.  
4. Outputs high-quality, trimmed reads for further analysis, discarding reads shorter than 70 bases.

---

### Step 3: Quantify Gene Expression with Kallisto

Kallisto is an alignment-free tool for RNA-seq data analysis. 
It quantifies gene expression by pseudoaligning reads to a reference transcriptome, 
enabling fast and accurate abundance estimation.  
Estimate transcript abundance without traditional alignment. 

### **1. Index the reference Transcriptome**  

- Kallisto requires an indexed reference transcriptome to process RNA-seq reads.
- **Reference transcriptome:** Generated using GFFread (FASTA format). (Step-1)
- Create a transcriptome index from the reference FASTA.
 
```bash
kallisto index -i <your_dir>/results/transcriptome.idx <your_dir>/results/<gffread_output.fa>

# `kallisto index -i`: Creates a transcriptome index file.
# `<gffread_output.fa>`: Input transcriptome in FASTA format.
```

#### **What It Does:**


### **2. Quantify Gene Expression** 

- Kallisto quantifies transcript abundances from paired-end reads.
- **Trimmed FASTQ files:** Preprocessed reads from RNA-seq experiments. (Step-2)
- Quantify expression (TPM) for each transcript using pseudoalignment
- Output: Abundance estimates for each condition.

```bash
kallisto quant -i <your_dir>/results/transcriptome.idx \
-o results/<condition> -b 100 \
<your_dir>/material/control/SRR9929273_1T.fq.gz \
<your_dir>/material/control/SRR9929273_2T.fq.gz

# `kallisto quant -i`: Uses the pre-built index for pseudoalignment.  
# `-o results/<condition>`: Saves results in a folder named after the condition.  
# `-b 100`: Generates 100 bootstraps to estimate abundance variability.  
# `<SRR9929273_1T.fq.gz>` and `<SRR9929273_2T.fq.gz>`: Paired-end trimmed FASTQ files.  
```

#### **What It Does:**

### **Repeat**  
Execute the `kallisto quant` command for each condition and replicate to analyze all samples.









 

 
























sbatch Tests_AK/methods/kallisto.sh

```bash

#!/bin/bash
# SBATCH --job-name=kallisto
# SBATCH --output=kallisto_%j.out
	
echo "starting script"

kallisto index -i <your_dir>/results/transcriptome.idx <your_dir>/results/<gffread_output.fa>

kallisto quant -i <your_dir>/results/transcriptome.idx \
-o results/<condition> -b 100 \
<your_dir>/material/control/SRR9929273_1T.fq.gz \
<your_dir>/material/control/SRR9929273_2T.fq.gz

echo "script finished!"

```














/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietM




```bash
#!/bin/bash

# Base dir where dietM is stored
base_dir="/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietM"

# Dir to adapters 
adapter_file="/usr/local/bin/trimmomatic/adapters/TruSeq3-SE.fa"

# Loop dmel and dsec
for species in dsec ; do
    species_dir="$base_dir/$species"

    # Check if dir exists
    if [ -d "$species_dir" ]; then
        # Loop for tissues
        for tissue in fat muscle whole_body gut; do
            tissue_dir="$species_dir/$tissue"

            # Check if tissue dir exists
            if [ -d "$tissue_dir" ]; then
                # Create "trimmed" directory inside tissue dir
                output_dir="$tissue_dir/trimmed"
                mkdir -p "$output_dir"

                # Check all fastq files
                for fastq in "$tissue_dir"/*.fastq.gz; do
                    # Get filename whitout extension
                    base_name=$(basename "$fastq" .fastq.gz)

                    # Execute trimmomatic for every file
                    java -jar /usr/local/bin/trimmomatic/trimmomatic-0.39.jar SE -threads 4 \
                    "$fastq" \
                    "$output_dir/${base_name}_trimmed.fastq" \
                    ILLUMINACLIP:$adapter_file:2:30:10:2 \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70
                done
            else
                echo "Error: $tissue_dir doesn't exist."
            fi
        done
    else
        echo "Error: $species_dir doesn't exist."
    fi
done
```



