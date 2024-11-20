
# EasyVectorOmics Analysis - Lab-book

In this project the goal is to reproduce the vector space analyses as done in
the _Cardamine hirsuta_ genome project [1]. The R-Code required to do so shall
be isolated and made usabel with _any_ data. Next, this code shall be used to
analyse public open access data of model species to evaluate whether the method
works for other data-sets as well. https://doi.org/10.1038/nplants.2016.167


# References

Watanabe, Kaori, Yasutetsu Kanaoka, Shoko Mizutani, Hironobu Uchiyama, Shunsuke Yajima, Masayoshi Watada, Tadashi Uemura, and Yukako Hattori. "Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species." Cell Reports 28, no. 10 (2019): 2594-2607.e7. [https://doi.org/10.1016/j.celrep.2019.08.011](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31064-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124719310642%3Fshowall%3Dtrue)

`Interspecies Comparative Analyses Reveal Distinct Carbohydrate-Responsive Systems among Drosophila Species`

This paper investigates dietary adaptability in *Drosophila* species, comparing generalists like *D. melanogaster* and *D. simulans* to specialists such as *D. sechellia*. Generalists can adapt to diets with varied nutrient balances due to carbohydrate-responsive gene regulation via TGF-β/Activin signaling, maintaining metabolic balance. Specialists lack this regulation, struggling with carbohydrate-rich diets, leading to metabolite accumulation and reduced adaptation. The findings highlight how gene-environment interactions shape nutritional flexibility, contributing to evolutionary differences between generalist and specialist species.


# ---------------------------------------------------------------------------------


# Data


- species under investigation: *D. melanogaster* and *D. simulans* to specialists such as *D. sechellia*
- Conditions: 3 Diets (M, P, ...)
- tissues: Muscle, Fat_body, whole body, ....

DATA AND CODE AVAILABILITY

All the RNA-sequencing data have been deposited and are available in the `DDBJ Sequence Read Archive.` 
BioProject accession number: PRJDB4481
The accession numbers for the data are DDBJ: DRA004295, DRA006831, and DRA007810 






### **Data Mining in R (DATR) exercises**

This workflow processes RNA-seq data to quantify gene expression, 
enabling the comparison of gene activity across different conditions. 
It ensures high-quality data and accurate expression quantification, 

- **Efficiency:** Streamlines RNA-seq data analysis with alignment-free quantification, saving time and computational resources.  
- **Accuracy:** Transforms raw sequencing reads into precise gene-level expression data, ensuring reliable results.  
- **Relevance:** Facilitates the identification of differentially expressed genes, aiding in the understanding of condition-specific gene regulation.


### Data 

https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&query_key=3&WebEnv=MCID_673684b340e26547c60a0f01&o=organism_s%3Aa%253Bacc_s%3Bacc_s%3Aa


#### -----------------------------------------------------------------------------------------

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

- dmel reference Transcriptome
for melanogaster we will use gff and fasta file from flybase, release FB2024_05, the current release by November 19 2024
https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/gff/dmel-all-r6.60.gff.gz
https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
 
they are stored in the server in /media/BioNAS/ag_hallab/EasyVectorOmics/material/references/dmel


- dsec reference Transcriptome




- which information do the gff and fasta gz exactly contain ?

#### -----------------------------------------------------------------------------------------

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

### **What It Does:**
1. Removes adapter sequences using the specified adapter file.  
2. Filters low-quality reads and bases.  
3. Ensures both reads in a pair are retained if possible.  
4. Outputs high-quality, trimmed reads for further analysis, discarding reads shorter than 70 bases.


#### -----------------------------------------------------------------------------------------

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

### **Repeat**  
Execute the `kallisto quant` command for each condition and replicate to analyze all samples.

#### -----------------------------------------------------------------------------------------







 

 
























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



