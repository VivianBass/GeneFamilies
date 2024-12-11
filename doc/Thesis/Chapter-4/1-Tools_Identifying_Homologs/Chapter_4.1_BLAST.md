
## **BLAST (Basic Local Alignment Search Tool) / Diamond** 

- **BLAST (Basic Local Alignment Search Tool)**

  - Performs local alignments by identifying matching regions between query sequences and database sequences
  - Calculates statistical significance and similarity scores for sequence alignments
  - Used to identify homologous sequences (orthologs, paralogs) through DNA, RNA, or protein similarity searches

    - _blastn_ (a tool from the BLAST suite for nucleotide sequence alignment) 

    - _blastp_ (a tool from the BLAST suite for protein sequence alignment) 


- **install BLAST:**

    ```R
    sudo apt-get install ncbi-blast+		(On Linux)
    ```



## **BLAST Usage Guide**

- Example for _all-vs-all alignment_
- Compares each sequence in a FASTA file against all other sequences in the same file
- Used to identify homologous relationships within a set of sequences 


### **1. Creating BLAST Database**

- **Input FASTA File**: DNA/protein sequences in FASTA format

    - FASTA files must be formatted with makeblastdb before use in local BLAST searches
    - You can directly use the available FASTA files (CDS) of the genes against each other. 
    - with the same FASTA file as both the query and the database. 


- **For Nucleotide Sequences**

```bash
# Format FASTA file into nucleotide database
makeblastdb -in family_genes.fasta -dbtype nucl -out family_db -parse_seqids
```

- **For Protein Sequences**

```bash
# Format FASTA file into protein database
makeblastdb -in protein_sequences.fasta -dbtype prot -out protein_db -parse_seqids
```

# --------------------------------------------------------------------------------

### **2. Running BLAST Search**

```bash
# Basic BLAST command structure
blastn -query <query_file> -db <database> -out <output_file> -outfmt 6

# Example for all-vs-all alignment
blastn -query family_genes.fasta -db family_db -out all_vs_all_results.txt -outfmt 6
```

- **Output Format** (`-outfmt 6`)


# --------------------------------------------------------------------------------

### **Key Parameters**

- **`-query <query_file>`**: Specifies the file containing nucleotide sequences in FASTA format
- **`-db <database>`**: nucleotide database for searching (can be pre-downloaded or custom-built)
- **`-out <output_file>`**: Names the output file for saving results
- **Output Format** (`-outfmt 6`)


### **Parameters**

```bash
# -evalue <value>: For more stringent matches, use a lower e-value,
blastn -query query.fasta -db nt -out output.txt -evalue 0.001

# -max_target_seqs: Limits the number of target sequences to show (default: `500`)
blastn -query query.fasta -db nt -out output.txt -max_target_seqs 10
```































