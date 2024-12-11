

## **Orthofinder**

### **Overview**

- Comprehensive tool for phylogenetic orthology inference
- Identifies orthogroups and gene families across multiple species
- Generates gene trees and species phylogeny


### **Process Flow**

1. **Input Preparation**

   - Protein FASTA files (one per species)
   - Each file contains all proteins for that species

2. **Analysis Steps**

   - All-vs-all sequence alignment (_BLAST_)
   - Gene tree inference
   - Species tree reconstruction
   - Ortholog inference

3. **Output Generation**

   - Orthogroups        (_Gene-families_)
   - Gene trees         (_Phylogenetic-trees_)
   - Species tree   
   - Ortholog relationships



### **Usage Example**

```bash
# Basic OrthoFinder run
orthofinder -f protein_sequences_directory

# Run with multiple threads
orthofinder -f protein_sequences_directory -t 8

# Using existing BLAST results
orthofinder -b previous_orthofinder_blast_results
````

### Parameters
- `-f`: Directory containing protein FASTA files
- `-t`: Number of threads (parallel processing)
- `-S`: Sequence search program (BLAST/DIAMOND)
- `-M`: MSA program (MAFFT/MUSCLE)

### Output Structure
- **Orthogroups/**
  - Orthogroups.tsv (main results table)
  - Orthogroups.GeneCount.tsv
  - Single_Copy_Orthogroups.txt
- **Species_Tree/**
  - SpeciesTree_rooted.txt
- **Gene_Trees/**
  - Individual gene family trees

Note: Results are organized in a directory named `Results_<date>` in your input directory
