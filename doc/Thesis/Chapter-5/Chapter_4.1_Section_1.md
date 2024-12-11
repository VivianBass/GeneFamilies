




#### **Section 1: Data Loading**
# ------------------------------------------------------------------------------------------------

- **Scripts**:
  - `load_gene_expression_data.R`: Handles RPKM data loading.
  - `load_gene_groups_data.R`: Imports ortholog and paralog relationships.
  - `load_gene_families_data.R`: Processes MCL output for gene family clustering.


- **Input Data**:


- **Expression Profiles**
  - RNA-seq data across multiple tissues
  - Normalized expression values (TPM/RPKM)

- **Gene Groups**
  - Conserved orthologs
  - In/out paralogs
  - Special in/out paralogs

- **Gene Families**
  - MCL clustering output
  - Family-wise gene organization into Gene-Family Clusters
  

  - Formats and files: Expression Profiles (regular/log2), Gene Groups, Gene Families.
  - Link to the outputs from Chapter 3.
  - load Information about the Orthologs and Tandems identified within the eight Brassicaceaen Genomes
- the txt Files contain Information about Pairwise Sequence Similarities (Orthologs & Paralogs)

  - orthologs (genes in different species with a common ancestor), They typically retain 
  the same function across species. genes from multiple species. These genes are found in the same cluster 
  because they share a high degree of sequence similarity, which suggests that they are orthologs.


### 1. `load_gene_expression_data.R`

Loads expression data for gene groups, typically in DNA FASTA format. Expression levels (e.g., RNA counts, RPKM) reflect gene activity: higher counts indicate active genes (ON), while lower counts indicate inactive genes (OFF). Profiles are multi-dimensional vectors representing gene activity across tissues. Normalization methods (e.g., TPM, FPKM) adjust for technical variations, ensuring accurate comparisons between samples.

USAGE: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

### 2. `load_gene_groups_data.R`

Loads ortholog and paralog relationships across species, with files categorized by group type. Each file follows a specific header convention:

- **Conserved Orthologs**: 
  One ortholog per species-gene pair; represents highly conserved, reliable relationships.
- **In Paralogs**: 
  Species-specific paralogs from recent duplications with no conserved ortholog.
- **Special In Paralogs**: 
  Species-specific paralogs with a conserved ortholog relationship, indicating potential conservation.
- **Out Paralogs**: 
  Paralogs from older duplications, not tied to a conserved ortholog.
- **Special Out Paralogs**: 
  Cross-species paralogs with a conserved ortholog relationship, suggesting retained functional ties.



  #### `con_orthologs.tsv`

  | Family      | Gene         | Gene_species | Ortholog     | Ortholog_species |
  |-------------|--------------|--------------|--------------|------------------|
  | OG0000000   | FBpp0117097  | dana         | FBpp0172663  | dmoj             |

  #### `in_paralogs.tsv`

  | Family      | Gene         | Gene_species | Paralog      | Paralog_species  |
  |-------------|--------------|--------------|--------------|------------------|
  | OG0000000   | FBpp0117097  | dana         | FBpp0172663  | dmoj             |

  #### `special_in_paralogs.tsv`

  | Family      | Gene         | Gene_species | Paralog      | Paralog_species  |
  |-------------|--------------|--------------|--------------|------------------|
  | OG0000000   | FBpp0117097  | dana         | FBpp0172663  | dmoj             |

  #### `out_paralogs.tsv`

  | Family      | Gene         | Gene_species | Paralog      | Paralog_species  |
  |-------------|--------------|--------------|--------------|------------------|
  | OG0000000   | FBpp0117097  | dana         | FBpp0172663  | dmoj             |

  #### `special_out_paralogs.tsv`

  | Family      | Gene         | Gene_species | Paralog      | Paralog_species  |
  |-------------|--------------|--------------|--------------|------------------|
  | OG0000000   | FBpp0117097  | dana         | FBpp0172663  | dmoj             |


  #### 3. `load_gene_families_data.R`

  Loads gene family clusters created by MCL. Each family clusters related genes based on sequence similarity, organizing them by shared evolutionary lineage.

- Clusters: In the context of gene families, a cluster refers to a group of genes that are related 
  by sequence similarity and are assumed to have evolved from a common ancestor. 
  Clusters represent groups of related genes, identified by their sequence similarity across different species

- Inputs: `mcl_output.txt` , `mcl_table.tsv` (from markov clustering tool mcl) 
- Outputs: `families.lst` & `families.genes.df, families.df` --> `families.RData`


  USAGE: Rscript exec/load_gene_families_data.R <families_file> 

  #### `<families_file>`

  | Family   | species1          	 | species2            | species3            |
  |----------|---------------------|---------------------|---------------------|
  | family_1 | gene1, gene2, gene3 | gene4, gene5, gene6 | gene7, gene8, gene9 |


  USAGE: Rscript exec/load_gene_expression_data.R <RPKM_counts_table.tsv>

#### `<RPKM_counts_table.tsv>`

  | id    | tissue  | expression |
  |-------|---------|------------|
  | gene1 | tissue1 | ####       |
  | gene2 | tissue1 | ####       |
  | gene1 | tissue2 | ####       |
  | gene2 | tissue2 | ####       |

<br>

- we need the mapping Information from the Fasta files to create rna.seq.exp.profiles, also to map the species name to the right gene we used the fasta file

#### `<rna_seq_exp_profiles>`

  | FBpp_ID     | Parent_FBgn | Species | fat (tissue)| gut         | muscle      | whole_body |
  |-------------|-------------|---------|-------------|-------------|-------------|------------|
  | FBpp0291548 | FBgn0085506 | dmel    | 0.000000000 | 0.000000000 | 0.000000000 | 1.0000000  |
  | FBpp0289382 | FBgn0259817 | dmel    | 0.000000000 | 0.003827632 | 0.006253107 | 0.9899193  |
  | FBpp0312442 | FBgn0085692 | dmel    | 0.007196562 | 0.000000000 | 0.000000000 | 0.9928034  |
  | FBpp0077828 | FBgn0002121 | dmel    | 0.111506516 | 0.266536495 | 0.283942387 | 0.3380146  |
  | FBpp0111921 | FBgn0031209 | dmel    | 0.009854926 | 0.756561172 | 0.010681603 | 0.2229023  |




# ------------------------------------------------------------------------------------------------




## **RNA-Seq Data Normalization**

### Purpose

- Corrects technical biases in expression data
- Enables fair comparison across samples
- Accounts for sequencing depth variations

### Common Methods

1. **RPKM** (Reads Per Kilobase Million)
   - Normalizes by gene length and library size
   - Formula: `(reads × 10⁹)/(total_reads × gene_length)`
   - Best for single-end sequencing

2. **FPKM** (Fragments Per Kilobase Million)
   - Similar to RPKM but for paired-end reads
   - Accounts for fragment rather than read counts

3. **TPM** (Transcripts Per Million)
   - More consistent across samples
   - Sums to same value (1 million) for each sample

### Normalization Process

1. **Length Normalization**
   - Divides counts by gene length
   - Corrects for longer genes generating more reads

2. **Depth Normalization**
   - Scales to million mapped reads
   - Accounts for sequencing depth differences

### Benefits

- Comparable expression values
- Reduced technical bias
- More accurate differential expression analysis








- **Normalization**: Adjusts raw expression data using RPKM methodology to ensure comparability.

- `Normalization of the Expression` (Transcriptome) data 

- Normalization Purpose: Adjusts for technical variations (e.g., sequencing depth, RNA quality) 
to ensure fair comparison of gene expression levels across different samples.

- Normalization Methods: Common methods include TPM (Transcripts Per Million), 
FPKM (Fragments Per Kilobase of transcript per Million mapped reads), 
and RPKM (Reads Per Kilobase of transcript per Million mapped reads).

- RPKM Overview: RPKM is used in RNA sequencing to normalize gene expression by accounting 
for both gene length and total read depth, enabling comparison between different samples.

- Gene Length Normalization: RPKM divides read counts by the gene length (in kilobases), 
preventing longer genes from appearing more highly expressed simply due to their length.

- Read Depth Adjustment: To account for variations in total read counts between samples, 
RPKM scales read counts to a million mapped reads, 
ensuring fair comparison of gene expression levels across samples.

