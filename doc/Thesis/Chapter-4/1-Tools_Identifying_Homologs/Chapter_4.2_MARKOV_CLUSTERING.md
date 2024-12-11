

## **Markov Clustering (MCL)**


### **Overview**

- Tool for clustering biological sequences based on similarity scores
- Version: MCL 14-137 (van Dongen)
- Purpose: Groups genes into families based on sequence similarities


### **Process Flow**

1. **Input Preparation**

   - Run all-vs-all BLAST search
   - Create similarity matrix from BLAST results

2. **Clustering**

   - Apply MCL algorithm to similarity matrix
   - Group sequences based on connection patterns

3. **Output Processing**

   - Convert results to table format
   - Assign cluster names to gene families



### **Usage Example**

```bash
# 1. Convert BLAST output to MCL input format
cut -f 1,2,11 blast_results.txt > mcl_input.txt

# 2. Run MCL clustering
mcl mcl_input.txt --abc -I 2.0 -o mcl_output.txt

# 3. Process results
mcl2orthoformat.pl mcl_output.txt > mcl_table.tsv
```


### **Parameters**

- `-I`: _Inflation value_ (default: 2.0)

  - Higher values = more granular clusters
  - Lower values = larger clusters

- `--abc`: Input format specification for similarity scores



### **Output Format**

  ```m
  | Cluster_ID | Species1 | Species2 | Species3 |
  |------------|----------|----------|----------|
  | Cluster1   | 2        | 1        | 3        |
  | Cluster2   | 1        | 2        | 1        |

  ```

`Note:` Numbers indicate gene count per species in each cluster
