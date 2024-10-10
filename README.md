# GeneFamilies

This is a branch with the latest tests of EasyVectorOmics.

See:

Gan, X. et al. The Cardamine hirsuta genome offers insight into the evolution of morphological diversity. Nature Plants 2, 16167 (2016).

## Install the package

### Requirements

**You need HyPhy version 2.2.1**

## Install

You need to install the required package `AHRD.on.gene.clusters` as follows:
```
require(devtools)
install_github("asishallab/AHRD_on_gene_clusters")
```

Now you can install GeneFamilies:
```
install_github("asishallab/GeneFamilies")
```

## Folders

- doc: general comunication of the project (last version is on the server)
- exec: R scripts to execute analysis
- experiments: data of flybase and orthofinder in their respective directory. RPKM_* directories are the experiments.
- inputs: old directory of inputs, now we use experiments folder.
- man: roxygen documentation
- R: helper functions for analysis
- vignettes: documentation about old and new analysis
- .env file: R scripts use this file to know where to put plot results and RData objects.

## Usage

1. In -env file you can find the data output directory (RData objects) and results directory (plots). Modify this file as your convenience.

2. Load gene expressions. In GeneFamilies directory execute:
    ```
    Rscript exec/loadRpkmRnaSeqCounts.R <expression_file>, for example:
    
    Rscript exec/loadRpkmRnaSeqCounts.R experiments/RPKM_flybase/RPKM.tsv
    ```

3. Load gene families. In GeneFamilies directory execute:
    ```
    Rscript exec/loadGeneFamilies.R <families_file> <counts_file>, for example:
    
    Rscript exec/loadGeneFamilies.R experiments/RPKM_flybase/families.tsv experiments/RPKM_flybase/counts.txt
    ```

4. Load orthologs and paralogs. In GeneFamilies directory execute:
    ```
    Rscript exec/loadOrthologsAndTandems.R <all_vs_all_file> <orthologs_file> <paralogs_file>, for example:
    
    Rscript exec/loadOrthologsAndTandems.R inputs/all_vs_all_results.txt experiments/RPKM_flybase/orthologs.csv experiments/RPKM_flybase/paralogs.csv
    ```

5. Compute expression profile distances. In GeneFamilies directory execute:
    ```
    Rscript exec/computeExpressionProfileDistances.R
    ```

6. Compute statistics of expression profile distances. In GeneFamilies directory execute:
    ```
    Rscript exec/investigateDistributionsOfExpressionProfileDistances.R <input_file>
    ```
