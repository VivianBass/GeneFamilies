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

