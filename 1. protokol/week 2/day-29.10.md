

**Date**: 29.10.2024  
**Git Branch**: `Gene-Families-tests-Andre`

**Tasks**:

- got new to do 

        - we need another way to load the data 
        - we need to load pairs of orthologs, pairs of paralogs (not only the list itself)

        - you have to modify the way the scripts are working, because, as I told you, we need to work with pairs of genes instead of a list of genes.
        - We only need you to measure distances as we are measure them now and also with logarithm transformed. 

        We will have 5 groups:
        conserved orthologs
        in paralogs with orthologs
        in paralogs without orthologs
        out paralogs with orthologs
        out paralogs without orthologs

        So you have to be able to load all 5 files, families and expression values, and measure the distances between the pairs

- adjusted exec/3.load_gene_groups_data.R for the new 5 gene groups

- thought to separate the all vs all (blast) from loading gene groups to a separate rscript just for blast etc. 

- adjusted messages in loading rscripts




**Doubts and Issues**:

**Next Steps**:

---

**Code**