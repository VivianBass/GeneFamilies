
## Dataset / Input files / Data / Output files


## Rscripts: 
-----------------------------------------------------------------------------------------

- generate




# `Section-1` - Loading

- Gene families
- Coding Sequences
- All vs all file (Blast / Diamond)
- Orthologs / Paralogs / Tandems

6. All vs all file , Use BLAST to generate this file.
    For Drosohpila we used aminoacid sequences to run Blast.


- 1.  `loadRpkmRnaSeqCounts.R`		-> Load gene expressions.	                                     
- 2.  `loadcodingSequences.R`       -> 
- 1.  `LoadOrthologsAndTandems.R`
- 3.  `loadGeneFamilies.R`

Input: RPKM_counts_table.tsv
Input: eight_brassicaceae_orthologs.txt 
Input: eight_brassicaceae_tandems.txt

Input: mcl_output.txt
Input: mcl_table.tsv
families.tsv
counts.txt




1.  Rscript path/2/GeneFamilies/exec/loadRpkmRnaSeqCounts.R RPKM_counts_table.tsv"
-   Rscript path/2/GeneFamilies/exec/loadRpkmRnaSeqCounts.R RPKM_counts_table.tsv"
- Rscript exec/loadRpkmRnaSeqCounts.R experiments/RPKM_flybase/RPKM.tsv

2.  Rscript exec/loadOrthologsAndTandems.R <all_vs_all_file> <orthologs_file> <paralogs_file>

3.  Rscript exec/loadGeneFamilies.R <families_file> <counts_file>,
-   Rscript path/2/GeneFamilies/exec/loadGeneFamilies.R path/2/GeneFamilies/data mcl_output.txt mcl_table.tsv"
-   Rscript exec/loadGeneFamilies.R experiments/RPKM_flybase/families.tsv experiments/RPKM_flybase/counts.txt





# `Section-2` - computing

- 1.  `generateRpkmExpressionProfiles.R`		                                       
- 4.  `computeExpressionProfileDistances.R` 	 
`investigateDistributionsOfExpressionProfileDistances.R`


6. Compute statistics of expression profile distances. In GeneFamilies directory execute:
    ```
    Rscript exec/investigateDistributionsOfExpressionProfileDistances.R <input_file>



5. Compute expression profile distances. In GeneFamilies directory execute:
    ```
    Rscript exec/computeExpressionProfileDistances.R
    ```

# `Section-3` - Plotting

	                                  
- 2.  `plotDistributionsOfExpressionProfileDistances.R`	
- 2.  `expressionAngles.R`	 















## Usage




- plots
- angles









- as for the CDS or Fasta file? for genome sequences, 
-> we only need coding sequences (CDS) if and only if we need to estimate expression levels from raw reads ourselves. 


- as for expression profiles
- two factor data (ST-Exp) (Experiment) -> min 4 family members (species) & multiple tissues
- to have enough axes (species) to form a viable expression vector space?
- Example:  D.sec_brain, D.mel_brain, D.sec_muscle, D.mel_muscle, etc .....
- therefor the expression data / tables should account for this 
- `what exactly should ST-exp contain ?`

- for now just ort
- Orthologs / Paralogs / Tandems / 
- as for orthologs: Ideally, we have several sets of data (Orthofinder/ Flybase (1 Hit))
- `How to distinguish between Orthologs and Paralogs? if a gene could be both?`
-> orthologs and paralogs are somewhat context specific terms
- Search for information about tandem duplicates.
- Additionally, you can write a short script that parses the GFF or GTF table for D. melanogaster Each gene identifier is provided its genomic coordinates in that file.
- If genes are in close neighborhood and are paralogs, you can conclude them to be tandem duplicates. Allow for a neighborhood of maybe five genes. Check the descriptions for GTF and GFF formats: https://www.ensembl.org/info/website/upload/gff.html
- 
- `how and where to find the tandems?`

-> could make mathematical set-intersection, to find the orthologs and paralogs of a gene family.
- `intersection? in context of mapping genes to the right gene-family or distinguish between orth / para ?`
-> see how orthologs and paralogs are cathegorized in Cardamine

- Use the FlyBase paralog file for verification; but based on OrthoFinder, all non orthologs should be considered paralogs.
- ???

- Additionally, df with the following information:  Combined-Factor-Abbrev  | Species               | Tissue
                                                    D.sec_brain             | Drosophila sechellia  | brain



- Gene-Family-File 
- gene identifier for _D. melanogaster_. (names of genes are basically the same for the Drosophila species)
- `naming convention of Drosophila genes?` different drosophila species but they have the same gene name idetifiers. ?
- Load ortholog and paralog information (genes) together with gene families. merge Information
- also provide information which genes are orthologs and which are paralogs. Are you a ortholog or paralog?
- ortholog and paralog groups ? and To which gene family do you belong?








    - Use all genes classified as orthologs and analyse our gene family data:
    - See, which gene families have more than one ortholog? Plot the number of
          orthologs per gene family in a boxplot or raincloud plot, i.e. the
          distribution of number of orthologs per gene family.

 Ideally, we have several sets of data:
    - A list of gene families, consisting of gene family ID, and genes
      belonging to the respective family
    - A list of orthogroups, same format, i.e. group-ID and genes belonging to
      the respective orthogroup

 What we need to do is provide _any means_ of obtaining, storing, and using
    this information for e.g. Drosophila species. 
    
    - The easiest I can think of is
    to have another table in the same format as our gene families file, in
    which the ortholog-groups and paralog-groups are defined, e.g.
    `Group-ID <TAB> Group-Type (ortholog|paralog|tandem cluster) <TAB> Gene-ID-1,Gene-ID-2,Gene-ID-3,...,Gene-ID-N`
    Use any format that works for you. You can also use separate tables for
    ortholog and paralog groups.



- as for Plots (Distribution plots / Boxplot)
- Include jitter points in all boxplots
- I need to know which plots used which distance measures:
- per gene family: distances between orthologs, distances between paralogs
- which plots exactly: ??  just the means, and or the medians etc. ??
- `should we use the ortholog/paralog groups or genefamily clusters? for printing the distribution plots?`

- as for angles
- - for the angles to the diagonal: measure the angle to the diag for each gene, not only the angle between the mean expression vector of e.g. orthologs (or paralogs) and the diag
(currently Angles_boxplot-with-pvalues.pdf)