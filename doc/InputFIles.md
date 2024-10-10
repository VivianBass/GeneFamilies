# Input files for Easy Vector Omics

This file contains all the descriptions for every input file that you need to run EasyVectorOmics package.

1. Gene expression table with genes of your XX species.

    #### Doubt:
    For Drosophila we only have the gene expression table for drosophila melanogaster. If we use this table, when we calculate the distances for orthologs, we got empty objects. This file was originally taken from Flybase.

    /media/BioNAS/ag_hallab/EasyVectorOmics/GeneFamilies/inputs/drosophila_file_RPKM_modify.tsv

    #### Response:

    - GeneFamilies can have various orthologs _of the same_ species. This is
      confusing, I know, but orthology and paralogy are context dependent
      terms. 
    - Use all genes classified as orthologs and analyse our gene family data:
        - See, which gene families have more than one ortholog? Plot the number of
          orthologs per gene family in a boxplot or raincloud plot, i.e. the
          distribution of number of orthologs per gene family.


2. Gene families

    This file contains genes belonging to a family per row. Every row are distinct families and has the structure:

    family \t gene_1 \t gene_2 \t gene_3

    #### Drosophila
    For Drosophila we used the file gene_group_data_fb_2024_03_curated.tsv file on the server and group by family to get the desired result.
    This file was originally taken from Flybase.

    /media/BioNAS/ag_hallab/prot-scriber/Flybase/material/gene_group_data_fb_2024_03_curated.tsv

    Here we have the family and all the associated genes, but we only have genes for malanogaster, so for the counting file should we relate somehow with the orthologs? Is that correct?

    Or should we use a tool to group the families and have the genes for the other species?

    #### Response:

    The family file is correct. The gene family file currently only uses gene
    identifier for _D. melanogaster_. The format is understood correctly. We
    have to associate this with another file (information source) about which
    genes are orthologs and which are paralogs. In the original `Gene_Families`
    project this was stored as binary RData (see [this file in
    Gene_Families](https://github.com/asishallab-group/GeneFamilies/blob/master/exec/loadOrthologsAndTandems.R)).
    The RData file `orthologsTandems.RData` contains the information which
    genes are orthologs and which are paralogs. Investigate:
    ```r
    load('./data/orthologsTandems.RData')
    ls()
    ```
    What we need to do is provide _any means_ of obtaining, storing, and using
    this information for e.g. Drosophila species. The easiest I can think of is
    to have another table in the same format as our gene families file, in
    which the ortholog-groups and paralog-groups are defined, e.g.
    `Group-ID <TAB> Group-Type (ortholog|paralog|tandem cluster) <TAB> Gene-ID-1,Gene-ID-2,Gene-ID-3,...,Gene-ID-N`
    Use any format that works for you. You can also use separate tables for
    ortholog and paralog groups.


3. Coding Sequences

    #### Drosophila
    For Drosophila we used gene coding sequences:
    >FBgn0012114 type=gene; loc=scaffold_13248:complement(3147334..3152273); ID=FBgn0012114; name=Dana\B-H1; dbxref=INTERPRO:IPR017970,FlyBase:FBgn0012122,GB:AY710106,GB:AY710107,GB:AY710108,GB:AY710109,GB:AY710110,GB:AY710111,GB:AY710112,GB:AY710113,GB:AY710114,GB:AY710115,GB:AY710116,GB:AY710117,GB:AY710118,GB:AY710119,GB:AY710120,GB:AY710121,GB:AY710122,GB:AY710123,GB:AY710124,GB:AY710125,GB:AY710126,GB:AY710127,GB:AY710128,GB:AY710129,GB:AY710130,GB:AY710131,GB:AY710132,GB:AY710133,GB:AY710134,GB:AY710135,GB:AY710136,GB:AY710137,GB:AY710138,GB:AY710139,GB:AY710140,GB:AY710141,GB:AY710142,GB:AY710143,GB:AY710144,GB:AY710145,GB:AY710146,GB:AY710147,GB:AY710148,GB:AY710149,GB:AY710150,GB:AY710151,GB:AY710152,GB:AY710153,GB:AY710154,GB:AY710155,GB:AY710156,GB:AY710157,GB:AY710158,GB:AY710159,GB:AY710160,GB:AY710161,GB:AY710162,GB:AY710163,GB:AY710164,GB:AY710165,GB:AY710166,GB:AY710167,GB:AY710168,GB:AY710169,GB:AY710170,GB:AY710171,GB:AY710172,GB:AY710173,GB:AY710174,GB:AY710175,GB:AY710176,GB:AY710177,GB:AY710178,GB:AY710179,GB:AY710180,GB:AY710181,GB:AY710182,GB:AY710183,GB:AY710184,GB:AY710185,GB:AY710186,GB:AY710187,GB:AY710188,GB:AY710189,GB:AY710190,GB:AY710191,GB:AY710192,GB:AY710193,GB:AY710194,GB:AY710195,GB:AY710196,GB:AY710197,GB:AY710198,GB:AY710199,GB:AY710200,GB:AY710201,GB:AY710202,GB:AY710203,GB:AY710204,GB:AY710205,GB:AY710206,GB:AY710207,GB:AY710208,GB:AY710209,GB:AY710210,GB:AY710211,GB:AY710212,GB:AY710213,GB:AY710214,GB:AY710215,GB:AY710216,GB:AY710217,GB:AY710218,GB:AY710219,GB:AY710220,GB:AY710221,GB:AY710222,GB:AY710223,GB:AY710224,GB:AY710225,GB:AY710226,GB:AY710227,GB:AY710228,GB:AY710229,GB:AY710230,GB:AY710231,GB:AY710232,GB:AY710233,GB:AY710234,GB:AY710235,GB:AY710236,GB:AY710237,GB:AY710238,GB:AY710239,GB:AY710240,GB:AY710241,GB:AY710242,GB:AY710243,GB:AY710244,GB:AY710245,GB:AY710246,GB:AY710247,GB:AY710248,GB:AY710249,GB:AY710250,GB:AY710251,GB:AY710252,GB:AY710253,GB:AY710254,GB:AY710255,GB:AY710256,GB:AY710257,GB:AY710258,GB:AY710259,GB:AY710260,GB:AY710261,GB:AY710262,GB:AY710263,GB:AY710264,GB:AY710265,GB:AY710266,GB:M59962,GB:M59963,GB_protein:AAA28381,GB:X56682,GB_protein:CAA40011,INTERPRO:IPR001356,INTERPRO:IPR009057,FlyBase_Annotation_IDs:GF22335,FlyBase:FBgn0012114,GLEANR:dana_GLEANR_6309,GNOMON:Dana_gnomon_100_gene.8394992,UniProt/TrEMBL:B3MW33,FlyBase_Annotation_IDs:B-H1,EntrezGene:6492915,GB_protein:EDV35178,UniProt/Swiss-Prot:P22544,GB_protein:EDV35178,INTERPRO:IPR020467,INTERPRO:IPR020479; MD5=2cc3a46a51101d0edb1b4cfa137449df; length=4940; release=r1.06; species=Dana; 

    Because it has the gene id and every file that we have has the gene ID, but we don't know if we should use CDS:

    >Dana\B-H1-PB type=CDS; loc=scaffold_13248:complement(join(3147441..3147703,3147867..3148204,3148325..3148553,3151062..3151989)); name=Dana\B-H1-RB; dbxref=FlyBase_Annotation_IDs:GF22335-PB,FlyBase:FBpp0348285,GNOMON:Dana_gnomon_100_7139227.p,REFSEQ:XP_001965663,GB_protein:EDV35178; MD5=e1cdbb8b9357b48ab1c2e0f9cbaacf60; length=1758; parent=FBgn0012114,FBtr0388591; release=r1.06; species=Dana; 

    If we use this type of file, we can't cross the information of our genes on script reconstructFamilyPhylogenes.R
    Maybe we should use this one but we have to map the parent genes somehow? Is that even correct?

    #### Response

    In the context of Easy Vector Omics, we only need coding sequences if and
    only if we need to estimate expression levels from raw reads ourselves. And
    in this case, we'd do that "outside" of the EasyVectorOmics R-Package;
    preparing our expression tables with e.g. Salmon or Bowtie2 ...

4. Orthologs

    A txt file that contains n columns according to your number of species and the best scored ortholog gene for the analized gene.

    #### Drosophila
    For Drosophila, we look for orthologs in Flybase and we found in previous releases a list of orthologs for melanogaster, we processed that file and used it, but we only have one ortholog gene for specie.

    We also ran Orthofinder and we got all ortholog genes belonging to an orthogroup.

    Should we use the one from flybase that only has one gene or use multiple genes for each specie?

    #### Response

    See 'Response' for point two, please. We need to merge the information of genes, i.e. for each gene we need data to answer the questions:
    - Are you a ortholog or paralog?
    - To which gene family to you belong?

    Ideally, we have several sets of data:
    - A list of gene families, consisting of gene family ID, and genes
      belonging to the respective family
    - A list of orthogroups, same format, i.e. group-ID and genes belonging to
      the respective orthogroup
    - The same for paralogs. Consequenctly by mathematical set-intersection, we
      find the orthologs and paralogs of a gene family.

    It is _very_ good, that we have the OrthoFinder results. If point two gives
    us not enough orthologs per family, we can check if OrthoFinder helps us
    out here.

5. Paralogs

    A txt file with two columns: Family and Gene.
    Each row should have just one gene belonging to that family.

    #### Drosophila
    For Drosophila, we used the dmel_paralogs file in Flybase (https://ftp.flybase.net/releases/current/precomputed_files/orthologs/dmel_paralogs_fb_2024_04.tsv.gz), change the name of the gene for paralog_cluster_XX, and deleted other columns we didn't need.

6. All vs all file
    
    Use BLAST to generate this file.

    #### Drosophila
    For Drosohpila we used aminoacid sequences to run Blast.

## Scripts

Ignore all scripts that deal with the detection of positive selection, like
FUBAR, MEME, pairwise Ka/Ks, etc. These look for signals of Darwinian selection
which we at this point do not want to include in our EasyVectorOmics analysis. 

- generateFubarBatchFiles.R 

    Maybe its a version issue but we had to modify the input file because hyphy was not working with the previous one.

    Also, we are getting a json output file from hyhpy and we have to convert that file into a table. We used this data in the json:

    ```
        "MLE":{
    "content":{
        "0":      [
    [11.34130318038266, 0.2232132768088211, -11.11808990357384, 0.9017435354148501, 0.07468873574934592, 0.2144993076964982, 0, 0],
        [4.274138816100098, 3.832016137967051, -0.4421226781330474, 0.4076560326842401, 0.5425947508839231, 3.152340193983333, 0, 0],
        [11.0575304431486, 0.1923876108068164, -10.86514283234179, 0.8884978450722705, 0.08654121116869828, 0.2517634295384992, 0, 0],
        [10.61462497726603, 6.8266128890981, -3.788012088167927, 0.5694354941296446, 0.3930587000395792, 1.720954382601847, 0, 0] ,
        ...
        ...
        ...
        ]
        },
    "headers":    [
    ["alpha", "Mean posterior synonymous substitution rate at a site"],
        ["beta", "Mean posterior non-synonymous substitution rate at a site"],
        ["beta-alpha", "Mean posterior beta-alpha"],
        ["Prob[alpha>beta]", "Posterior probability of negative selection at a site"],
        ["Prob[alpha<beta]", "Posterior probability of positive selection at a site"],
        ["BayesFactor[alpha<beta]", "Empiricial Bayes Factor for positive selection at a site"] 
        ]
    }
    ```
    Is that correct?

    When we executed generateFubarBatchFiles.R sometimes we don't have trees results for some clusters, is that normal?



- loadFubarResults.R

    When we executed loadFubarResults.R, the FUBAR table has too little data, maybe we are not using the correct fasta files?

    Also, Fubar fams decisive evidence is empty because we don't have Bayes Factor > 100


- General
    Step 4 plots everything that we created, but we need the geneSets, and for that we need fubar, busted, etc.