


Function Breakdown
Input Arguments:

gene.accessions: This is a vector of gene identifiers (orthologs from various species).
expression.profiles: A data frame containing gene expression profiles across multiple tissues (default is rna.seq.exp.profils).
expr.prof.gene.col: The column name that holds the main gene accession IDs (default is "Gene").
expr.prof.gene.variant.col: The column name that holds gene variants (default is "gene.exp.var").
tissues: A vector of tissue names to limit the distance calculations to specific tissues (default is NULL, meaning all tissues).
dist.method: The method for calculating distances (default is "euclidean").
per.tissue: A logical flag indicating whether to calculate distances per tissue (if TRUE) or across all tissues at once (default is FALSE).
Function Flow:

Identifying Tissues: If tissues is not provided, the function automatically extracts tissue names by excluding the columns for gene and gene variant information from the data frame. These are treated as the columns with expression data.

Matching Genes: The function locates the rows in expression.profiles that match any of the gene accessions provided in gene.accessions. These accessions are expected to match either the expr.prof.gene.col (main gene accessions) or the expr.prof.gene.variant.col (gene variants).

Extracting Expression Profiles: Once the relevant rows (genes) are identified, the function extracts their expression profiles across the specified tissues.

Distance Calculation:

If per.tissue = TRUE, it calculates the distance between expression profiles for each tissue separately. This means that for each tissue, a distance matrix will be computed comparing the expression profiles of all orthologous genes.
If per.tissue = FALSE, it calculates a single distance matrix across all tissues combined.
Output:

The function returns a distance matrix (or a list of distance matrices if computed per tissue). These matrices contain pairwise distances between the expression profiles of the orthologous genes, which are useful for clustering or comparing the genes based on their expression similarity.
Error Handling:

The function includes several tryCatch blocks to gracefully handle errors during row name setting or distance calculation.
If only one matching gene is found or an error occurs, the function returns NA.
Input Data Example:
The input data, as you provided, represents clusters of orthologous genes from different species. Here is an example of one cluster:

r
Code kopieren
ortholog_cluster_1: [
"311263",       # Gene from species 1
"AT1G01790",    # Gene from species 2
"Carubv10011770m", # Gene from species 3
"CARHR000530",  # Gene from species 4
"aet_11800",    # Gene from species 5
"Bra032651",    # Gene from species 6
"Thhalv10006595m", # Gene from species 7
"Tp1g00490"     # Gene from species 8
]
Each cluster contains genes from different species, all of which are orthologs (genes that have evolved from a common ancestor). The function is applied to each cluster individually to compute the expression profile distances between these genes.

Example Process:
For ortholog_cluster_1:

The function retrieves the expression profiles for the genes in this cluster (311263, AT1G01790, etc.) from the expression.profiles data frame.
It calculates the Euclidean distances between the expression profiles of these genes, either for each tissue (if per.tissue = TRUE) or across all tissues combined (if per.tissue = FALSE).
The result is a distance matrix showing how similar or different the expression patterns of these orthologous genes are.
mclapply Call:
The two mclapply calls at the end of the script show that the function is applied to all ortholog clusters (like ortholog_cluster_1, ortholog_cluster_2, etc.):

orthologs.exp.prof.dists <- mclapply(orthologs.lst, expressionProfilesDists_Viv) applies the function across all ortholog clusters without calculating distances per tissue.
orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, expressionProfilesDists_Viv, per.tissue = TRUE) applies the function but calculates distances separately for each tissue.
Summary:
This function helps compute how similar or different the expression profiles of orthologous genes are across various tissues (or within individual tissues), producing pairwise distance matrices for further analysis, such as clustering.