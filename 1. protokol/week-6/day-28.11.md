
**Date**:  28.11.2024
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:








- maybe also both log2 for rpkm counts not normalized and normalized / any difference?
- So please, get the logarithm of the original tpm values and after that normalize the values between 0 and 1, and then you can calculate the distances
yes , so basically rpkm values (counts) -> log2 transform -> normalize -> distances etc
- you can use a count so you don't get negative values, for example 1, as you said before. 
Please check if you have negative values, if so:
Get the log2(tpm+1), then normalize the values between 0 and 1 and then calculate the distances
- exactly, if you get negative values you can add 1 to the counts so you get log2(value+1), don't add 1 before the log2 calculation
- if i log2 transform normilezed values -> negative

- log2 transform values before distance calculation -> rna.seq.exp.profils_log2


- see if any negative values afterwards, and if this causes an issue

- separate plots for t-test and wilcoxon

- calculate statistics for everything and make plots for everything




- adding angles distance calculation scripts
your vectors, as I understand your code, are the genes in all tissues, so for gene1 if you have 4 tissues you will have something like
 
gene_id | tissue1 | tissue2 | tissue3 | tissue4
fbPP            0.9          3.4              5.4           0.2
 
so your vector will be [0.9,3.4,5.4,0.2]
 
this is how you are calculating the euclidiean distance now
You can do something like this I think
# Expression values in vector form
gene1 <- c(23.19, 12.34, 8.45)
gene2 <- c(2.50, 9.87, 11.23)
 
# Dot product
dot_product <- sum(gene1 * gene2)
 
# Magnitudes
magnitude_gene1 <- sqrt(sum(gene1^2))
magnitude_gene2 <- sqrt(sum(gene2^2))
 
# Cosine
cosine_angle <- dot_product / (magnitude_gene1 * magnitude_gene2)





**Doubts and Issues**:

- should i make separate scripts for cosine angles distnaces and the euclidean distances? or all in one script?

- how to handle log transformed negative values? 3 cases of log transformed values:
- log transform before of after distance calculation?? i did afterwards 

- for cosine angles distances, they are already per tissue, somehow, do i need per tissue and overall like in euclidean and how ?
- there are cosine distances and also cosine angles , Which should i use, or basically just the formula i got 


- need to clarify which statistics i need to calculate , 

- what about genefamilies

- log transform rna.seq.exp.profils first ??


**Next Steps**:


---

**Code:**