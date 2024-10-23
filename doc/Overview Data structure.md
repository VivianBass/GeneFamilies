

-----------------------------------------------------------------------------------------
# `Section - 1: Load tissue specific expression data`
-----------------------------------------------------------------------------------------

> str(`rpkm.rna.seq.counts`)

  # 'data.frame':   1977070 obs. of  5 variables:
  # $ id        : chr  "dmel-FBgn0000003" "dmel-FBgn0000008" "dmel-FBgn0000014" "dmel-FBgn0000015" ...
  # $ tissue    : chr  "mE_mRNA_A_MateF_1d_head" "mE_mRNA_A_MateF_1d_head" "mE_mRNA_A_MateF_1d_head" "mE_mRNA_A_MateF_1d_head" ...
  # $ rank      : chr  "NA" "NA" "NA" "NA" ...
  # $ expression: num  0 10 0 0 8 1 0 58 22 2 ...
  # $ variance  : num  NA NA NA NA NA NA NA NA NA NA ...


- Output: > str(`RNA_Seq_RPKM_and_profiles.tsv`) (`rna.seq.exp.profils`)

# 'data.frame':   17681 obs. of  113 variables:

# $ mE_mRNA_A_MateF_1d_head        : num  0 0.00979 0 0 0.00701 ...
# $ mE_mRNA_A_MateF_4d_ovary       : num  0 0.0196 0 0 0.0158 ...
# $ mE_mRNA_A_MateM_1d_head        : num  0 0.0137 0 0 0.0114 ...
# $ mE_mRNA_A_VirF_1d_head         : num  0 0.01077 0 0 0.00701 ...
# $ mE_mRNA_A_VirF_4d_head         : num  0 0.00686 0 0 0.00788 ...
# $ mE_mRNA_A_MateF_20d_head       : num  0 0.00686 0 0 0.00525 ...
# $ mE_mRNA_A_MateF_4d_head        : num  0 0.00979 0 0 0.00613 ...







-----------------------------------------------------------------------------------------
# `2. Compute pairwise gene expression profile distances within gene groups`
-----------------------------------------------------------------------------------------


> dim(orthologs)
  [1] 16077    13

> summary(orthologs)

  #  Orthogroup            dana               dere               dgri
  Length:16077       Length:16077       Length:16077       Length:16077
  Class :character   Class :character   Class :character   Class :character
  Mode  :character   Mode  :character   Mode  :character   Mode  :character
  #    dmel               dmoj               dper               dpse
  Length:16077       Length:16077       Length:16077       Length:16077
  Class :character   Class :character   Class :character   Class :character
  Mode  :character   Mode  :character   Mode  :character   Mode  :character
  #    dsec               dsim               dvir               dwil
  Length:16077       Length:16077       Length:16077       Length:16077
  Class :character   Class :character   Class :character   Class :character
  Mode  :character   Mode  :character   Mode  :character   Mode  :character
  #    dyak
  Length:16077
  Class :character
  Mode  :character

r$> head(`paralogs`) (Dataframe, 223055 Objects/Elements)

  #               Family                Gene
  # 1 paralogs_cluster_0    dmel-FBgn0029830
  # 2 paralogs_cluster_0    dmel-FBgn0029835
  # 3 paralogs_cluster_0    dmel-FBgn0001263
  # 4 paralogs_cluster_0    dmel-FBgn0067864
  # 5 paralogs_cluster_0    dmel-FBgn0000163
  # 6 paralogs_cluster_0    dmel-FBgn0026313


r$> str(orthologs.lst) List of 16077 objects (ortholog clusters)

# $ ortholog_cluster_1    : Named chr [1:13] "OG0000000" "FBgn0090932, FBgn0092744, FBgn0092910, FBgn0092945, FBgn0095994, FBgn0096182, FBgn0096292, FBgn0097172, FBgn009"| #__truncated__ "FBgn0103302, FBgn0103306, FBgn0103309, FBgn0103312, FBgn0104863, FBgn0104867, FBgn0105228, FBgn0105263, FBgn010"| __truncated__ "FBgn0117630, FBgn0117652, #FBgn0118969, FBgn0130789, FBgn0130791, FBgn0130796, FBgn0130801, FBgn0130926, FBgn013"| __truncated__ ...
#  ..- attr(*, "names")= chr [1:13] "Orthogroup" "dana" "dere" "dgri" ...
# $ ortholog_cluster_2    : Named chr [1:13] "OG0000001" "FBgn0090950, FBgn0091434, FBgn0092911, FBgn0092946, FBgn0093040, FBgn0095967, FBgn0095992, FBgn0096019, FBgn009"| #__truncated__ "FBgn0103253, FBgn0103288, FBgn0103300, FBgn0103304, FBgn0103307, FBgn0103310, FBgn0103313, FBgn0103354, FBgn010"| __truncated__ "FBgn0117663, FBgn0121413, #FBgn0129951, FBgn0130790, FBgn0130802, FBgn0130873, FBgn0130947, FBgn0130972, FBgn013"| __truncated__ ...
#  ..- attr(*, "names")= chr [1:13] "Orthogroup" "dana" "dere" "dgri" ...
# $ ortholog_cluster_3    : Named chr [1:13] "OG0000002" "FBgn0088059, FBgn0092909, FBgn0092912, FBgn0092944, FBgn0095878, FBgn0095993, FBgn0096291, FBgn0096319, FBgn009"| #__truncated__ "FBgn0103254, FBgn0103267, FBgn0103305, FBgn0103308, FBgn0103311, FBgn0103315, FBgn0104862, FBgn0104865, FBgn010"| __truncated__ "FBgn0129952, FBgn0130787, # # #FBgn0130792, FBgn0130794, FBgn0130797, FBgn0130800, FBgn0130948, FBgn0130960, FBgn013"| __truncated__ ...
#  ..- attr(*, "names")= chr [1:13] "Orthogroup" "dana" "dere" "dgri" ...

r$> str(paralogs.lst) List of 10772 objects (tandem_clusters)

# List of 10772
# $ paralogs_cluster_0    : chr [1:28] "dmel-FBgn0029830" "dmel-FBgn0029835" "dmel-FBgn0001263" "dmel-FBgn0067864" ...
# $ paralogs_cluster_1    : chr [1:82] "dmel-FBgn0030058" "dmel-FBgn0026411" "dmel-FBgn0020617" "dmel-FBgn0008636" ...
# $ paralogs_cluster_2    : chr [1:86] "dmel-FBgn0029697" "dmel-FBgn0030058" "dmel-FBgn0026411" "dmel-FBgn0020617" ...
# $ paralogs_cluster_3    : chr [1:60] "dmel-FBgn0000382" "dmel-FBgn0003079" "dmel-FBgn0003731" "dmel-FBgn0028484" ...
# $ paralogs_cluster_4    : chr [1:24] "dmel-FBgn0004170" "dmel-FBgn0002561" "dmel-FBgn0000137" "dmel-FBgn0011276" ...
# $ paralogs_cluster_5    : chr [1:33] "dmel-FBgn0029690" "dmel-FBgn0000326" "dmel-FBgn0034736" "dmel-FBgn0034972" ...
...


r$> str(families.lst) List of List of 1244 objects (clusters)

 # $ cluster_1   :List of 47
 # ..$ : chr "dmel-FBgn0034997"
 # ..$ : chr "dana-FBgn0144047"
 # ..$ : chr "dvir-FBgn0088427"
 # ..$ : chr "dpse-FBgn0208706"
 # ..$ : chr "dsec-FBgn0077426"
 # ..$ : chr "dwil-FBgn0166790"
 # ..$ : chr "dper-FBgn0224185"
 # ..$ : chr "dyak-FBgn0157602"
 # ..$ : chr "dsim-FBgn0232001"
 # ..$ : chr "dere-FBgn0183587"
 # ...
 # $ cluster_2   :List of 72
 # ..$ : chr "dmel-FBgn0263598"
 # ..$ : chr "dana-FBgn0140828"
 # ..$ : chr "dvir-FBgn0091241"
 # ..$ : chr "dpse-FBgn0250311"
 # ..$ : chr "dsec-FBgn0165036"
 # ..$ : chr "dwil-FBgn0217312"
 # ..$ : chr "dper-FBgn0236016"
 # ..$ : chr "dyak-FBgn0195230"

> str(families.df)

# 'data.frame':   1244 obs. of  14 variables:
# $ id  : chr  "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
# $ dmel: num  4 6 1 3 12 15 1 5 19 4 ...
# $ dana: num  4 6 1 3 12 14 1 5 18 4 ...
# $ dvir: num  4 6 1 3 12 14 1 5 18 4 ...
# $ dpse: num  4 6 1 3 12 13 1 5 18 4 ...
# $ dsec: num  4 6 1 3 12 13 1 5 18 4 ...
# $ dwil: num  4 6 1 2 10 12 1 5 18 4 ...
# $ dper: num  4 6 1 2 10 11 1 5 18 4 ...
# $ dyak: num  4 6 1 3 12 13 1 5 18 4 ...
# $ dsim: num  4 6 1 2 12 14 1 5 18 4 ...
# $ dere: num  4 6 1 2 12 13 1 5 18 4 ...
# $ dgri: num  3 6 1 2 10 10 1 5 15 4 ...
# $ dmoj: num  4 6 1 2 10 11 1 5 16 4 ...
# $ size: num  43 66 11 28 126 142 11 55 196 44 ...




  str(`families.exp.prof.dists`) List of 1244 (Key: Value pairs)

# $ cluster_1   : 'dist' num [1:6] 0.185 0.304 0.306 0.239 0.245 ...
#  ..- attr(*, "Size")= int 4
#  ..- attr(*, "Labels")= chr [1:4] "dmel-FBgn0034997" "dmel-FBgn0035421" "dmel-FBgn0039768" "dmel-FBgn0039769"
#  ..- attr(*, "Diag")= logi FALSE
#  ..- attr(*, "Upper")= logi FALSE
#  ..- attr(*, "method")= chr "euclidean"
#  ..- attr(*, "call")= language dist(x = exp.profs[, tissues], method = dist.method)
# $ cluster_2   : 'dist' num [1:3] 0.403 0.395 0.144
#  ..- attr(*, "Size")= int 3
#  ..- attr(*, "Labels")= chr [1:3] "dmel-FBgn0032464" "dmel-FBgn0263598" "dmel-FBgn0265262"
#  ..- attr(*, "Diag")= logi FALSE
#  ..- attr(*, "Upper")= logi FALSE
#  ..- attr(*, "method")= chr "euclidean"
#  ..- attr(*, "call")= language dist(x = exp.profs[, tissues], method = dist.method)

  str(`families.exp.prof.dists.tissue`) [1] 1244

#  ..$ mE_mRNA_L3_Wand_saliv          : 'dist' num [1:6] 0.000814 0.0005 0.0005 0.001315 0.001315 ...
#  .. ..- attr(*, "Size")= int 4
#  .. ..- attr(*, "Labels")= chr [1:4] "dmel-FBgn0034997" "dmel-FBgn0035421" "dmel-FBgn0039768" "dmel-FBgn0039769"
#  .. ..- attr(*, "Diag")= logi FALSE
#  .. ..- attr(*, "Upper")= logi FALSE
#  .. ..- attr(*, "method")= chr "euclidean"
#  .. ..- attr(*, "call")= language dist(x = setNames(exp.profs[, tissue], rownames(exp.profs)), method = dist.method)
#  ..$ mE_mRNA_A_VirF_20d_head        : 'dist' num [1:6] 0.00212 0.00475 0.00475 0.00263 0.00263 ...
#  .. ..- attr(*, "Diag")= logi FALSE
#  .. ..- attr(*, "Upper")= logi FALSE
#  .. ..- attr(*, "method")= chr "euclidean"
#  .. ..- attr(*, "call")= language dist(x = setNames(exp.profs[, tissue], rownames(exp.profs)), method = dist.method)
#  ..$ mE_mRNA_A_VirF_20d_head        : 'dist' num [1:6] 0.00212 0.00475 0.00475 0.00263 0.00263 ...
#  .. ..- attr(*, "method")= chr "euclidean"
#  .. ..- attr(*, "call")= language dist(x = setNames(exp.profs[, tissue], rownames(exp.profs)), method = dist.method)



> str(`paralogs.exp.prof.dists`) List of [1] 10772

# $ paralogs_cluster_96   : 'dist' num [1:351] 0.0578 0.0721 0.0593 0.0413 0.0722 ...
#  ..- attr(*, "Size")= int 27
#  ..- attr(*, "Labels")= chr [1:27] "dmel-FBgn0005322" "dmel-FBgn0013563" "dmel-FBgn0016983" "dmel-FBgn0020369" ...
#  ..- attr(*, "Diag")= logi FALSE
#  ..- attr(*, "Upper")= logi FALSE
#  ..- attr(*, "method")= chr "euclidean"
#  ..- attr(*, "call")= language dist(x = exp.profs[, tissues], method = dist.method)
# $ paralogs_cluster_97   : 'dist' num [1:6] 0.189 0.306 0.12 0.273 0.109 ...
#  ..- attr(*, "Size")= int 4
#  ..- attr(*, "Labels")= chr [1:4] "dmel-FBgn0001086" "dmel-FBgn0030011" "dmel-FBgn0034937" "dmel-FBgn0262699"
#  ..- attr(*, "Diag")= logi FALSE
#  ..- attr(*, "Upper")= logi FALSE
#  ..- attr(*, "method")= chr "euclidean"
#  ..- attr(*, "call")= language dist(x = exp.profs[, tissues], method = dist.method)

> str(`paralogs.exp.prof.dists.tissue`) List of [1] 10772

#  ..$ mE_mRNA_WPP_fat                : 'dist' num [1:378] 0.003125 0.001315 0.002061 0.000271 0.003383 ...
#  .. ..- attr(*, "Size")= int 28
#  .. ..- attr(*, "Labels")= chr [1:28] "dmel-FBgn0000163" "dmel-FBgn0001263" "dmel-FBgn0001624" "dmel-FBgn0010620" ...
#  .. ..- attr(*, "Diag")= logi FALSE
#  .. ..- attr(*, "Upper")= logi FALSE
#  .. ..- attr(*, "method")= chr "euclidean"
#  .. ..- attr(*, "call")= language dist(x = setNames(exp.profs[, tissue], rownames(exp.profs)), method = dist.method)
#  ..$ mE_mRNA_WPP_saliv              : 'dist' num [1:378] 0 0.00197 0.00355 0.01273 0.01374 ...
#  .. ..- attr(*, "Size")= int 28
#  .. ..- attr(*, "Labels")= chr [1:28] "dmel-FBgn0000163" "dmel-FBgn0001263" "dmel-FBgn0001624" "dmel-FBgn0010620" ...
#  .. ..- attr(*, "Diag")= logi FALSE
#  .. ..- attr(*, "Upper")= logi FALSE
#  .. ..- attr(*, "method")= chr "euclidean"
#  .. ..- attr(*, "call")= language dist(x = setNames(exp.profs[, tissue], rownames(exp.profs)), method = dist.method)






