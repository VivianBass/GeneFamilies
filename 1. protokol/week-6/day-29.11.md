
**Date**:  28.11.2024
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:


- figure out , how to plot complete, create some plot scripts for complete
- überlegen wie man t tests und wilcox für die completen machen kann

load(file.path(output_data_dir, "exp.prof.dists.RData"))

 [4] "con_orthologs_v.lst_dists"
 [5] "con_orthologs_v.lst_dists_tissue"

 [9] "in_paralogs_v.lst_dists"
[10] "in_paralogs_v.lst_dists_tissue"

[11] "out_paralogs_v.lst_dists"
[12] "out_paralogs_v.lst_dists_tissue"

[17] "special_in_paralogs_v.lst_dists"
[18] "special_in_paralogs_v.lst_dists_tissue"

[19] "special_out_paralogs_v.lst_dists"
[20] "special_out_paralogs_v.lst_dists_tissue"


load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))

 [4] "con_orthologs_v.lst_dists_log2"
 [5] "con_orthologs_v.lst_dists_tissue_log2"

 [9] "in_paralogs_v.lst_dists_log2"
[10] "in_paralogs_v.lst_dists_tissue_log2"

[11] "out_paralogs_v.lst_dists_log2"
[12] "out_paralogs_v.lst_dists_tissue_log2"

[17] "special_in_paralogs_v.lst_dists_log2"
[18] "special_in_paralogs_v.lst_dists_tissue_log2"

[19] "special_out_paralogs_v.lst_dists_log2"
[20] "special_out_paralogs_v.lst_dists_tissue_log2"


load(file.path(output_data_dir, "exp.prof.angles.RData"))

 [4] "con_orthologs_v.lst_cos_angles_dists"
                     
 [8] "in_paralogs_v.lst_cos_angles_dists"

 [9] "out_paralogs_v.lst_cos_angles_dists"

[14] "special_in_paralogs_v.lst_cos_angles_dists"

[15] "special_out_paralogs_v.lst_cos_angles_dists"


load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))

 [4] "con_orthologs_v.lst_cos_angles_dists_log2"

 [8] "in_paralogs_v.lst_cos_angles_dists_log2"

 [9] "out_paralogs_v.lst_cos_angles_dists_log2"

[14] "special_in_paralogs_v.lst_cos_angles_dists_log2"

[15] "special_out_paralogs_v.lst_cos_angles_dists_log2"


- t-test, wilcox test for the complete expression distances

- new experiment with filtered orthogroups ( with at least 6 gene pairs)

- combine all plot scripts, redundant code etc. 

- documentation 



**Doubts and Issues**:

**Next Steps**:


---

**Code:**