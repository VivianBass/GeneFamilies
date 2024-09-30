

# calculate expression profile distances for orthologs and paralogs


# Inputs:
# rna.seq.exp.profils
# paralogs.lst, orthologs.lst, families.lst

# Input files used in this script:
# vv_df3_ortho.lst as (orthologs.lst)
# vv_df3_para.lst as (paralogs.lst)
# df14_P_dmel as (rna.seq.exp.profils) 

# function used in this script:
# expressionProfilesDists (from calculate_expression_profile_dists_function.R)

library(parallel)

# orthologs exp.prof.dists

orthologs.exp.prof.dists <- mclapply(orthologs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = FALSE))

orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))


# paralogs exp.prof.dists

paralogs.exp.prof.dists <- mclapply(paralogs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = FALSE))

paralogs.exp.prof.dists.tissue <- mclapply(paralogs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))


save(paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
    orthologs.exp.prof.dists.tissue, file = file.path("data", "ExpressionProfileDistances.RData"))

# families exp.prof.dists (not used in further analysis)

#non.singleton.fams <- families.df$id[which(families.df$size > 1)]

#families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], 
#function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))

#families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], 
#function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))

#save(paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
#    orthologs.exp.prof.dists.tissue, families.exp.prof.dists, families.exp.prof.dists.tissue, 
#    file = file.path("data", "ExpressionProfileDistances.RData"))