


# ------------------------------------------------------------------------------------
# expressionAngles.R
# ------------------------------------------------------------------------------------

library(parallel)
# for mclapply

load("RData-Files/GeneGroups.RData")
load("RData-Files/RNA_Seq_RPKM_and_profiles.RData")
load("RData-Files/rpkmExpressionProfiles.RData")

source("R/expression_funks.R")

# tands.w.orths --> in GeneGroups.RData
# rpkm.expr.profiles.df --> in rpkmExpressionProfiles.RData
# exprVecSpaceEvolutionAfterDupl (function) --> from expression_funks.R

tands.w.orths.angles.df <- Reduce(rbind, mclapply(names(tands.w.orths), 
    function(fam.name) {
        genes <- intersect(tands.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, tand.classifier)
    }))

print(head(tands.w.orths.angles.df))