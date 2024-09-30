
# t-tests and plots based on:
# https://www.datanovia.com/en/fr/blog/comment-effectuer-un-test-t-multiple-dans-r-pour-differentes-variables/
# located in folder: plots/t-test-plots

library(tidyverse)
library(rstatix)
library(ggpubr)
library(tibble)


# --------------------------------------------------------------------------------

# Angles

# Inputs:
# orths.expr.angle.diag.df
# paralog.expr.angle.diag.df

p.lst <- list(paralog = orths.expr.angle.diag.df, ortholog = paralog.expr.angle.diag.df )
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

----------------------------------------------------------------------------------

# Angles

mydata <- p.df %>% as_tibble()

mydata.long <- mydata %>% pivot_longer(-gene.type, names_to = "variables", values_to = "value")
mydata.long <- mydata.long %>% filter(!is.na(value) & !is.infinite(value))

stat.test <- mydata.long %>%
  group_by(variables) %>%  
  t_test(value ~ gene.type) %>%  
  adjust_pvalue(method = "BH") %>%  
  add_significance()

# Create the boxplot again
myplot <- ggboxplot(
  mydata.long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
facet_wrap(~variables) 

# Add p-values to the plot
stat.test <- stat.test %>% add_xy_position(x = "gene.type") 
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Save the plot as a PDF
ggsave("plots/t-test-plots/Angles_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 10, height = 7)
write_tsv(stat.test, "plots/t-test-plots/Angles_statistical_results.tsv")

# --------------------------------------------------------

# exp.prof.dists

paralogs.exp.prof.dists_dataframe 
orthologs.exp.prof.dists_dataframe 

colnames(orthologs.exp.prof.dists_dataframe) <- c("Family", "exp.prof.dists")

# calculate means function 
paralogs.exp.prof.dists_dataframe <- calculate_means(paralogs.exp.prof.dists_dataframe , column_name = "exp.prof.dists")

orthologs.exp.prof.dists_dataframe <- calculate_means(orthologs.exp.prof.dists_dataframe , column_name = "exp.prof.dists")

p.lst <- list(paralog = paralogs.exp.prof.dists_dataframe , ortholog = orthologs.exp.prof.dists_dataframe)
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, exp.prof.dist = p.lst[[gene.type]]$exp.prof.dists, 
        stringsAsFactors = FALSE)
}))

mydata <- p.df %>% as_tibble()

mydata.long <- mydata %>% pivot_longer(-gene.type, names_to = "variables", values_to = "value")
mydata.long <- mydata.long %>% filter(!is.na(value) & !is.infinite(value))

stat.test <- mydata.long %>%
  group_by(variables) %>%  
  t_test(value ~ gene.type) %>%  
  adjust_pvalue(method = "BH") %>%  
  add_significance()

# Create the boxplot again
myplot <- ggboxplot(
  mydata.long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
facet_wrap(~variables) 

# Add p-values to the plot
stat.test <- stat.test %>% add_xy_position(x = "gene.type") 
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Save the plot as a PDF
ggsave("plots/t-test-plots/exp.prof.dists_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 10, height = 7)
write_tsv(stat.test, "plots/t-test-plots/exp.prof.dists_statistical_results.tsv")


# --------------------------------------------------------------------------------

"df_median_mean_paralogs.exp.prof.dists.tissue"
"df_median_mean_orthologs.exp.prof.dists.tissue"

"df_median_mean_paralogs.exp.prof.dists.tissue_vertical"
"df_median_mean_orthologs.exp.prof.dists.tissue_vertical"