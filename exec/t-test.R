
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


# exp.prof.dists.tissue

df1 <- df_median_mean_paralogs.exp.prof.dists.tissue
df2 <- df_median_mean_orthologs.exp.prof.dists.tissue

> names(df1)
[1] "Cluster"           "Mean_Fat_Body"     "Mean_Gut"
[4] "Mean_Muscle"       "Mean_Whole_Body"   "Median_Fat_Body"
[7] "Median_Gut"        "Median_Muscle"     "Median_Whole_Body"

> names(df2)
[1] "Cluster"           "Mean_Fat_Body"     "Mean_Gut"
[4] "Mean_Muscle"       "Mean_Whole_Body"   "Median_Fat_Body"
[7] "Median_Gut"        "Median_Muscle"     "Median_Whole_Body"

c("Cluster" , "Median_Fat_Body", "Median_Gut","Median_Muscle_ortho", "Median_Whole_Body_ortho")

c("Cluster", "Mean_Fat_Body_ortho", "Mean_Gut_ortho", "Mean_Muscle_ortho", "Mean_Whole_Body_ortho")

df_ortho_median <- df1[, c("Cluster" , "Median_Fat_Body", "Median_Gut","Median_Muscle", "Median_Whole_Body")]

df_ortho_mean <- df1[, c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]

df_paralog_median <- df2[, c("Cluster" , "Median_Fat_Body", "Median_Gut","Median_Muscle", "Median_Whole_Body")]

df_paralog_mean <- df2[, c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]




df_ortho_para_mean 
df_ortho_para_median

library(parallel)

p.lst <- list(paralog = df_paralog_mean, ortholog = df_ortho_mean)


p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
  data.frame(
    gene.type = gene.type,
    Mean_Fat_Body = p.lst[[gene.type]]$Mean_Fat_Body, 
    Mean_Gut = p.lst[[gene.type]]$Mean_Gut,
    Mean_Muscle = p.lst[[gene.type]]$Mean_Muscle,
    Mean_Whole_Body = p.lst[[gene.type]]$Mean_Whole_Body,
    stringsAsFactors = FALSE
  )
}))


library(tidyverse)
library(rstatix)
library(ggpubr)


mydata <- p.df  
mydata_long <- mydata %>% pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value")
mydata_long <- mydata_long %>% filter(!is.na(value) & !is.infinite(value))

stat.test <- mydata_long %>%
  group_by(variables) %>%
  pairwise_t_test(value ~ gene.type, p.adjust.method = "BH")

myplot <- ggboxplot(
  mydata_long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",  
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~variables)  

stat.test <- stat.test %>% add_xy_position(x = "gene.type")
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

ggsave("plots/t-test-plots/tissues_mean_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 12, height = 8)

write_tsv(stat.test, "plots/t-test-plots/tissues_mean_statistical_results.tsv")

# --------------------------------------------------------------------------------

df_ortho_para_median

library(parallel)

p.lst <- list(paralog = df_paralog_median, ortholog = df_ortho_median)


p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
  data.frame(
    gene.type = gene.type,
    Median_Fat_Body = p.lst[[gene.type]]$Median_Fat_Body, 
    Median_Gut = p.lst[[gene.type]]$Median_Gut,
    Median_Muscle = p.lst[[gene.type]]$Median_Muscle,
    Median_Whole_Body = p.lst[[gene.type]]$Median_Whole_Body,
    stringsAsFactors = FALSE
  )
}))


mydata <- p.df  
mydata_long <- mydata %>% pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value")
mydata_long <- mydata_long %>% filter(!is.na(value) & !is.infinite(value))

stat.test <- mydata_long %>%
  group_by(variables) %>%
  pairwise_t_test(value ~ gene.type, p.adjust.method = "BH")

myplot <- ggboxplot(
  mydata_long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",  
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~variables)  

stat.test <- stat.test %>% add_xy_position(x = "gene.type")
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

ggsave("plots/t-test-plots/tissues_median_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 12, height = 8)

write_tsv(stat.test, "plots/t-test-plots/tissues_median_statistical_results.tsv")