
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggpubr)  
  library(tidyverse)
  library(tibble)

# Wilcoxon Signed-Rank Test: Used to compare two related groups (paired data), similar to the paired t-test.

# R Libraries for Wilcoxon tests, rstatix and ggpubr

perform_wilcox_test <- function(data_long, output_file) {
  
  # Perform the Wilcoxon test
  wilcox_result <- data_long %>%
    wilcox_test(value ~ group, paired = TRUE) %>%  # If your data is paired
    adjust_pvalue(method = "BH") %>%               # Adjust p-values using Benjamini-Hochberg (BH)
    add_significance()                             # Add significance stars
  
  # Save the results to a .tsv file
  write.table(wilcox_result, 
              file = file.path("plots", paste0(output_file, ".tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Return results for further analysis or plotting
  return(list(results = wilcox_result))
}






mydata <- p.df  
mydata_long <- mydata %>% 
  pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>% 
  filter(!is.na(value) & !is.infinite(value))

# Perform the Wilcoxon test and adjust p-values using Benjamini-Hochberg (BH) method
stat.test <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% 
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "gene.type")

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

# Add p-values to the plot
bxp_with_pvals <- myplot + 
  stat_pvalue_manual(stat.test, tip.length = 0) +  
  labs(subtitle = get_test_label(stat.test, detailed = TRUE))  

print(stat.test)

effect_size <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_effsize(value ~ gene.type)

print(effect_size)

ggsave("plots/wilcox-test-boxplot-with-pvalues.pdf", plot = bxp_with_pvals, width = 8, height = 6)
write_tsv(stat.test, "plots/wilcox-test-statistical-results.tsv")


# -----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------





# exp.prof.dists_mean

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

#-----------------------------------------------------------------------------


paralogs.exp.prof.dists_stats_df
orthologs.exp.prof.dists_stats_df

df1 <- paralogs.exp.prof.dists_stats_df
df2 <- orthologs.exp.prof.dists_stats_df


calculate_means <- function(df, column_name) {
  library(dplyr)

  df_mean <- df %>%
    rowwise() %>%  
    mutate(
      mean_value = mean(as.numeric(unlist(strsplit(as.character(get(column_name)), ","))), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Family, mean_value)
  
  colnames(df_mean) <- c("Family", column_name)

  return(df_mean)
}

df1 <- calculate_means(df1 , column_name = "mean")
df2 <- calculate_means(df2 , column_name = "mean")



para_stats_median <- df1[ , c("Family", "median")]
ortho_stats_median <- df2[ , c("Family", "median")]

para_stats_mean <- df1[ , c("Family", "mean")]
ortho_stats_mean <- df2[ , c("Family", "mean")]


p.lst <- list(paralog = para_stats_mean , ortholog = ortho_stats_mean)
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, median = p.lst[[gene.type]]$mean, 
        stringsAsFactors = FALSE)
}))



mydata <- p.df  
mydata_long <- mydata %>% 
  pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>% 
  filter(!is.na(value) & !is.infinite(value))

# Perform the Wilcoxon test and adjust p-values using Benjamini-Hochberg (BH) method
stat.test <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% 
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "gene.type")

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

# Add p-values to the plot
bxp_with_pvals <- myplot + 
  stat_pvalue_manual(stat.test, tip.length = 0) +  
  labs(subtitle = get_test_label(stat.test, detailed = TRUE))  

print(stat.test)

effect_size <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_effsize(value ~ gene.type)

print(effect_size)

ggsave("plots/wilcox-test-plots/mean_exp.prof.dists_boxplot-with-pvalues.pdf", plot = bxp_with_pvals, width = 8, height = 6)
write_tsv(stat.test, "plots/wilcox-test-plots/mean_exp.prof.dists_statistical-results.tsv")


# -----------------------------------------------------------------------------
# Angles

# orths.expr.angle.diag.df
# paralog.expr.angle.diag.df

p.lst <- list(paralog = orths.expr.angle.diag.df, ortholog = paralog.expr.angle.diag.df )
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

mydata <- p.df  
mydata_long <- mydata %>% 
  pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>% 
  filter(!is.na(value) & !is.infinite(value))

# Perform the Wilcoxon test and adjust p-values using Benjamini-Hochberg (BH) method
stat.test <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% 
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "gene.type")

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

# Add p-values to the plot
bxp_with_pvals <- myplot + 
  stat_pvalue_manual(stat.test, tip.length = 0) +  
  labs(subtitle = get_test_label(stat.test, detailed = TRUE))  

print(stat.test)

effect_size <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_effsize(value ~ gene.type)

print(effect_size)

ggsave("plots/wilcox-test-plots/Angles_boxplot-with-pvalues.pdf", plot = bxp_with_pvals, width = 8, height = 6)
write_tsv(stat.test, "plots/wilcox-test-plots/Angles_statistical-results.tsv")


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

library(parallel)


p.lst <- list(paralog = df_ortho_median, ortholog = df_paralog_median)


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




p.lst <- list(paralog = df_ortho_mean, ortholog = df_paralog_mean)

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

mydata <- p.df  
mydata_long <- mydata %>% 
  pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>% 
  filter(!is.na(value) & !is.infinite(value))

# Perform the Wilcoxon test and adjust p-values using Benjamini-Hochberg (BH) method
stat.test <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% 
  add_significance()

# Add XY position for p-value annotations
stat.test <- stat.test %>% add_xy_position(x = "gene.type")

# Create the boxplot for four tissues
myplot <- ggboxplot(
  mydata_long,                
  x = "gene.type",           
  y = "value",               
  fill = "gene.type",        
  palette = "npg",           
  legend = "none",           
  ggtheme = theme_pubr(border = TRUE)  
) +
  facet_wrap(~variables)  # This allows separate facets for each tissue


bxp_with_pvals <- myplot + 
  stat_pvalue_manual(stat.test, tip.length = 0) +  
  labs(subtitle = stat.test %>% 
                   group_by(variables) %>%
                   summarise(label = paste0("Wilcoxon test, W = ", round(statistic, 2), 
                                            ", p = ", round(p, 3), 
                                            ", n = ", n1 + n2)) %>%
                   pull(label) %>%
                   paste(collapse = "\n"))

print(stat.test)

# Compute effect size for the Wilcoxon test
effect_size <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_effsize(value ~ gene.type)

# Print effect sizes
print(effect_size)

# Save the plot and statistical results
ggsave("plots/wilcox-test-plots/tissues_median_boxplot_with_pvalues.pdf", plot = bxp_with_pvals, width = 12, height = 8)
write_tsv(stat.test, "plots/wilcox-test-plots/tissues_median_statistical_results.tsv")



