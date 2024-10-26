
require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)
library(RColorBrewer)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/plot_expression_angles.R")
                
load(file.path(output_data_dir,"expr.angle.diag.RData"))

# plotting the expression angle data
plot.df <- expr.angle.diag.df[!is.nan(expr.angle.diag.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = c("ortholog", "paralog"))
pdf(file.path(output_data_dir, "expressionAngleToDiagonalBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(angle.diag ~ gene.type, data = plot.df, xlab = "Type of Gene", 
    ylab = "relative tissue specificity", border = pushed.colors, col = addAlpha(pushed.colors))
dev.off()

plot.df$rel.vers <- 1 - plot.df$angle.diag
pdf(file.path(output_data_dir, "relativeExpressionVersatilityBoxplot.pdf"))
boxplot(rel.vers ~ gene.type, data = plot.df, xlab = "Type of Gene", ylab = "relative tissue versatility", 
    border = pushed.colors, col = addAlpha(pushed.colors), outline = FALSE)
dev.off()

message("DONE")

# -------------------------------------------------------------

# t-tests and plots based on:
# https://www.datanovia.com/en/fr/blog/comment-effectuer-un-test-t-multiple-dans-r-pour-differentes-variables/
# located in folder: plots/t-test-plots

library(tidyverse)
library(rstatix)
library(ggpubr)
library(tibble)



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

# -------------------------------------------------------------------------

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