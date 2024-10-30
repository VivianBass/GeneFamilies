require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(ggplot2)
library(ggsignif)
library(gridExtra)

library(dplyr)
library(purrr)

library(tidyverse)
library(rstatix)
library(ggpubr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(file.path(output_data_dir, "exp.prof.sists_statistics.RData.RData"))
load(file.path(output_data_dir, "exp.prof.dists_tissue.RData"))
load(file.path(output_data_dir, "exp.prof.sists.RData"))

significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# mean_p.df_long
# median_p.df_long

# Perform t-tests for median and mean per tissue between orthologs and paralogs group 
# and adjust p-values, then apply the significance_level function
# group_by() -> defines the tissue column

t_test_median_tissue <- median_p.df_long %>%
  group_by(variables) %>%                             
  t_test(value ~ Type, alternative = "greater") %>%  
  adjust_pvalue(method = "BH") %>%                   
  mutate(significance = sapply(p, significance_level),  
         analysis = "Median")  

t_test_mean_tissue <- mean_p.df_long %>%
  group_by(variables) %>%                             
  t_test(value ~ Type, alternative = "greater") %>%  
  adjust_pvalue(method = "BH") %>%                   
  mutate(significance = sapply(p, significance_level),  
         analysis = "Mean")  

t_test_summary_tissue <- bind_rows(t_test_median_tissue, t_test_mean_tissue)

write.csv(t_test_summary_tissue, file.path(results_dir, "t_test_summary_tissue.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------

# geom_signif(): This function is used to add significance annotations directly to the plot 
# so plots already have the significance level !!

# Creating a Boxplot for Mean Distances
# Calculate the appropriate height for the PDF based on the number of facets
n_facets <- length(unique(mean_p.df_long$variables))
plot_height <- 8 * n_facets  # 8 inches per facet

boxplot_mean <- ggplot(mean_p.df_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances Tissue", y = "Distance") +
  theme_pubr(border = TRUE) +
  geom_signif(comparisons = list(c("ortholog", "paralog")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))) +
  facet_grid(variables ~ ., scales = "free_y", space = "free_y")

ggsave(file.path(results_dir, "tissues_mean_boxplot12.pdf"), plot = boxplot_mean, width = 12, height = plot_height)




# Creating a Boxplot for Median Distances
# Calculate the appropriate height for the PDF based on the number of facets
n_facets <- length(unique(median_p.df_long$variables))
plot_height <- 8 * n_facets  # 8 inches per facet

boxplot_median <- ggplot(median_p.df_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances Tissue", y = "Distance") +
  theme_pubr(border = TRUE) +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))) +
  facet_grid(variables ~ ., scales = "free_y", space = "free_y")

ggsave(file.path(results_dir, "tissues_median_boxplot12.pdf"), plot = boxplot_median, width = 12, height = plot_height)

