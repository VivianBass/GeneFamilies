
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggpubr)  

  library(tidyverse)
  library(tibble)

# -----------------------------------------------------------------------------

# t-test and wilcox test on long format data , if it isnt already
# Reshape the data to long format
mydata_long <- p.df %>% pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>%
  filter(!is.na(value) & !is.infinite(value))




# Function for performing t-tests and Wilcoxon tests on pairwise distribution data

perform_analysis <- function(data_long, output_file) {

  # Perform t-test 
  ttest_result <- t_test(value ~ group, data = data_long, paired = TRUE)


  # Perform Wilcoxon test with Benjamini-Hochberg (BH) p-value adjustment
  wilcox_result <- data_long %>% group_by(variables) %>%
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% add_significance()

  wilcox_result <- data_long %>% wilcox_test(value ~ group, paired = TRUE) %>%  
  adjust_pvalue(method = "BH") %>% add_significance()
  

  write.table(ttest_result, file = file.path("plots", paste0(output_file, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(wilcox_result, file = file.path("plots", paste0(output_file, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

}













