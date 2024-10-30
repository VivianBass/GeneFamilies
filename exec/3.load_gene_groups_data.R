require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/load_gene_groups_data.R  ...")

message("We will have 5 groups:")
message("Orthologs_Header: Family / Gene / Gene_species / Ortholog / Ortholog_species")
message("input.args[[1]]: conserved orthologs file")

message("Paralogs_Header: Family / Gene / Gene_species / Paralog / Paralog_species ")
message("input.args[[2]]: in paralogs with orthologs")
message("input.args[[3]]: in paralogs without orthologs")
message("input.args[[4]]: out paralogs with orthologs")
message("input.args[[5]]: out paralogs without orthologs")

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[1]] <- "experiments/test/conserved_orthologs_test.txt"
input.args[[2]] <- "experiments/test/in_paralogs_test.tsv"
input.args[[3]] <- "experiments/test/special_in_paralogs_test.tsv"
input.args[[4]] <- "experiments/test/out_paralogs_test.tsv"
input.args[[5]] <- "experiments/test/special_out_paralogs_test.tsv"

# Load data:
# we need to load pairs of orthologs, pairs of paralogs (not only the list itself)
# we need to work with pairs of genes instead of a list of genes.


# Orthologs
# conserved orthologs
con_orthologs <- read.table(input.args[[1]], header = TRUE, sep = "\t", 
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

con_orthologs.lst <- con_orthologs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
               mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
               select(cluster_name, Gene) %>% deframe()

# Paralogs
# in paralogs with orthologs
in_paralogs <- read.table(input.args[[2]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

in_paralogs.lst <- in_paralogs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
               mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
               select(cluster_name, Gene) %>% deframe()

# in paralogs without orthologs
out_paralogs <- read.table(input.args[[3]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

out_paralogs.lst <- in_paralogs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
               mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
               select(cluster_name, Gene) %>% deframe()

# out paralogs with orthologs
special_in_paralogs <- read.table(input.args[[4]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

special_in_paralogs.lst <- in_paralogs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
               mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
               select(cluster_name, Gene) %>% deframe()

# out paralogs without orthologs
special_out_paralogs <- read.table(input.args[[5]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

special_out_paralogs.lst <- in_paralogs %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>%
               mutate(cluster_name = paste("Orthogroup_", row_number(), sep = "")) %>%
               select(cluster_name, Gene) %>% deframe()


save(con_orthologs, con_orthologs.lst,
     in_paralogs, in_paralogs.lst,
     out_paralogs, out_paralogs.lst,
     special_in_paralogs, special_in_paralogs.lst,
     special_out_paralogs, special_out_paralogs.lst,
     file = file.path(output_data_dir, "gene_groups.RData"))

message("DONE")

