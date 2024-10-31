require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/load_gene_groups_data.R  input.args ...")

message("PURPOSE: This R script loads five different gene groups of pairwise orthologs and paralogs data, including one conserved ortholog file and four different paralog files. 
All four paralog files should follow the same header naming convention.")

message("Orthologs_Header: Family | Gene | Gene_species | Ortholog | Ortholog_species")
message("input.args[[1]]:  path/2/<conserved_orthologs.tsv>")

message("Paralogs_Header:  Family | Gene | Gene_species | Paralog | Paralog_species ")
message("input.args[[2]]:  path/2/<in_paralogs.tsv>")
message("input.args[[3]]:  path/2/<special_in_paralogs.tsv> ")
message("input.args[[4]]:  path/2/<out_paralogs.tsv>")
message("input.args[[5]]:  path/2/<special_out_paralogs.tsv>")

input.args <- commandArgs(trailingOnly = TRUE)

# orthologs
# conserved orthologs
con_orthologs <- read.table(input.args[[1]], header = TRUE, sep = "\t", 
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

con_orthologs.lst <- con_orthologs %>% 
               group_by(Family, Gene_species, Gene, Ortholog_species) %>%
               summarise(Ortholog = list(Ortholog), .groups = "drop") %>%
               group_by(Family, Gene_species, Gene) %>%
               summarise(orthologs = list(setNames(Ortholog, Ortholog_species)), .groups = "drop") %>%
               # Create the nested list by Family
               group_by(Family) %>%
               summarise(gene_info = list(setNames(orthologs, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>% deframe()

# paralogs
# in paralogs with orthologs
in_paralogs <- read.table(input.args[[2]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

in_paralogs.lst <- in_paralogs %>%
               group_by(Family, Gene_species, Gene, Paralog_species) %>%
               summarise(Paralog = list(Paralog), .groups = "drop") %>%
               group_by(Family, Gene_species, Gene) %>%
               summarise(paralogs = list(setNames(Paralog, Paralog_species)), .groups = "drop") %>%
               group_by(Family) %>%
               summarise(gene_info = list(setNames(paralogs, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>% 
               deframe()


# in paralogs without orthologs
out_paralogs <- read.table(input.args[[3]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

out_paralogs.lst <- out_paralogs %>%
               group_by(Family, Gene_species, Gene, Paralog_species) %>%
               summarise(Paralog = list(Paralog), .groups = "drop") %>%
               group_by(Family, Gene_species, Gene) %>%
               summarise(paralogs = list(setNames(Paralog, Paralog_species)), .groups = "drop") %>%
               group_by(Family) %>%
               summarise(gene_info = list(setNames(paralogs, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>% 
               deframe()

# out paralogs with orthologs
special_in_paralogs <- read.table(input.args[[4]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

special_in_paralogs.lst <- special_in_paralogs %>%
               group_by(Family, Gene_species, Gene, Paralog_species) %>%
               summarise(Paralog = list(Paralog), .groups = "drop") %>%
               group_by(Family, Gene_species, Gene) %>%
               summarise(paralogs = list(setNames(Paralog, Paralog_species)), .groups = "drop") %>%
               group_by(Family) %>%
               summarise(gene_info = list(setNames(paralogs, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>% 
               deframe()

# out paralogs without orthologs
special_out_paralogs <- read.table(input.args[[5]], header = TRUE, sep = "\t",          
               comment.char = "", quote = "", na.strings = "", 
               colClasses = rep("character", 5))

special_out_paralogs.lst <- special_out_paralogs %>%
               group_by(Family, Gene_species, Gene, Paralog_species) %>%
               summarise(Paralog = list(Paralog), .groups = "drop") %>%
               group_by(Family, Gene_species, Gene) %>%
               summarise(paralogs = list(setNames(Paralog, Paralog_species)), .groups = "drop") %>%
               group_by(Family) %>%
               summarise(gene_info = list(setNames(paralogs, paste0("(", Gene_species, ", ", Gene, ")"))), .groups = "drop") %>% 
               deframe()


save(con_orthologs, con_orthologs.lst,
     in_paralogs, in_paralogs.lst,
     out_paralogs, out_paralogs.lst,
     special_in_paralogs, special_in_paralogs.lst,
     special_out_paralogs, special_out_paralogs.lst,
     file = file.path(output_data_dir, "gene_groups.RData"))

message("DONE")



