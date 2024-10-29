require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("input.args[[1]]: <all_vs_all_file>")
message("<all_vs_all_file> is a tabular Blast Output txt file, generated in an 'all vs all' approach") 

input.args <- commandArgs(trailingOnly = TRUE)

# colClasses should be adjusted for different data
all.vs.all.sim <- fread(input.args[[1]], data.table = FALSE, header = FALSE, 
    stringsAsFactors = FALSE, sep = "\t", na.strings = "", 
    colClasses = c(rep("character", 2), rep("numeric", 10)))

save(all.vs.all.sim, file = file.path(output_data_dir, "blast.RData"))