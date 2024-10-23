

## `Dataframes`
--------------------------------------------------------------------------------------
- View() / str() / summary()/ colSums() / dim(df) / head() / class() / which() 
--------------------------------------------------------------------------------------

- df <- read.table("path/2/df.tsv", sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

- write.table(df, "orthologs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

- df <- tibble(Family = character(), col1 = c(1, 2, 3), stringsAsFactors = FALSE)
- df <- data.frame(Family = character(), col1 = c(1, 2, 3), stringsAsFactors = FALSE)

--------------------------------------------------------------------------------------

- `names`:
- df[1:5, "column_name"] 
- colnames(df) <- c("1", "2", "3")
- df$column_name  <- NA
- rownames(df)    <- NULL
- rownames(df)    <- paste0("cluster_", 1:nrow(df))
- df$column_name  <- paste0("cluster_", 1:nrow(df))
- df$column_name  <- sub("^dmel-", "", df$column_name)
- df$column_name  <- sub("^[^-]+-", "", df$column_name)
- df$column_name  <- as.numeric(df$column_name)


## `dplyr Package Dataframes` + `tidyr Package Dataframes` + `tibble Package Dataframes`
--------------------------------------------------------------------------------------

- library(dplyr)
- library(tibble)
- library(tidyr)


- `select`:
--------------------------------------------------------------------------------------

- new_df <- df %>% select(column1, column2, column3)


- `filter/remove rows & columns`:
--------------------------------------------------------------------------------------

- df <- df %>% distinct(Gene, .keep_all = TRUE)             (Filter df)
- df <- df %>% filter(column1 == "value1")                  (Filter df)
- df <- df %>% rowwise() %>% filter()                       (Filter df)

- df <- df[!(df$column_name %in% c("value1", "value2")), ]  (remove rows)
- df <- df[, !(names(df) %in% c("col1", "col2"))]           (remove columns)
- df <- subset(df, select = -c(col1, col2))                 (remove columns)


- `mutate`:
--------------------------------------------------------------------------------------

- df  <- df %>% mutate(new_column = NA)
- df  <- df %>% rowwise() %>% mutate(function(Gene)) %>% ungroup()


- `merge`: 
--------------------------------------------------------------------------------------

- result <- df1 %>% left_join(df2, by = 'Family')
- result <- df1 %>% inner_join(df2, by = 'Family')    
- result <- df1 %>% full_join(df2, by = 'Family')
- result <- df1 %>% semi_join(df2, by = 'Family')
- result <- df1 %>% anti_join(df2, by = 'Family')

- Summary:

- left_join: Keeps all rows from the left data frame and adds matching rows from the right.
- inner_join: Keeps only rows with matching values in both data frames.
- full_join: Keeps all rows from both data frames, with NA where no match is found.
- semi_join: Returns rows from the left data frame that have matches in the right data frame.
- anti_join: Returns rows from the left data frame that do not have matches in the right data frame.

- or standart:

- df <- rbind(df1, df2)
- df <- cbind(df, new_column = my_vector)


- `df --> list` (using dplyr):
--------------------------------------------------------------------------------------

- list <- df %>% group_by(Family) %>% summarise(Gene = list(Gene)) %>% deframe() 

- list <- split(df$Gene, df$Family)
- list <- lapply(list, function(x) list(x))

- row_values <- unlist(df[n, !is.na(df[n,])])           (rows to list)
- row_values <- unlist(df11[2, !is.na(df11[2,])])       (rows to list)
- list <- as.list(df2_new_data[2:4, ])                  (rows to list)


## `Lists`
--------------------------------------------------------------------------------------

- save/load:
- saveRDS(families.lst, "families.lst.rds")
- loaded_list <- readRDS("my_list.rds")

- list(vectors, variables) / unlist(list) 

- Indexing Data in List:
- my_list$a       -> Access element by name
- my_list[["b"]]	-> Access element by name
- my_list[[1]]    -> Access first element            

- add prefix to list elements:
- names(list) <- paste0("cluster_", names(list))

- names       <- paste0("cluster_", seq_along(list))
- names(list) <- names


- apply a function to each element of a list:
- sapply(my_list, function(x))  
- lapply(my_list, function(x))  
- mclapply(list, function(x))

- intersect() -> takes x and y (2 vector Arguments), vector(Schnittmenge)
- setdiff()   -> takes x and y (2 vector Arguments), vector(elements that are present in x but not in y) 
- union()     -> takes x and y (2 vector Arguments), vector(containing all the unique elements from input vectors) 
- Reduce ()   -> takes 2 arguments (a binary function() and list/vector), applies function to all elements in List/Vector
- match()     -> takes x and y (2 vector Arguments)

- sorting:
- cluster_names <- names(your_list)
- sorted_names <- sort(cluster_names)
- sorted_names <- sort(cluster_names, function(x) as.numeric(sub("cluster_", "", x)))
- sorted_list <- your_list[sorted_names]

- remove list nesting level:
- list <- lapply(list, function(cluster) {cluster[-1]})
- list <- lapply(list, function(x) x[[1]])


-----------------------------------------------------------------------------------
# Server  
-----------------------------------------------------------------------------------

- user: kosanenko
- Password: pwd_Ko@n_2024
- IP: 143.93.91.124

- ssh remote_username@remote_host
- ssh kosanenko@143.93.91.124

- /media/BioNAS/ag_hallab/EasyVectorOmics/GeneFamilies 
- /media/BioNAS/ag_hallab/GeneFamilies

- Download Files:
- scp kosanenko@143.93.91.124:/media/BioNAS/ag_hallab/GeneFamilies/data/  /path/2/local
- scp -r user@remote.host:/path/to/remote_folder /home/localuser/destination  (whole folder)

- Upload Files:
- scp /path/2/local kosanenko@143.93.91.124:/media/BioNAS/ag_hallab/EasyVectorOmics


-----------------------------------------------------------------------------------
# Git-basic-Commands 
-----------------------------------------------------------------------------------

- 1. Navigate to Your local Git Repository, via Ubuntu Terminal
- cd /mnt/c/Users/andre/Desktop/GeneFamilies-tests-Andre1   

- create Repository on GitHub: EasyVectorOmics

…or create a new repository on the command line
echo "# EasyVectorOmics" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/andre3322/EasyVectorOmics.git
git push -u origin main


…or push an existing repository from the command line:
git remote add origin https://github.com/andre3322/EasyVectorOmics.git
git branch -M main
git push -u origin main
git checkout GeneFamilies-tests-Andrej

git remote add origin https://github.com/VivianBass/GeneFamilies.git


https://github.com/VivianBass/GeneFamilies.git


- Branches:
- git branch -a 		(To see all branches, including remote branches)
- git branch 		(Check Current Branch you are using)
- git checkout GeneFamilies-tests-Andre  (switching from the master branch to the GeneFamilies-Andre)


- pull


- Upload: 
- git add .
- git commit -m "msg"
- git push

- GitHub/ username:	andre3322
- Token as Pass: 	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 






## `Slurm-Script`
-----------------------------------------------------------------------------------

#!/bin/bash

	# SBATCH --job-name=my_job
	# SBATCH --output=my_job%j.out
	
	echo "Starting test.R script"

	Rscript exec/test.R

	echo "test.R script finished!"


- sbatch test_slurm.sh (execute)

- squeue (check status)

- scancel job_id



