
# Extract IDs using regular expressions
extract_ids <- function(header) {
  fbpp_match <- regexpr("FBpp\\d+", header)
  fbgn_match <- regexpr("FBgn\\d+", header)
  
  fbpp_id <- regmatches(header, fbpp_match)
  fbgn_id <- regmatches(header, fbgn_match)
  
  data.frame(
    FBpp_ID = fbpp_id,
    Parent_ID = fbgn_id
  )
}

# Read FASTA file and process headers
process_fasta <- function(fasta_file) {
  # Read FASTA headers
  headers <- readLines(fasta_file)
  headers <- headers[grep("^>", headers)]
  headers <- sub("^>", "", headers)
  
  # Process each header and combine results
  result_df <- do.call(rbind, lapply(headers, extract_ids))
  return(result_df)
}


mapping_df <- process_fasta("dsim-all-CDS-r2.02.fasta")

save(mapping_df, "mapping_df.RData")

# -----------------------------------------------------------------------------

extract_ids <- function(header) {
  # Extract parent FBgn ID
  parent_match <- regexpr("parent=([^,]+)", header)
  parent_str <- regmatches(header, parent_match)
  fbgn_id <- sub("parent=", "", strsplit(parent_str, ",")[[1]][1])
  
  # Extract FlyBase Annotation ID
  annot_match <- regexpr("FlyBase_Annotation_IDs:[^,]+", header)
  annot_str <- regmatches(header, annot_match)
  flybase_id <- sub("FlyBase_Annotation_IDs:", "", annot_str)
  
  data.frame(
    FlyBase_Annot_ID = flybase_id,
    Parent_FBgn = fbgn_id,
    stringsAsFactors = FALSE
  )
}

process_fasta <- function(fasta_file) {
  headers <- readLines(fasta_file)
  headers <- headers[grep("^>", headers)]
  headers <- sub("^>", "", headers)
  
  result_df <- do.call(rbind, lapply(headers, extract_ids))
  return(result_df)
}

mapping_df <- process_fasta("dsim-all-CDS-r2.02.fasta")
save(mapping_df_fasta_cds_dsim, file = "mapping_df_fasta_cds_dsim.RData")

mapping_df_fasta_cds_dsim <- mapping_df

# --------------------------------------------------------------------

extract_ids <- function(header) {
  # Extract FBpp ID from start of header
  fbpp_id <- sub("^>*(FBpp\\d+).*", "\\1", header)
  
  # Extract species
  species_match <- regexpr("species=([^;]+)", header)
  species_str <- regmatches(header, species_match)
  species_id <- sub("species=", "", species_str)
  
  # Extract parent FBgn ID
  parent_match <- regexpr("parent=([^,]+)", header)
  parent_str <- regmatches(header, parent_match)
  fbgn_id <- sub("parent=", "", strsplit(parent_str, ",")[[1]][1])
  
  # Extract FlyBase Annotation ID
  annot_match <- regexpr("FlyBase_Annotation_IDs:[^,]+", header)
  annot_str <- regmatches(header, annot_match)
  flybase_id <- sub("FlyBase_Annotation_IDs:", "", annot_str)
  
  data.frame(
    FBpp_ID = fbpp_id,
    Species = species_id,
    Parent_FBgn = fbgn_id,
    FlyBase_Annot_ID = flybase_id,
    stringsAsFactors = FALSE
  )
}

process_fasta <- function(fasta_file) {
  headers <- readLines(fasta_file)
  headers <- headers[grep("^>", headers)]
  headers <- sub("^>", "", headers)
  
  result_df <- do.call(rbind, lapply(headers, extract_ids))
  return(result_df)
}

mapping_df <- process_fasta("data/combined_larger_seq.fasta")
save(mapping_df, file = "mapping_df_fasta_protein_dana.RData")
