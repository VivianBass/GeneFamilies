
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


mapping_df <- process_fasta("combined_larger_seq.fasta")

save(mapping_df, "mapping_df.RData")