

# >FBpp0113192 type=protein; loc=scaffold_13337:join(17236694..17236810,17241571..17241816,17241904..17242057,17245092..17245726,17245789..17245947,17246015..17246324,17246412..17246706,17247147..17247429,17247619..17247861); ID=FBpp0113192; name=Dana\Om(2D)-PA; parent=FBgn0010400,FBtr0114700; dbxref=FlyBase:FBpp0113192,FlyBase_Annotation_IDs:GF10000-PA,GB_protein:EDV40849,REFSEQ:XP_001958043,FlyMine:FBpp0113192; MD5=f03f12164a7494c47fe0763bff201e8a; length=813; release=r1.3; species=Dana;

# provide the fasta files you need, 

extract_ids <- function(header) {
  # Extract FBpp ID from start of header
  fbpp_id <- sub("^>*(FBpp\\d+).*", "\\1", header)
  
  # Extract species with NA handling
  species_match <- regexpr("species=([^;]+)", header)
  species_id <- if(species_match > 0) {
    species_str <- regmatches(header, species_match)
    sub("species=", "", species_str)
  } else {
    NA
  }
  
  # Extract parent FBgn with NA handling
  parent_match <- regexpr("parent=([^;]+)", header)
  fbgn_id <- if(parent_match > 0) {
    parent_str <- regmatches(header, parent_match)
    parent_parts <- strsplit(sub("parent=", "", parent_str), ",")[[1]]
    if(length(parent_parts) > 0) parent_parts[1] else NA
  } else {
    NA
  }
  
  data.frame(
    FBpp_ID = fbpp_id,
    Parent_FBgn = fbgn_id, 
    Species = species_id,
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
save(mapping_df, mapping_df_unique, file = "mapping_df.RData")