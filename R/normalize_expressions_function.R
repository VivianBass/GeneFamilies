
# Assuming df is your dataframe

# df <- data.frame(
#  Gene = c("FBgn0000008", "FBgn0000014", "FBgn0000015", "FBgn0000017", "FBgn0000018"),
#  Tissue1 = c("4.9", "17.7", "1.9", NA, "7.4"),
#  Tissue2 = c("3.8", "15.3", "2.3", "7.7", "7.1"),
#  Tissue3 = c("3.2", "13.3", "1.9", "7.0", "8.3"),
#  Tissue4 = c("4.8", "18.7", "2.0", "8.6", "7.7"),
#  stringsAsFactors = FALSE)

# Convert expression columns to numeric
df[, -1] <- lapply(df[, -1], function(x) as.numeric(gsub(",", ".", x)))

# Function to normalize expression values for a single gene
normalize_expression <- function(expressions) {
  # Replace NA with 0
  expressions[is.na(expressions)] <- 0
  return(expressions / sum(expressions, na.rm = TRUE))
}

normalized_df <- df
normalized_df[, -1] <- t(apply(df[, -1], 1, normalize_expression))








