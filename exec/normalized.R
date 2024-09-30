# Assuming df12 is your dataframe
# Assuming df12 is your dataframe
df12 <- data.frame(
  Gene = c("FBgn0000008", "FBgn0000014", "FBgn0000015", "FBgn0000017", "FBgn0000018"),
  Tissue1 = c("4.9", "17.7", "1.9", NA, "7.4"),
  Tissue2 = c("3.8", "15.3", "2.3", "7.7", "7.1"),
  Tissue3 = c("3.2", "13.3", "1.9", "7.0", "8.3"),
  Tissue4 = c("4.8", "18.7", "2.0", "8.6", "7.7"),
  stringsAsFactors = FALSE
)

# Convert expression columns to numeric
df14[, -1] <- lapply(df14[, -1], function(x) as.numeric(gsub(",", ".", x)))

# Function to normalize expression values for a single gene
normalize_expression <- function(expressions) {
  # Replace NA with 0
  expressions[is.na(expressions)] <- 0
  return(expressions / sum(expressions, na.rm = TRUE))
}

# Apply normalization to each row (gene) in the dataframe
normalized_df <- df14
normalized_df[, -1] <- t(apply(df14[, -1], 1, normalize_expression))

# Print the normalized dataframe
print(normalized_df)













# Lade dplyr Paket
library(dplyr)

# Beispiel Dataframe mit deinen Headernamen
df <- data.frame(
  Gene = c("Gene1", "Gene2", "Gene3"),
  dmel_M_whole_body_3rd = c(1, 2, 3),
  dmel_M_whole_body_3rd_2 = c(4, 5, 6), # Beispiel für zweite Spalte
  dmel_M_whole_body_3rd_3 = c(7, 8, 9),
  dmel_P_whole_body_3rd = c(10, 11, 12),
  dmel_P_whole_body_3rd_2 = c(13, 14, 15),
  dmel_P_whole_body_3rd_3 = c(16, 17, 18)
  # Weitere Spalten hier hinzufügen...
)

library(dplyr)

# Angenommen, df13 ist dein ursprüngliches DataFrame
# Eindeutige Spaltennamen erstellen
colnames(df13) <- make.names(colnames(df13), unique = TRUE)

# Berechnung der Mittelwerte mit den neuen, eindeutigen Spaltennamen
df_mean <- df13 %>%
  rowwise() %>%
  mutate(
    Mean_dmel_M_Whole_Body = mean(c_across(starts_with("dmel_M_whole_body_3rd")), na.rm = TRUE),
    Mean_dmel_P_Whole_Body = mean(c_across(starts_with("dmel_P_whole_body_3rd")), na.rm = TRUE),
    Mean_dsec_M_Whole_Body = mean(c_across(starts_with("dsec_M_whole_body_3rd")), na.rm = TRUE),
    Mean_dsec_P_Whole_Body = mean(c_across(starts_with("dsec_P_whole_body_3rd")), na.rm = TRUE),
    Mean_dmel_M_Muscle = mean(c_across(starts_with("dmel_M_muscle")), na.rm = TRUE),
    Mean_dmel_P_Muscle = mean(c_across(starts_with("dmel_P_muscle")), na.rm = TRUE),
    Mean_dsec_M_Muscle = mean(c_across(starts_with("dsec_M_muscle")), na.rm = TRUE),
    Mean_dsec_P_Muscle = mean(c_across(starts_with("dsec_P_muscle")), na.rm = TRUE),
    Mean_dmel_M_Gut = mean(c_across(starts_with("dmel_M_gut")), na.rm = TRUE),
    Mean_dmel_P_Gut = mean(c_across(starts_with("dmel_P_gut")), na.rm = TRUE),
    Mean_dsec_M_Gut = mean(c_across(starts_with("dsec_M_gut")), na.rm = TRUE),
    Mean_dsec_P_Gut = mean(c_across(starts_with("dsec_P_gut")), na.rm = TRUE),
    Mean_dmel_M_Fat_Body = mean(c_across(starts_with("dmel_M_fat_body")), na.rm = TRUE),
    Mean_dmel_P_Fat_Body = mean(c_across(starts_with("dmel_P_fat_body")), na.rm = TRUE),
    Mean_dsec_M_Fat_Body = mean(c_across(starts_with("dsec_M_fat_body")), na.rm = TRUE),
    Mean_dsec_P_Fat_Body = mean(c_across(starts_with("dsec_P_fat_body")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(Gene, 
         Mean_dmel_M_Whole_Body, 
         Mean_dmel_P_Whole_Body,
         Mean_dsec_M_Whole_Body,
         Mean_dsec_P_Whole_Body,
         Mean_dmel_M_Muscle,
         Mean_dmel_P_Muscle,
         Mean_dsec_M_Muscle,
         Mean_dsec_P_Muscle,
         Mean_dmel_M_Gut,
         Mean_dmel_P_Gut,
         Mean_dsec_M_Gut,
         Mean_dsec_P_Gut,
         Mean_dmel_M_Fat_Body,
         Mean_dmel_P_Fat_Body,
         Mean_dsec_M_Fat_Body,
         Mean_dsec_P_Fat_Body)