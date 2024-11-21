#!/bin/bash

# Variables principales
BASE_DIR="/projects/EasyVectorOmics/material/dmel"         
INDEX_FILE="/projects/EasyVectorOmics/material/references/dmel_transcriptome.idx"  
RESULTS_DIR="/projects/EasyVectorOmics/results"   

# Loop for every tissue
for TISSUE_DIR in "$BASE_DIR"/*/trimmed; do
  TISSUE=$(basename $(dirname "$TISSUE_DIR")) 
  
  echo "Processing tissue: $TISSUE"

  # Get all fastq files 
  for FILE in "$TISSUE_DIR"/*.fastq.gz; do
    FILENAME=$(basename "$FILE" .fastq.gz) 
    echo "Processing file: $FILENAME"
    
    # Temporal file with read lenghts
    zcat "$FILE" | awk '(NR%4==2) {print length($1)}' > "read_lengths_$FILENAME.txt"
    
    # Calculate mean and sd for this file
    MEAN_SD=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}' "read_lengths_$FILENAME.txt")
    MEAN_LENGTH=$(echo $MEAN_SD | cut -d' ' -f1)
    STD_DEV=$(echo $MEAN_SD | cut -d' ' -f2)
    
    echo "Mean length: $MEAN_LENGTH, Std Dev: $STD_DEV"
    
    # Execute kallisto
    OUTPUT_DIR="$RESULTS_DIR/$TISSUE/$FILENAME"
    mkdir -p "$OUTPUT_DIR"
    kallisto quant -i "$INDEX_FILE" \
      -o "$OUTPUT_DIR" -b 100 \
      --single -l "$MEAN_LENGTH" -s "$STD_DEV" \
      "$FILE"
    
    # Delete temporal read lenghts file
    rm "read_lengths_$FILENAME.txt"
    
    echo "$FILENAME completed."
  done
done

echo "DONE for all tissues"

