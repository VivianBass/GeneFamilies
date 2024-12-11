#!/bin/bash

# Base dir where dietM is stored
base_dir="/media/BioNAS/ag_hallab/EasyVectorOmics/material/dietM"

# Dir to adapters 
adapter_file="/usr/local/bin/trimmomatic/adapters/TruSeq3-SE.fa"

# Loop dmel and dsec
for species in dsec ; do
    species_dir="$base_dir/$species"

    # Check if dir exists
    if [ -d "$species_dir" ]; then
        # Loop for tissues
        for tissue in fat muscle whole_body gut; do
            tissue_dir="$species_dir/$tissue"

            # Check if tissue dir exists
            if [ -d "$tissue_dir" ]; then
                # Create "trimmed" directory inside tissue dir
                output_dir="$tissue_dir/trimmed"
                mkdir -p "$output_dir"

                # Check all fastq files
                for fastq in "$tissue_dir"/*.fastq.gz; do
                    # Get filename whitout extension
                    base_name=$(basename "$fastq" .fastq.gz)

                    # Execute trimmomatic for every file
                    java -jar /usr/local/bin/trimmomatic/trimmomatic-0.39.jar SE -threads 4 \
                    "$fastq" \
                    "$output_dir/${base_name}_trimmed.fastq" \
                    ILLUMINACLIP:$adapter_file:2:30:10:2 \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70
                done
            else
                echo "Error: $tissue_dir doesn't exist."
            fi
        done
    else
        echo "Error: $species_dir doesn't exist."
    fi
done