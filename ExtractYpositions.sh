#!/bin/bash

# Ensure that samtools is installed
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found. Please install it first."
    exit
fi

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/bam/directory /path/to/output/file"
    exit 1
fi

# Input BAM directory and output file provided as arguments
BAM_DIR="$1"
OUTPUT_FILE="$2"

# Check if BAM directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: Directory $BAM_DIR does not exist."
    exit 1
fi

# Initialize output file (empty or create a new one)
> $OUTPUT_FILE

# Loop through all BAM files in the specified directory
for bam_file in "$BAM_DIR"/*.bam
do
   # Extract the sample name from the BAM filename (remove directory and extension)
    sample_name=$(basename "$bam_file")

    # Ensure it's a valid BAM file
   # if ! samtools quickcheck "$bam_file" 2>/dev/null; then
    #    echo "Skipping non-BAM file: $sample_name"
     #   continue
   # fi

    echo "Processing $bam_file..."
   
   # Check if BAM index file exists (.bai)
    if [ ! -f "$bam_file.bai" ]; then
        echo "Index file for $bam_file not found. Indexing now..."
        samtools index "$bam_file"
    fi

    # Extract chromosome Y reads and their positions, append to the output file
    samtools view "$bam_file" Y | awk -v sample="$sample_name" '{print sample"\t"$4}' >> $OUTPUT_FILE
done

# Display completion message
echo "All positions have been extracted and saved to $OUTPUT_FILE."
