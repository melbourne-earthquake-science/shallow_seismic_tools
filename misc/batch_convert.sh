#!/bin/bash

# Define the directory containing MiniSEED files
MINISEED_DIR="/mnt/data/" # This is the path inside the Docker container

# Check if an argument was passed
if [ ! -z "$1" ]; then
  MINISEED_DIR="$1"
fi

# Continue using $MINISEED_DIR in your script
echo "MINISEED_DIR is set to: $MINISEED_DIR"

# Loop through all .mseed files in the directory
for mseed_file in "$MINISEED_DIR"/*.mseed; do
    # Get the base filename without extension
    base_name=$(basename "$mseed_file" .mseed)

    # Define the output SEG-2 filename
    seg2_file="$MINISEED_DIR/${base_name}.seg2"

    # Convert the MiniSEED file to SEG-2 format using Geopsy
    geopsy "$mseed_file" -export "$seg2_file" -export-format seg2

    # Check if the SEG-2 file was created successfully
    if [ -f "$seg2_file" ]; then
        # Rename the .seg2 file to .dat
        mv "$seg2_file" "$MINISEED_DIR/${base_name}.dat"
        echo "Converted and renamed $mseed_file to ${base_name}.dat"
    else
        echo "Failed to convert $mseed_file to SEG-2 format"
    fi
done
