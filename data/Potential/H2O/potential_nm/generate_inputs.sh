#!/bin/bash

# Check if filenames are provided as arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <xyz_filename> <orca_template> <submit_script>"
    exit 1
fi

xyz_filename="$1"
orca_template="$2"
submit_script="$3"
counter=0
block=""
first_structure=true

# Read the XYZ file
while IFS= read -r line; do
    # Check if the line indicates the start of a new block
    if [[ "$line" =~ ^[0-9]+$ ]]; then
        if [ "$first_structure" = true ]; then
            first_structure=false
        else
            counter=$((counter + 1))
            dir_name="P$counter"
            mkdir -p "$dir_name"
            geom_file="$dir_name/geometry.xyz"
            echo "$block" > "$geom_file"  # Write the entire block to the new file

            # Create additional file using the template
            orca_file="$dir_name/orca_calculation.inp"
            cp "$orca_template" "$orca_file"

            # Copy the submit script to the directory
            cp "$submit_script" "$dir_name/"

            block=""
        fi
    fi

    block="$block$line"$'\n'  # Append the line to the current block
done < "$xyz_filename"

# Write the last block after the loop
counter=$((counter + 1))
dir_name="P$counter"
mkdir -p "$dir_name"
geom_file="$dir_name/geometry.xyz"
echo "$block" > "$geom_file"

# Create additional file using the template for the last structure
orca_file="$dir_name/orca_calculation.inp"
cp "$orca_template" "$orca_file"

# Copy the submit script to the last directory
cp "$submit_script" "$dir_name/"
