#!/bin/bash

# Path to the first folder
first_folder="32"
first_file="$first_folder/sim.c"

# Open the first sim.c in VS Code
code "$first_file"

# Loop through all folders in the parent directory except the first
for folder in */; do
    # Skip the first folder
    if [[ "$folder" == "$first_folder/" ]]; then
        continue
    fi
    
    # Path to the sim.c in the current folder
    other_file="${folder}sim.c"
    
    # Check if sim.c exists in the other folder
    if [[ -f "$other_file" ]]; then
        # Open the diff view in VS Code between the first sim.c and the current one
        code -d "$first_file" "$other_file"
    else
        echo "No sim.c found in $folder"
    fi
done