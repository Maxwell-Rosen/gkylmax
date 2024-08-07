#!/bin/bash

# Initialize variables to keep track of the highest number and the corresponding file
max_number=-1
matching_files=$(ls gk_wham-elc_M0_*)
for file in $matching_files; do
  # # Extract the number between the prefix and suffix
  number=$(echo "$file" | sed 's/gk_wham-elc_M0_//; s/\.gkyl//')
  # If the number is greater than the current max_number, update max_number and max_file
  if [ "$number" -gt "$max_number" ]; then
    max_number="$number"
  fi
done

# Print the file with the highest number
if [ "$max_number" -ne -1 ]; then
  echo "The highest number is $max_number"
fi