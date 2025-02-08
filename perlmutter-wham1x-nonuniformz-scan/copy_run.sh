#!/bin/bash

module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.0
module load nccl/2.18.3-cu12

# Define arrays
# cell_numbers=(48 64 80 96)
map_fractions=(0.500 0.600 0.700 0.800 0.900 0.999)
cell_numbers=(48)
# map_fractions=(0.000 0.999)
# Loop over all combinations of cell_numbers and map_fractions
rm "core/sim"

for cell_number in "${cell_numbers[@]}"; do
  for map_fraction in "${map_fractions[@]}"; do

    # Create the folder structure
    folder_name="${cell_number}/${map_fraction}"

    mkdir -p "$folder_name"

    # Copy core files into the folder
    cp core/* "$folder_name/"

    # Change into the folder
    cd "$folder_name" || exit

    # Update line 797 and line 554 in sim.c
    sed -i "803s/.*/        .map_strength = $map_fraction,/" sim.c
    sed -i "564s/.*/  int Nz = $cell_number;/" sim.c
    sed -i "4/.*/#SBATCH -J w-$map_fraction-$cell_number" jobscript-gkyl-perlmutter

    # Build the simulation
    make sim

    # Submit the job
    sbatch jobscript-gkyl-perlmutter

    # Print confirmation
    echo "submitted job for cell_number = $cell_number and mapping fraction = $map_fraction"

    # Change back to the root directory
    cd - || exit
  done
done