#!/bin/bash

module load nvhpc/25.5
module load cudatoolkit/12.9
module load openmpi/cuda-12.9/nvhpc-25.5/4.1.8
module load openblas/0.3.x

# Define arrays
cell_numbers=(200 250 325 400)
vpar_cell_numbers=(24 32 48 64)
mu_cell_numbers=(8 16 24 32)

rm core/sim
rm core/sim.d

# Submit jobs for z scans
for cell_number in "${cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="z-scan/${cell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i "160s/.*/  int Nz = $cell_number;/" sim.c
  sed -i "5s/.*/#SBATCH -J poa-z-$cell_number/" jobscript_gkyl-stellar-amd
  
  # Build the simulation
  make

  # Submit the job
  sbatch jobscript_gkyl-stellar-amd

  # Print confirmation
  echo "submitted job for cell_number = $cell_number"

  # Change back to the root directory
  cd - || exit
done

# Submit jobs for vpar scans
for vcell_number in "${vpar_cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="vpar-scan/${vcell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i "161s/.*/  int Nvpar = $vcell_number;/" sim.c
  sed -i "5s/.*/#SBATCH -J poa-vpar-$vcell_number/" jobscript_gkyl-stellar-amd

  # Build the simulation
  make

  # Submit the job
  sbatch jobscript_gkyl-stellar-amd

  # Print confirmation
  echo "submitted job for vpar = $vcell_number"

  # Change back to the root directory
  cd - || exit
done

# Submit jobs for mu scans
for vcell_number in "${mu_cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="mu-scan/${vcell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i "162s/.*/  int Nmu = $vcell_number;/" sim.c
  sed -i "5s/.*/#SBATCH -J poa-mu-${vcell_number}/" jobscript_gkyl-stellar-amd

  # Build the simulation
  make

  # Submit the job
  sbatch jobscript_gkyl-stellar-amd

  # Print confirmation
  echo "submitted job for mu = $vcell_number"

  # Change back to the root directory
  cd - || exit
done