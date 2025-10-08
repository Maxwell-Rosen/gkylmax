#!/bin/bash

module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.4
module load nccl/2.18.3-cu12

# Define arrays
cell_numbers=(144 216 288 432 576)
vcell_numbers=(16 48 64 96)

rm core/sim

# Submit jobs for z scans
for cell_number in "${cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="z-scan/${cell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i "440s/.*/  int Nz = $cell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-z-$cell_number/" jobscript-gkyl-perlmutter
  
  # Build the simulation
  make

  ./sim -s1

  # Submit the job
  # bash submit-restarts.sh

  # Print confirmation
  echo "submitted job for cell_number = $cell_number"

  # Change back to the root directory
  cd - || exit
done

# # Submit jobs for vpar scans
# for vcell_number in "${vcell_numbers[@]}"; do
#   # Create the folder structure
#   folder_name="vpar-scan/${vcell_number}"

#   mkdir -p "$folder_name"

#   # Copy core files into the folder
#   cp core/* "$folder_name/"

#   # Change into the folder
#   cd "$folder_name" || exit

#   sed -i "441s/.*/  int Nvpar = $vcell_number;/" sim.c
#   sed -i "4s/.*/#SBATCH -J poa-vpar-$vcell_number/" jobscript-gkyl-perlmutter

#   # Build the simulation
#   make

#   # Submit the job
#   # bash submit-restarts.sh

#   # Print confirmation
#   echo "submitted job for vpar = $vcell_number"

#   # Change back to the root directory
#   cd - || exit
# done

# # Submit jobs for mu scans
# for vcell_number in "${vcell_numbers[@]}"; do
#   # Create the folder structure
#   folder_name="mu-scan/${vcell_number}"

#   mkdir -p "$folder_name"

#   # Copy core files into the folder
#   cp core/* "$folder_name/"

#   # Change into the folder
#   cd "$folder_name" || exit

#   sed -i "442s/.*/  int Nmu = $vcell_number;/" sim.c
#   sed -i "4s/.*/#SBATCH -J poa-both-$cell_number-$vcell_number/" jobscript-gkyl-perlmutter

#   # Build the simulation
#   make

#   # Submit the job
#   # bash submit-restarts.sh

#   # Print confirmation
#   echo "submitted job for mu = $vcell_number"

#   # Change back to the root directory
#   cd - || exit
# done