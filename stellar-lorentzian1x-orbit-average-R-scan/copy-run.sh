#!/bin/bash

module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.4
module load nccl/2.18.3-cu12

# Define arrays for magnetic field parameters
mcB_values=(2.130115 2.665626 3.691260 4.490901 5.168764)
gamma_values=(0.451454 0.331696 0.226381 0.182792 0.157441)
R_values=(3 5 10 15 20)

rm core/sim

# Submit jobs for paired mcB and gamma scans
for i in "${!mcB_values[@]}"; do
  mcB="${mcB_values[$i]}"
  gamma="${gamma_values[$i]}"
  R="${R_values[$i]}"
  
  # Create the folder structure
  folder_name="R-scan/R-${R}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i "372s/.*/  double mcB = $mcB;/" sim.c
  sed -i "373s/.*/  double gamma = $gamma;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-R-${R}/" jobscript-gkyl-stellar-amd
  
  # Build the simulation
  make

  ./sim -s1

  # Submit the job
  # sbatch jobscript-gkyl-stellar-amd

  # Print confirmation
  echo "submitted job for R = $R (mcB = $mcB, gamma = $gamma)"

  # Change back to the root directory
  cd - || exit
done