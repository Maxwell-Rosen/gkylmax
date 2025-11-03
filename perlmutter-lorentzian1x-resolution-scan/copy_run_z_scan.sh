#!/bin/bash
module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.4
module load nccl/2.18.3-cu12

#.Disable CUDA-ware MPI, since it causes problems on Perlmutter and we use NCCL alone.
export MPICH_GPU_SUPPORT_ENABLED=0

#.On Perlmutter some jobs get warnings about DVS_MAXNODES (used in file stripping).
#.We set it to 24 for now, but really this depends on the amount/size of I/O being performed.
#.See online NERSC docs and the intro_mpi man page.
export DVS_MAXNODES=24_
export MPICH_MPIIO_DVS_MAXNODES=24

# Define arrays
cell_numbers=(144 192 216 250)
vpar_cell_numbers=(32 48 64 96)
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

  sed -i "367s/.*/  int Nz = $cell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-z-$cell_number/" jobscript-gkyl-perlmutter
  
  # Build the simulation
  make

  # Submit the job
  bash submit-restarts.sh

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

  sed -i "368s/.*/  int Nvpar = $vcell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-vpar-$vcell_number/" jobscript-gkyl-perlmutter

  # Build the simulation
  make

  # Submit the job
  bash submit-restarts.sh

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

  sed -i "369s/.*/  int Nmu = $vcell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-both-$cell_number-$vcell_number/" jobscript-gkyl-perlmutter

  # Build the simulation
  make

  # Submit the job
  bash submit-restarts.sh

  # Print confirmation
  echo "submitted job for mu = $vcell_number"

  # Change back to the root directory
  cd - || exit
done