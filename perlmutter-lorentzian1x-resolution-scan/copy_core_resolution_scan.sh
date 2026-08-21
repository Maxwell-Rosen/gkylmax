#!/bin/bash
set -euo pipefail

# Set to true to use interactive allocations, or false to submit jobs with sbatch.
use_interactive_sessions=true

if [[ $use_interactive_sessions != true && $use_interactive_sessions != false ]]; then
  echo "use_interactive_sessions must be true or false, not '$use_interactive_sessions'." >&2
  exit 1
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
cd "$script_dir"

module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/9.0.1
module load cudatoolkit/13.0
module load nccl/2.29.2-cu13
module load cray-libsci/25.09.0  #.Disable CUDA-ware MPI, since it causes problems on Perlmutter and we use NCCL alone.
export MPICH_GPU_SUPPORT_ENABLED=0

#.On Perlmutter some jobs get warnings about DVS_MAXNODES (used in file stripping).
#.We set it to 24 for now, but really this depends on the amount/size of I/O being performed.
#.See online NERSC docs and the intro_mpi man page.
export DVS_MAXNODES=24_
export MPICH_MPIIO_DVS_MAXNODES=24

# Safely route GPUDirect RDMA over the Host Bridge
export NCCL_NET_GDR_LEVEL=PHB

# Tell NCCL to use the Libfabric plugin
export NCCL_NET="AWS Libfabric"
export NCCL_CROSS_NIC=1

# Disable host registration and eager messages to prevent Slingshot 11 hangs
export FI_CXI_DISABLE_HOST_REGISTER=1
export FI_CXI_RDZV_GET_MIN=0
export FI_CXI_RDZV_THRESHOLD=0
export FI_CXI_RDZV_EAGER_SIZE=0

# Define arrays
cell_numbers=(144 192 216 250)
vpar_cell_numbers=(32 48 64 96)
mu_cell_numbers=(8 16 24 32)

run_job()
{
  if [[ $use_interactive_sessions == true ]]; then
    bash run-interactive-until-complete.sh
  else
    sbatch jobscript-gkyl-perlmutter
  fi
}

rm -f core/sim core/sim.d

# Start runs for z scans
for cell_number in "${cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="z-scan/${cell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i -E "s/^[[:space:]]*int Nz = [0-9]+;/  int Nz = $cell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-z-$cell_number/" jobscript-gkyl-perlmutter
  
  # Build the simulation
  make

  # Start the run using the selected submission mode
  run_job

  # Print confirmation
  echo "started run for cell_number = $cell_number"

  # Change back to the root directory
  cd - || exit
done

# Start runs for vpar scans
for vcell_number in "${vpar_cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="vpar-scan/${vcell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i -E "s/^[[:space:]]*int Nvpar = [0-9]+;/  int Nvpar = $vcell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-vpar-$vcell_number/" jobscript-gkyl-perlmutter

  # Build the simulation
  make

  # Start the run using the selected submission mode
  run_job

  # Print confirmation
  echo "started run for vpar = $vcell_number"

  # Change back to the root directory
  cd - || exit
done

# Start runs for mu scans
for vcell_number in "${mu_cell_numbers[@]}"; do
  # Create the folder structure
  folder_name="mu-scan/${vcell_number}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

  sed -i -E "s/^[[:space:]]*int Nmu = [0-9]+;/  int Nmu = $vcell_number;/" sim.c
  sed -i "4s/.*/#SBATCH -J poa-mu-$vcell_number/" jobscript-gkyl-perlmutter

  # Build the simulation
  make

  # Start the run using the selected submission mode
  run_job

  # Print confirmation
  echo "started run for mu = $vcell_number"

  # Change back to the root directory
  cd - || exit
done
