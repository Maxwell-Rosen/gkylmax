#!/bin/bash
set -euo pipefail

# Set to true to use interactive allocations, or false to submit jobs with sbatch.
use_interactive_sessions=false

if [[ $use_interactive_sessions != true && $use_interactive_sessions != false ]]; then
  echo "use_interactive_sessions must be true or false, not '$use_interactive_sessions'." >&2
  exit 1
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
cd "$script_dir"

: "${PS1:=}"
module purge
module load nvhpc/26.5 # For NCCL.
module load cudatoolkit/13.3
module load openmpi/nvhpc-26.5/5.0.10
module load openblas/0.3.x
export LD_LIBRARY_PATH="${NVHPC_ROOT:?}/comm_libs/nccl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export NCCL_DEBUG=WARN UCX_LOG_LEVEL=warn UCX_TLS=self,sm,cuda

# Define arrays
cell_numbers=(128 192 256 320)
# vpar_cell_numbers=(32 48 56)
# mu_cell_numbers=(24 40)

run_job()
{
  if [[ $use_interactive_sessions == true ]]; then
    bash run-interactive-until-complete.sh
  else
    sbatch jobscript-gkyl-stellar-ai
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
  sed -i "4s/.*/#SBATCH -J poa-z-$cell_number/" jobscript-gkyl-stellar-ai
  
  # Build the simulation
  make

  # Start the run using the selected submission mode
  run_job

  # Print confirmation
  echo "started run for cell_number = $cell_number"

  # Change back to the root directory
  cd - || exit
done

# # Start runs for vpar scans
# for vcell_number in "${vpar_cell_numbers[@]}"; do
#   # Create the folder structure
#   folder_name="vpar-scan/${vcell_number}"

#   mkdir -p "$folder_name"

#   # Copy core files into the folder
#   cp core/* "$folder_name/"

#   # Change into the folder
#   cd "$folder_name" || exit

#   sed -i -E "s/^[[:space:]]*int Nvpar = [0-9]+;/  int Nvpar = $vcell_number;/" sim.c
#   sed -i "4s/.*/#SBATCH -J poa-vpar-$vcell_number/" jobscript-gkyl-stellar-ai

#   # Build the simulation
#   make

#   # Start the run using the selected submission mode
#   run_job

#   # Print confirmation
#   echo "started run for vpar = $vcell_number"

#   # Change back to the root directory
#   cd - || exit
# done

# # Start runs for mu scans
# for vcell_number in "${mu_cell_numbers[@]}"; do
#   # Create the folder structure
#   folder_name="mu-scan/${vcell_number}"

#   mkdir -p "$folder_name"

#   # Copy core files into the folder
#   cp core/* "$folder_name/"

#   # Change into the folder
#   cd "$folder_name" || exit

#   sed -i -E "s/^[[:space:]]*int Nmu = [0-9]+;/  int Nmu = $vcell_number;/" sim.c
#   sed -i "4s/.*/#SBATCH -J poa-mu-$vcell_number/" jobscript-gkyl-stellar-ai

#   # Build the simulation
#   make

#   # Start the run using the selected submission mode
#   run_job

#   # Print confirmation
#   echo "started run for mu = $vcell_number"

#   # Change back to the root directory
#   cd - || exit
# done
