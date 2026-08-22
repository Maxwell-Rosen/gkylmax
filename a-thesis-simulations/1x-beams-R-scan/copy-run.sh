#!/bin/bash

# Perlmutter-specific modules (commented out for stellar-amd)
# module load PrgEnv-gnu/8.5.0
# module load craype-accel-nvidia80
# module load cray-mpich/8.1.28
# module load cudatoolkit/12.4
# module load nccl/2.18.3-cu12

# Define arrays for magnetic field parameters
mcB_values=(2.130115 2.665626 3.691260 4.490901 5.416264 8.125522)
gamma_values=(0.451454 0.331696 0.226381 0.182792 0.149893 0.098619)
R_values=(3 5 10 15 22 50)
src_amp=(133.489434322 167.08296763 187.529468419 193.893947813 197.854501103 202.174079918)
src_temp=(25141.1582175 25127.2323251 25120.9096778 25117.3014755 25114.8563657 25126.3674234)

mkdir -p R-scan

rm -f core/sim

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

	if [ "$R" -eq 3 ]; then
		sed -i 's/^[[:space:]]*int Nz = .*/  int Nz = 192;/' sim.c
	elif [ "$R" -eq 5 ]; then
		sed -i 's/^[[:space:]]*int Nz = .*/  int Nz = 256;/' sim.c
	elif [ "$R" -eq 50 ]; then
		sed -i 's/^[[:space:]]*int Nz = .*/  int Nz = 800;/' sim.c
		sed -i 's/^[[:space:]]*bool is_positivity_enabled_fdp = .*/  bool is_positivity_enabled_fdp = true;/' sim.c
		sed -i 's|^[[:space:]]*double cfl_factor_times_omega_max = .*|  double cfl_factor_times_omega_max = 1/20.0; // CFL factor for fixed factor times omega max multiplier.|' sim.c
	fi

  sed -i "s/^[[:space:]]*double gamma0 = .*/  double gamma0 = ${src_amp[$i]}; \/\/ Beam intM0 = 3.1099679832479321e+20/" sim.c
  sed -i "s/^[[:space:]]*double E_beam = .*/  double E_beam = ${src_temp[$i]} * GKYL_ELEMENTARY_CHARGE; \/\/ Beam intM2 = 1.2940166363523933e+06/" sim.c
  sed -i "s/^[[:space:]]*double mcB = .*/  double mcB = $mcB;/" sim.c
  sed -i "s/^[[:space:]]*double gamma = .*/  double gamma = $gamma;/" sim.c
  sed -i "s|^[[:space:]]*\.filename_psi = .*|    .filename_psi = \"/home/mr1884/scratch/gkylmax/generate_efit/lorentzian_R${R}.geqdsk_psi.gkyl\", // psi file to use|" sim.c

  sed -i "s/^#SBATCH -J .*/#SBATCH -J poa-bem-R-${R}/" jobscript-gkyl-perlmutter

  # Build the simulation
	make clean
  make


  # Fine-tune gamma0/E_beam against the Hamiltonian-moments integrated
  # diagnostic before submitting.
  bash optimize_source_params.sh

  # Submit the job
  sbatch jobscript-gkyl-perlmutter
	# bash submit-restarts.sh

	wait

  # Print confirmation
  echo "submitted job for R = $R (mcB = $mcB, gamma = $gamma)"

  # Change back to the root directory
  cd - || exit
done
