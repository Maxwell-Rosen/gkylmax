#!/bin/bash

# Perlmutter-specific modules (commented out for stellar-amd)
# module load PrgEnv-gnu/8.5.0
# module load craype-accel-nvidia80
# module load cray-mpich/8.1.28
# module load cudatoolkit/12.4
# module load nccl/2.18.3-cu12

# Define arrays for magnetic field parameters
mcB_values=(2.130115 2.665626 3.691260 4.490901 5.416264)
gamma_values=(0.451454 0.331696 0.226381 0.182792 0.149893)
R_values=(3 5 10 15 22)
E_beam_midmap_values=(25488.0 25592.9 25678.6 25708.5 25727.8)
T_beam_midmap_values=(204.6 211.6 216.0 217.6 218.6)
gamma_midmap_values=(7.457e+02 8.753e+02 9.519e+02 9.745e+02 9.883e+02)

mkdir -p R-scan

rm -f core/sim

# Submit jobs for paired mcB and gamma scans
for i in "${!mcB_values[@]}"; do
  mcB="${mcB_values[$i]}"
  gamma="${gamma_values[$i]}"
  R="${R_values[$i]}"
	E_beam_midmap="${E_beam_midmap_values[$i]}"
	T_beam_midmap="${T_beam_midmap_values[$i]}"
	gamma_midmap="${gamma_midmap_values[$i]}"
	
  # Create the folder structure
  folder_name="R-scan/R-${R}"

  mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

	if [ "$R" -eq 3 ]; then
		sed -i "353s/.*/  int Nz = 200;/" sim.c
	elif [ "$R" -eq 5 ]; then
		sed -i "353s/.*/  int Nz = 300;/" sim.c
	fi
	sed -i "s/double gamma0 = .*/double gamma0 = $gamma_midmap; /" sim.c
	sed -i "s/double T_beam = .*/double T_beam = $T_beam_midmap * GKYL_ELEMENTARY_CHARGE; /" sim.c
	sed -i "s/double E_beam = .*/double E_beam = $E_beam_midmap * GKYL_ELEMENTARY_CHARGE; /" sim.c

  sed -i "362s/.*/  double mcB = $mcB;/" sim.c
  sed -i "363s/.*/  double gamma = $gamma;/" sim.c
  sed -i "880s|.*|    .filename_psi = \"/home/mr1884/scratch/gkylmax/generate_efit/lorentzian_R${R}.geqdsk_psi.gkyl\",|" sim.c

  sed -i "5s/.*/#SBATCH -J poa-bem-R-${R}/" jobscript-gkyl-stellar-amd

  # Build the simulation
	make clean
  make


  # Submit the job
  # ./sim -s1
  # sbatch jobscript-gkyl-stellar-amd
	bash submit-restarts.sh
	wait

  # Print confirmation
  echo "submitted job for R = $R (mcB = $mcB, gamma = $gamma)"

  # Change back to the root directory
  cd - || exit
done