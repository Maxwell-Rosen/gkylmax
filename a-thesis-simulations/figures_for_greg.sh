pgkyl 1x-beams-extra-dilated/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl 1x-beams-extra-dilated/zzim-ion_source_integrated_moms.gkyl ev "f0 f1 / 2 *" sel -c0 pl -g --ymax 2 --ymin 0 --title "Boundary Flux / Source from extra dilated" &

pgkyl 1x-beams-anneal/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl 1x-beams-anneal/zzim-ion_source_integrated_moms.gkyl ev "f0 f1 / 2 *" sel -c0 pl -g --ymax 2 --ymin 0 --title "Boundary Flux / Source from anneal" &

pgkyl 1x-beams-ecycles/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl 1x-beams-ecycles/zzim-ion_source_integrated_moms.gkyl ev "f0 f1 / 2 *" sel -c0 pl -g --ymax 2 --ymin 0 --title "Boundary Flux / Source from ecycles" &

pgkyl 1x-beams-higher-alpha/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl 1x-beams-higher-alpha/zzim-ion_source_integrated_moms.gkyl ev "f0 f1 / 2 *" sel -c0 pl -g --ymax 2 --ymin 0 --title "Boundary Flux / Source from higher alpha" &

pgkyl 1x-beams-resolution-scan/mu-scan/40/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl 1x-beams-resolution-scan/mu-scan/40/zzim-ion_source_integrated_moms.gkyl ev "f0 f1 / 2 *" sel -c0 pl -g --ymax 2 --ymin 0 --title "Boundary Flux / Source N_mu = 40" &

pgkyl 1x-beams-anneal/zzim-ion_integrated_moms.gkyl  -l "Anneal" 1x-beams-extra-dilated/zzim-ion_integrated_moms.gkyl -l "Extra Dilated" sel -c0 pl -g --title "Anneal vs Extra Dilated" --ylabel "Integrated Density" -f0 &
