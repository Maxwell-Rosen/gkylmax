maxwellian_folder="/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average"

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_BiMaxwellianMoments_100.gkyl" gk_lorentzian_mirror-ion_BiMaxwellianMoments_200.gkyl interp pl -f0 --legend "maxwellian,beams" --subplot-titles "density, Upar, Vtpar, Vtperp" &

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl" gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl interp pl -f0 --legend "maxwellian,beams" --subplot-titles "Source density, Source Upar, Source E" &

# pgkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_integrated_moms.gkyl" gk_lorentzian_mirror-ion_integrated_moms.gkyl pl -f0 --legend "maxwellian,beams" --subplot-titles "Integrated M0, integrated M1, integrated M2par, integrated M2perp" &

# pgkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl" gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl pl -f0 --legend "maxwellian,beams" --subplot-titles "bflux M0, bflux M1, bflux M2par, bflux M2perp" --logy &

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl" gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl interp sel -c2 --z0 -0.98:0.98 integ 0 info

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_M0_0.gkyl" gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_M2_0.gkyl" gk_lorentzian_mirror-ion_source_M2_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info

max0=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_M0_0.gkyl" interp sel --z0 -0.98:0.98 integ 0 info | grep -m1 'Maximum:' | awk '{print $3}')
echo "Maxwellian integ M0: $max0"

max0=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "gk_lorentzian_mirror-ion_source_M0_0.gkyl" interp sel --z0 -0.98:0.98 integ 0 info | grep -m1 'Maximum:' | awk '{print $3}')
echo "Beam integ M0: $max0"

max2=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "$maxwellian_folder/gk_lorentzian_mirror-ion_source_M2_0.gkyl" interp sel --z0 -0.98:0.98 integ 0 info | grep -m1 'Maximum:' | awk '{print $3}')
echo "Maxwellian integ M2: $max2"

max2=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "gk_lorentzian_mirror-ion_source_M2_0.gkyl" interp sel --z0 -0.98:0.98 integ 0 info | grep -m1 'Maximum:' | awk '{print $3}')
echo "Beam integ M2: $max2"