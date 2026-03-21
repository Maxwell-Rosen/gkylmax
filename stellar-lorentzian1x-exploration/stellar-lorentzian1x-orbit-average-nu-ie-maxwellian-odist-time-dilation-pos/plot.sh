mkdir -p python-plots

# Frame settings
last_oap_start_frame=40
end_oap_frame=45
end_fdp_frame=65

# # For some reason, simulations with very weak initial conditions (small density, small f) seem to crash rapidly

# # Strange oscillations
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f : -m -v sel --z0 0.0 anim --float --xlabel "Parallel Velocity $v_\parallel$ (m/s)" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_z_eq_0.mp4" --no-show &

# # Examining the beam in mu
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f : -m -v sel --z2 46 anim --float --xlabel "Z, m" --ylabel "Parallel Velocity $v_\parallel$ (m/s)" --saveas "python-plots/ion_distf_mu_46.mp4" --no-show &
# # Study vpar = 0
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f : -m -v sel --z1 0.0 anim --float --xlabel "Z, m" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_vpar_eq_0.mp4" --no-show&


# # Seems to go away during a long FDP# Strange oscillations
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${last_oap_start_frame}": -m -v sel --z0 0.0 anim --float --xlabel "Parallel Velocity $v_\parallel$ (m/s)" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_z_eq_0_last_phase.mp4" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${last_oap_start_frame}": -m -v sel --z2 46 anim --float --xlabel "Z, m" --ylabel "Parallel Velocity $v_\parallel$ (m/s)" --saveas "python-plots/ion_distf_mu_46_last_phase.mp4" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${last_oap_start_frame}": -m -v sel --z1 0.0 anim --float --xlabel "Z, m" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_vpar_eq_0_last_phase.mp4" --no-show &

# # There is a big difference between the OAP equilibrium and FDP equilibrium

# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_oap_frame}" -m -v sel --z0 0.0 pl --title "Ion Distribution Function at z=0.0 at end of last OAP" --xlabel "Parallel Velocity $v_\parallel$ (m/s)" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_z_eq_0_end_OAP.png" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_fdp_frame}" -m -v sel --z0 0.0 pl --title "Ion Distribution Function at z=0.0 at end of last FDP" --xlabel "Parallel Velocity $v_\parallel$ (m/s)" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_z_eq_0_end_FDP.png" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_oap_frame}" -m -v sel --z2 46 pl --title "Ion Distribution Function mu near beam at end of last OAP" --xlabel "Z, m" --ylabel "Parallel Velocity $v_\parallel$ (m/s)" --saveas "python-plots/ion_distf_mu_46_end_OAP.png" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_fdp_frame}" -m -v sel --z2 46 pl --title "Ion Distribution Function mu near beam at end of last FDP" --xlabel "Z, m" --ylabel "Parallel Velocity $v_\parallel$ (m/s)" --saveas "python-plots/ion_distf_mu_46_end_FDP.png" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_oap_frame}" -m -v sel --z1 0.0 pl --title "Ion Distribution Function at vpar=0.0 at end of last OAP" --xlabel "Z, m" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_vpar_eq_0_end_OAP.png" --no-show &
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f "${end_fdp_frame}" -m -v sel --z1 0.0 pl --title "Ion Distribution Function at vpar=0.0 at end of last FDP" --xlabel "Z, m" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_vpar_eq_0_end_FDP.png" --no-show &

# # Final BiMaxwellian moments
# frame="${end_fdp_frame}"
# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_${frame}.gkyl" interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --xlabel 'Axial Position z (m)' --subplot-ylabels "Density, U_||, T_||, T_perp" --saveas "python-plots/bimaxwellian_moments_frame_${frame}.png" --no-show &

# # movie of BiMaxwellian Moments
# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' anim  --float --saveas "python-plots/bimaxwellian_moments.mp4" --no-show &

# # Movie of integrated mu
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f : -m -v integ 2 anim --float --xlabel "Z, m" --ylabel "Parallel Velocity $v_\parallel$ (m/s)" --saveas "python-plots/ion_distf_integ_mu.mp4" --no-show &
# # Movie of integrated vpar
# pgkyl gk-distf -n gk_lorentzian_mirror -s ion -f : -m -v integ 1 anim --float --xlabel "Z, m" --ylabel "Magnetic Moment $\mu$ (J/T)" --saveas "python-plots/ion_distf_integ_vpar.mp4" --no-show&

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_${end_fdp_frame}.gkyl ../stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil10/gk_lorentzian_mirror-ion_BiMaxwellianMoments_${end_fdp_frame}.gkyl ../stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil15/gk_lorentzian_mirror-ion_BiMaxwellianMoments_${end_fdp_frame}.gkyl ../stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil20/gk_lorentzian_mirror-ion_BiMaxwellianMoments_${end_fdp_frame}.gkyl interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --xlabel 'Axial Position z (m)' --subplot-ylabels "Density, U_||, T_||, T_perp" --legend "1/6,1/10,1/15,1/20" -f0  &


pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl ../stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil20/gk_lorentzian_mirror-ion_BiMaxwellianMoments_90.gkyl interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --xlabel 'Axial Position z (m)' --subplot-ylabels "Density, U_||, T_||, T_perp" --legend "1/6,1/20" -f0&

stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil10
wait