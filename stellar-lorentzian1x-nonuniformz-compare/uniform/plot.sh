mkdir -p python-plots

pgkyl "gk_lorentzian_mirror-field_[0-9]*.gkyl" interp sel --z0 0.0 col pl --title 'Midplane electrostatic potential' --xlabel 'Time (s)' --ylabel 'Electric Potential Æ (V)' --saveas python-plots/field_potential_z0_vs_time.png --no-show &

pgkyl gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl pl --title 'Upper Boundary Flux' --xlabel 'Time (s)' --logy --scatter --saveas python-plots/ion_bflux_xupper_M0_vs_time.png --no-show --no-legend &

pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl pl --title 'Integrated Ion Moments' --xlabel 'Time (s)' --scatter --saveas python-plots/ion_integrated_moments_vs_time.png --no-show --no-legend &

pgkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp col ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --xlabel 'Time (s)' --ylabel "Z, m" --saveas python-plots/ion_BiMaxwellian_moments_vs_time.png --no-show --no-legend &

frame=$(ls gk_lorentzian_mirror-ion_cflrate_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 0.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.0 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z0.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 0.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.5 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z0.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.0 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.5 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.9 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.9 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.9_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp integ 2 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --logz --zmin 1e-35 --title "Ion Distribution Integrated over ¼ (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_integrated_mu_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp integ 1 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-10 --title "Ion Distribution Integrated over $v_\parallel$ (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_integrated_vpar_frame${frame}.png --no-show &

pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_source_0.gkyl gk_lorentzian_mirror-ion_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Source vs Distribution at $z_2=0.5$ (Frame 0) - WHAM1X' --logz --zmin 1e-25 --saveas python-plots/ion_source_vs_distf_z0.5_frame0.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_source_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Source at $z_2=0.5$ (Frame 0) - WHAM1X' --logz --zmin 1e-25 --saveas python-plots/ion_source_z0.5_frame0.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Distribution at $z_2=0.5$ (Frame 0) - WHAM1X' --logz --zmin 1e-25 --saveas python-plots/ion_distf_z0.5_frame0.png --no-show &

pgkyl gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl interp pl --xlabel 'Axial Position z (m)' --title 'Source Hamiltonian Moments' --logz --zmin 1e-25 --saveas python-plots/ion_source_HamiltonianMoments_frame0.png --no-show &

# Plot the cfl rate at the last frame at z=0.0
pgkyl gk_lorentzian_mirror-ion_cflrate_${frame}.gkyl sel --z0 0.0 pl --title 'CFL Rate at z=0.0' --xlabel '$v_\parallel$ computational' --ylabel '$\mu$ computational' --logz --saveas python-plots/cfl_rate_z0.0_frame${frame}.png --no-show & 
# Plot cfl rate 5 frames ago
pgkyl gk_lorentzian_mirror-ion_cflrate_$(($frame - 5)).gkyl sel --z0 0.0 pl --title 'CFL Rate at z=0.0 (5 frames ago)' --xlabel '$v_\parallel$ computational' --ylabel '$\mu$ computational' --logz --saveas python-plots/cfl_rate_z0.0_frame$(($frame - 5)).png --no-show &

# Save data for the time trace of the field at z= 0.0, 0.5, 1.0, 1.5, and 1.9
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 0.0 col write -s -f "field_time_trace_z0_eq_0p0.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 0.5 col write -s -f "field_time_trace_z0_eq_0p5.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 0.98 col write -s -f "field_time_trace_z0_eq_0p98.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.0 col write -s -f "field_time_trace_z0_eq_1p0.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.5 col write -s -f "field_time_trace_z0_eq_1p5.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.9 col write -s -f "field_time_trace_z0_eq_1p9.gkyl" &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 2.5 col write -s -f "field_time_trace_z0_eq_2p5.gkyl" &

# Plot bimaxwellian moments at the final frame
pgkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_${frame}.gkyl interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --title 'Final Bi-Maxwellian Moments' --saveas python-plots/ion_BiMaxwellianMoments_frame${frame}.png --no-show --no-legend &
