mkdir -p python-plots

pgkyl "gk_wham-field_[0-9]*.gkyl" interp sel --z0 0.0 col pl --title 'Midplane electrostatic potential' --xlabel 'Time (s)' --ylabel 'Electric Potential φ (V)' --saveas python-plots/field_potential_z0_vs_time.png --no-show &

pgkyl gk_wham-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl sel -c0 pl --title 'Ion Particle Flux M0 Moment - Upper Boundary Loss Cone (Modified)' --xlabel 'Time (s)' --ylabel 'Particle Flux M0 Moment (particles/s)' --logy --scatter --saveas python-plots/ion_bflux_xupper_M0_vs_time.png --no-show &

pgkyl gk_wham-ion_integrated_moms.gkyl pl --title 'Integrated Ion Moments vs Time - WHAM1X Orbit Averaging' --xlabel 'Time (s)' --ylabel 'Integrated Moments' --scatter --saveas python-plots/ion_integrated_moments_vs_time.png --no-show &

pgkyl "gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp col pl --title 'Ion BiMaxwellian Moments Evolution - WHAM1X Simulation' --xlabel 'Time (s)' --ylabel 'BiMaxwellian Moments' --saveas python-plots/ion_BiMaxwellian_moments_vs_time.png --no-show &

frame=$(ls gk_wham-ion-cflrate_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp sel --z0 0.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.0 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z0.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp sel --z0 0.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.5 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z0.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.0 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.5 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.9 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.9 (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_z1.9_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp integ 2 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --logz --zmin 1e-35 --title "Ion Distribution Integrated over μ (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_integrated_mu_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_wham-ion_mapc2p_vel.gkyl gk_wham-ion_${frame}.gkyl interp integ 1 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-10 --title "Ion Distribution Integrated over $v_\parallel$ (Frame ${frame}) - WHAM1X" --saveas python-plots/ion_distf_integrated_vpar_frame${frame}.png --no-show &


# pgkyl gk_wham-ion_BiMaxwellianMoments_12.gkyl gk_wham-ion_BiMaxwellianMoments_0.gkyl ../initial-conditions/boltz-elc-288z-nu2000/gk_wham-ion_BiMaxwellianMoments_1500.gkyl interp pl -f0

# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_275.gkyl integ 1 pl --title 'cfl OAP integ 1 old' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_275.gkyl integ 2 pl --title 'cfl OAP integ 2 old' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_285.gkyl integ 1 pl --title 'cfl RDP integ 1' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_285.gkyl integ 2 pl --title 'cfl RDP integ 2' --logz &


# pgkyl gk_wham-ion-cflrate_7.gkyl integ 1 pl --title 'cfl OAP integ 1 relaxed' --logz &
# pgkyl gk_wham-ion-cflrate_7.gkyl integ 2 pl --title 'cfl OAP integ 2 relaxed' --logz &

# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_285.gkyl sel --z0 0.95 pl --title 'cfl OAP integ 1 old' --logz &
# pgkyl gk_wham-ion-cflrate_1.gkyl sel --z0 0.95 pl --title 'cfl OAP integ 2 relaxed' --logz &

# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_275.gkyl gk_wham-ion-cflrate_1.gkyl sel --z0 0.0 ev 'f[0] f[1] -' pl --title 'difference' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_wham-ion-cflrate_275.gkyl gk_wham-ion-cflrate_1.gkyl sel --z0 0.0 pl &

# pgkyl "gk_wham-field_[0-9]*.gkyl" interp anim &
# pgkyl good-long-RDP-few-cycles/gk_wham-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl sel -c0 pl --title 'bflux M0 moment xupper true loss cone' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy --scatter &
# pgkyl "../stellar-wham1x-288z-restart-true-collisions/gk_wham-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl" sel -c0 pl --title 'bflux M0 moment xupper true loss cone' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy --scatter --xmin '463e-6'&

# pgkyl gk_wham-ion_BiMaxwellianMoments_0.gkyl gk_wham-ion_BiMaxwellianMoments_790.gkyl "../stellar-wham1x-288z-restart-true-collisions/gk_wham-ion_BiMaxwellianMoments_80.gkyl" interp pl --title 'BiMaxwellianMoments' -f0 --logy &

# pgkyl gk_wham-ion_9.gkyl interp sel --z0 0.0 pl --title 'f at t=0' --logz --zmin 1e-25 &
# pgkyl gk_wham-ion_19.gkyl interp sel --z0 0.0 pl --title 'f at t=19' --logz --zmin 1e-25 &
# pgkyl gk_wham-ion_29.gkyl interp sel --z0 0.0 pl --title 'f at t=29' --logz --zmin 1e-25 &


# pgkyl gk_wham-ion_11.gkyl gk_wham-ion_19.gkyl interp sel --z2 0.5 ev 'f[0] f[1] -' pl --title 'df after OAP' &
# pgkyl gk_wham-ion_20.gkyl gk_wham-ion_21.gkyl interp sel --z2 0.2 ev 'f[0] f[1] -' pl --title 'df after RDP' &


# pgkyl gk_wham-ion_20.gkyl gk_wham-ion_21.gkyl interp integ 2 ev 'f[0] f[1] - abs' pl --title 'df after RDP' --logz --zmin 1e-25&

# pgkyl "gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp sel -c0 anim --logy &