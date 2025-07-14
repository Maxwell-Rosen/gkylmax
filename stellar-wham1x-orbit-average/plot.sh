pgkyl "gk_wham-field_[0-9]*.gkyl" interp sel --z0 0.0 col pl --title 'field z = 0.0 POA' --xlabel 't, s' --ylabel 'phi, V' &
pgkyl gk_wham-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl sel -c0 pl --title 'bflux M0 moment xupper modified loss cone POA' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy --scatter &
# pgkyl gk_wham-ion-cflrate_1.gkyl sel --z0 0.0 pl &
# pgkyl gk_wham-ion-cflrate_33.gkyl sel --z0 0.0 pl &

frame=6
# pgkyl gk_wham-ion-cflrate_${frame}.gkyl sel --z0 0.0 pl --title 'cfl z = 0.0'&
# pgkyl gk_wham-ion_${frame}.gkyl interp sel --z0 0.0 ev 'f abs' pl --logz --zmin 1e-25 --title 'z = 0.0 relaxed' &
# pgkyl gk_wham-ion_15.gkyl interp sel --z0 0.0 ev 'f abs' pl --logz --zmin 1e-25 --title 'z = 0.0 old' &
# pgkyl gk_wham-ion_${frame}.gkyl interp sel --z0 0.5 pl --logz --zmin 1e-25 --title 'z = 0.5' &
# pgkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.0 pl --logz --zmin 1e-25 --title 'z = 1.0' &
# pgkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.5 pl --logz --zmin 1e-25 --title 'z = 1.5' &
# pgkyl gk_wham-ion_${frame}.gkyl interp sel --z0 1.9 pl --logz --zmin 1e-25 --title 'z = 1.9' &
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

# pgkyl "gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp col pl --title 'BiMaxwellianMoments' &
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