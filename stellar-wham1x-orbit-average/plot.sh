# pgkyl gk_wham-ion-cflrate_1.gkyl sel --z0 0.0 pl &
# pgkyl gk_wham-ion-cflrate_33.gkyl sel --z0 0.0 pl &

pgkyl gk_wham-ion_32.gkyl interp sel --z0 0.0 pl --logz --zmin 1e-25 --title 'z = 0.0' &
pgkyl gk_wham-ion_32.gkyl interp sel --z0 0.5 pl --logz --zmin 1e-25 --title 'z = 0.5' &
pgkyl gk_wham-ion_32.gkyl interp sel --z0 1.0 pl --logz --zmin 1e-25 --title 'z = 1.0' &
pgkyl gk_wham-ion_32.gkyl interp sel --z0 1.5 pl --logz --zmin 1e-25 --title 'z = 1.5' &
pgkyl gk_wham-ion_32.gkyl interp sel --z0 1.9 pl --logz --zmin 1e-25 --title 'z = 1.9' &
# pgkyl gk_wham-ion_BiMaxwellianMoments_12.gkyl gk_wham-ion_BiMaxwellianMoments_0.gkyl ../initial-conditions/boltz-elc-288z-nu2000/gk_wham-ion_BiMaxwellianMoments_1500.gkyl interp pl -f0