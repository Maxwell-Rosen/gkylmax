# pgkyl --c2p-vel 32/misc/gk_wham-elc_mapc2p_vel.gkyl 32/Distributions/gk_wham-elc_30.gkyl interp sel --z0 2 pl --title "32z" --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 64/misc/gk_wham-elc_mapc2p_vel.gkyl 64/Distributions/gk_wham-elc_30.gkyl interp sel --z0 4 pl --title "64z" --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 96/misc/gk_wham-elc_mapc2p_vel.gkyl 96/Distributions/gk_wham-elc_30.gkyl interp sel --z0 6 pl --title "96z" --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 128/misc/gk_wham-elc_mapc2p_vel.gkyl 128/Distributions/gk_wham-elc_30.gkyl interp sel --z0 8 pl --title "128z" --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 192/misc/gk_wham-elc_mapc2p_vel.gkyl 192/Distributions/gk_wham-elc_30.gkyl interp sel --z0 12 pl --title "192z" --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 288/misc/gk_wham-elc_mapc2p_vel.gkyl 288/Distributions/gk_wham-elc_30.gkyl interp sel --z0 18 pl --title "288z" --xlabel "v_par" --ylabel "mu" &

# pgkyl --c2p-vel 32/misc/gk_wham-elc_mapc2p_vel.gkyl 32/Distributions/gk_wham-elc_30.gkyl interp sel --z0 2 ev "f[:] abs" pl --title "32z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 64/misc/gk_wham-elc_mapc2p_vel.gkyl 64/Distributions/gk_wham-elc_30.gkyl interp sel --z0 4 ev "f[:] abs" pl --title "64z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 96/misc/gk_wham-elc_mapc2p_vel.gkyl 96/Distributions/gk_wham-elc_30.gkyl interp sel --z0 6 ev "f[:] abs" pl --title "96z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 128/misc/gk_wham-elc_mapc2p_vel.gkyl 128/Distributions/gk_wham-elc_30.gkyl interp sel --z0 8 ev "f[:] abs" pl --title "128z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 192/misc/gk_wham-elc_mapc2p_vel.gkyl 192/Distributions/gk_wham-elc_30.gkyl interp sel --z0 12 ev "f[:] abs" pl --title "192z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &
# pgkyl --c2p-vel 288/misc/gk_wham-elc_mapc2p_vel.gkyl 288/Distributions/gk_wham-elc_30.gkyl interp sel --z0 18 ev "f[:] abs" pl --title "288z" --logz --zmin 1e-20 --xlabel "v_par" --ylabel "mu" &

pgkyl 32/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "32z" &
pgkyl 64/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "64z" &
pgkyl 96/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "96z" &
pgkyl 128/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "128z" &
pgkyl 192/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "192z" &
pgkyl 288/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_30.gkyl interp pl --title "288z" &