# pgkyl sim-nonunif/Geometry/gk_wham-mc2nu_pos.gkyl interp sel -c2 pl --saveas "plots/nunif-mc2nu.png" &
# pgkyl sim-nonunif/M/gk_wham-ion_M0_0.gkyl sim-unif/M/gk_wham-ion_M0_0.gkyl interp pl -f0 --title "M0_0" --saveas "plots/M0_frame_0.png" &
# pgkyl sim-nonunif/M/gk_wham-ion_M0_1.gkyl sim-unif/M/gk_wham-ion_M0_1.gkyl interp pl -f0 --title "M0_1" --saveas "plots/M0_frame_1.png" &
pgkyl sim-nonunif/Geometry/gk_wham-jacobgeo.gkyl interp pl -f0 --title "non-uniform jacobgeo" --saveas "plots/jacobgeo.png" --xlim -0.1 0.1 &