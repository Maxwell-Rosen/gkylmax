pgkyl gk_wham-elc_M0_0.gkyl interp pl --title "elc M0" --logy --saveas "python-plots/gk_wham-elc_M0_0 logy.png" --no-show &
pgkyl gk_wham-elc_M0_0.gkyl interp pl --title "elc M0" --saveas "python-plots/gk_wham-elc_M0_0.png" --no-show &
pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-elc_M0_144.gkyl interp pl --title "old IC elc M0" --logy --saveas "python-plots/gk_wham-elc_M0_144 logy.png" --no-show &
pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-elc_M0_144.gkyl interp pl --title "old IC elc M0" --saveas "python-plots/gk_wham-elc_M0_144.png" --no-show &

pgkyl gk_wham-field_0.gkyl interp pl --title "field" --logy --saveas "python-plots/gk_wham-field_0 logy.png" --no-show &
pgkyl gk_wham-field_0.gkyl interp pl --title "field" --saveas "python-plots/gk_wham-field_0.png" --no-show &
pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-field_144.gkyl interp pl --title "old IC field" --logy --saveas "python-plots/gk_wham-field_144 logy.png" --no-show &
pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-field_144.gkyl interp pl --title "old IC field" --saveas "python-plots/gk_wham-field_144.png" --no-show &



pgkyl gk_wham-ion_M0_0.gkyl interp pl --title "Ion M0" --logy --saveas "python-plots/gk_wham-ion_M0_0 logy.png" --no-show &
pgkyl gk_wham-ion_M0_0.gkyl interp pl --title "Ion M0" --saveas "python-plots/gk_wham-ion_M0_0.png" --no-show &

pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-ion_M0_144.gkyl interp pl --title "old IC Ion M0" --logy --saveas "python-plots/gk_wham-ion_M0_144 logy.png" --no-show &
pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-ion_M0_144.gkyl interp pl --title "old IC I M0" --saveas "python-plots/gk_wham-ion_M0_144.png" --no-show &
# pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-elc_M0_144.gkyl interp pl --title "old IC Electron M0" &

# pgkyl gk_wham-jacobgeo.gkyl interp pl --title "Jacobgeo" &
# pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-jacobgeo.gkyl interp pl --title "old IC Jacobgeo" &


# pgkyl gk_wham-ion_0.gkyl  interp sel --z0 0.0 pl --title "Ion Distribution at z=0.0" &
# pgkyl gk_wham-elc_0.gkyl interp sel --z0 0.0 pl --title "Electron Distribution at z=0.0" &
# pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-ion_144.gkyl interp sel --z0 0.0 pl --title "old IC Ion Distribution at z=0.0" &

# pgkyl gk_wham-ion_0.gkyl  interp sel --z0 1.0 pl --title "Ion Distribution at z=1.0" &
# pgkyl gk_wham-elc_0.gkyl interp sel --z0 0.0 pl --title "Electron Distribution at z=0.0" &
# pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-ion_144.gkyl interp sel --z0 1.5707963268 pl --title "old IC Ion Distribution at z=1.5707" &


# pgkyl ../initial-conditions/kinet-elc-288z-nu2000/old-geom/gk_wham-elc_144.gkyl interp sel --z0 0.0 pl --title "old IC Electron Distribution at z=0.0" &

# pgkyl gk_wham-ion_0.gkyl interp sel --z0 1.0 pl --title "Ion Distribution at z=1.0" &
# pgkyl gk_wham-elc_0.gkyl interp sel --z0 1.0 pl --title "Electron Distribution at z=1.0" &
# pgkyl gk_wham-ion_0.gkyl interp sel --z0 1.8 pl --title "Ion Distribution at z=1.8" &
# pgkyl gk_wham-elc_0.gkyl interp sel --z0 1.8 pl --title "Electron Distribution at z=1.8" &