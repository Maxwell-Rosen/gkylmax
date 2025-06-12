saveLoc="python-plots/gk_wham"

# pgkyl Field/gk_wham-field_[0-9]*.gkyl interp col sel --z1 0.0 pl --saveas "$saveLoc-field-midplane.png" --title "Field at midplane" --xlabel "time, s" --ylabel "Potential, V" --no-show &
# pgkyl Field/gk_wham-field_0.gkyl Field/gk_wham-field_74.gkyl interp pl --saveas "$saveLoc-field-0-74.png" --title "Field 0-74" --xlabel "z" --ylabel "Potential, V" -f0 --no-show &
# pgkyl Field/gk_wham-field_0.gkyl Field/gk_wham-field_74.gkyl interp ev "f[:] grad abs" pl --saveas "$saveLoc-field-grad-0-74.png" --title "Grad abs Field 0-74" --xlabel "z" --ylabel "abs(E), V/m" -f0 --logy --no-show &
# pgkyl misc/gk_wham-ion_integrated_moms.gkyl misc/gk_wham-elc_integrated_moms.gkyl sel -c0 pl -f0 --saveas "$saveLoc-integrated-moms.png" --title "Integrated density" --xlabel "time, s" --ylabel "Density, m^-2" --no-show &

# pgkyl BiMaxwellianMoments/gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl interp col pl --saveas "$saveLoc-ion-BiMaxwellianMoments.png" --title "Ion BiMaxwellian moments" --xlabel "time, s" --no-show &
# pgkyl BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_[0-9]*.gkyl interp col pl --saveas "$saveLoc-elc-BiMaxwellianMoments.png" --title "Electron BiMaxwellian moments" --xlabel "time, s" --no-show &

# pgkyl BiMaxwellianMoments/gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl interp anim --title "Ion BiMaxwellian moments" --float &
# pgkyl BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_[0-9]*.gkyl interp anim --title "Electron BiMaxwellian moments" --float &
# pgkyl "M/gk_wham-ion_source_M0_[0-9]*.gkyl" interp anim --title "Ion source M0 moments" --float &
# pgkyl "M/gk_wham-elc_source_M0_[0-9]*.gkyl" interp anim --title "Electron source M0 moments" --float &

pgkyl misc/gk_wham-ion_source_integrated_moms.gkyl misc/gk_wham-elc_source_integrated_moms.gkyl sel -c0 pl -f0 --saveas "$saveLoc-source-integrated-moms.png" --title "Source integrated density" --xlabel "time, s" -s &

# pgkyl gk_wham-field_[0-9]*.gkyl interp col sel --z1 0.0 pl --title "Field" --xlabel "time, s" --ylabel "Potential, V" --saveas "$saveLoc-field-midplane.png" --no-show &

