
pgkyl "gk_wham-field_[0-9]*.gkyl" interp sel --z0 0.0 col pl --title 'field z = 0.0' --xlabel 't, s' --ylabel 'phi, V' &
pgkyl "gk_wham-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp col pl --title 'BiMaxwellianMoments' &
# pgkyl "gk_wham-field_[0-9]*.gkyl" interp anim &
 pgkyl gk_wham-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl sel -c0 pl --title 'bflux M0 moment xupper' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy &