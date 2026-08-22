pgkyl "z-scan/288/gk_lorentzian_mirror-field_0.gkyl" -l "IC" \
      "z-scan/144/gk_lorentzian_mirror-field_20.gkyl" -l "$ N_z = $ 144" \
      "z-scan/216/gk_lorentzian_mirror-field_20.gkyl" -l "$ N_z = $ 216" \
      "z-scan/288/gk_lorentzian_mirror-field_20.gkyl" -l "$ N_z = $ 288" \
      "z-scan/432/gk_lorentzian_mirror-field_20.gkyl" -l "$ N_z = $ 432" \
      interp pl --legend --xlabel 'z [m]' --ylabel 'Potential (V)' -f0 --saveas compare-plots/lorentzian1x-z-resolution-scan-field.png --figsize "12,10" --no-show \
      --title "Resolution scan in the direction of $z$ - Field" &

pgkyl "z-scan/288/gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl" -l "IC" \
      "z-scan/144/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_z = $ 144" \
      "z-scan/216/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_z = $ 216" \
      "z-scan/288/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_z = $ 288" \
      "z-scan/432/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_z = $ 432" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --legend --xlabel 'z [m]' -f3 --saveas compare-plots/lorentzian1x-z-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $ z $" &

      # "z-scan/576/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_z = $ 576" \

      
pgkyl "vpar-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl" -l "IC" \
      "vpar-scan/16/gk_lorentzian_mirror-ion_BiMaxwellianMoments_19.gkyl" -l "$ N_{v_\parallel} = $ 16" \
      "vpar-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_19.gkyl" -l "$ N_{v_\parallel} = $ 32" \
      "vpar-scan/48/gk_lorentzian_mirror-ion_BiMaxwellianMoments_19.gkyl" -l "$ N_{v_\parallel} = $ 48" \
      "vpar-scan/64/gk_lorentzian_mirror-ion_BiMaxwellianMoments_19.gkyl" -l "$ N_{v_\parallel} = $ 64" \
      "vpar-scan/96/gk_lorentzian_mirror-ion_BiMaxwellianMoments_19.gkyl" -l "$ N_{v_\parallel} = $ 96" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --legend --xlabel 'z [m]' -f3 --saveas compare-plots/lorentzian1x-vpar-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $ v_\parallel $" &
      
pgkyl "mu-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl" -l "IC" \
      "mu-scan/16/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_\mu = $ 16" \
      "mu-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_\mu = $ 32" \
      "mu-scan/48/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_\mu = $ 48" \
      "mu-scan/64/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_\mu = $ 64" \
      "mu-scan/96/gk_lorentzian_mirror-ion_BiMaxwellianMoments_20.gkyl" -l "$ N_\mu = $ 96" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --legend --xlabel 'z [m]' -f3 --saveas compare-plots/lorentzian1x-mu-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $\mu$" & 
