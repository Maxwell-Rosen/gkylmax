pgkyl "z-scan/200/gk_lorentzian_mirror-field_65.gkyl" -l "$ N_z = $ 200" \
      "z-scan/250/gk_lorentzian_mirror-field_65.gkyl" -l "$ N_z = $ 250" \
      "z-scan/326/gk_lorentzian_mirror-field_65.gkyl" -l "$ N_z = $ 326" \
      "z-scan/400/gk_lorentzian_mirror-field_65.gkyl" -l "$ N_z = $ 400" \
      interp pl --xlabel 'z [m]' --ylabel 'Potential (V)' -f0 --saveas compare-plots/lorentzian1x-z-resolution-scan-field.png --figsize "12,10" --no-show \
      --title "Resolution scan in the direction of $z$ - Field" --legend "N_z = 200, N_z = 250, N_z = 326, N_z = 400"

pgkyl "z-scan/200/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_z = $ 200" \
      "z-scan/250/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_z = $ 250" \
      "z-scan/326/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_z = $ 326" \
      "z-scan/400/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_z = $ 400" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --xlabel 'z [m]' -f3 --saveas compare-plots/lorentzian1x-z-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $ z $" --legend "N_z = 200, N_z = 250, N_z = 326, N_z = 400"

pgkyl "vpar-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_{v_\parallel} = $ 32" \
      "vpar-scan/48/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_{v_\parallel} = $ 48" \
      "vpar-scan/64/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_{v_\parallel} = $ 64" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --xlabel "z [m]" -f3 --saveas compare-plots/lorentzian1x-vpar-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $ v_\parallel $" --legend "N_{v_\parallel} = 32, N_{v_\parallel} = 48, N_{v_\parallel} = 64"

pgkyl "mu-scan/8/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_\mu = $ 8" \
      "mu-scan/16/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_\mu = $ 16" \
      "mu-scan/24/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_\mu = $ 24" \
      "mu-scan/32/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" -l "$ N_\mu = $ 32" \
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --xlabel 'z [m]' -f3 --saveas compare-plots/lorentzian1x-mu-resolution-scan-ion-bimaxwellian-moments.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $\mu$" --legend "N_\mu = 8, N_\mu = 16, N_\mu = 24, N_\mu = 32" --ymin 1e-2
